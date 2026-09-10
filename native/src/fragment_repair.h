// Reuse complete conditional fragment results across changed earlier histories.
// This is NOT witness-only repair: the whole compressed candidate result must
// have identical constraints. Different histories can pass that dependency test.
class FragmentRepair {
    Source base, active_source;
    const Target& target;
    std::vector<Pair> active_cuts;
    size_t budget, resident=0, peak_resident=0;
    std::unique_ptr<ExtensionCache> extensions;
    struct Entry {
        std::string key;
        GrowResult result;
        GrowthDependencies dependencies;
        std::vector<int> mapping;
        std::vector<Pair> islands, cuts, deferred, pair_reads;
        std::vector<uint8_t> row_reads;
        size_t bytes=0;
        Entry(int n) : dependencies(n) {}
    };
    using Entries=std::list<Entry>;
    Entries lru;
    std::unordered_map<std::string,std::vector<Entries::iterator>> cache;
    long calls=0, hits=0, history_hits=0, topology_hits=0, boundary_rebases=0;
    long checks=0, evictions=0, mapping_misses=0, island_misses=0, topology_misses=0, boundary_misses=0;
    long logical_extensions=0, reused_extensions=0, logical_certificates=0, reused_certificates=0;

    static std::vector<Pair> normalized(std::vector<Pair> pairs) {
        for (auto& [a,b]:pairs) if (a>b) std::swap(a,b);
        std::sort(pairs.begin(),pairs.end());
        pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
        return pairs;
    }
    static std::vector<Pair> difference(const std::vector<Pair>& a,const std::vector<Pair>& b) {
        std::vector<Pair> out;
        std::set_symmetric_difference(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(out));
        return out;
    }
    template<class Value> static void append(std::string& key,const Value& value) {
        key.append(reinterpret_cast<const char*>(&value),sizeof(value));
    }
    static size_t entry_bytes(const Entry& e) {
        size_t out=sizeof(e)+e.key.capacity()+256+e.mapping.capacity()*sizeof(int)
            +e.row_reads.capacity()+e.dependencies.mapping_reads.capacity()+e.dependencies.boundary_atoms.capacity();
        for (const auto* v:{&e.islands,&e.cuts,&e.deferred,&e.pair_reads,&e.dependencies.local_deferred})
            out+=v->capacity()*sizeof(Pair);
        for (const auto& [r,members]:e.dependencies.island_reads)
            out+=sizeof(r)+sizeof(members)+members.capacity()*sizeof(int);
        for (const auto& iso:e.result.isos) {
            out+=sizeof(iso)+(iso.mapping.capacity()+iso.deferred_edges.capacity())*sizeof(Pair)
                +iso.fragment.capacity()*sizeof(int);
            const auto& c=iso.cand;
            out+=(c.img.capacity()+c.mapped.capacity()+c.exact_fixed.capacity())*sizeof(int);
            for (const auto* blocks:{&c.blocks,&c.automorph}) {
                out+=blocks->capacity()*sizeof(Block);
                for (const auto& b:*blocks) out+=(b.r.capacity()+b.p.capacity())*sizeof(int);
            }
        }
        return out;
    }
    void insert(Entry entry) {
        entry.bytes=entry_bytes(entry);
        if (entry.bytes>budget) return;
        while (!lru.empty() && resident+entry.bytes>budget) {
            auto old=lru.begin();
            auto bucket=cache.find(old->key);
            auto& entries=bucket->second;
            entries.erase(std::remove(entries.begin(),entries.end(),old),entries.end());
            if (entries.empty()) cache.erase(bucket);
            resident-=old->bytes;lru.erase(old);++evictions;
        }
        resident+=entry.bytes;
        peak_resident=std::max(peak_resident,resident);
        lru.push_back(std::move(entry));
        auto added=std::prev(lru.end());cache[added->key].push_back(added);
    }
    bool valid(const Entry& entry,const std::vector<int>& mapping,
               const std::vector<Pair>& islands,const std::vector<Pair>& cuts,
               const std::vector<Pair>& deferred) {
        ++checks;
        for (auto [a,b]:difference(entry.cuts,cuts)) {
            if (entry.row_reads[a] || entry.row_reads[b] ||
                std::binary_search(entry.pair_reads.begin(),entry.pair_reads.end(),Pair{a,b})) {
                ++topology_misses;return false;
            }
        }
        for (int r=0;r<base.n;++r) if (entry.dependencies.mapping_reads[r] && entry.mapping[r]!=mapping[r]) {
            ++mapping_misses;return false;
        }
        for (const auto& [r,old_members]:entry.dependencies.island_reads) {
            int label=-1;
            for (auto [atom,iid]:islands) if (atom==r) {label=iid;break;}
            std::vector<int> members;
            if (label>=0) for (auto [atom,iid]:islands) if (iid==label) members.push_back(atom);
            if (members!=old_members) {++island_misses;return false;}
        }
        for (auto [a,b]:difference(entry.deferred,deferred))
            if (entry.dependencies.boundary_atoms[a] || entry.dependencies.boundary_atoms[b]) {
                ++boundary_misses;return false;
            }
        return true;
    }
public:
    FragmentRepair(const PySource& source,const PyTarget& product,size_t cache_bytes,size_t extension_bytes=0)
        :base(source.g),active_source(source.g),target(product.g),budget(cache_bytes) {
        base.reads=nullptr;active_source.reads=nullptr;
        if (extension_bytes) extensions=std::make_unique<ExtensionCache>(extension_bytes);
    }
    std::pair<GrowResult,long> compute(const std::vector<Pair>& cuts_in,int seed,const std::vector<int>& mapping,
                    double graph_floor,double iso_tol,int min_lock_size,long max_branches,
                    const std::vector<Pair>* islands_ptr,const std::vector<Pair>& prior_deferred,bool allow_mapped_seed) {
        const std::vector<Pair> empty_islands;
        const auto& islands=islands_ptr ? *islands_ptr : empty_islands;
        GrowResult result;
        long saved_extensions=0;
        {
            ++calls;
            auto cuts=normalized(cuts_in),deferred=normalized(prior_deferred);
            if (cuts!=active_cuts) {
                for (auto [a,b]:cuts)
                    if (a<0 || b<0 || a>=base.n || b>=base.n || !base.has_edge(a,b))
                        throw std::invalid_argument("cut is not an active edge of the repair source");
                active_source=base;active_cuts=cuts;
                for (auto [a,b]:cuts) {
                    active_source.edge[(size_t)a*base.n+b]=active_source.edge[(size_t)b*base.n+a]=0;
                    for (auto [x,y]:{Pair{a,b},Pair{b,a}}) {
                        auto& neighbors=active_source.adj[x];
                        neighbors.erase(std::remove(neighbors.begin(),neighbors.end(),y),neighbors.end());
                    }
                }
            }
            std::string key;
            append(key,seed);append(key,graph_floor);append(key,iso_tol);
            append(key,min_lock_size);append(key,max_branches);append(key,allow_mapped_seed);
            std::vector<uint8_t> occupied(target.n,0);
            for (int p:mapping) if (p>=0) occupied[p]=1;
            key.append(reinterpret_cast<const char*>(occupied.data()),occupied.size());
            // Occupancy fixes the same pointwise target stabilizer. Relabeling
            // untouched locked R atoms changes color names, not the equivalence
            // relation between candidates. Touched mapping values are checked.
            auto bucket=cache.find(key);
            auto chosen=lru.end();
            if (bucket!=cache.end()) for (auto item:bucket->second)
                if (valid(*item,mapping,islands,cuts,deferred)) {chosen=item;break;}
            if (chosen!=lru.end()) {
                ++hits;
                history_hits+=(chosen->mapping!=mapping || chosen->islands!=islands);
                topology_hits+=(chosen->cuts!=cuts);
                result=chosen->result;
                if (chosen->deferred!=deferred) {
                    ++boundary_rebases;
                    auto rebased=deferred;
                    rebased.insert(rebased.end(),chosen->dependencies.local_deferred.begin(),chosen->dependencies.local_deferred.end());
                    rebased=normalized(std::move(rebased));
                    for (auto& iso:result.isos) iso.deferred_edges=rebased;
                }
                saved_extensions=result.profile.extend_calls;
                reused_extensions+=saved_extensions;
                reused_certificates+=result.profile.certificate_calls;
                lru.splice(lru.end(),lru,chosen);
            } else {
                Entry entry(base.n);
                GrowthTrace trace(active_source);trace.store_checkpoints=false;
                result=grow_island(active_source,target,seed,mapping,graph_floor,iso_tol,min_lock_size,
                    max_branches,islands_ptr,prior_deferred,allow_mapped_seed,&trace,nullptr,-1,&entry.dependencies,extensions.get());
                entry.key=std::move(key);entry.result=result;entry.mapping=mapping;
                entry.islands=islands;entry.cuts=cuts;entry.deferred=deferred;
                entry.row_reads.resize(base.n);
                for (int r=0;r<base.n;++r) entry.row_reads[r]=(trace.reads.rows[r]!=-2);
                for (int a=0;a<base.n;++a) for (int b=a+1;b<base.n;++b)
                    if (!entry.row_reads[a] && !entry.row_reads[b] && trace.reads.pairs[(size_t)a*base.n+b]!=-2)
                        entry.pair_reads.push_back({a,b});
                insert(std::move(entry));
            }
            logical_extensions+=result.profile.extend_calls;
            logical_certificates+=result.profile.certificate_calls;
        }
        return {std::move(result),saved_extensions};
    }
    py::object grow(const std::vector<Pair>& cuts,int seed,const std::vector<int>& mapping,
                    double graph_floor,double iso_tol,int min_lock_size,long max_branches,
                    py::object islands_obj,const std::vector<Pair>& deferred,bool allow_mapped_seed) {
        std::vector<Pair> islands;
        const std::vector<Pair>* islands_ptr=nullptr;
        if (!islands_obj.is_none()) {islands=islands_obj.cast<std::vector<Pair>>();islands_ptr=&islands;}
        std::pair<GrowResult,long> computed;
        {
            py::gil_scoped_release release;
            computed=compute(cuts,seed,mapping,graph_floor,iso_tol,min_lock_size,max_branches,
                             islands_ptr,deferred,allow_mapped_seed);
        }
        auto out=grow_result_dict(computed.first).cast<py::dict>();
        out["reused_extensions"]=computed.second;
        return out;
    }
    py::dict stats() const {
        py::dict out;
        out["calls"]=calls;out["full_hits"]=hits;out["prefix_hits"]=0;
        out["history_hits"]=history_hits;out["topology_hits"]=topology_hits;out["boundary_rebases"]=boundary_rebases;
        out["dependency_checks"]=checks;out["mapping_misses"]=mapping_misses;out["island_misses"]=island_misses;
        out["topology_misses"]=topology_misses;out["boundary_misses"]=boundary_misses;
        out["logical_extensions"]=logical_extensions;out["reused_extensions"]=reused_extensions;
        out["logical_certificates"]=logical_certificates;out["reused_certificates"]=reused_certificates;
        out["entries"]=lru.size();out["resident_bytes"]=resident;out["peak_resident_bytes"]=peak_resident;
        out["cache_budget_bytes"]=budget;out["evictions"]=evictions;
        if (extensions) out.attr("update")(extensions->stats());
        return out;
    }
};

void register_fragment_repair(py::module_& mod) {
    py::class_<FragmentRepair>(mod,"FragmentRepair")
        .def(py::init<const PySource&,const PyTarget&,size_t,size_t>(),
            py::arg("source"),py::arg("target"),py::arg("cache_bytes"),py::arg("extension_cache_bytes")=0,py::keep_alive<1,3>())
        .def("grow",&FragmentRepair::grow,py::arg("cuts"),py::arg("seed"),py::arg("mapping"),
             py::arg("graph_floor"),py::arg("iso_tol"),py::arg("min_lock_size"),py::arg("max_branches"),
             py::arg("islands"),py::arg("prior_deferred_edges"),py::arg("allow_mapped_seed"))
        .def("stats",&FragmentRepair::stats);
}
