// Experimental exact source-cut reuse. No growth/branching decisions change.
// The session owns a fixed source template and uses one immutable target.
class GrowthReplay {
    Source base;
    Source active_source;
    std::vector<Pair> active_cuts;
    const Target& target;
    size_t budget, resident = 0;
    bool checkpoints, reference_only;
    struct Entry {
        std::string key;
        std::shared_ptr<GrowthTrace> trace;
        size_t bytes;
    };
    using Entries = std::list<Entry>;
    Entries lru;
    std::unordered_map<std::string, std::vector<Entries::iterator>> cache;
    long calls=0, full_hits=0, prefix_hits=0, evictions=0;
    long source_rebuilds=0, recorded_calls=0;
    long skipped_extensions=0, total_extensions=0, skipped_certificates=0, total_certificates=0;
    size_t peak_resident=0;

    static void append(std::string& key, int64_t value) {
        key.append(reinterpret_cast<const char*>(&value), sizeof(value));
    }
    static void append_real(std::string& key, double value) {
        key.append(reinterpret_cast<const char*>(&value), sizeof(value));
    }
    void insert(std::string key, std::shared_ptr<GrowthTrace> trace) {
        size_t bytes = trace->bytes()+key.capacity()+sizeof(Entry)+128;
        if (bytes>budget) return; // cache residency only; never limits the search
        while (!lru.empty() && resident+bytes>budget) {
            auto item=lru.begin();
            auto found=cache.find(item->key);
            auto& group=found->second;
            group.erase(std::remove(group.begin(),group.end(),item),group.end());
            if (group.empty()) cache.erase(found);
            resident-=item->bytes;
            lru.erase(item);
            ++evictions;
        }
        lru.push_back({std::move(key),std::move(trace),bytes});
        auto item=std::prev(lru.end());
        cache[item->key].push_back(item);
        resident+=bytes;
        peak_resident=std::max(peak_resident,resident);
    }
public:
    GrowthReplay(const PySource& source, const PyTarget& product, size_t cache_bytes,
                 bool checkpoint_replay, bool uncut_reference_only)
        : base(source.g),active_source(source.g),target(product.g),budget(cache_bytes),
          checkpoints(checkpoint_replay),reference_only(uncut_reference_only) {
        base.reads=nullptr;
        active_source.reads=nullptr;
    }

    py::object grow(const std::vector<Pair>& cuts, int seed, const std::vector<int>& mapping,
                    double graph_floor, double iso_tol, int min_lock_size, long max_branches,
                    py::object islands_obj, const std::vector<Pair>& prior_deferred,
                    bool allow_mapped_seed) {
        std::vector<Pair> islands;
        const std::vector<Pair>* islands_ptr=nullptr;
        if (!islands_obj.is_none()) { islands=islands_obj.cast<std::vector<Pair>>(); islands_ptr=&islands; }
        GrowResult result;
        long reused_extensions=0, reused_certificates=0;
        {
            py::gil_scoped_release release;
            ++calls;
            auto normalized_cuts=cuts;
            for (auto& e:normalized_cuts) if (e.first>e.second) std::swap(e.first,e.second);
            std::sort(normalized_cuts.begin(),normalized_cuts.end());
            normalized_cuts.erase(std::unique(normalized_cuts.begin(),normalized_cuts.end()),normalized_cuts.end());
            if (normalized_cuts!=active_cuts) {
                for (auto [a,b] : normalized_cuts)
                    if (a<0 || b<0 || a>=base.n || b>=base.n || !base.has_edge(a,b))
                        throw std::invalid_argument("cut is not an active edge of the replay source");
                ++source_rebuilds;
                active_source=base;
                active_cuts=normalized_cuts;
                for (auto [a,b] : active_cuts) {
                    active_source.edge[(size_t)a*base.n+b]=active_source.edge[(size_t)b*base.n+a]=0;
                    for (auto [x,y] : {Pair{a,b},Pair{b,a}}) {
                        auto& neighbors=active_source.adj[x];
                        neighbors.erase(std::remove(neighbors.begin(),neighbors.end(),y),neighbors.end());
                    }
                }
            }
            const Source& source=active_source;
            std::string key;
            append(key,seed); append_real(key,graph_floor); append_real(key,iso_tol);
            append(key,min_lock_size); append(key,max_branches); append(key,allow_mapped_seed);
            append(key,mapping.size()); for (int image:mapping) append(key,image);
            append(key,islands_ptr ? (int64_t)islands.size() : -1);
            for (auto [a,b]:islands) { append(key,a);append(key,b); }
            auto deferred=prior_deferred;
            for (auto& e:deferred) if (e.first>e.second) std::swap(e.first,e.second);
            std::sort(deferred.begin(),deferred.end());
            deferred.erase(std::unique(deferred.begin(),deferred.end()),deferred.end());
            append(key,deferred.size()); for (auto [a,b]:deferred) {append(key,a);append(key,b);}

            std::shared_ptr<GrowthTrace> previous;
            int resume=-1;
            bool full=false;
            auto found=cache.find(key);
            auto chosen=lru.end();
            if (found!=cache.end()) for (auto entry:found->second) {
                const auto& candidate=entry->trace;
                int prefix=candidate->reusable_prefix(source);
                if (prefix==(int)candidate->checkpoints.size()) {
                    previous=candidate; chosen=entry; full=true; break;
                }
                if (checkpoints && prefix>=0 &&
                    (!previous || candidate->checkpoints[prefix]->profile.extend_calls>reused_extensions)) {
                    previous=candidate; chosen=entry; resume=prefix;
                    reused_extensions=candidate->checkpoints[prefix]->profile.extend_calls;
                    reused_certificates=candidate->checkpoints[prefix]->certificate_calls;
                }
            }
            if (chosen!=lru.end()) lru.splice(lru.end(),lru,chosen);
            if (full) {
                ++full_hits;
                result=previous->result;
                reused_extensions=result.profile.extend_calls;
                reused_certificates=result.profile.certificate_calls;
            } else {
                if (previous) ++prefix_hits;
                std::shared_ptr<GrowthTrace> trace;
                if (!reference_only || cuts.empty()) {
                    ++recorded_calls;
                    trace=std::make_shared<GrowthTrace>(source);
                    trace->store_checkpoints=checkpoints;
                }
                result=grow_island(source,target,seed,mapping,graph_floor,iso_tol,min_lock_size,
                                   max_branches,islands_ptr,prior_deferred,allow_mapped_seed,
                                   trace.get(),previous.get(),resume);
                if (trace) {
                    trace->result=result;
                    insert(std::move(key),std::move(trace));
                }
            }
            total_extensions+=result.profile.extend_calls;
            total_certificates+=result.profile.certificate_calls;
            skipped_extensions+=reused_extensions;
            skipped_certificates+=reused_certificates;
        }
        auto out=grow_result_dict(result).cast<py::dict>();
        out["reused_extensions"]=reused_extensions;
        return out;
    }
    py::dict stats() const {
        py::dict out;
        out["calls"]=calls; out["full_hits"]=full_hits; out["prefix_hits"]=prefix_hits;
        out["logical_extensions"]=total_extensions; out["reused_extensions"]=skipped_extensions;
        out["logical_certificates"]=total_certificates; out["reused_certificates"]=skipped_certificates;
        out["entries"]=lru.size(); out["resident_bytes"]=resident; out["peak_resident_bytes"]=peak_resident;
        out["cache_budget_bytes"]=budget; out["evictions"]=evictions;
        out["source_rebuilds"]=source_rebuilds; out["recorded_calls"]=recorded_calls;
        return out;
    }
};

void register_growth_replay(py::module_& mod) {
    py::class_<GrowthReplay>(mod,"GrowthReplay")
        .def(py::init<const PySource&,const PyTarget&,size_t,bool,bool>(),
             py::arg("source"),py::arg("target"),py::arg("cache_bytes"),
             py::arg("checkpoint_replay")=true,py::arg("uncut_reference_only")=false,
             py::keep_alive<1,3>())
        .def("grow",&GrowthReplay::grow,py::arg("cuts"),py::arg("seed"),py::arg("mapping"),
             py::arg("graph_floor"),py::arg("iso_tol"),py::arg("min_lock_size"),
             py::arg("max_branches"),py::arg("islands"),py::arg("prior_deferred_edges"),
             py::arg("allow_mapped_seed"))
        .def("stats",&GrowthReplay::stats);
}
