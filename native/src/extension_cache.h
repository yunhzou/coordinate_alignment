// Exact per-candidate extension reuse. Scoped to one fixed source-WBO/target
// pair; source cuts may change. Dedupe, admission and caps still execute normally.
class ExtensionCache {
    struct Key {
        std::shared_ptr<const CandidateCacheKey> candidate;
        std::string context;
        size_t hash;
        bool operator==(const Key& other) const {
            return context==other.context && (candidate==other.candidate || candidate->value==other.candidate->value);
        }
    };
    struct KeyHash { size_t operator()(const Key& key) const {return key.hash;} };
    struct Entry { Key key; std::vector<Cand> children; size_t bytes; };
    std::list<Entry> lru;
    std::unordered_map<Key, std::list<Entry>::iterator, KeyHash> entries;
    size_t budget, resident=0, peak=0;
    long calls=0, hits=0, negative_hits=0, evictions=0;

    template<class T> static void scalar(std::string& key,const T& value) {
        key.append(reinterpret_cast<const char*>(&value),sizeof(value));
    }
    template<class T> static void vector(std::string& key,const std::vector<T>& values) {
        scalar(key,values.size());
        key.append(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(T));
    }
    static void blocks(std::string& key,const std::vector<Block>& values) {
        scalar(key,values.size());
        for (const auto& b:values) {vector(key,b.r);vector(key,b.p);scalar(key,b.extendable);}
    }
    static size_t bytes(const Cand& c) {
        size_t n=sizeof(c)+(c.img.capacity()+c.mapped.capacity()+c.exact_fixed.capacity())*sizeof(int);
        for (const auto* bs:{&c.blocks,&c.automorph}) {
            n+=bs->capacity()*sizeof(Block);
            for (const auto& b:*bs) n+=(b.r.capacity()+b.p.capacity())*sizeof(int);
        }
        return n;
    }
    static Key key_for(const Cand& c,const Context& ctx) {
        if (!c.extension_key) {
            std::string value;
            vector(value,c.img);blocks(value,c.blocks);vector(value,c.exact_fixed);
            scalar(value,c.mult);blocks(value,c.automorph);
            size_t hash=std::hash<std::string>()(value);
            c.extension_key=std::make_shared<CandidateCacheKey>(CandidateCacheKey{std::move(value),hash});
        }
        std::string key;
        scalar(key,ctx.n);scalar(key,ctx.iso_tol);vector(key,ctx.fragment_old);
        vector(key,ctx.bonded_in_frag);vector(key,ctx.r_wbos);
        scalar(key,ctx.strict_r);scalar(key,ctx.strict_w);
        const bool merge=ctx.is_merge();scalar(key,merge);
        if (merge) {
            vector(key,ctx.island_atoms);
            scalar(key,(*ctx.mapping)[ctx.n]);
            for (int r:ctx.island_atoms) scalar(key,(*ctx.mapping)[r]);
            // Every source topology test in island_merge_wbo_consistent.
            for (int r:ctx.island_atoms) for (int a:c.mapped)
                scalar(key,ctx.R->has_edge(r,a));
            // Added island atoms are also included in the merged child's mapping.
            for (int r:ctx.island_atoms) for (int a:ctx.island_atoms)
                scalar(key,ctx.R->has_edge(r,a));
        } else {
            // Only targets of n's element are read by the free extension.
            int code=ctx.R->ecode[ctx.n];scalar(key,code);
            if (code>=0 && code<(int)ctx.P->same_element.size())
                for (int p:ctx.P->same_element[code]) scalar(key,ctx.locked_p[p]);
        }
        // All topology inputs of r_compatible_with_block, including failed tests.
        // Calling has_edge while constructing a hit key also records dependencies
        // for the enclosing whole-fragment repair cache.
        for (const auto& b:c.blocks) if (b.open()) {
            for (int a=0;a<ctx.R->n;++a) if (ctx.fragment_old[a] && c.has(a) && !contains(b.r,a)) {
                scalar(key,ctx.R->has_edge(ctx.n,a));
                for (int r:b.r) scalar(key,ctx.R->has_edge(r,a));
            }
        }
        size_t hash=c.extension_key->hash ^ (std::hash<std::string>()(key)<<1);
        return {c.extension_key,std::move(key),hash};
    }
public:
    explicit ExtensionCache(size_t cache_bytes):budget(cache_bytes) {}
    void extend(const Cand& c,const Context& ctx,std::vector<Cand>& out) {
        ++calls;
        auto key=key_for(c,ctx);
        auto found=entries.find(key);
        if (found!=entries.end()) {
            ++hits;negative_hits+=found->second->children.empty();
            out.insert(out.end(),found->second->children.begin(),found->second->children.end());
            lru.splice(lru.end(),lru,found->second);return;
        }
        std::vector<Cand> children;
        if (ctx.is_merge()) {
            Cand child;
            if (extend_locked_merge(c,ctx,child)) children.push_back(std::move(child));
        } else extend_free_atom(c,ctx,children);
        size_t size=sizeof(Entry)+sizeof(Key)+2*key.context.capacity()+key.candidate->value.capacity()
            +sizeof(CandidateCacheKey)+192+children.capacity()*sizeof(Cand);
        for (const auto& child:children) size+=bytes(child);
        if (size<=budget) {
            while (!lru.empty() && resident+size>budget) {
                resident-=lru.front().bytes;entries.erase(lru.front().key);lru.pop_front();++evictions;
            }
            resident+=size;peak=std::max(peak,resident);
            lru.push_back({key,children,size});entries.emplace(std::move(key),std::prev(lru.end()));
        }
        out.insert(out.end(),std::make_move_iterator(children.begin()),std::make_move_iterator(children.end()));
    }
    py::dict stats() const {
        py::dict d;
        d["extension_calls"]=calls;d["extension_hits"]=hits;d["extension_negative_hits"]=negative_hits;
        d["extension_evictions"]=evictions;d["extension_resident_bytes"]=resident;
        d["extension_peak_bytes"]=peak;d["extension_budget_bytes"]=budget;
        return d;
    }
};
