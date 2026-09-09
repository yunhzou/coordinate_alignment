// Included inside engine.cpp's implementation namespace.
struct GrowthCheckpoint {
    std::vector<Cand> cands;
    std::vector<uint8_t> fragment, used_edge;
    int fragment_size;
    std::vector<Pair> deferred;
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> heap;
    GrowProfile profile;
    long certificate_calls;

    size_t bytes() const {
        size_t size = sizeof(*this) + fragment.capacity() + used_edge.capacity()
            + deferred.capacity()*sizeof(Pair) + heap.size()*sizeof(HeapItem)
            + cands.capacity()*sizeof(Cand);
        for (const auto& c : cands) {
            size += (c.img.capacity()+c.mapped.capacity()+c.exact_fixed.capacity())*sizeof(int);
            for (const auto* blocks : {&c.blocks, &c.automorph}) {
                size += blocks->capacity()*sizeof(Block);
                for (const auto& b : *blocks) size += (b.r.capacity()+b.p.capacity())*sizeof(int);
            }
        }
        return size;
    }
};

struct GrowthTrace {
    SourceReads reads;
    std::vector<uint8_t> source_edges;
    std::vector<std::shared_ptr<const GrowthCheckpoint>> checkpoints;
    GrowResult result;
    bool store_checkpoints = true;
    explicit GrowthTrace(const Source& source) : reads(source.n), source_edges(source.edge) {}

    // -1: initialization already differs; size(): the whole result is reusable.
    int reusable_prefix(const Source& source) const {
        int first = (int)checkpoints.size();
        for (int a=0; a<source.n; ++a) for (int b=a+1; b<source.n; ++b) {
            size_t pair = (size_t)a*source.n+b;
            if (source_edges[pair] == source.edge[pair]) continue;
            for (int step : {reads.rows[a], reads.rows[b], reads.pairs[pair]})
                if (step != -2) first = std::min(first, step);
        }
        return first;
    }
    void inherit_prefix(const GrowthTrace& other, int step) {
        checkpoints.assign(other.checkpoints.begin(), other.checkpoints.begin()+step);
        for (size_t i=0; i<reads.rows.size(); ++i)
            reads.rows[i] = other.reads.rows[i] < step ? other.reads.rows[i] : -2;
        for (size_t i=0; i<reads.pairs.size(); ++i)
            reads.pairs[i] = other.reads.pairs[i] < step ? other.reads.pairs[i] : -2;
    }
    size_t bytes() const {
        size_t size = sizeof(*this)+source_edges.capacity()
            +(reads.rows.capacity()+reads.pairs.capacity())*sizeof(int)
            +checkpoints.capacity()*sizeof(void*);
        // Shared checkpoints are conservatively charged to each cache entry.
        for (const auto& state : checkpoints) size += state->bytes();
        for (const auto& iso : result.isos) {
            size += sizeof(iso)+(iso.mapping.capacity()+iso.deferred_edges.capacity())*sizeof(Pair)
                +iso.fragment.capacity()*sizeof(int);
            const auto& c=iso.cand;
            size += (c.img.capacity()+c.mapped.capacity()+c.exact_fixed.capacity())*sizeof(int);
            for (const auto* blocks : {&c.blocks,&c.automorph}) {
                size += blocks->capacity()*sizeof(Block);
                for (const auto& b : *blocks) size += (b.r.capacity()+b.p.capacity())*sizeof(int);
            }
        }
        return size;
    }
};
