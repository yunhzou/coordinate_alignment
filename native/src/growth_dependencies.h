// Conditional dependencies of one complete fragment-growth operation.
// The caller separately fixes target occupancy (hence its point stabilizer).
struct GrowthDependencies {
    std::vector<uint8_t> mapping_reads, boundary_atoms;
    std::vector<std::pair<int, std::vector<int>>> island_reads;
    std::vector<Pair> local_deferred;

    explicit GrowthDependencies(int n) : mapping_reads(n, 0), boundary_atoms(n, 0) {}
    void mapping_at(int atom) { mapping_reads[atom] = 1; }
    void island_at(int atom, const std::vector<Pair>* islands, const std::vector<int>& island_of) {
        std::vector<int> members;
        if (islands && island_of[atom] >= 0)
            for (auto [r, label] : *islands)
                if (label == island_of[atom]) { members.push_back(r); mapping_at(r); }
        for (const auto& item : island_reads) if (item.first == atom) return;
        island_reads.push_back({atom, std::move(members)});
    }
    void boundary(const std::vector<uint8_t>& fragment) {
        for (size_t r=0; r<fragment.size(); ++r) boundary_atoms[r] |= fragment[r];
    }
};
