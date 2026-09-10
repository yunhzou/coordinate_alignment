// One conditional greedy growth, plus lazily finalized earlier frontiers.
class FragmentChoices {
    const Source& source;
    const Target& target;
    std::vector<int> mapping;
    double tolerance;
    long cap;
    GrowResult greedy;
    std::vector<FragmentFrontier> frontiers;
public:
    FragmentChoices(const PySource& r,const PyTarget& p,int seed,std::vector<int> locked,
                    double floor,double tol,long limit,const std::vector<Pair>& islands,
                    const std::vector<Pair>& deferred)
        :source(r.g),target(p.g),mapping(std::move(locked)),tolerance(tol),cap(limit) {
        py::gil_scoped_release release;
        greedy=grow_island(source,target,seed,mapping,floor,tol,1,cap,&islands,deferred,false,
                           nullptr,nullptr,-1,nullptr,nullptr,&frontiers);
    }
    py::object normal() const {return grow_result_dict(greedy);}
    py::list options() const {
        py::list out;
        for (size_t i=0;i<frontiers.size();++i) {
            const auto& f=frontiers[i];
            out.append(py::make_tuple(i,std::count(f.fragment.begin(),f.fragment.end(),uint8_t{1}),
                                      f.next_atom,f.cands.size()));
        }
        return out;
    }
    py::object close(size_t index) const {
        GrowResult result;
        {
            py::gil_scoped_release release;
            const auto& f=frontiers.at(index);
            // Deferred boundary evidence prevents quotienting away placements
            // with different interactions with the not-yet-matched atoms.
            auto boundary=f.deferred;
            for (int a=0;a<source.n;++a) if (f.fragment[a])
                for (int b:source.adj[a]) if (!f.fragment[b])
                    boundary.push_back({std::min(a,b),std::max(a,b)});
            std::sort(boundary.begin(),boundary.end());
            boundary.erase(std::unique(boundary.begin(),boundary.end()),boundary.end());
            Canonicalizer canon(&target,mapping);long certificates=0;
            finish_fragment(source,target,mapping,f.cands,f.fragment,boundary,tolerance,1,cap,
                            canon,certificates,result);
        }
        return grow_result_dict(result);
    }
};

void register_fragment_choices(py::module_& mod) {
    py::class_<FragmentChoices>(mod,"FragmentChoices")
        .def(py::init<const PySource&,const PyTarget&,int,std::vector<int>,double,double,long,
                      const std::vector<Pair>&,const std::vector<Pair>&>(),
             py::keep_alive<1,2>(),py::keep_alive<1,3>())
        .def("normal",&FragmentChoices::normal)
        .def("options",&FragmentChoices::options)
        .def("close",&FragmentChoices::close);
}
