"""Direction selection outside the directed, symmetry-compressed AAM engine."""
from dataclasses import dataclass, replace

from .domain import AAMProblem, AAMSearchConfig


@dataclass(frozen=True)
class AAMSearchPlan:
    """Keep input R/P identities separate from search source/target identities.

    The raw AAM result and its symmetry actions retain search orientation.
    Convert only a concrete mapping after selecting/realizing its search path.
    Do not invert generators or reinterpret a compressed graph in place.
    """
    input_problem: AAMProblem
    problem: AAMProblem
    config: AAMSearchConfig
    reversed: bool

    def to_input_mapping(self, mapping):
        pairs = dict(mapping)
        return {r:p for p,r in pairs.items()} if self.reversed else pairs

    def to_search_mapping(self, mapping):
        return self.to_input_mapping(mapping)

    @property
    def direction(self):
        return 'P_to_R' if self.reversed else 'R_to_P'


def plan_aam_search(problem, config=None):
    """Search smaller explicit-atom endpoint first; retain R→P on ties."""
    config = config or AAMSearchConfig()
    reverse = problem.source_atom_count > problem.target_atom_count
    search_problem = AAMProblem(problem.product, problem.reactant, problem.name) if reverse else problem
    search_config = replace(config, anchors=tuple((p,r) for r,p in config.anchors)) if reverse else config
    return AAMSearchPlan(problem, search_problem, search_config, reverse)
