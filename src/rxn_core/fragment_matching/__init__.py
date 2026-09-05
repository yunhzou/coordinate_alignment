"""Reusable augmented fragment detection."""

from .detection import detect_fragments, prepare_fragment_target
from .models import (
    FragmentCandidate,
    FragmentDerivation,
    FragmentDetectionConfig,
    FragmentDetectionResult,
    FragmentTargetContext,
)
from .parallel import FragmentDetectionExecution, detect_fragments_parallel
from .progressive import (
    ProgressiveFragmentMatchingResult,
    ProgressiveFragmentPlacement,
    ProgressiveSelection,
    progressive_fragment_matching,
)
from .symmetry import (
    FragmentOrbitLimitExceeded,
    materialize_target_coverage_orbit,
)

__all__ = [
    "FragmentCandidate",
    "FragmentDerivation",
    "FragmentDetectionConfig",
    "FragmentDetectionResult",
    "FragmentDetectionExecution",
    "FragmentTargetContext",
    "FragmentOrbitLimitExceeded",
    "detect_fragments",
    "detect_fragments_parallel",
    "materialize_target_coverage_orbit",
    "prepare_fragment_target",
    "ProgressiveFragmentMatchingResult",
    "ProgressiveFragmentPlacement",
    "ProgressiveSelection",
    "progressive_fragment_matching",
]
