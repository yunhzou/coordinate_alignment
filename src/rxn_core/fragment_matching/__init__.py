"""Reusable augmented fragment detection."""

from .detection import detect_fragments, prepare_fragment_target
from .models import (
    FragmentCandidate,
    FragmentDetectionConfig,
    FragmentDetectionResult,
    FragmentTargetContext,
)
from .symmetry import (
    FragmentOrbitLimitExceeded,
    materialize_target_coverage_orbit,
)

__all__ = [
    "FragmentCandidate",
    "FragmentDetectionConfig",
    "FragmentDetectionResult",
    "FragmentTargetContext",
    "FragmentOrbitLimitExceeded",
    "detect_fragments",
    "materialize_target_coverage_orbit",
    "prepare_fragment_target",
]
