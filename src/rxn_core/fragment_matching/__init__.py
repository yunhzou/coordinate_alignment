"""Reusable augmented fragment detection."""

from .detection import detect_fragments, prepare_fragment_target
from .models import (
    FragmentCandidate,
    FragmentDetectionConfig,
    FragmentDetectionResult,
    FragmentTargetContext,
)

__all__ = [
    "FragmentCandidate",
    "FragmentDetectionConfig",
    "FragmentDetectionResult",
    "FragmentTargetContext",
    "detect_fragments",
    "prepare_fragment_target",
]
