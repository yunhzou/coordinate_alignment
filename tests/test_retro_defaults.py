"""Recommendation bond flexibility must not change core AAM defaults."""
import inspect
from pathlib import Path
import sys

from rxn_core.fragment import FragmentMatchConfig
from rxn_core.fragment_matching import FragmentDetectionConfig
from rxn_core.retrosynthesis import assemble_fragment_cover
from rxn_core.retrosynthesis.catalog_index import CandidateIndexConfig
from rxn_core.retrosynthesis.config import DEFAULT_ISO_TOLERANCE


def test_retro_tolerance_is_separate_from_core_matching():
    assert DEFAULT_ISO_TOLERANCE == 1.0
    assert CandidateIndexConfig().iso_tolerance == 1.0
    assert inspect.signature(assemble_fragment_cover).parameters["iso_tolerance"].default == 1.0
    assert FragmentDetectionConfig().iso_tolerance == 0.5
    assert FragmentMatchConfig().iso_tolerance == 0.5


def test_bank_scanner_uses_retro_default_and_allows_explicit_override():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
    from search_mcule_retro import _parser
    required = ["--target-smiles", "C", "--catalog", "bank.csv.gz", "--output", "part.gz"]
    assert _parser().parse_args(required).iso_tolerance == 1.0
    assert _parser().parse_args(required + ["--iso-tolerance", "0.5"]).iso_tolerance == 0.5
