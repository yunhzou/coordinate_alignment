"""Validate one AAM benchmark record against versioned regression metadata."""
from __future__ import annotations

import json
from pathlib import Path


def load_contract(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _bond_set(value):
    return {tuple(sorted(map(int, pair))) for pair in value or ()}


def _mechanism_events(record):
    return {
        (
            frozenset(_bond_set(item.get("broken_bonds_R"))),
            frozenset(_bond_set(item.get("formed_bonds_R"))),
        )
        for item in record.get("mechanisms") or ()
    }


def evaluate_record(record, contract, performance_profile=None):
    """Return hard failures and improvement advisories for one record."""
    case_name = str(record.get("case"))
    case = (contract.get("cases") or {}).get(case_name)
    if case is None:
        return {
            "case": case_name,
            "failures": [f"case is absent from regression contract: {case_name}"],
            "advisories": [],
        }

    failures = []
    advisories = []

    def equal(field, actual, expected):
        if actual != expected:
            failures.append(f"{field}: expected {expected!r}, got {actual!r}")

    equal("atom_count", int(record.get("atom_count", -1)),
          int(case["atom_count"]))
    equal("selected_mechanism_count",
          int(record.get("selected_mechanism_count", -1)),
          int(case.get("selected_mechanism_count", 1)))

    expected_event_certificates = set(
        map(str, case.get("mechanism_event_certificate_digests") or ()))
    observed_event_certificates = {
        str(item.get("event_certificate_digest"))
        for item in record.get("mechanisms") or ()
        if item.get("event_certificate_digest") is not None
    }
    expected_events = {
        (
            frozenset(_bond_set(item.get("broken_bonds_R"))),
            frozenset(_bond_set(item.get("formed_bonds_R"))),
        )
        for item in case.get("mechanism_events") or ()
    }
    if (expected_event_certificates
            and observed_event_certificates != expected_event_certificates):
        failures.append(
            "symmetry-canonical mechanism event certificates differ from "
            "ground truth")
    elif (not expected_event_certificates and expected_events
          and _mechanism_events(record) != expected_events):
        failures.append("mechanism bond-change events differ from ground truth")

    for mechanism in record.get("mechanisms") or ():
        violations = mechanism.get("index_chirality_violations")
        if violations != int(case.get("max_chirality_violations", 0)):
            failures.append(
                "index chirality violations: expected "
                f"{case.get('max_chirality_violations', 0)}, got {violations!r}")

    observed_degeneracy = {
        (group.get("center_R"),
         tuple(sorted(map(int, group.get("r_atoms") or ()))))
        for mechanism in record.get("mechanisms") or ()
        for group in mechanism.get("degeneracy_groups") or ()
    }
    required_degeneracy = {
        (group.get("center_R"),
         tuple(sorted(map(int, group.get("r_atoms") or ()))))
        for group in case.get("required_degeneracy_groups") or ()
    }
    missing_degeneracy = required_degeneracy - observed_degeneracy
    if missing_degeneracy:
        failures.append(
            "required degeneracy groups are missing: "
            f"{sorted(missing_degeneracy, key=repr)}")

    metrics = record.get("regression_metrics") or {}
    configured_cap = metrics.get("configured_max_branches")
    expected_cap = int(case.get("max_branches", 100))
    if configured_cap is not None and int(configured_cap) != expected_cap:
        failures.append(
            f"configured_max_branches: expected {expected_cap}, "
            f"got {configured_cap}")
    max_live = metrics.get("max_live_branches")
    if max_live is not None and int(max_live) > expected_cap:
        failures.append(
            f"max_live_branches exceeded cap: {max_live}>{expected_cap}")
    candidate_limit = case.get("max_growth_candidates")
    max_candidates = metrics.get("max_growth_candidates")
    if (candidate_limit is not None and max_candidates is not None
            and int(max_candidates) > int(candidate_limit)):
        failures.append(
            "max_growth_candidates exceeded case ceiling: "
            f"{max_candidates}>{candidate_limit}")

    calls = metrics.get("post_aam_call_counts") or {}
    for label, maximum in (case.get("hard_call_limits") or {}).items():
        actual = int(calls.get(label, 0))
        if actual > int(maximum):
            failures.append(f"{label} call count exceeded: {actual}>{maximum}")
    for label, maximum in (case.get("target_call_limits") or {}).items():
        actual = int(calls.get(label, 0))
        if actual > int(maximum):
            advisories.append(f"{label} remains redundant: {actual}>{maximum}")

    post_metrics = metrics.get("pipeline_post_aam_metrics") or {}
    requests = post_metrics.get("completed_candidate_group_requests")
    calculations = post_metrics.get(
        "completed_candidate_group_calculations")
    cache_hits = post_metrics.get("completed_candidate_group_cache_hits")
    if requests is not None:
        if calculations is None or cache_hits is None:
            failures.append("completed candidate group cache metrics incomplete")
        elif int(calculations) + int(cache_hits) != int(requests):
            failures.append(
                "completed candidate group cache accounting mismatch: "
                f"{calculations}+{cache_hits}!={requests}")

    if performance_profile:
        profile = ((case.get("performance_profiles") or {})
                   .get(performance_profile))
        if profile is None:
            failures.append(
                f"case has no performance profile {performance_profile!r}")
        else:
            for field, maximum in profile.items():
                actual = record.get(field)
                if actual is None:
                    failures.append(f"missing performance metric: {field}")
                elif float(actual) > float(maximum):
                    failures.append(
                        f"{field} exceeded {performance_profile}: "
                        f"{actual}>{maximum}")

    return {
        "case": case_name,
        "performance_profile": performance_profile,
        "failures": failures,
        "advisories": advisories,
        "passed": not failures,
    }
