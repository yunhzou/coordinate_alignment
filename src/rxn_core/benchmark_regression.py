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

    expected_events = {
        (
            frozenset(_bond_set(item.get("broken_bonds_R"))),
            frozenset(_bond_set(item.get("formed_bonds_R"))),
        )
        for item in case.get("mechanism_events") or ()
    }
    if expected_events and _mechanism_events(record) != expected_events:
        failures.append("mechanism bond-change events differ from ground truth")

    for mechanism in record.get("mechanisms") or ():
        violations = mechanism.get("index_chirality_violations")
        if violations != int(case.get("max_chirality_violations", 0)):
            failures.append(
                "index chirality violations: expected "
                f"{case.get('max_chirality_violations', 0)}, got {violations!r}")

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

    calls = metrics.get("post_aam_call_counts") or {}
    for label, maximum in (case.get("hard_call_limits") or {}).items():
        actual = int(calls.get(label, 0))
        if actual > int(maximum):
            failures.append(f"{label} call count exceeded: {actual}>{maximum}")
    for label, maximum in (case.get("target_call_limits") or {}).items():
        actual = int(calls.get(label, 0))
        if actual > int(maximum):
            advisories.append(f"{label} remains redundant: {actual}>{maximum}")

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
