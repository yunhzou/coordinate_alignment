"""Node-level compatibility policy for weighted graph matching."""
from __future__ import annotations

from collections.abc import Callable, Hashable, Sequence
from dataclasses import dataclass
from typing import Any


def _node_attr(g, node, name, default=None):
    data = g.nodes[node]
    if name in data:
        return data.get(name, default)
    features = data.get("features")
    if isinstance(features, dict):
        return features.get(name, default)
    return default


class ElementNodeMatchPolicy:
    """Default chemistry policy: nodes are compatible iff elements match."""

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return (
            _node_attr(g_query, query_node, "element")
            == _node_attr(g_target, target_node, "element")
        )

    def key(self, g, node) -> Hashable:
        return _node_attr(g, node, "element")


@dataclass(frozen=True)
class AttributeNodeMatchPolicy:
    """Match nodes by one or more node attributes or feature entries."""

    fields: tuple[str, ...]

    def __init__(self, fields: str | Sequence[str]):
        if isinstance(fields, str):
            fields = (fields,)
        object.__setattr__(self, "fields", tuple(fields))
        if not self.fields:
            raise ValueError("AttributeNodeMatchPolicy requires at least one field")

    def key(self, g, node) -> Hashable:
        values = tuple(_node_attr(g, node, field) for field in self.fields)
        return values[0] if len(values) == 1 else values

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return self.key(g_query, query_node) == self.key(g_target, target_node)


@dataclass(frozen=True)
class CallableNodeMatchPolicy:
    """Policy wrapper for user-supplied node compatibility functions."""

    node_match: Callable[[Any, int, Any, int], bool]
    node_key: Callable[[Any, int], Hashable] | None = None

    def compatible(self, g_query, query_node, g_target, target_node) -> bool:
        return bool(self.node_match(g_query, query_node, g_target, target_node))

    def key(self, g, node) -> Hashable:
        if self.node_key is not None:
            return self.node_key(g, node)
        return DEFAULT_NODE_POLICY.key(g, node)


DEFAULT_NODE_POLICY = ElementNodeMatchPolicy()


def as_node_match_policy(policy=None, *, node_match=None, node_key=None):
    """Normalize policy/callback input to a policy object.

    Callers may pass:
    - ``None`` for the default same-element policy
    - an object with ``compatible(...)`` and ``key(...)``
    - a callable ``(g_query, i, g_target, j) -> bool``
    - explicit ``node_match=...`` and optional ``node_key=...``
    """
    if policy is not None and node_match is not None:
        raise ValueError("pass either policy or node_match, not both")
    if policy is None and node_match is None:
        if node_key is not None:
            raise ValueError("node_key requires a custom node_match or policy")
        return DEFAULT_NODE_POLICY
    if policy is None:
        return CallableNodeMatchPolicy(node_match, node_key)
    if isinstance(policy, str):
        if node_key is not None:
            raise ValueError("node_key is not used with string attribute policy")
        return AttributeNodeMatchPolicy(policy)
    if isinstance(policy, (tuple, list)):
        if node_key is not None:
            raise ValueError("node_key is not used with attribute policy")
        return AttributeNodeMatchPolicy(tuple(policy))
    if callable(policy) and not hasattr(policy, "compatible"):
        return CallableNodeMatchPolicy(policy, node_key)
    if not hasattr(policy, "compatible") or not hasattr(policy, "key"):
        raise TypeError("node policy must define compatible(...) and key(...)")
    if node_key is not None:
        raise ValueError("node_key is not used when policy object is supplied")
    return policy
