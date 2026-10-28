"""Exact maximum cliques with edge-weight tie-breaking and search pruning."""

import logging
import math
from numbers import Integral

import networkx as nx

LOGGER = logging.getLogger(__name__)


def max_edge_weight_clique(graph, weight="coor"):
    """Return the best clique by size, edge sum, then sorted feature IDs.

    Two bitset branch-and-bound searches first establish the maximum size,
    then optimize weight at that fixed size. Greedy coloring bounds clique
    size; optimistic edge sums bound weight. No maximal-clique enumeration
    is used. The graph is not modified, and self-loops are ignored.

    Edges must have the requested attribute, containing a finite integer or
    float weight. Negative weights are supported. Float weights are scaled
    exactly to integers internally so pruning cannot lose near-equal ties.
    Equal-size, equal-weight cliques are resolved deterministically by their
    lexicographically sorted feature IDs. Empty graphs return ([], 0).
    Worst-case running time is exponential.
    """
    if graph.is_directed() or graph.is_multigraph():
        raise nx.NetworkXNotImplemented("expected a simple undirected graph")
    if not graph:
        return [], 0

    try:
        sorted_features = sorted(graph)
    except TypeError:
        sorted_features = sorted(
            graph,
            key=lambda node: (
                type(node).__module__,
                type(node).__qualname__,
                repr(node),
            ),
        )
    feature_rank = {feature: rank for rank, feature in enumerate(sorted_features)}
    nodes = sorted(
        graph,
        key=lambda node: (-graph.degree[node], feature_rank[node]),
    )
    node_index = {node: index for index, node in enumerate(nodes)}
    adjacency = [0] * len(nodes)
    edge_ratios = []
    scale = 1
    for source, target, data in graph.edges(data=True):
        if source == target:
            continue
        source, target = node_index[source], node_index[target]
        adjacency[source] |= 1 << target
        adjacency[target] |= 1 << source
        value = data[weight]
        if isinstance(value, Integral):
            numerator, denominator = int(value), 1
        else:
            value = float(value)
            if not math.isfinite(value):
                raise ValueError("clique edge weights must be finite")
            numerator, denominator = value.as_integer_ratio()
        scale = max(scale, denominator)
        edge_ratios.append((source, target, numerator, denominator))

    # Binary floats have power-of-two denominators, so the largest is a
    # common denominator. These integer scores preserve their exact sums.
    weights = [{} for _ in nodes]
    for source, target, numerator, denominator in edge_ratios:
        value = numerator * (scale // denominator)
        weights[source][target] = value
        weights[target][source] = value
    del edge_ratios
    if not any(adjacency):
        return [sorted_features[0]], 0

    weighted_degrees = [sum(row.values()) for row in weights]
    all_nodes = (1 << len(nodes)) - 1
    score_scale = 1 << len(nodes)
    feature_scores = [1 << (len(nodes) - feature_rank[node] - 1) for node in nodes]

    def color_sort(candidates, gains=None):
        # Each color class is independent: a clique can take at most one
        # node from it. Prefix color numbers are upper bounds on clique size.
        order = []
        bounds = []
        color = 0
        remaining = candidates
        while remaining:
            color += 1
            available = remaining
            color_class = []
            while available:
                bit = available & -available
                vertex = bit.bit_length() - 1
                color_class.append(vertex)
                remaining ^= bit
                available ^= bit
                available &= ~adjacency[vertex]
            if gains is not None:
                # Traversal is reversed, so promising nodes are placed last.
                color_class.sort(key=lambda v: (gains[v], weighted_degrees[v]))
            order.extend(color_class)
            bounds.extend([color] * len(color_class))
        return order, bounds

    # Seed the size search with a greedy clique to make early bounds useful.
    best_clique = []
    candidates = all_nodes
    while candidates:
        choices = []
        remaining = candidates
        while remaining:
            bit = remaining & -remaining
            choices.append(bit.bit_length() - 1)
            remaining ^= bit
        vertex = max(
            choices,
            key=lambda v: (
                (adjacency[v] & candidates).bit_count(),
                weighted_degrees[v],
            ),
        )
        best_clique.append(vertex)
        candidates &= adjacency[vertex]

    def find_size(clique, candidates):
        nonlocal best_clique
        if len(clique) + candidates.bit_count() <= len(best_clique):
            return
        order, bounds = color_sort(candidates)
        for index in range(len(order) - 1, -1, -1):
            if len(clique) + bounds[index] <= len(best_clique):
                return
            vertex = order[index]
            candidates &= ~(1 << vertex)
            clique.append(vertex)
            following = candidates & adjacency[vertex]
            if following:
                find_size(clique, following)
            elif len(clique) > len(best_clique):
                best_clique = clique[:]
            clique.pop()

    find_size([], all_nodes)
    target_size = len(best_clique)
    best_edge_weight = sum(
        weights[vertex][other]
        for index, vertex in enumerate(best_clique)
        for other in best_clique[:index]
    )
    best_score = best_edge_weight * score_scale + sum(
        feature_scores[vertex] for vertex in best_clique
    )
    # LOGGER.info("Maximum clique size is %s; optimizing edge weight.", target_size)

    max_edge = max(max(row.values()) for row in weights if row) * score_scale
    # Rows ordered once let the tighter bound collect the strongest possible
    # edges to remaining candidates without sorting at each search node.
    ranked_edges = [
        sorted(
            ((other, value * score_scale) for other, value in row.items()),
            key=lambda edge: edge[1],
            reverse=True,
        )
        for row in weights
    ]

    def find_weight(clique, candidates, gains, current_score):
        nonlocal best_clique, best_score
        needed = target_size - len(clique)
        if not needed:
            if current_score > best_score:
                best_clique = clique[:]
                best_score = current_score
            return
        if candidates.bit_count() < needed:
            return
        if needed == 1:
            # Every candidate already connects to the whole partial clique.
            # Only its gain matters when selecting the final vertex.
            remaining = candidates
            while remaining:
                bit = remaining & -remaining
                vertex = bit.bit_length() - 1
                candidate_score = current_score + gains[vertex]
                if candidate_score > best_score:
                    best_clique = clique + [vertex]
                    best_score = candidate_score
                remaining ^= bit
            return
        order, bounds = color_sort(candidates, gains)
        if bounds[-1] < needed:
            return

        class_gains = [None] * bounds[-1]
        for vertex, color in zip(order, bounds):
            previous = class_gains[color - 1]
            if previous is None or gains[vertex] > previous:
                class_gains[color - 1] = gains[vertex]
        upper_score = (
            current_score
            + sum(sorted(class_gains, reverse=True)[:needed])
            + needed * (needed - 1) // 2 * max_edge
        )
        if upper_score <= best_score:
            return

        # Each selected vertex contributes its gain to the current clique
        # plus half its incident edges within the completion. Taking its
        # strongest possible (needed - 1) edges overestimates that contribution.
        # Maximize per color class, then take the best `needed` classes.
        class_scores = [None] * bounds[-1]
        for vertex, color in zip(order, bounds):
            missing = needed - 1
            contribution = 2 * gains[vertex]
            if missing:
                for other, value in ranked_edges[vertex]:
                    if candidates & (1 << other):
                        contribution += value
                        missing -= 1
                        if not missing:
                            break
            if missing:
                continue
            previous = class_scores[color - 1]
            if previous is None or contribution > previous:
                class_scores[color - 1] = contribution
        scores = [score for score in class_scores if score is not None]
        if len(scores) < needed:
            return
        doubled_upper = 2 * current_score + sum(sorted(scores, reverse=True)[:needed])
        if doubled_upper <= 2 * best_score:
            return

        for index in range(len(order) - 1, -1, -1):
            # Equality must remain searchable to resolve edge-weight ties.
            if bounds[index] < needed:
                return
            vertex = order[index]
            candidates &= ~(1 << vertex)
            following = candidates & adjacency[vertex]
            if following.bit_count() < needed - 1:
                continue
            next_gains = gains.copy()
            remaining = following
            while remaining:
                bit = remaining & -remaining
                other = bit.bit_length() - 1
                next_gains[other] += weights[vertex][other] * score_scale
                remaining ^= bit
            clique.append(vertex)
            find_weight(clique, following, next_gains, current_score + gains[vertex])
            clique.pop()

    find_weight([], all_nodes, feature_scores[:], 0)
    best_edge_weight = best_score // score_scale
    total_weight = best_edge_weight if scale == 1 else best_edge_weight / scale
    clique = sorted(
        (nodes[vertex] for vertex in best_clique),
        key=feature_rank.__getitem__,
    )
    return clique, total_weight
