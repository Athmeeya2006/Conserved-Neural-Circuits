"""
Mutual-isomorphism pruning (core algorithm).
"""


def has_edge(out_adj, a, b):
    """True if the directed edge a -> b exists in this adjacency map."""
    s = out_adj.get(a)
    return s is not None and b in s


def pair_consistent(ti, tj, outs, ds_list):
    """
    Under the cell-type bijection, the directed edges ti->tj and tj->ti must be
    present in all datasets or in none. Otherwise the induced subgraphs differ.
    """
    for a, b in ((ti, tj), (tj, ti)):
        flags = [has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list]
        if any(flags) and not all(flags):
            return False
    return True


def consistent_degree(idx, cur, outs, ds_list):
    """Count agreeing (present-in-all) edges from node idx to the kept set."""
    a = cur[idx]
    deg = 0
    for k, b in enumerate(cur):
        if k == idx:
            continue
        if all(has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list) or \
           all(has_edge(outs[ds], b[ds], a[ds]) for ds in ds_list):
            deg += 1
    return deg


def largest_isomorphic_subset(anchors, outs, ds_list, verbose=True):
    """
    Prune the anchor set down to the largest mutually-isomorphic induced
    subgraph under the fixed cross-dataset bijection.

    Two anchors conflict when their connectivity disagrees across datasets
    (see pair_consistent). The kept set must be conflict-free, i.e. an
    independent set in the conflict graph; the largest such set is the maximum
    independent set. That is NP-hard, so this uses the standard greedy
    approximation: repeatedly remove the single most-conflicting anchor, with
    ties broken toward the FEWEST consistent edges so peripheral/padding nodes
    go before real circuit nodes. Removing one node at a time (rather than
    batch-removing every conflicting node) is essential: in a dense conflict
    graph almost every node has some conflict, so a batch pass collapses the
    whole set to nothing.

    The conflict graph and per-node consistent-edge counts are precomputed
    once, then nodes are peeled with incremental degree updates.
    """
    n = len(anchors)
    if n < 2:
        return list(anchors)

    conflict = [set() for _ in range(n)]
    cdeg = [0] * n  # consistent (present-in-all) edges per node
    for i in range(n):
        a = anchors[i]
        for j in range(i + 1, n):
            b = anchors[j]
            if not pair_consistent(a, b, outs, ds_list):
                conflict[i].add(j)
                conflict[j].add(i)
            elif all(has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list) or \
                    all(has_edge(outs[ds], b[ds], a[ds]) for ds in ds_list):
                cdeg[i] += 1
                cdeg[j] += 1

    alive = [True] * n
    deg = [len(conflict[i]) for i in range(n)]
    removed = 0
    while True:
        best = -1
        best_key = None
        for i in range(n):
            if not alive[i] or deg[i] == 0:
                continue
            key = (deg[i], -cdeg[i])  # most conflicts, then fewest real edges
            if best_key is None or key > best_key:
                best_key = key
                best = i
        if best == -1:
            break  # conflict-free: remaining set is mutually isomorphic
        alive[best] = False
        removed += 1
        for j in conflict[best]:
            if alive[j]:
                deg[j] -= 1
        if verbose and removed % 100 == 0:
            print(f"      removed {removed}, kept {n - removed}", flush=True)

    return [anchors[i] for i in range(n) if alive[i]]
