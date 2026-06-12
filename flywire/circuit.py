"""
Connected-circuit extraction (the biologically meaningful part) and output.
"""

from collections import defaultdict

from .isomorphism import has_edge, pair_consistent


def largest_connected_component(subset, outs, ds_list):
    """
    Within the consistent set, return the subset of tuples forming the largest
    weakly-connected component over edges present in all datasets. This strips
    the edgeless padding and exposes the real shared circuit.
    """
    adj = defaultdict(set)
    n = len(subset)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = subset[i], subset[j]
            fwd = all(has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list)
            rev = all(has_edge(outs[ds], b[ds], a[ds]) for ds in ds_list)
            if fwd or rev:
                adj[i].add(j)
                adj[j].add(i)
    seen = set()
    best = []
    for start in range(n):
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack.extend(y for y in adj[x] if y not in seen)
        if len(comp) > len(best):
            best = comp
    return [subset[i] for i in best]


def _agreement_adjacency(anchors, outs, ds_list):
    """
    Directed and undirected adjacency over anchor indices for the *agreement
    graph*: an edge i -> j exists iff the directed edge is present in ALL
    datasets. Returns (succ, undirected).
    """
    n = len(anchors)
    succ = [set() for _ in range(n)]
    und = [set() for _ in range(n)]
    for i in range(n):
        a = anchors[i]
        for j in range(n):
            if i == j:
                continue
            b = anchors[j]
            if all(has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list):
                succ[i].add(j)
                und[i].add(j)
                und[j].add(i)
    return succ, und


def _weak_components(nodes, und):
    """Weakly-connected components (over agreement edges) within `nodes`."""
    nodes = set(nodes)
    seen = set()
    comps = []
    for start in nodes:
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack.extend(y for y in und[x] if y in nodes and y not in seen)
        comps.append(comp)
    return comps


def _largest_weak_piece(nodes, und):
    comps = _weak_components(nodes, und)
    return set(max(comps, key=len)) if comps else set()


def _prune_component_to_conflict_free(comp, anchors, outs, ds_list, und):
    """
    Within one agreement-graph component, remove conflicting nodes until no pair
    is inconsistent, then return the largest connected piece that survives.

    The conflict graph is precomputed once; nodes are peeled highest-conflict
    first, ties broken toward the FEWEST agreement edges so densely-wired
    circuit nodes are kept and peripheral conflicting nodes go first.
    Connectivity is only re-evaluated at the end (collapsing to the largest
    piece after every single removal discards good structure prematurely).
    """
    comp = list(comp)
    m = len(comp)
    comp_set = set(comp)
    conf = [set() for _ in range(m)]
    for x in range(m):
        for y in range(x + 1, m):
            if not pair_consistent(anchors[comp[x]], anchors[comp[y]],
                                   outs, ds_list):
                conf[x].add(y)
                conf[y].add(x)
    agdeg = [len(und[comp[x]] & comp_set) for x in range(m)]
    alive = [True] * m
    deg = [len(conf[x]) for x in range(m)]
    while True:
        best = -1
        best_key = None
        for x in range(m):
            if not alive[x] or deg[x] == 0:
                continue
            key = (deg[x], -agdeg[x])  # most conflicts, then least-wired
            if best_key is None or key > best_key:
                best_key = key
                best = x
        if best == -1:
            break
        alive[best] = False
        for y in conf[best]:
            if alive[y]:
                deg[y] -= 1
    survivors = {comp[x] for x in range(m) if alive[x]}
    return _largest_weak_piece(survivors, und)


def largest_conserved_circuit(anchors, outs, ds_list, verbose=True):
    """
    Connectivity-first search for the largest *connected, conflict-free*
    conserved circuit (the directed induced subgraph that is identical across
    all datasets).

    Unlike maximising the node count (which is biased to keep isolated padding
    and discard wired neurons), this searches directly for a circuit:

      1. Build the agreement graph: directed edge u -> v iff present in ALL
         datasets. Its weakly-connected components are edge-rich circuit
         candidates by construction.
      2. For each component (largest first), prune the minimum nodes to remove
         any internal inconsistency while keeping it connected.
      3. Return the largest surviving connected, conflict-free component.

    The result is the biologically meaningful conserved circuit; it does not
    contain edgeless padding.
    """
    n = len(anchors)
    if n < 2:
        return list(anchors)
    _, und = _agreement_adjacency(anchors, outs, ds_list)
    comps = sorted((c for c in _weak_components(range(n), und) if len(c) >= 2),
                   key=len, reverse=True)
    if verbose:
        print(f"      agreement graph: {len(comps)} components, "
              f"largest {len(comps[0]) if comps else 0}", flush=True)
    best = set()
    for comp in comps:
        if len(comp) <= len(best):
            break  # pruning only shrinks; no smaller component can win
        pruned = _prune_component_to_conflict_free(comp, anchors, outs,
                                                   ds_list, und)
        if len(pruned) > len(best):
            best = pruned
    return [anchors[i] for i in sorted(best)]


def internal_edges(subset, outs, ds_list):
    """Directed edges present in all datasets among the subset."""
    edges = []
    n = len(subset)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = subset[i], subset[j]
            if all(has_edge(outs[ds], a[ds], b[ds]) for ds in ds_list):
                edges.append((i, j))
    return edges


def write_submission(subset, ds_list, path):
    cols = ds_list[:3]
    with open(path, 'w') as fh:
        fh.write(','.join(cols) + '\n')
        for t in subset:
            fh.write(','.join(str(t[ds]) for ds in cols) + '\n')
    print(f"  wrote {path}: {len(subset)} rows x {len(cols)} cols")


def visualise(subset, ds_list, outs, annots, out_png):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError:
        print("  (matplotlib/networkx not installed; skipping figure)")
        return
    ds = ds_list[0]
    ann = annots[ds]
    graph = nx.DiGraph()
    for t in subset:
        nid = t[ds]
        graph.add_node(nid, label=ann.get(nid, {}).get('type', str(nid)))
    ids = {t[ds] for t in subset}
    for t in subset:
        for u in outs[ds].get(t[ds], ()):
            if u in ids:
                graph.add_edge(t[ds], u)
    if graph.number_of_nodes() == 0:
        return
    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(
        graph, seed=7, k=1.5 / max(1, graph.number_of_nodes() ** 0.5))
    nx.draw_networkx_nodes(graph, pos, node_size=420, node_color='#3b6ea5',
                           alpha=0.9, ax=ax)
    nx.draw_networkx_edges(graph, pos, arrows=True, arrowsize=10,
                           edge_color='#555', width=0.8, ax=ax)
    if graph.number_of_nodes() <= 60:
        nx.draw_networkx_labels(
            graph, pos, {n: graph.nodes[n]['label'][:12] for n in graph.nodes()},
            font_size=6, ax=ax)
    ax.set_title(f"Shared circuit in {ds}: "
                 f"{graph.number_of_nodes()} neurons, "
                 f"{graph.number_of_edges()} edges")
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  wrote {out_png}")
