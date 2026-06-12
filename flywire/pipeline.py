"""
End-to-end solve: correspondence -> conserved circuit -> deliverables.
"""

import random

from .loaders import load_edges
from .anchors import anchors_from_correspondence
from .isomorphism import largest_isomorphic_subset, has_edge
from .circuit import (
    largest_conserved_circuit, write_submission, visualise,
)
from .report import write_circuit_report

DS_LIST = ['BANC', 'FAFB', 'MANC']
EDGE_FILES = {
    'BANC': 'banc_626_edge_list.csv',
    'FAFB': 'fafb_783_edge_list.csv',
    'MANC': 'manc_1.2.1_edge_list.csv',
}


def induced_edges(subset, out_adj, ds):
    """Directed edges (i, j) of the induced subgraph in one dataset."""
    ids = {i: t[ds] for i, t in enumerate(subset)}
    edges = set()
    for i, a in ids.items():
        succ = out_adj.get(a)
        if not succ:
            continue
        for j, b in ids.items():
            if i != j and b in succ:
                edges.add((i, j))
    return edges


def verify_isomorphic(subset, outs, ds_list):
    """True iff the induced subgraph is identical across all datasets."""
    per_ds = {ds: induced_edges(subset, outs[ds], ds) for ds in ds_list}
    ref = per_ds[ds_list[0]]
    return all(per_ds[ds] == ref for ds in ds_list), per_ds


def random_pair_check(subset, outs, ds_list, k=20, seed=0):
    """Re-verify edge consistency on k random pairs (independent check)."""
    random.seed(seed)
    n = len(subset)
    if n < 2:
        return 0, 0
    ok = checked = 0
    for _ in range(min(k, n * (n - 1) // 2)):
        i, j = random.sample(range(n), 2)
        consistent = True
        for u, v in ((subset[i], subset[j]), (subset[j], subset[i])):
            flags = [has_edge(outs[ds], u[ds], v[ds]) for ds in ds_list]
            if any(flags) and not all(flags):
                consistent = False
        checked += 1
        ok += int(consistent)
    return ok, checked


def solve(correspondence, ds_list=None, edge_files=None, write=True):
    """
    Run the full pipeline. ``correspondence`` is a DataFrame with one column per
    dataset plus ``cell_type`` and ``neurotransmitter``. Returns a results dict.
    """
    ds_list = ds_list or DS_LIST
    edge_files = edge_files or EDGE_FILES
    type_of = dict(zip(correspondence[ds_list[0]], correspondence.cell_type))
    nt_of = dict(zip(correspondence[ds_list[0]],
                     correspondence.neurotransmitter.fillna('')))
    nt_source_of = dict(zip(correspondence[ds_list[0]],
                            correspondence.nt_source.fillna('')
                            if 'nt_source' in correspondence.columns
                            else [''] * len(correspondence)))

    outs = {ds: load_edges(edge_files[ds])[0] for ds in ds_list}
    anchors = anchors_from_correspondence(correspondence, ds_list)
    print(f"\nanchors (curated 1:1 matches): {len(anchors)}")

    # Headline: connectivity-first conserved circuit.
    circuit = largest_conserved_circuit(anchors, outs, ds_list, verbose=True)
    iso_ok, per_ds_edges = verify_isomorphic(circuit, outs, ds_list)
    n_edges = len(per_ds_edges[ds_list[0]])
    print(f"conserved circuit: {len(circuit)} neurons, {n_edges} edges, "
          f"isomorphic across all datasets: {iso_ok}")

    # Secondary: max node-count consistent set (reported for the N tension).
    star = largest_isomorphic_subset(anchors, outs, ds_list, verbose=False)
    print(f"max-N consistent set (mostly edgeless padding): {len(star)}")
    ok, checked = random_pair_check(star, outs, ds_list)
    print(f"random pair re-check on consistent set: {ok}/{checked} consistent")

    results = {
        'anchors': anchors, 'circuit': circuit, 'star': star,
        'per_ds_edges': per_ds_edges, 'iso_ok': iso_ok, 'n_edges': n_edges,
        'type_of': type_of, 'nt_of': nt_of, 'nt_source_of': nt_source_of,
        'outs': outs, 'ds_list': ds_list,
    }
    if write:
        write_outputs(results)
    return results


def write_outputs(r):
    ds_list = r['ds_list']
    write_submission(r['circuit'], ds_list, 'network.csv')
    write_submission(r['star'], ds_list, 'consistent_set_maxN.csv')
    write_circuit_report(r['circuit'], ds_list, r['per_ds_edges'],
                         r['type_of'], r['nt_of'], r['nt_source_of'],
                         len(r['star']), r['iso_ok'])
    annots = {ds_list[0]: {t[ds_list[0]]: {'type': r['type_of'].get(
        t[ds_list[0]], '?')} for t in r['circuit']}}
    visualise(r['circuit'], ds_list, r['outs'], annots,
              'circuit_visualization.png')
