"""
Command-line entry point.

Subcommands:
  selftest   run the synthetic algorithm check (no data or network needed)
  fetch      pull the Codex cross-dataset correspondence via CAVE (needs token)
  solve      build the conserved circuit + all deliverables from a correspondence
  search     original singleton cell-type search over --annot files (alternative)
"""

import argparse
import os
import sys
from collections import defaultdict
from itertools import combinations

from .loaders import (
    DEFAULT_EDGE_FILES,
    load_edges,
    load_annotations,
    fetch_cave_annotations,
)
from .anchors import singleton_anchors
from .isomorphism import largest_isomorphic_subset
from .circuit import (
    largest_connected_component,
    largest_conserved_circuit,
    write_submission,
    visualise,
)


# --------------------------------------------------------------------------- #
# selftest
# --------------------------------------------------------------------------- #
def selftest():
    """Synthetic check: prune one inconsistent neuron, recover the circuit."""
    print("SELF-TEST")
    maps = {
        'A': {f't{i}': 10 + i for i in range(1, 7)},
        'B': {f't{i}': 20 + i for i in range(1, 7)},
        'C': {f't{i}': 30 + i for i in range(1, 7)},
    }
    ds_list = ['A', 'B', 'C']

    def mk(edges):
        out = defaultdict(set)
        for a, b in edges:
            out[a].add(b)
        return dict(out)

    chain = [('t1', 't2'), ('t2', 't3')]

    def inst(m, et):
        return [(m[a], m[b]) for a, b in et]

    outs = {
        'A': mk(inst(maps['A'], chain)),
        'B': mk(inst(maps['B'], chain)),
        'C': mk(inst(maps['C'], chain) + inst(maps['C'], [('t6', 't1')])),
    }
    annots = {d: {maps[d][ct]: {'type': ct, 'nt': None} for ct in maps[d]}
              for d in ds_list}

    anchors = singleton_anchors(annots, ds_list)
    star = largest_isomorphic_subset(anchors, outs, ds_list, verbose=False)
    kept = sorted(t['_type'] for t in star)
    assert kept == ['t1', 't2', 't3', 't4', 't5'], kept

    cc = largest_connected_component(star, outs, ds_list)
    assert sorted(t['_type'] for t in cc) == ['t1', 't2', 't3'], cc

    circuit = largest_conserved_circuit(anchors, outs, ds_list, verbose=False)
    assert sorted(t['_type'] for t in circuit) == ['t1', 't2', 't3'], circuit

    print(f"  max-N consistent set : {kept}")
    print(f"  conserved circuit    : {sorted(t['_type'] for t in circuit)}")
    print("  PASSED: inconsistent neuron pruned, edgeless padding separated, "
          "feedforward circuit recovered")


# --------------------------------------------------------------------------- #
# fetch / solve  (Codex curated-match correspondence)
# --------------------------------------------------------------------------- #
def cmd_fetch(args):
    from .cave import fetch_and_build
    fetch_and_build(out_csv=args.out)


def cmd_solve(args):
    import pandas as pd
    from .pipeline import solve, DS_LIST
    corr = pd.read_csv(args.correspondence,
                       dtype={ds: int for ds in DS_LIST})
    solve(corr, write=True)


# --------------------------------------------------------------------------- #
# search  (original singleton cell-type methodology over --annot files)
# --------------------------------------------------------------------------- #
def parse_kv(items):
    d = {}
    for it in items or []:
        if '=' not in it:
            raise ValueError(f"expected NAME=path, got {it}")
        k, v = it.split('=', 1)
        d[k.strip().upper()] = v
    return d


def load_datasets(edge_paths, annot_paths):
    graphs, annots = {}, {}
    for ds, p in edge_paths.items():
        if not os.path.exists(p):
            print(f"  missing edge file for {ds}: {p}")
            continue
        out_adj, in_adj, nodes = load_edges(p)
        graphs[ds] = {'out': out_adj, 'in': in_adj, 'nodes': nodes}
    for ds in graphs:
        if ds in annot_paths and os.path.exists(annot_paths[ds]):
            annots[ds] = load_annotations(annot_paths[ds])
        else:
            print(f"  no annotation file for {ds}; trying CAVE ...")
            annots[ds] = fetch_cave_annotations(ds)
        if not annots.get(ds):
            print(f"  WARNING: no cell types for {ds}; only the degenerate "
                  f"edgeless solution exists. Provide --annot {ds}=...")
    return graphs, annots


def search_best(graphs, annots, available, require_nt=False):
    best = None
    print("\nSearching dataset combinations:")
    for k in range(3, len(available) + 1):
        for combo in combinations(available, k):
            ds_list = list(combo)
            outs = {ds: graphs[ds]['out'] for ds in ds_list}
            anchors = singleton_anchors(annots, ds_list, require_nt=require_nt)
            if len(anchors) < 2:
                print(f"  {combo}: {len(anchors)} anchors (skip)")
                continue
            star = largest_isomorphic_subset(anchors, outs, ds_list,
                                             verbose=False)
            circuit = largest_conserved_circuit(anchors, outs, ds_list,
                                                verbose=False)
            print(f"  {combo}: anchors={len(anchors)}  consistent "
                  f"N={len(star)}  conserved circuit={len(circuit)}")
            score = (len(circuit), len(star))
            if best is None or score > best['score']:
                best = {'score': score, 'combo': ds_list, 'star': star,
                        'circuit': circuit, 'outs': outs}
    return best


def cmd_search(args):
    edge_paths = parse_kv(args.edges) or dict(DEFAULT_EDGE_FILES)
    annot_paths = parse_kv(args.annot)
    graphs, annots = load_datasets(edge_paths, annot_paths)
    available = [d for d in graphs if annots.get(d)]
    if len(available) < 3:
        print(f"\nNeed >=3 datasets with edges and annotations. "
              f"Have: {available}")
        sys.exit(1)
    best = search_best(graphs, annots, available, require_nt=args.require_nt)
    if not best:
        print("No viable correspondence found.")
        sys.exit(1)
    ds_list, star, circuit = best['combo'], best['star'], best['circuit']
    print(f"\nBEST: datasets={tuple(ds_list)}  consistent N={len(star)}  "
          f"conserved circuit={len(circuit)}")
    write_submission(circuit, ds_list, args.out)
    write_submission(star, ds_list, 'consistent_set_maxN.csv')
    visualise(circuit, ds_list, best['outs'], annots,
              'circuit_visualization.png')


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(prog='flywire',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd')

    sub.add_parser('selftest', help='run the synthetic algorithm check')

    p_fetch = sub.add_parser('fetch', help='pull Codex correspondence (CAVE)')
    p_fetch.add_argument('--out',
                         default='codex_correspondence_banc_fafb_manc.csv')

    p_solve = sub.add_parser('solve', help='build conserved circuit + outputs')
    p_solve.add_argument('--correspondence',
                         default='codex_correspondence_banc_fafb_manc.csv')

    p_search = sub.add_parser('search',
                              help='singleton cell-type search over --annot')
    p_search.add_argument('--edges', nargs='*')
    p_search.add_argument('--annot', nargs='*')
    p_search.add_argument('--require-nt', action='store_true')
    p_search.add_argument('--out', default='submission.csv')

    # Back-compat: `--selftest` flag with no subcommand.
    ap.add_argument('--selftest', action='store_true', help=argparse.SUPPRESS)

    args = ap.parse_args(argv)

    if args.selftest or args.cmd == 'selftest':
        selftest()
    elif args.cmd == 'fetch':
        cmd_fetch(args)
    elif args.cmd == 'solve':
        cmd_solve(args)
    elif args.cmd == 'search':
        cmd_search(args)
    else:
        ap.print_help()
