"""
Synthetic tests for the conserved-circuit algorithm (no data or network).

Three datasets share a feedforward chain t1 -> t2 -> t3; t4, t5 are edgeless
padding; t6 carries one inconsistent edge (t6 -> t1) present in only one
dataset. The algorithm must drop t6, separate the padding, and recover the
t1 -> t2 -> t3 circuit.
"""

from collections import defaultdict

from flywire.anchors import singleton_anchors
from flywire.isomorphism import largest_isomorphic_subset
from flywire.circuit import (
    largest_connected_component, largest_conserved_circuit,
)

DS_LIST = ['A', 'B', 'C']


def _scenario():
    maps = {d: {f't{i}': base + i for i in range(1, 7)}
            for d, base in (('A', 10), ('B', 20), ('C', 30))}
    chain = [('t1', 't2'), ('t2', 't3')]

    def mk(edges):
        out = defaultdict(set)
        for a, b in edges:
            out[a].add(b)
        return dict(out)

    def inst(m, et):
        return [(m[a], m[b]) for a, b in et]

    outs = {
        'A': mk(inst(maps['A'], chain)),
        'B': mk(inst(maps['B'], chain)),
        'C': mk(inst(maps['C'], chain) + inst(maps['C'], [('t6', 't1')])),
    }
    annots = {d: {maps[d][ct]: {'type': ct, 'nt': None} for ct in maps[d]}
              for d in DS_LIST}
    return outs, singleton_anchors(annots, DS_LIST)


def test_max_n_consistent_set_drops_inconsistent_neuron():
    outs, anchors = _scenario()
    star = largest_isomorphic_subset(anchors, outs, DS_LIST, verbose=False)
    assert sorted(t['_type'] for t in star) == ['t1', 't2', 't3', 't4', 't5']


def test_connected_component_separates_padding():
    outs, anchors = _scenario()
    star = largest_isomorphic_subset(anchors, outs, DS_LIST, verbose=False)
    cc = largest_connected_component(star, outs, DS_LIST)
    assert sorted(t['_type'] for t in cc) == ['t1', 't2', 't3']


def test_connectivity_first_recovers_circuit():
    outs, anchors = _scenario()
    circuit = largest_conserved_circuit(anchors, outs, DS_LIST, verbose=False)
    assert sorted(t['_type'] for t in circuit) == ['t1', 't2', 't3']
