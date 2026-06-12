"""
FlyWire Codex qualification challenge.

Find the largest neuronal circuit (directed induced subgraph) that is mutually
isomorphic across at least three of the five connectomes: BANC, FAFB, MANC,
MAOL, MCNS.

The naive reading, "maximise N subject to mutually isomorphic induced
subgraphs", has a degenerate answer. Any set of neurons with no edges among
them induces the empty graph in every dataset, which is trivially isomorphic,
so arbitrary pairing yields a huge but meaningless edgeless result. The
challenge values rigour over raw size, so this package avoids that trap by:

  1. Fixing the cross-dataset correspondence biologically, using cell-type
     annotations (homologous neurons) rather than arbitrary pairing. Singleton
     cell types, those with exactly one neuron of that type per dataset, give
     an unambiguous bijection that serves as the isomorphism's node mapping.
  2. Finding the largest subset of those anchors whose induced directed
     connectivity is identical across the datasets (true mutual isomorphism).
  3. Reporting both the full consistent set and its largest connected
     component, keeping edgeless padding separate from the real wiring.

The primary correspondence uses FlyWire Codex's own curated neuron-level
cross-dataset match ids (fafb_783_match_id, manc_121_match_id), a real
homology bijection rather than string-matched cell types. The headline result
is the largest *connected* conserved circuit (connectivity-first search), not
the max node count, which would just reproduce the edgeless trap.

Module layout:
  loaders      edge-list and annotation loading
  cave         Codex/CAVE fetch and curated-match correspondence
  anchors      singleton and correspondence-based anchor construction
  isomorphism  mutual-isomorphism pruning
  circuit      connectivity-first conserved-circuit extraction and output
  report       circuit_report.txt generation
  pipeline     end-to-end solve
  cli          command-line entry point and self-test

Author: Athmeeya M Kashyap, IIIT Hyderabad
"""

from .loaders import (
    DEFAULT_EDGE_FILES,
    detect_columns,
    load_edges,
    load_annotations,
    fetch_cave_annotations,
)
from .anchors import singleton_anchors, anchors_from_correspondence
from .isomorphism import (
    has_edge,
    pair_consistent,
    consistent_degree,
    largest_isomorphic_subset,
)
from .circuit import (
    largest_connected_component,
    largest_conserved_circuit,
    internal_edges,
    write_submission,
    visualise,
)

__all__ = [
    'DEFAULT_EDGE_FILES',
    'detect_columns',
    'load_edges',
    'load_annotations',
    'fetch_cave_annotations',
    'singleton_anchors',
    'anchors_from_correspondence',
    'has_edge',
    'pair_consistent',
    'consistent_degree',
    'largest_isomorphic_subset',
    'largest_connected_component',
    'largest_conserved_circuit',
    'internal_edges',
    'write_submission',
    'visualise',
]
