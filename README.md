<div align="center">

# Conserved Neural Circuits

**Finding the largest neuronal circuit with identical wiring across multiple *Drosophila* connectomes.**

[![CI](https://github.com/Athmeeya2006/Conserved-Neural-Circuits/actions/workflows/ci.yml/badge.svg)](https://github.com/Athmeeya2006/Conserved-Neural-Circuits/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![Lint: ruff](https://img.shields.io/badge/lint-ruff-261230.svg)](https://github.com/astral-sh/ruff)

</div>

## Overview

Five whole-animal *Drosophila* connectomes (BANC, FAFB, MANC, MAOL, MCNS) have each been reconstructed independently. This project finds the largest directed circuit that is mutually isomorphic, meaning it has identical wiring with direction preserved, across at least three of them. Such a circuit is a neural motif that is conserved between individuals. The approach and its handling of the degenerate trivial solution are described in [Approach](#approach).

> **Author:** Athmeeya M Kashyap, IIIT Hyderabad

### Submission files

| File | Purpose |
| --- | --- |
| [`network.csv`](network.csv) | Matched neuron table: the three chosen datasets as columns, one matched neuron per row |
| [`README.md`](README.md) | Technical approach, assumptions, and reproduction steps (this file) |
| [`science.md`](science.md) | One-page scientific summary of the conserved circuit |

## Table of contents

- [Result](#result)
- [Approach](#approach)
- [Installation](#installation)
- [Usage](#usage)
- [Repository structure](#repository-structure)
- [Datasets](#datasets)
- [Outputs](#outputs)
- [Reproducibility and testing](#reproducibility-and-testing)
- [Scope and assumptions](#scope-and-assumptions)

## Result

A 4-neuron divergent motif, conserved with identical connectivity across BANC, FAFB, and MANC:

```
              ┌──▶ DNge039   acetylcholine
  DNge076 ────┼──▶ DNge068   glutamate
   (GABA)     └──▶ DNge019   acetylcholine
```

A single GABAergic descending neuron drives three other descending neurons. The induced subgraph and its degree sequence `[(0,3), (1,0), (1,0), (1,0)]` are identical in all three connectomes: a female brain (FAFB), a female brain plus nerve cord (BANC), and a male nerve cord (MANC).

The search found no conserved circuit larger than 4 neurons across 30 connected components of the cross-dataset agreement graph. The circuit is the submission table, [`network.csv`](#outputs). Its biological analysis is in [`science.md`](science.md).

## Approach

**A note on the degenerate solution.** The objective (maximise N subject to mutually isomorphic induced subgraphs) has a trivial solution: any set of neurons with no edges between them induces the empty graph in every dataset, which is trivially isomorphic. Arbitrary pairing then produces a large but meaningless disconnected set. This project avoids that in two ways: neurons are matched by Codex-curated homology rather than arbitrary pairing, and the search targets the largest *connected* conserved circuit rather than raw node count.

```
 FlyWire Codex, BANC codex_annotations (CAVE, materialization v626)
 |   fafb_783_match_id, manc_121_match_id, cell_type, neurotransmitter
 v
 Curated 1:1 correspondence            BANC to FAFB to MANC   (1,218 matched triples)
 v
 Agreement graph                       directed edge u to v iff present in ALL datasets
 v
 Per connected component               remove conflicting nodes (keep the well-wired core),
                                        then take the largest connected piece
 v
 Largest connected, conflict-free conserved circuit   --->  network.csv
```

**1. Biological correspondence.** Neurons are matched across datasets using Codex's curated `fafb_783_match_id` and `manc_121_match_id` annotations, which are versioned to the FAFB v783 and MANC v1.2.1 edge lists. This is a real homology mapping rather than a match on cell-type strings. After verifying every ID against the edge lists and enforcing a strict 1:1 mapping, 1,218 matched triples remain.

**2. Connectivity-first circuit search.** Instead of maximising node count, the algorithm builds an agreement graph, where an edge exists only if it is present in all three datasets. Within each connected component it removes the fewest conflicting nodes, highest-conflict first, keeping densely wired neurons, and returns the largest connected, conflict-free component. The implementation is in `flywire/circuit.py` (`largest_conserved_circuit`), and `circuit_report.txt` records the resulting neurons, edges, and per-dataset degree sequences.

## Installation

Requires Python 3.9 or newer.

```bash
git clone https://github.com/Athmeeya2006/Conserved-Neural-Circuits.git
cd Conserved-Neural-Circuits

python3 -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

| Dependency | Role |
| --- | --- |
| `matplotlib`, `networkx` | Circuit visualization |
| `pandas`, `pyarrow` | Correspondence tables and cached data |
| `caveclient` | FlyWire Codex / CAVE access |

## Usage

The package is driven through four subcommands.

```bash
python -m flywire selftest        # synthetic algorithm check, no data or network
python -m flywire fetch           # pull the Codex correspondence, requires a CAVE token
python -m flywire solve           # build the conserved circuit and all deliverables
python -m flywire search --annot FAFB=fafb.csv BANC=banc.csv MANC=manc.csv
```

<details>
<summary><b>Authenticating <code>fetch</code> with a CAVE token</b></summary>

`fetch` reads BANC annotations from the FlyWire CAVE backend, which needs a one-time token:

```python
from caveclient import CAVEclient
CAVEclient().auth.save_token(token="<your token>", overwrite=True)
```

Create a token (logged in to the same Google account as Codex, with each datastack's terms accepted) at https://global.daf-apis.com/auth/api/v1/create_token. The fetched correspondence is cached, so `fetch` only needs to run once.

</details>

<details>
<summary><b>The <code>search</code> subcommand (alternative methodology)</b></summary>

`search` implements the original approach. It anchors the correspondence on cell types that are singletons in every dataset, then searches all combinations of three or more datasets for the best conserved circuit. It relies on per-dataset annotation files (`--annot NAME=path`) and acts as a fallback when curated match IDs are not available. The `solve` path is preferred.

</details>

## Repository structure

```
Conserved-Neural-Circuits/
├── flywire/                     Python package
│   ├── cli.py                   subcommands: selftest, fetch, solve, search
│   ├── loaders.py               stream edge lists and annotation files
│   ├── cave.py                  Codex/CAVE fetch and curated-match correspondence
│   ├── anchors.py               singleton and correspondence-based anchors
│   ├── isomorphism.py           pairwise consistency and conflict-graph pruning
│   ├── circuit.py               connectivity-first conserved-circuit extraction
│   ├── report.py                circuit_report.txt generation
│   └── pipeline.py              end-to-end solve()
├── tests/test_pipeline.py       synthetic self-test (pytest)
├── network.csv                  submission: matched neurons of the conserved circuit
├── science.md                   one-page scientific summary
├── *_edge_list.csv              five input connectomes (provided separately)
├── requirements.txt
└── .github/workflows/ci.yml     self-test, pytest, ruff
```

## Datasets

Each edge list has the header `source neuron id,target neuron id`. IDs are integers, edges are directed and unweighted.

| Key  | File                       | Nodes   | Edges     | Description |
| ---- | -------------------------- | ------- | --------- | --- |
| BANC | `banc_626_edge_list.csv`   | 112,885 | 2,676,592 | Brain and nerve cord (v626) |
| FAFB | `fafb_783_edge_list.csv`   | 138,584 | 3,732,460 | Full adult female brain, FlyWire (v783) |
| MANC | `manc_1.2.1_edge_list.csv` |  23,641 | 5,305,638 | Male adult nerve cord (v1.2.1) |
| MAOL | `maol_1.1_edge_list.csv`   |  51,669 | 6,484,936 | Male adult optic lobe (v1.1) |
| MCNS | `mcns_0.9_edge_list.csv`   | 165,820 | 6,239,112 | Male central nervous system (v0.9) |

## Outputs

Produced by `python -m flywire solve`:

| File | Committed | Contents |
| --- | :---: | --- |
| `network.csv` | yes | The submission table: the conserved circuit, one column per dataset, one neuron per row |
| `consistent_set_maxN.csv` | yes | Max node-count consistent set, largely edgeless, included for transparency |
| `circuit_report.txt` | yes | Neurons, edges, per-dataset degree sequences, motif, provenance |
| `circuit_visualization.png` | yes | Rendered circuit network graph |
| `codex_correspondence_*.csv`, `codex_annotations_626.parquet` | no | Cached intermediates, regenerate via `fetch` |

`network.csv` holds the connected conserved circuit, not the larger node count. The max-N consistent set is reported separately in `consistent_set_maxN.csv` and is never presented as the result.

## Reproducibility and testing

```bash
python -m flywire selftest        # synthetic end-to-end check
python -m pytest tests/ -q         # unit tests
ruff check flywire/ tests/         # lint, matches CI
```

The self-test builds three synthetic datasets that share a `t1 -> t2 -> t3` chain, with edgeless padding and one inconsistent edge, then asserts that the inconsistent neuron is dropped, the padding is separated, and the circuit is recovered. CI runs the self-test, the unit tests, and linting on Python 3.9, 3.11, and 3.12.

## Scope and assumptions

- **Three datasets.** Codex publishes curated cross-dataset match IDs from BANC to FAFB and MANC only. MAOL and MCNS have no match columns and are not on the FlyWire CAVE server, so the defensible trio is `{BANC, FAFB, MANC}`. No matches were fabricated to reach more datasets.
- **Unweighted edges.** Edge presence, not synapse count, defines isomorphism, following the problem statement. Thresholding by synapse weight could enlarge the circuit but would change the stated problem, so it is left as future work.
- **Curated matches as ground truth.** Codex `*_match_id` entries are treated as homology. Matches that are not 1:1 are discarded rather than resolved heuristically.

