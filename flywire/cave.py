"""
FlyWire Codex (CAVE) access: fetch BANC annotations and build the curated
cross-dataset neuron correspondence.

The BANC ``codex_annotations`` table (datastack ``brain_and_nerve_cord_public``,
materialization v626, matching ``banc_626_edge_list.csv``) is entity-attribute:
each row is ``(pt_root_id, classification_system, value)`` where ``value`` is
stored in the ``cell_type`` column. The classification systems used here:

  ``cell_type``                   BANC cell type (e.g. DNp32)
  ``neurotransmitter_verified``   curated neurotransmitter (preferred)
  ``neurotransmitter_predicted``  predicted neurotransmitter (fallback)
  ``fafb_783_match_id``           Codex-curated match to a FAFB v783 root id
  ``manc_121_match_id``           Codex-curated match to a MANC v1.2.1 body id

``fafb_783_match_id`` / ``manc_121_match_id`` are Codex's own curated
neuron-level homology assignments, versioned exactly to the FAFB v783 and MANC
v1.2.1 edge lists. Using them as the isomorphism's node mapping fixes the
cross-dataset correspondence by biology rather than arbitrary pairing.

A CAVE token is required::

    from caveclient import CAVEclient
    CAVEclient().auth.save_token(token="<token>", overwrite=True)

Get one (logged in to the same Google account as Codex, datastack terms
accepted) at https://global.daf-apis.com/auth/api/v1/create_token
"""

import collections
import os

from .loaders import load_edges

DATASTACK = 'brain_and_nerve_cord_public'
MAT_VERSION = 626
ANNOT_CACHE = 'codex_annotations_626.parquet'
CORRESPONDENCE_CSV = 'codex_correspondence_banc_fafb_manc.csv'

# Edge lists the match ids index into.
MATCH_EDGE_FILES = {
    'BANC': 'banc_626_edge_list.csv',
    'FAFB': 'fafb_783_edge_list.csv',
    'MANC': 'manc_1.2.1_edge_list.csv',
}


def fetch_codex_annotations(cache=ANNOT_CACHE):
    """Return the BANC ``codex_annotations`` table, cached to parquet."""
    import pandas as pd
    if cache and os.path.exists(cache):
        print(f"  using cached {cache}")
        return pd.read_parquet(cache)
    from caveclient import CAVEclient
    client = CAVEclient(DATASTACK)
    client.materialize.version = MAT_VERSION
    print(f"  querying {DATASTACK} codex_annotations @ v{MAT_VERSION} ...")
    df = client.materialize.query_table('codex_annotations')
    if cache:
        df.to_parquet(cache)
        print(f"  cached {len(df):,} rows to {cache}")
    return df


def build_correspondence(df, edge_files=None, out_csv=CORRESPONDENCE_CSV):
    """
    Turn the raw annotation table into a strict 1:1 BANC/FAFB/MANC
    correspondence. Returns a DataFrame with columns
    ``BANC, FAFB, MANC, cell_type, neurotransmitter`` and writes ``out_csv``.
    """
    import pandas as pd
    edge_files = edge_files or MATCH_EDGE_FILES

    def values(system):
        s = df[df['classification_system'] == system][['pt_root_id', 'cell_type']]
        return dict(zip(s['pt_root_id'], s['cell_type']))

    cell_type = values('cell_type')
    fafb_match = values('fafb_783_match_id')
    manc_match = values('manc_121_match_id')
    nt_pred = values('neurotransmitter_predicted')
    nt_ver = values('neurotransmitter_verified')

    print("  verifying matched ids against the edge lists ...")
    nodes = {ds: load_edges(edge_files[ds])[2] for ds in ('BANC', 'FAFB', 'MANC')}

    rows = []
    for banc_id, fafb_raw in fafb_match.items():
        if banc_id not in manc_match:
            continue
        try:
            bid, fid, mid = int(banc_id), int(fafb_raw), int(manc_match[banc_id])
        except (TypeError, ValueError):
            continue
        if (bid not in nodes['BANC'] or fid not in nodes['FAFB']
                or mid not in nodes['MANC']):
            continue
        nt = (nt_ver.get(banc_id) or nt_pred.get(banc_id) or '').strip()
        rows.append((bid, fid, mid, str(cell_type.get(banc_id, '')).strip(), nt))
    print(f"  raw valid triples: {len(rows)}")

    for col in range(3):  # strict 1:1 in BANC, FAFB, MANC
        counts = collections.Counter(r[col] for r in rows)
        dup = {k for k, v in counts.items() if v > 1}
        rows = [r for r in rows if r[col] not in dup]
    print(f"  strict 1:1 bijection: {len(rows)}")

    out = pd.DataFrame(rows, columns=['BANC', 'FAFB', 'MANC',
                                      'cell_type', 'neurotransmitter'])
    if out_csv:
        out.to_csv(out_csv, index=False)
        print(f"  wrote {out_csv}: {len(out)} matched neuron triples")
    return out


def fetch_and_build(out_csv=CORRESPONDENCE_CSV):
    """Convenience: fetch the table and build the correspondence in one call."""
    return build_correspondence(fetch_codex_annotations(), out_csv=out_csv)
