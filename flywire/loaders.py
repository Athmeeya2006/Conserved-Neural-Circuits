"""
Edge-list and annotation loading.

Integer IDs and dict-of-sets adjacency keep memory and set operations fast even
on the 100 MB+ edge lists.
"""

import csv
from collections import defaultdict


# Local edge lists shipped with this repo. Used when --edges is not given.
DEFAULT_EDGE_FILES = {
    'BANC': 'banc_626_edge_list.csv',
    'FAFB': 'fafb_783_edge_list.csv',
    'MANC': 'manc_1.2.1_edge_list.csv',
    'MAOL': 'maol_1.1_edge_list.csv',
    'MCNS': 'mcns_0.9_edge_list.csv',
}


def detect_columns(header_fields):
    """Return (src_idx, tgt_idx) by matching common naming patterns."""
    src = tgt = None
    for i, field in enumerate(header_fields):
        name = field.strip().lower()
        if src is None and any(k in name for k in ('source', 'pre', 'src', 'from')):
            src = i
        if tgt is None and any(k in name for k in ('target', 'post', 'tgt', 'to')):
            tgt = i
    if src is None or tgt is None:
        src, tgt = 0, 1  # fallback: first two columns
    return src, tgt


def load_edges(path):
    """
    Stream a possibly large edge-list CSV into directed adjacency.

    Returns (out_adj, in_adj, nodes) where out_adj and in_adj are
    dict[int, set[int]] and nodes is a set[int].
    """
    print(f"  loading edges: {path}", flush=True)
    out_adj = defaultdict(set)
    in_adj = defaultdict(set)
    nodes = set()

    with open(path, 'r', newline='') as fh:
        first = fh.readline()
        if not first:
            raise ValueError(f"empty file: {path}")
        fields = first.rstrip('\n').split(',')
        header_like = any(c.isalpha() for c in first)
        if header_like:
            si, ti = detect_columns(fields)
            print(f"    columns: source[{si}]='{fields[si].strip()}', "
                  f"target[{ti}]='{fields[ti].strip()}'")
            data_iter = fh
        else:
            si, ti = 0, 1
            data_iter = _chain_first(first, fh)

        n_edges = 0
        for line in data_iter:
            if not line:
                continue
            parts = line.rstrip('\n').split(',')
            if len(parts) <= max(si, ti):
                continue
            try:
                s = int(parts[si])
                t = int(parts[ti])
            except ValueError:
                continue
            out_adj[s].add(t)
            in_adj[t].add(s)
            nodes.add(s)
            nodes.add(t)
            n_edges += 1
            if n_edges % 1_000_000 == 0:
                print(f"    ... {n_edges:,} edges", flush=True)

    print(f"    {len(nodes):,} nodes, {n_edges:,} edges")
    return dict(out_adj), dict(in_adj), nodes


def _chain_first(first_line, rest):
    yield first_line
    for line in rest:
        yield line


def load_annotations(path):
    """
    Return {neuron_id: {'type': str, 'nt': str|None}} from a metadata CSV with
    auto-detected columns.
    """
    print(f"  loading annotations: {path}", flush=True)
    with open(path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        hl = [h.strip().lower() for h in header]

        def find(keys):
            for i, h in enumerate(hl):
                if any(k in h for k in keys):
                    return i
            return None

        id_i = find(('root_id', 'neuron_id', 'root id', 'neuron id', 'id'))
        ty_i = find(('cell_type', 'cell type', 'type', 'label', 'class'))
        nt_i = find(('neurotransmitter', 'nt', 'transmitter'))

        if id_i is None or ty_i is None:
            raise ValueError(
                f"could not find id/type columns in {path}; header={header}")

        ann = {}
        for row in reader:
            if len(row) <= max(id_i, ty_i):
                continue
            try:
                nid = int(row[id_i])
            except ValueError:
                continue
            ctype = row[ty_i].strip()
            if not ctype or ctype.lower() in ('nan', 'none', 'unknown', 'na'):
                continue
            nt = None
            if nt_i is not None and len(row) > nt_i:
                nt = row[nt_i].strip() or None
            ann[nid] = {'type': ctype, 'nt': nt}

    print(f"    {len(ann):,} annotated neurons")
    return ann


def fetch_cave_annotations(dataset):
    """
    Optional fallback via CAVEclient. Datastacks and table schemas vary by
    dataset and require an auth token (set up with client.auth.save_token).
    Returns {} on any failure so the run can continue with whatever file-based
    annotations exist.
    """
    datastacks = {
        'FAFB': 'flywire_fafb_public',
        'BANC': 'brain_and_nerve_cord',
        'MANC': 'manc_public_v1_0',
        'MCNS': 'flywire_male_public',
        'MAOL': 'flywire_fafb_public',
    }
    name = datastacks.get(dataset)
    if not name:
        return {}
    try:
        from caveclient import CAVEclient
        client = CAVEclient(name)
        tables = client.materialize.get_tables()
        for tbl in ('cell_info', 'cell_type_local', 'neuron_information',
                    'codex_annotations', 'classification'):
            if tbl in tables:
                df = client.materialize.query_table(tbl)
                ann = {}
                cols = {c.lower(): c for c in df.columns}
                idc = next((cols[c] for c in cols
                            if 'root_id' in c or c == 'id'), None)
                tyc = next((cols[c] for c in cols
                            if 'cell_type' in c or c == 'type'), None)
                if idc and tyc:
                    for _, r in df.iterrows():
                        try:
                            ann[int(r[idc])] = {'type': str(r[tyc]), 'nt': None}
                        except Exception:
                            pass
                if ann:
                    return ann
    except Exception as e:
        print(f"    CAVE fetch failed for {dataset}: {e}")
    return {}
