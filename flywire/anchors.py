"""
Correspondence construction (biological anchoring).
"""

from collections import defaultdict


def singleton_anchors(annots, ds_list, require_nt=False):
    """
    Build unambiguous cross-dataset correspondences from cell types that occur
    exactly once in each of the chosen datasets (singletons). Optionally also
    require the neurotransmitter to agree across datasets.

    Returns a list of dicts: {ds: neuron_id, '_type': cell_type}.
    """
    # cell type -> ds -> list of (neuron_id, nt)
    by_type = defaultdict(lambda: defaultdict(list))
    for ds in ds_list:
        for nid, meta in annots[ds].items():
            by_type[meta['type']][ds].append((nid, meta.get('nt')))

    anchors = []
    for ctype, perds in by_type.items():
        if not all(ds in perds for ds in ds_list):
            continue
        if not all(len(perds[ds]) == 1 for ds in ds_list):
            continue  # ambiguous: skip to keep the bijection rigorous
        nts = [perds[ds][0][1] for ds in ds_list]
        if require_nt and len({x for x in nts if x}) > 1:
            continue  # neurotransmitters disagree, reject as homolog
        anchor = {ds: perds[ds][0][0] for ds in ds_list}
        anchor['_type'] = ctype
        anchors.append(anchor)
    return anchors


def anchors_from_correspondence(df, ds_list):
    """
    Build anchor dicts from a curated correspondence table (e.g. Codex
    cross-dataset match ids). ``df`` has one integer id column per dataset plus
    a ``cell_type`` column. Returns a list of {ds: neuron_id, '_type': type}.
    """
    anchors = []
    for _, row in df.iterrows():
        anchor = {ds: int(row[ds]) for ds in ds_list}
        anchor['_type'] = str(row['cell_type'])
        anchors.append(anchor)
    return anchors
