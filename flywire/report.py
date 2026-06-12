"""
Human-readable deliverables: circuit_report.txt and METHODS.md.
"""


def write_circuit_report(circuit, ds_list, per_ds_edges, type_of, nt_of,
                         nt_source_of, n_star, iso_ok,
                         path='circuit_report.txt'):
    """Write the per-circuit report (neurons, edges, degree sequences, motif)."""
    edges = sorted(per_ds_edges[ds_list[0]])

    def degree_seq(ds):
        outd = [0] * len(circuit)
        ind = [0] * len(circuit)
        for i, j in per_ds_edges[ds]:
            outd[i] += 1
            ind[j] += 1
        return sorted(zip(ind, outd))

    degseqs = {ds: degree_seq(ds) for ds in ds_list}
    recurrent = any((j, i) in set(edges) for i, j in edges)
    motif = ("recurrent (contains a reciprocal/loop edge)" if recurrent
             else "feedforward (no reciprocal edges)")

    with open(path, 'w') as fh:
        fh.write(f"CONSERVED CIRCUIT REPORT ({' / '.join(ds_list)})\n")
        fh.write("=" * 52 + "\n\n")
        fh.write(f"Consistent mutually-isomorphic set N : {n_star}\n")
        fh.write(f"Largest connected conserved circuit  : {len(circuit)} "
                 f"neurons, {len(edges)} directed edges\n")
        fh.write(f"Induced edges identical across all   : {iso_ok}\n\n")

        fh.write("Neurons in the circuit\n" + "-" * 52 + "\n")
        hdr = f"{'idx':<4}{'cell_type':<12}{'neurotransmitter':<18}{'nt_source':<12}"
        hdr += "".join(f"{ds:<20}" for ds in ds_list)
        fh.write(hdr.rstrip() + "\n")
        for i, t in enumerate(circuit):
            key = t[ds_list[0]]
            line = (f"{i:<4}{type_of.get(key, '?'):<12}"
                    f"{(nt_of.get(key) or '?'):<18}"
                    f"{(nt_source_of.get(key) or '?'):<12}")
            line += "".join(f"{t[ds]:<20}" for ds in ds_list)
            fh.write(line.rstrip() + "\n")

        fh.write("\nDirected edges (cell-type -> cell-type, present in all)\n")
        fh.write("-" * 52 + "\n")
        for i, j in edges:
            ki, kj = circuit[i][ds_list[0]], circuit[j][ds_list[0]]
            fh.write(f"  {type_of.get(ki, '?')} ({i}) -> "
                     f"{type_of.get(kj, '?')} ({j})\n")

        fh.write("\nDegree sequence (in,out) sorted, per dataset "
                 "(identical => isomorphic)\n")
        fh.write("-" * 52 + "\n")
        for ds in ds_list:
            fh.write(f"  {ds}: {degseqs[ds]}\n")
        fh.write(f"  all identical: "
                 f"{all(degseqs[d] == degseqs[ds_list[0]] for d in ds_list)}\n")
        fh.write(f"\nMotif: {motif}; {len(edges)} edges over "
                 f"{len(circuit)} neurons.\n")

        fh.write("\nProvenance\n" + "-" * 52 + "\n")
        fh.write(
            "Edges:  banc_626_edge_list.csv, fafb_783_edge_list.csv,\n"
            "        manc_1.2.1_edge_list.csv (directed, unweighted).\n"
            "Types/NT/matches: FlyWire Codex CAVE datastack\n"
            "        brain_and_nerve_cord_public, materialization v626,\n"
            "        table codex_annotations (cell_type,\n"
            "        neurotransmitter_verified/_predicted,\n"
            "        fafb_783_match_id, manc_121_match_id). See flywire/cave.py.\n")
    print(f"  wrote {path}")
