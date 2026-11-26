#!/usr/bin/env python3
"""Label each AMRFinder hit with where it sits: chromosome or a MOB-suite plasmid.

usage: amr_plasmid_location.py <sample_id> <amrfinder.tsv> <contig_report.txt> <mobtyper_results.txt> <out.tsv>

Joins on contig id. Hits on contigs MOB-recon could not place are reported as
"unassigned" rather than guessed.
"""
import csv
import sys

sample_id, amr_path, contig_path, mobtyper_path, out_path = sys.argv[1:6]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# contig id -> (molecule_type, primary_cluster_id, rep_type)
contigs = {}
for row in read_tsv(contig_path):
    contigs[row["contig_id"]] = (
        row["molecule_type"],
        row["primary_cluster_id"],
        row.get("rep_type(s)", "-"),
    )

# plasmid cluster -> (predicted_mobility, rep_type for the whole plasmid)
# mobtyper sample_id looks like "<prefix>:<primary_cluster_id>"
plasmids = {}
for row in read_tsv(mobtyper_path):  # empty file when no plasmids were found
    cluster = row["sample_id"].rsplit(":", 1)[-1]
    plasmids[cluster] = (row["predicted_mobility"], row["rep_type(s)"])

amr_rows = read_tsv(amr_path)
extra = ["location", "plasmid_cluster", "plasmid_rep_types", "plasmid_mobility"]
fields = (list(amr_rows[0].keys()) if amr_rows else ["Name", "Contig id"]) + extra

with open(out_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for hit in amr_rows:
        mol, cluster, _ = contigs.get(hit["Contig id"], ("unassigned", "-", "-"))
        if mol == "plasmid":
            mobility, reps = plasmids.get(cluster, ("-", "-"))
        else:
            cluster, mobility, reps = "-", "-", "-"
        hit.update(
            location=mol,
            plasmid_cluster=cluster,
            plasmid_rep_types=reps,
            plasmid_mobility=mobility,
        )
        writer.writerow(hit)
