#!/usr/bin/env python3
"""BED of K-12 resistance-relevant loci, for per-gene coverage and variant counts.

usage: make_locus_bed.py <genomic.gff> <out.bed>

Only loci that exist in BOTH K-12 and the isolates belong here: comparing an
acquired gene to a reference that lacks it says nothing. Name column is
"<group>::<gene>" so downstream steps can keep the grouping.
"""
import re
import sys

gff_path, out_path = sys.argv[1:3]

GROUPS = {
    # porin loss is the classic carbapenem-resistance route with no carbapenemase
    "Porins & envelope": ["ompF", "ompC", "ompA", "phoE", "lamB", "envZ", "ompR"],
    # mutations here de-repress AcrAB-TolC
    "Efflux regulators": ["marR", "marA", "acrR", "soxR", "soxS", "rob", "envR"],
    "Efflux pumps": ["acrA", "acrB", "tolC", "acrE", "acrF", "mdfA", "emrE", "emrD"],
    "Fluoroquinolone targets": ["gyrA", "gyrB", "parC", "parE"],
    "Penicillin-binding proteins": ["ftsI", "mrdA", "mrcA", "mrcB"],
    "Colistin / LPS": ["basR", "basS", "phoP", "phoQ", "mgrB", "eptA", "arnA"],
    # nitrofurantoin resistance follows loss of either nitroreductase
    "Nitrofurantoin": ["nfsA", "nfsB", "ribE"],
    "Fosfomycin": ["glpT", "uhpT", "murA"],
    "Other targets": ["rpoB", "folP", "cyaA", "cirA"],
}
WANTED = {gene: group for group, genes in GROUPS.items() for gene in genes}

rows = {}
with open(gff_path) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"(?:^|;)gene=([^;]+)", f[8])
        if not m or m.group(1) not in WANTED:
            continue
        gene = m.group(1)
        # GFF is 1-based inclusive; BED is 0-based half-open
        rows[gene] = (f[0], int(f[3]) - 1, int(f[4]), f"{WANTED[gene]}::{gene}")

missing = sorted(set(WANTED) - set(rows))
if missing:
    print(f"not found in GFF, skipped: {', '.join(missing)}", file=sys.stderr)

with open(out_path, "w") as out:
    for chrom, start, end, name in sorted(rows.values(), key=lambda r: r[1]):
        out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
print(f"wrote {len(rows)} loci", file=sys.stderr)
