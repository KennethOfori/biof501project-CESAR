#!/usr/bin/env python3
"""Per-locus depth and variant count for one sample, against K-12.

usage: locus_table.py <sample_id> <loci.bed> <mosdepth.regions.bed.gz> \
                      <mosdepth.summary.txt> <sample.vcf.gz> <out.tsv>

depth_ratio normalises each locus against the sample's genome-wide mean, so the
value is comparable across samples: ~1 means present at normal copy number,
near 0 means the locus is deleted or unmappable.
"""
import gzip
import sys

sample_id, bed_path, regions_path, summary_path, vcf_path, out_path = sys.argv[1:7]

loci = []  # (chrom, start, end, group, gene)
with open(bed_path) as fh:
    for line in fh:
        chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
        group, gene = name.split("::", 1)
        loci.append((chrom, int(start), int(end), group, gene))

# mosdepth mean per region, keyed on coordinates (its name column is not relied on)
depth = {}
with gzip.open(regions_path, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        depth[(f[0], int(f[1]), int(f[2]))] = float(f[-1])

genome_mean = None
with open(summary_path) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] == "total":
            genome_mean = float(f[3])
if not genome_mean:
    sys.exit("no genome-wide mean depth in mosdepth summary")

# variant positions from the VCF; counted into loci below
variants = []
with gzip.open(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 2)
        variants.append((f[0], int(f[1])))

with open(out_path, "w") as out:
    out.write("sample\tgroup\tgene\tdepth_ratio\tn_variants\n")
    for chrom, start, end, group, gene in loci:
        mean = depth.get((chrom, start, end))
        if mean is None:
            continue
        n = sum(1 for c, pos in variants if c == chrom and start < pos <= end)
        out.write(f"{sample_id}\t{group}\t{gene}\t{mean / genome_mean:.3f}\t{n}\n")
