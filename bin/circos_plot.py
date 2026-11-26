#!/usr/bin/env python3
"""Usage: circos_plot.py <sample_id> <mosdepth.regions.bed.gz> <variants.vcf.gz>"""
import sys
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from pycirclize import Circos

sample, cov_file, vcf_file = sys.argv[1:4]

cov = pd.read_csv(cov_file, sep="\t", header=None, names=["chrom", "start", "end", "depth"])
vcf = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, usecols=[0, 1, 3, 4],
                  names=["chrom", "pos", "ref", "alt"])
vcf["is_snp"] = (vcf["ref"].str.len() == 1) & (vcf["alt"].str.len() == 1)

sizes = cov.groupby("chrom")["end"].max().to_dict()
# leave a 30 degree gap at the top of the circle to hold the ring labels
GAP = 30
circos = Circos(sizes, start=GAP / 2, end=360 - GAP / 2, space=0 if len(sizes) == 1 else 3)
circos.text(sample, size=14)

# ring labels, written in the gap (r = middle of each ring's radius range)
BIN = 10_000  # bin size for the SNP and indel density histograms
for label, r in [("Genome", 96.5), ("Read depth", 77.5), ("SNPs / 10 kb", 55.5), ("Indels / 10 kb", 38.5)]:
    circos.text(label, r=r, deg=0, size=7, orientation="horizontal")


def density_track(sector, positions, r_lim, color):
    """Histogram ring: variants per BIN, with its own 0..max scale."""
    edges = np.append(np.arange(0, sector.size, BIN), sector.size)  # last bin may be short
    counts, _ = np.histogram(positions, bins=edges)
    vmax = max(int(counts.max()), 1)
    track = sector.add_track(r_lim, r_pad_ratio=0.05)
    track.axis()
    track.yticks([0, vmax], labels=["0", str(vmax)], vmin=0, vmax=vmax, side="left", label_size=6)
    track.bar((edges[:-1] + edges[1:]) / 2, counts, width=np.diff(edges), vmin=0, vmax=vmax,
              color=color, lw=0)

for sector in circos.sectors:
    c = cov[cov["chrom"] == sector.name]
    v = vcf[vcf["chrom"] == sector.name]
    mid = ((c["start"] + c["end"]) / 2).to_numpy()
    depth = c["depth"].to_numpy()
    cap = np.percentile(depth, 99)            # stop single spikes squashing the track

    # outer ring: genome axis with Mb ticks
    outer = sector.add_track((95, 98))
    outer.axis(fc="black")
    outer.xticks_by_interval(500_000, label_formatter=lambda x: f"{x/1e6:.1f} Mb",
                             label_orientation="vertical")

    # ring 2: read depth (red bars = windows with ~no coverage -> likely deletions)
    cov_track = sector.add_track((65, 90), r_pad_ratio=0.05)
    cov_track.axis()
    cov_track.yticks([0, cap / 2, cap], labels=["0", f"{cap/2:.0f}x", f"{cap:.0f}x"],
                     vmin=0, vmax=cap, side="left", label_size=6)
    cov_track.fill_between(mid, np.minimum(depth, cap), vmin=0, vmax=cap, color="steelblue")
    # shade each low-coverage window across the full height of the ring. Below 10% of the
    # sample's median depth is treated as "effectively absent": loose enough to tolerate
    # normal coverage dips, strict enough that only likely deletions are flagged.
    low =c[c["depth"] < 0.1 * np.median(depth)]
    for start, end in zip(low["start"], low["end"]):
        cov_track.rect(start, end, color="red", alpha=0.6, lw=0)

    snps, indels = v[v["is_snp"]], v[~v["is_snp"]]

    # rings 3 and 4: SNP and indel density (per 10 kb bin) -- shows hotspots that a dot ring
    # hides. Each ring has its own scale: indels are far rarer than SNPs.
    density_track(sector, snps["pos"], (48, 63), "darkorange")
    density_track(sector, indels["pos"], (31, 46), "purple")

fig = circos.plotfig()
legend = [
    Patch(color="steelblue", label="Read depth (5 kb windows)"),
    Patch(color="red", alpha=0.6, label="Low coverage (<10% of median)"),
    Patch(color="darkorange", label="SNP density (SNPs per 10 kb)"),
    Patch(color="purple", label="Indel density (indels per 10 kb)"),
]
circos.ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)
fig.savefig(f"{sample}.circos.png", dpi=300, bbox_inches="tight")
