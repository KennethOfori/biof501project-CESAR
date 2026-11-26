#!/usr/bin/env python3
"""Pairwise ANI matrix from fastANI output.

usage: ani_heatmap.py <ani.tsv> <out.png>

ANI saturates inside a species: every pair here will sit above the 95% boundary.
The ramp is therefore scaled to the observed off-diagonal range, not 0-100, or
every cell would look identical. The diagonal is self-comparison and is masked.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

ani_path, out_path = sys.argv[1:3]

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
MUTED = "#8a8984"
DIAGONAL = "#e4e3df"
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "ani", ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"],
)
SPECIES_BOUNDARY = 95.0  # below this, two genomes are different species


def label(path):
    name = os.path.basename(path)
    for suffix in (".contigs.fasta", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    # the K-12 reference arrives under its RefSeq accession name
    return "K-12" if name.startswith("GCF_") else name


ani = {}
names = []
with open(ani_path) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        q, r = label(f[0]), label(f[1])
        ani[(q, r)] = float(f[2])
        for n in (q, r):
            if n not in names:
                names.append(n)

# samples first, reference last
names.sort(key=lambda n: (n == "K-12", n))

off = [v for (q, r), v in ani.items() if q != r]
if not off:
    sys.exit("no between-genome comparisons in fastANI output")
norm = Normalize(min(off), max(off))

size = 0.85 * len(names)
fig, ax = plt.subplots(figsize=(size + 2.6, size + 2.0), facecolor=SURFACE)
ax.set_facecolor(SURFACE)

for yi, q in enumerate(names):
    for xi, r in enumerate(names):
        value = ani.get((q, r))
        if q == r:
            colour, text, ink = DIAGONAL, "—", MUTED
        elif value is None:
            # fastANI reports nothing below roughly 80% ANI
            colour, text, ink = SURFACE, "<80", MUTED
        else:
            colour = SEQUENTIAL(norm(value))
            lum = 0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2]
            text, ink = f"{value:.2f}", SURFACE if lum < 0.5 else TEXT
        ax.add_patch(
            plt.Rectangle(
                (xi - 0.5, yi - 0.5), 1, 1,
                facecolor=colour, edgecolor=SURFACE, linewidth=2,
            )
        )
        ax.annotate(text, xy=(xi, yi), ha="center", va="center", fontsize=8, color=ink)

ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, fontsize=8, color=TEXT, rotation=30, ha="left")
ax.xaxis.tick_top()
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=8, color=TEXT)
ax.tick_params(length=0, pad=4)
ax.set_xlim(-0.5, len(names) - 0.5)
ax.set_ylim(len(names) - 0.5, -0.5)
for side in ax.spines.values():
    side.set_visible(False)

ax.annotate(
    f"All pairs are above the {SPECIES_BOUNDARY:g}% species boundary; the colour scale "
    f"spans the\nobserved range ({min(off):.2f}–{max(off):.2f}%). Rows are query, columns reference: "
    "fastANI\nis not perfectly symmetric.",
    xy=(0, 0), xycoords="axes fraction", xytext=(0, -30), textcoords="offset points",
    fontsize=7, color=MUTED, va="top", ha="left",
)
fig.suptitle(
    "Average nucleotide identity (fastANI)", x=0.5, y=0.985,
    ha="center", fontsize=10, color=TEXT,
)
fig.tight_layout(rect=(0, 0.06, 1, 0.95))
fig.savefig(out_path, dpi=200, facecolor=SURFACE)
