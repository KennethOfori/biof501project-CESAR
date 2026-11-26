#!/usr/bin/env python3
"""Gene x sample presence grid from amr_location_summary.tsv.

usage: amr_location_heatmap.py <amr_location_summary.tsv> <out.png>

Each cell marks a gene/mutation found in a sample; colour and shape say whether it
sits on the chromosome or on a MOB-suite plasmid. Rows are grouped by AMRFinder type.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

summary_path, out_path = sys.argv[1:3]

# palette: categorical slots 1-2 (validated pair) + text/grid tokens
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
STYLE = {  # location -> (colour, marker)
    "chromosome": ("#2a78d6", "s"),
    "plasmid": ("#eb6834", "D"),
}

GROUPS = [  # (label, test on a row) in display order
    ("Resistance genes", lambda r: r.Type == "AMR" and r.Subtype == "AMR"),
    ("Resistance\npoint mutations", lambda r: r.Type == "AMR" and r.Subtype.startswith("POINT")),
    ("Metal resistance", lambda r: r.Type == "STRESS" and r.Subtype == "METAL"),
    ("Biocide / acid\nresistance", lambda r: r.Type == "STRESS" and r.Subtype != "METAL"),
    ("Virulence", lambda r: r.Type == "VIRULENCE"),
]

df = pd.read_csv(summary_path, sep="\t", dtype=str).fillna("-")
df["group"] = None
for label, test in GROUPS:
    df.loc[df.apply(test, axis=1), "group"] = label
df = df.dropna(subset=["group"])

# one cell per gene x sample; a gene with any copy on a plasmid counts as plasmid
cells = (
    df.groupby(["group", "Element symbol", "Name"])["location"]
    .agg(lambda s: "plasmid" if (s == "plasmid").any() else "chromosome")
    .reset_index()
)

samples = sorted(df["Name"].unique())
group_order = [g for g, _ in GROUPS if g in set(cells["group"])]

# row order: by group, then most widespread first, then name
rows = []
for g in group_order:
    sub = cells[cells["group"] == g]
    counts = sub.groupby("Element symbol")["Name"].nunique()
    genes = sorted(counts.index, key=lambda gene: (-counts[gene], gene.lower()))
    rows.append((g, genes))

# y positions with a one-row gap between groups
y_of, y, bands = {}, 0, []
for g, genes in rows:
    start = y
    for gene in genes:
        y_of[gene] = y
        y += 1
    bands.append((g, start, y - 1))
    y += 1
n_rows = y - 1

fig_h = 1.6 + 0.19 * n_rows
fig, ax = plt.subplots(figsize=(6.4, fig_h), facecolor=SURFACE)
ax.set_facecolor(SURFACE)

# recessive guide lines for each column and row
for x in range(len(samples)):
    ax.axvline(x, color=GRID, lw=0.6, zorder=0)
for gene, yy in y_of.items():
    ax.axhline(yy, color=GRID, lw=0.4, zorder=0)

for loc, (colour, marker) in STYLE.items():
    sub = cells[cells["location"] == loc]
    ax.scatter(
        [samples.index(s) for s in sub["Name"]],
        [y_of[g] for g in sub["Element symbol"]],
        s=58, c=colour, marker=marker, edgecolors=SURFACE, linewidths=1.2, zorder=3,
    )

ax.set_yticks(list(y_of.values()))
ax.set_yticklabels(list(y_of.keys()), fontsize=7, color=TEXT)
ax.set_xticks(range(len(samples)))
ax.set_xticklabels(samples, fontsize=8, color=TEXT, rotation=30, ha="left")
ax.xaxis.tick_top()
ax.tick_params(length=0, pad=3)
ax.set_xlim(-0.6, len(samples) - 0.4)
ax.set_ylim(n_rows + 0.8, -0.8)
for side in ax.spines.values():
    side.set_visible(False)

# group labels and separators in the left margin
for g, start, end in bands:
    ax.annotate(
        g, xy=(0, (start + end) / 2), xycoords=("figure fraction", "data"),
        xytext=(8, 0), textcoords="offset points",
        fontsize=8, fontweight="bold", color=TEXT_2, va="center", ha="left",
    )
    if start > 0:
        ax.axhline(start - 1, color=TEXT_2, lw=0.8, zorder=1)

handles = [
    Line2D([], [], ls="", marker=m, ms=7, mfc=c, mec=SURFACE, label=loc.capitalize())
    for loc, (c, m) in STYLE.items()
]
fig.legend(
    handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=2,
    frameon=False, fontsize=8, labelcolor=TEXT, handletextpad=0.3, columnspacing=1.2,
)
fig.suptitle(
    "AMRFinderPlus hits by sample and location", x=0.99, y=0.995, ha="right",
    fontsize=10, color=TEXT,
)

fig.subplots_adjust(left=0.42, right=0.88, top=1 - 1.05 / fig_h, bottom=0.25 / fig_h)
fig.savefig(out_path, dpi=200, facecolor=SURFACE)
