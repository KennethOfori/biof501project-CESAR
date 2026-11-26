#!/usr/bin/env python3
"""Linear maps of the contigs that carry clustered resistance genes.

usage: resistance_region_plot.py <amr_location_summary.tsv> <out.png>

Short-read contigs are fragments, not closed replicons, so each is drawn as its
own linear bar to scale. Only contigs carrying >= MIN_GENES resistance genes are
shown: those are the integron / IS-flanked multidrug regions.
"""
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrow, Rectangle

summary_path, out_path = sys.argv[1:3]

MIN_GENES = 2  # contigs with fewer resistance genes are not "regions"
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e4e3df"
COLOUR = {"chromosome": "#2a78d6", "plasmid": "#eb6834"}

NODE_RE = re.compile(r"length_(\d+)_cov_([\d.]+)")

df = pd.read_csv(summary_path, sep="\t", dtype=str).fillna("-")
# acquired resistance genes, plus qacEdelta1: the 3' marker of a class 1 integron.
# Metal resistance (the chromosomal pco/sil island) is a separate story and its
# 102 kb contig would flatten the scale here.
df = df[(df["Subtype"] == "AMR") | (df["Element symbol"] == "qacEdelta1")]
df[["Start", "Stop"]] = df[["Start", "Stop"]].astype(int)

regions = []
for (sample, contig), hits in df.groupby(["Name", "Contig id"]):
    if len(hits) < MIN_GENES:
        continue
    m = NODE_RE.search(contig)
    if not m:
        continue
    first = hits.iloc[0]
    regions.append(
        {
            "sample": sample,
            "contig": contig,
            "length": int(m.group(1)),
            "cov": float(m.group(2)),
            "location": first["location"],
            "cluster": first["plasmid_cluster"],
            "reps": first["plasmid_rep_types"],
            "genes": hits.sort_values("Start")[
                ["Element symbol", "Start", "Stop", "Strand"]
            ].values.tolist(),
        }
    )

# tallest, gene-richest regions first within each sample
regions.sort(key=lambda r: (r["sample"], -len(r["genes"]), -r["length"]))

max_len = max(r["length"] for r in regions)
row_h = 0.78
fig_h = 1.5 + row_h * len(regions)
fig, ax = plt.subplots(figsize=(9.5, fig_h), facecolor=SURFACE)
ax.set_facecolor(SURFACE)

for i, reg in enumerate(regions):
    y = -i
    colour = COLOUR[reg["location"]]
    # the contig itself, to scale
    ax.add_patch(
        Rectangle((0, y - 0.1), reg["length"], 0.2, facecolor=GRID, edgecolor="none", zorder=1)
    )
    for j, (gene, start, stop, strand) in enumerate(reg["genes"]):
        x, dx = (start, stop - start) if strand != "-" else (stop, start - stop)
        head = min(abs(dx) * 0.4, max_len * 0.006)
        ax.add_patch(
            FancyArrow(
                x, y, dx, 0, width=0.26, head_width=0.26, head_length=head,
                length_includes_head=True, facecolor=colour, edgecolor=SURFACE,
                linewidth=0.8, zorder=3,
            )
        )
        # stagger labels over two levels so adjacent genes do not collide
        level = y + (0.19 if j % 2 == 0 else 0.42)
        ax.plot(
            [(start + stop) / 2] * 2, [y + 0.13, level - 0.02],
            color=GRID, lw=0.6, zorder=2,
        )
        ax.annotate(
            gene, xy=((start + stop) / 2, level), ha="center", va="bottom",
            fontsize=7, color=TEXT, zorder=4,
        )

    label = f"{reg['sample']}  {reg['contig'].split('_cov')[0]}"
    detail = f"{reg['length'] / 1000:.1f} kb · {reg['cov']:.0f}×"
    if reg["location"] == "plasmid":
        reps = reg["reps"] if reg["reps"] != "-" else "no replicon"
        detail += f" · plasmid {reg['cluster']} ({reps.split(',')[0]})"
    else:
        detail += " · unplaced / chromosome"
    ax.annotate(
        label, xy=(0, y + 0.06), xycoords=("figure fraction", "data"), xytext=(8, 0),
        textcoords="offset points", fontsize=7.5, color=TEXT, va="bottom", ha="left",
    )
    ax.annotate(
        detail, xy=(0, y - 0.06), xycoords=("figure fraction", "data"), xytext=(8, 0),
        textcoords="offset points", fontsize=6.5, color=MUTED, va="top", ha="left",
    )

ax.set_xlim(-max_len * 0.01, max_len * 1.02)
ax.set_ylim(-len(regions) + 0.35, 0.95)
ax.set_yticks([])
ax.set_xlabel("position on contig (kb)", fontsize=8, color=TEXT_2)
ax.set_xticks([t for t in range(0, max_len + 1, 1000)])
ax.set_xticklabels([f"{t // 1000}" for t in range(0, max_len + 1, 1000)], fontsize=7.5, color=TEXT_2)
ax.tick_params(length=0, pad=2)
for name, side in ax.spines.items():
    side.set_visible(False)
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
ax.set_axisbelow(True)

handles = [
    Rectangle((0, 0), 1, 1, facecolor=COLOUR[k], edgecolor="none")
    for k in ("chromosome", "plasmid")
]
fig.legend(
    handles=handles, labels=["Chromosome / unplaced contig", "Plasmid contig"],
    loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=2, frameon=False,
    fontsize=8, labelcolor=TEXT, handlelength=1.1, handletextpad=0.5, columnspacing=1.4,
)
fig.suptitle(
    "Clustered resistance genes on assembled contigs", x=0.99, y=0.995, ha="right",
    fontsize=10, color=TEXT,
)
fig.subplots_adjust(left=0.26, right=0.98, top=1 - 0.75 / fig_h, bottom=0.75 / fig_h)
fig.savefig(out_path, dpi=200, facecolor=SURFACE)
