#!/usr/bin/env python3
"""Resistance-relevant K-12 loci: coverage and variant burden per sample.

usage: locus_heatmap.py <locus_summary.tsv> <out.png>

Left panel diverges around 1.0 (the sample's own genome-wide mean depth), so a
deleted or disrupted locus reads red and a duplicated one blue. Right panel is
the count of bcftools calls inside the locus, a sequential magnitude.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

summary_path, out_path = sys.argv[1:3]

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
# diverging blue <-> red with a neutral gray midpoint; sequential blue for counts
DIVERGING = LinearSegmentedColormap.from_list(
    "loss_gain", ["#9e2020", "#d03b3b", "#dc7272", "#e8a9a9", "#f0efec",
                  "#9ec5f4", "#5598e7", "#2a78d6", "#184f95"],
)
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "count", ["#f5f5f2", "#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"],
)

df = pd.read_csv(summary_path, sep="\t")
samples = sorted(df["sample"].unique())
group_order = list(dict.fromkeys(df["group"]))

# y positions, grouped, with a gap between groups
y_of, labels, y, bands = {}, [], 0, []
for group in group_order:
    genes = sorted(df.loc[df["group"] == group, "gene"].unique())
    start = y
    for gene in genes:
        y_of[gene] = y
        labels.append(gene)
        y += 1
    bands.append((group, start, y - 1))
    y += 1
n_rows = y - 1

depth = df.pivot(index="gene", columns="sample", values="depth_ratio")
counts = df.pivot(index="gene", columns="sample", values="n_variants")

fig_h = 1.5 + 0.21 * n_rows
fig, axes = plt.subplots(
    1, 2, figsize=(8.6, fig_h), facecolor=SURFACE, sharey=True,
    gridspec_kw={"wspace": 0.12},
)

panels = [
    (axes[0], depth, DIVERGING, Normalize(0, 2), "Depth relative to genome mean", "{:.2f}"),
    (axes[1], counts, SEQUENTIAL, Normalize(0, max(1, counts.to_numpy().max())),
     "Variants vs K-12", "{:.0f}"),
]

for ax, table, cmap, norm, title, fmt in panels:
    ax.set_facecolor(SURFACE)
    for gene, yy in y_of.items():
        for xi, sample in enumerate(samples):
            value = table.at[gene, sample]
            colour = cmap(norm(value))
            ax.add_patch(
                plt.Rectangle(
                    (xi - 0.46, yy - 0.42), 0.92, 0.84,
                    facecolor=colour, edgecolor=SURFACE, linewidth=1.2,
                )
            )
            # white ink only where the cell is dark enough to need it
            lum = 0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2]
            ax.annotate(
                fmt.format(value), xy=(xi, yy), ha="center", va="center",
                fontsize=6, color=SURFACE if lum < 0.5 else TEXT_2,
            )
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, fontsize=7.5, color=TEXT, rotation=30, ha="left")
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=3)
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_ylim(n_rows + 0.7, -0.7)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.annotate(
        title, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 34),
        textcoords="offset points", ha="center", fontsize=8.5, color=TEXT,
    )

axes[0].set_yticks(list(y_of.values()))
axes[0].set_yticklabels(labels, fontsize=7, color=TEXT)

for group, start, end in bands:
    axes[0].annotate(
        group, xy=(0, (start + end) / 2), xycoords=("figure fraction", "data"),
        xytext=(8, 0), textcoords="offset points",
        fontsize=7.5, fontweight="bold", color=TEXT_2, va="center", ha="left",
    )
    if start > 0:
        for ax in axes:
            ax.axhline(start - 1, color=GRID, lw=0.9, zorder=0)

fig.suptitle(
    "K-12 resistance loci: coverage and variant burden", x=0.99, y=0.995,
    ha="right", fontsize=10, color=TEXT,
)
fig.subplots_adjust(left=0.30, right=0.97, top=1 - 1.15 / fig_h, bottom=0.2 / fig_h)
fig.savefig(out_path, dpi=200, facecolor=SURFACE)
