"""Add vertical space between the title and the diagram in an nf-metro SVG.

nf-metro draws the title at a fixed offset with no option to change it, so this
grows the canvas upward by GAP pixels and moves only the title into the new
space. Then it rasterizes the result to PNG.

Usage: python add_title_space.py metro_map.svg metro_map.png [gap_px] [scale]
Needs resvg_py (pip install resvg_py), which nf-metro already installs.
"""
import re
import sys

import resvg_py

svg_path, png_path = sys.argv[1], sys.argv[2]
gap = float(sys.argv[3]) if len(sys.argv) > 3 else 50
scale = float(sys.argv[4]) if len(sys.argv) > 4 else 2

svg = open(svg_path).read()

m = re.search(r'width="([\d.]+)" height="([\d.]+)" viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
if not m:
    sys.exit("SVG root size/viewBox not found (already processed?)")
w, h = float(m.group(1)), float(m.group(2))
top_pad = 20  # extra margin above the title so it doesn't touch the image edge
grow = gap + top_pad
new_h = h + grow
svg = svg.replace(
    m.group(0),
    f'width="{w:g}" height="{new_h:g}" viewBox="0 {-grow:g} {w:g} {new_h:g}"',
    1,
)

# stretch the background so it also covers the new strip at the top
svg, n_bg = re.subn(
    r'(<rect x="0") y="0"( width="[\d.]+") height="[\d.]+"( [^>]*class="nf-metro-bg")',
    rf'\1 y="{-grow:g}"\2 height="{new_h:g}"\3',
    svg,
    count=1,
)

# move the title up into the new space
def shift_title(t):
    return re.sub(r'y="([\d.]+)"', lambda y: f'y="{float(y.group(1)) - gap:g}"', t.group(0), count=1)

svg, n_title = re.subn(r'<text[^>]*class="nf-metro-title"[^>]*>', shift_title, svg, count=1)
if not (n_bg and n_title):
    sys.exit(f"background found: {bool(n_bg)}, title found: {bool(n_title)}")

open(svg_path, "w").write(svg)
png = resvg_py.svg_to_bytes(svg_string=svg, zoom=int(scale))
open(png_path, "wb").write(bytes(png))
print(f"added {gap:g}px above the diagram -> {svg_path}, {png_path}")
