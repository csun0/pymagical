"""Plot unique features and circuits in the reported MAGICAL circuit tables."""
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
import pandas as pd
import seaborn as sns

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

sys.path.insert(0, "/mnt/home/csun1/palette")
from plot_ex import hf, style_axes

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
MM_TO_IN = 1 / 25.4

husl_12 = sns.color_palette("husl", 12)
color_dict = {
    "Excitatory Neurons": husl_12[-4],
    "Inhibitory Neurons": husl_12[0],
    "Astrocytes": husl_12[4],
    "Microglia": husl_12[1],
    "Oligodendrocytes": husl_12[-3],
    "OPCs": husl_12[-2],
    "Endothelial": "goldenrod",
    "Unassigned": "gray",
}
names = {
    "excitatory_neurons": "Excitatory Neurons",
    "inhibitory_neurons": "Inhibitory Neurons",
    "astrocytes": "Astrocytes",
    "microglia": "Microglia",
    "oligodendrocytes": "Oligodendrocytes",
    "opcs": "OPCs",
}
labels = list(names.values())
colors = [color_dict[label] for label in labels]

# Count distinct reported pairs and triads directly from the circuit tables.
tf_pattern = re.compile(r"([A-Za-z0-9_.-]+) \(\d+\.\d+, [+-] \[")
circuit_counts = []
for cell_type in names:
    table = pd.read_csv(OUT / cell_type / f"{cell_type}_py_2000.txt", sep="\t")
    pairs, triads = set(), set()
    for row in table.itertuples(index=False, name=None):
        pair = (row[0], row[3], int(row[4]), int(row[5]))
        pairs.add(pair)
        tfs = tf_pattern.findall(row[7])
        if not tfs:
            raise ValueError(f"Unable to parse TFs for {cell_type}: {row[7]}")
        triads.update((*pair, tf) for tf in tfs)
    circuit_counts.append({
        "cell_type": cell_type,
        "reported_rows": len(table),
        "unique_gene_peak_pairs": len(pairs),
        "unique_gene_peak_tf_triads": len(triads),
        "unique_peaks": len({pair[1:] for pair in pairs}),
        "unique_tfs": len({triad[-1] for triad in triads}),
        "unique_target_genes": len({pair[0] for pair in pairs}),
    })
counts = pd.DataFrame(circuit_counts).set_index("cell_type")
if not (
    (counts["unique_gene_peak_pairs"] >= counts["unique_peaks"]).all()
    and (counts["unique_gene_peak_pairs"] >= counts["unique_target_genes"]).all()
    and (counts["unique_gene_peak_tf_triads"] >= counts["unique_gene_peak_pairs"]).all()
):
    raise ValueError("Reported circuit counts are inconsistent")
counts.to_csv(OUT / "summary_circuit_counts.csv")
df = counts

axis_font = hf.blf.copy()
axis_font.set_size(16)
value_font = hf.blf.copy()
value_font.set_size(12)

fig, ax = plt.subplots(
    figsize=(183 * MM_TO_IN, 137.25 * MM_TO_IN),
    dpi=600, layout="constrained",
)
# Keep the requested black background as an editable vector rectangle even
# when exporting with transparent=True.
fig.add_artist(Rectangle(
    (0, 0), 1, 1, transform=fig.transFigure, facecolor="black",
    edgecolor="none", zorder=-100, clip_on=False,
))
fig.patch.set_facecolor("black")

metrics = [
    ("unique_peaks", "Peaks", "//", -0.22),
    ("unique_tfs", "TFs", "", 0),
    ("unique_target_genes", "Target genes", "..", 0.22),
]
feature_max = df[[metric[0] for metric in metrics]].to_numpy().max()
plt.rcParams["hatch.linewidth"] = 1.0
ax.set_facecolor("black")
for column, label, hatch, offset in metrics:
    values = df[column].to_numpy()
    positions = [x + offset for x in range(len(df))]
    ax.bar(positions, values, width=0.16, color=colors,
           edgecolor="white", linewidth=1.0, hatch=hatch)
    for x, value in enumerate(values):
        ax.text(x + offset, value + feature_max * 0.016, f"{value:,}", rotation=90,
                ha="center", va="bottom", color="white", fontproperties=value_font)

ax.set_xticks(range(len(df)), labels)
ax.set_xlim(-0.6, len(df) - 0.4)
ax.set_ylim(0, feature_max * 1.45)
ax.set_ylabel("# in output circuits", fontproperties=axis_font, color="white", labelpad=10)
ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
style_axes(ax, lw=1.5)
ax.tick_params(colors="white", width=1.5, length=5)
for spine in ax.spines.values():
    spine.set_color("white")
for tick in ax.get_xticklabels():
    tick.set_fontproperties(axis_font)
    tick.set_color("white")
    tick.set_rotation(45)
    tick.set_rotation_mode("anchor")
    tick.set_ha("right")
for tick in ax.get_yticklabels():
    tick.set_fontproperties(axis_font)
    tick.set_color("white")
handles = [Patch(facecolor="#666666", edgecolor="white", linewidth=1.0,
                 hatch=hatch, label=label) for _, label, hatch, _ in metrics]
ax.legend(handles=handles, loc="upper right", ncol=3, frameon=False,
          labelcolor="white", prop=axis_font, handlelength=1.1,
          columnspacing=0.7, handletextpad=0.35)

sns.despine(fig=fig, top=True, right=True)
for stem, include_triads in [("summary_counts", False), ("summary_counts_triads", True)]:
    annotations = []
    for x, cell_type in enumerate(df.index):
        pair_count = counts.loc[cell_type, "unique_gene_peak_pairs"]
        triad_count = counts.loc[cell_type, "unique_gene_peak_tf_triads"]
        text = f"{pair_count:,} pairs"
        if include_triads:
            text += f"\n{triad_count:,} triads"
        tallest = df.loc[cell_type, [metric[0] for metric in metrics]].max()
        annotations.append(ax.text(
            x, tallest + feature_max * 0.19, text, ha="center", va="bottom",
            color="white", fontproperties=value_font, linespacing=1.3,
        ))
    for extension in ("pdf", "svg", "png"):
        path = OUT / f"{stem}.{extension}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0, transparent=True)
        print(path)
    for annotation in annotations:
        annotation.remove()
plt.close(fig)
