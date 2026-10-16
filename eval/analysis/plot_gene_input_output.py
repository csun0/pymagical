"""Compare supplied candidate genes with unique genes in reported circuits.

Run from the repository with: uv run eval/analysis/plot_gene_input_output.py
The peak script reuses this plotting function to keep both figures identical.
"""
from pathlib import Path
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
}
names = {
    "excitatory_neurons": "Excitatory Neurons",
    "inhibitory_neurons": "Inhibitory Neurons",
    "astrocytes": "Astrocytes",
    "microglia": "Microglia",
    "oligodendrocytes": "Oligodendrocytes",
    "opcs": "OPCs",
}


def plot_input_output(feature):
    """Count raw candidate inputs and distinct reported output features."""
    if feature not in {"genes", "peaks"}:
        raise ValueError("feature must be 'genes' or 'peaks'")
    rows = []
    for cell_type in names:
        candidates = pd.read_csv(
            ROOT / "inputs" / cell_type / f"sig_cr_{feature}.txt",
            sep="\t", header=None,
        )
        circuits = pd.read_csv(OUT / cell_type / f"{cell_type}_py_2000.txt", sep="\t")
        if feature == "genes":
            supplied = set(candidates.iloc[:, 0])
            reported = set(circuits["Gene_symbol"])
            supplied_for_validation = supplied
        else:
            supplied = set(candidates.itertuples(index=False, name=None))
            reported = set(circuits[["Peak_chr", "Peak_start", "Peak_end"]]
                           .itertuples(index=False, name=None))
            # Match the loader's chromosome encoding for validation only:
            # it preserves chr1..chr22 and writes other chromosomes as chr0.
            autosomes = {f"chr{i}" for i in range(1, 23)}
            supplied_for_validation = {
                (chrom if chrom in autosomes else "chr0", start, end)
                for chrom, start, end in supplied
            }
        if not reported.issubset(supplied_for_validation):
            raise ValueError(f"Reported {feature} absent from inputs for {cell_type}")
        rows.append({
            "cell_type": cell_type,
            f"input_{feature}": len(supplied),
            f"output_{feature}": len(reported),
        })
    counts = pd.DataFrame(rows).set_index("cell_type")
    stem = f"summary_{feature}_input_output"
    counts.to_csv(OUT / f"{stem}.csv")
    print(counts.to_string())

    axis_font = hf.blf.copy()
    axis_font.set_size(16)
    value_font = hf.blf.copy()
    value_font.set_size(12)
    fig, ax = plt.subplots(
        figsize=(183 * MM_TO_IN, 137.25 * MM_TO_IN),
        dpi=600, layout="constrained",
    )
    fig.add_artist(Rectangle(
        (0, 0), 1, 1, transform=fig.transFigure, facecolor="black",
        edgecolor="none", zorder=-100, clip_on=False,
    ))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    plt.rcParams["hatch.linewidth"] = 1.0
    colors = [color_dict[label] for label in names.values()]
    maximum = counts.to_numpy().max()
    metrics = [(f"input_{feature}", "Input", "//", -0.13),
               (f"output_{feature}", "Output", "", 0.13)]
    for column, label, hatch, offset in metrics:
        values = counts[column].to_numpy()
        ax.bar([x + offset for x in range(len(counts))], values,
               width=0.16, color=colors, edgecolor="white", linewidth=1.0,
               hatch=hatch)
        for x, value in enumerate(values):
            ax.text(x + offset, value + maximum * 0.016, f"{value:,}",
                    rotation=90, ha="center", va="bottom", color="white",
                    fontproperties=value_font)
    ax.set_xticks(range(len(counts)), list(names.values()))
    ax.set_xlim(-0.6, len(counts) - 0.4)
    ax.set_ylim(0, maximum * 1.28)
    ax.set_ylabel(f"# of {feature}", fontproperties=axis_font,
                  color="white", labelpad=10)
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
    ax.legend(handles=handles, loc="upper right", ncol=2, frameon=False,
              labelcolor="white", prop=axis_font, handlelength=1.1,
              columnspacing=0.7, handletextpad=0.35)
    sns.despine(fig=fig, top=True, right=True)
    for extension in ("pdf", "svg", "png"):
        path = OUT / f"{stem}.{extension}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0, transparent=True)
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    plot_input_output("genes")
