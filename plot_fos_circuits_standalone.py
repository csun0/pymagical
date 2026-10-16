import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MPath
import seaborn as sns

# --- Embedded plotting style (no local helper imports) -----------------------
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

colors_theme = {
    "dark_blue": "#071D49", "medium_blue": "#A6B5E0", "light_blue": "#EDF0FF",
    "red_dark": "#CF451C", "red_light": "#F7634F",
    "cobalt_dark": "#0066F5", "cobalt_light": "#00A1FF",
    "medium_gray": "#B9B4B4", "dark_gray": "#4B4C4E",
    "black": "#000000", "beige": "#F2D9D0",
}
STRUCT = colors_theme["dark_blue"]

# Use Matplotlib's default font family and weight; retain text sizes.
hf = SimpleNamespace()
for prefix in ("", "b", "i"):
    for size_code, points in {"l": 12, "m": 9, "s": 6, "xs": 3}.items():
        properties = FontProperties(size=points)
        if prefix == "b":
            properties.set_weight("bold")
        elif prefix == "i":
            properties.set_style("italic")
        setattr(hf, f"{prefix}{size_code}f", properties)

MM_TO_IN = 1 / 25.4

THEMES = {
    "light": dict(
        suffix="", bg=None, struct=colors_theme["dark_blue"],
        peak=colors_theme["red_dark"], gene=colors_theme["cobalt_dark"], other=colors_theme["medium_gray"],
        loop=colors_theme["dark_gray"], band=colors_theme["beige"], band_alpha=0.8,
        atac={"CTRL": colors_theme["medium_gray"], "CR": colors_theme["cobalt_dark"]},
    ),
    "dark": dict(
        suffix="_dark", bg=colors_theme["black"], struct=colors_theme["light_blue"],
        peak=colors_theme["red_light"], gene=colors_theme["cobalt_light"], other=colors_theme["medium_gray"],
        loop=colors_theme["medium_blue"], band=colors_theme["red_dark"], band_alpha=0.3,
        atac={"CTRL": colors_theme["medium_gray"], "CR": colors_theme["cobalt_light"]},
    ),
}
T = SimpleNamespace(**THEMES["light"])


def set_theme(name):
    T.__dict__.update(THEMES[name])


def style_axes(ax, lw=0.5):
    """plot_ex.style_axes, but colored by the active theme."""
    for spine in ax.spines.values():
        spine.set_linewidth(lw)
        spine.set_color(T.struct)
    ax.tick_params(width=lw, length=2, colors=T.struct)
    for lbl in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        lbl.set_fontproperties(hf.mf)


def blank(ax):
    ax.set_yticks([])
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def save(fig, path):
    """Write `path` (pdfs/<kind>/<name>.pdf) plus a 600-dpi PNG copy at pngs/<kind>/<name>.png.

    Light: transparent background. Dark: filled with the theme background.
    """
    pdf_root = next(p for p in path.parents if p.name == "pdfs")
    png = pdf_root.parent / "pngs" / path.relative_to(pdf_root).with_suffix(".png")
    for d in (path.parent, png.parent):
        d.mkdir(parents=True, exist_ok=True)
    bg = dict(transparent=True) if T.bg is None else dict(transparent=True, facecolor=T.bg)
    for out in (path, png):
        fig.savefig(out, dpi=600, **bg)
        print(f"Wrote {out}")

# --- Parameters -------------------------------------------------------------
INPUTS = Path("/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/inputs")
CELLTYPE_DIR = INPUTS / "excitatory_neurons"
REFSEQ = INPUTS / "rhemac10_refseq.txt"
FRAG_DIR = Path("/mnt/home/agebrain/ceph/anderson/snRNA_ATAC_BA46/raw")
GTF = "/mnt/home/agebrain/ceph/resources/genomes/rheMac10/ucsc/rheMac10.ensGene.gtf"
OUT_DIR = Path(__file__).resolve().parent
CACHE_DIR = OUT_DIR / "cache"
LIB_CACHE = CACHE_DIR / "exc_lib.tsv"
PDF_DIR = OUT_DIR / "pdfs"

CIRCUITS = {
    # gene: chrom, TSS, peak, peak-gene prob, peak TFs (ap1_circuits.tsv, by prob), zoom inset
    "TBK1":    dict(chrom="chr11", tss=64001781, peak=(64010485, 64010986), prob=0.866, zoom=False,
                    tfs=["ZNF169", "ZNF557", "FOSL2", "ARID3A", "FOS"]),
    "PRKCA":   dict(chrom="chr16", tss=62679586, peak=(62808762, 62809263), prob=0.8738,
                    tfs=["ETS2", "NFE2L1", "ZNF22", "FOS"]),
    "LRRTM4":  dict(chrom="chr13", tss=31165391, peak=(31568555, 31569056), prob=0.8908,
                    tfs=["ZNF77", "FOS", "FOSL2"]),
    "SLC17A5": dict(chrom="chr4",  tss=99282133, peak=(99356278, 99356779), prob=0.814,
                    tfs=["ZNF613", "PRDM6", "ZNF557", "ZNF384", "ZNF182", "FOSL2"]),
}
CONDITIONS = ["CTRL", "CR"]
BIN_BP = 25
PAD_FRAC = 0.12          # wide window = peak..TSS span padded by this fraction per side
MIN_PAD = 2000

FIG_W_MM = 183
FONT, IFONT = hf.mf, hf.mf  # one text size everywhere (9 pt)
ATAC_ALPHA = 0.65
ZOOM_FLANK_BP = 1000     # inset shows peak +/- this
ZOOM_W, ZOOM_H = 0.2, 0.8  # inset size, as a fraction of the ATAC axes
ZOOM_Y = 0.2             # inset bottom, as a fraction of the ATAC axes


# --- Data -------------------------------------------------------------------
def read_meta():
    meta = pd.read_csv(CELLTYPE_DIR / "atac_meta.txt", sep="\t", header=None,
                       names=["idx", "barcode", "celltype", "sample", "condition"])
    meta["bc"] = meta.barcode.str.split(":", n=1).str[1]
    return meta


def load_lib():
    """Total in-peak counts per condition, plus pseudobulk CPM of each circuit peak."""
    if LIB_CACHE.exists():
        print(f"Loaded cached library sizes from {LIB_CACHE}")
        return pd.read_csv(LIB_CACHE, sep="\t")

    import pyarrow.csv as pv

    peaks = pd.read_csv(CELLTYPE_DIR / "atac_peaks.txt", sep="\t", header=None,
                        names=["idx", "chr", "start", "end"])
    key = peaks.set_index(["chr", "start"]).idx
    peak_idx = {g: key.loc[(c["chrom"], c["peak"][0])] for g, c in CIRCUITS.items()}
    meta = read_meta()
    n_cells = meta.idx.max() + 1
    cell_tot = np.zeros(n_cells)
    peak_cell = {g: np.zeros(n_cells) for g in CIRCUITS}

    print("Streaming atac_counts.txt ...")
    reader = pv.open_csv(CELLTYPE_DIR / "atac_counts.txt",
                         read_options=pv.ReadOptions(column_names=["peak", "cell", "count"],
                                                     block_size=1 << 28),
                         parse_options=pv.ParseOptions(delimiter="\t"))
    for i, batch in enumerate(reader):
        p = batch.column("peak").to_numpy()
        c = batch.column("cell").to_numpy()
        n = batch.column("count").to_numpy().astype(float)
        cell_tot += np.bincount(c, weights=n, minlength=n_cells)
        for g, pi in peak_idx.items():
            m = p == pi
            peak_cell[g] += np.bincount(c[m], weights=n[m], minlength=n_cells)
        if i % 20 == 0:
            print(f"  batch {i}", flush=True)

    cond = meta.set_index("idx").condition
    rows = []
    for cnd in CONDITIONS:
        cells = cond.index[cond == cnd].values
        lib = cell_tot[cells].sum()
        row = {"condition": cnd, "lib": lib}
        row.update({g: peak_cell[g][cells].sum() / lib * 1e6 for g in CIRCUITS})
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(LIB_CACHE, sep="\t", index=False)
    print(f"Wrote {LIB_CACHE}")
    return out


def load_coverage(gene, chrom, window, lib):
    """Binned fragment-coverage CPM per condition across `window`."""
    cache = CACHE_DIR / f"{gene}_coverage.tsv"
    if cache.exists():
        cov = pd.read_csv(cache, sep="\t")
        if cov.bin_start.iloc[0] == window[0] and cov.bin_start.iloc[1] - cov.bin_start.iloc[0] == BIN_BP:
            print(f"Loaded cached coverage from {cache}")
            return cov
        print(f"Cached coverage for {gene} does not match window; recomputing")

    import pysam

    meta = read_meta()
    width = window[1] - window[0]
    diff = {c: np.zeros(width + 1) for c in CONDITIONS}
    n_frag = {c: 0 for c in CONDITIONS}
    for sample, cells in meta.groupby("sample"):
        cond_of_bc = dict(zip(cells.bc, cells.condition))
        tbx = pysam.TabixFile(str(FRAG_DIR / sample / "atac_fragments.tsv.gz"))
        for line in tbx.fetch(chrom, window[0], window[1]):
            _, start, end, bc, _ = line.split("\t")
            cond = cond_of_bc.get(bc)
            if cond is None:
                continue
            s0, e0 = max(int(start), window[0]) - window[0], min(int(end), window[1]) - window[0]
            diff[cond][s0] += 1
            diff[cond][e0] -= 1
            n_frag[cond] += 1
    print(f"{gene}: excitatory-neuron fragments in window:", n_frag)

    n_bins = width // BIN_BP
    cov = pd.DataFrame({"bin_start": window[0] + np.arange(n_bins) * BIN_BP})
    for cond in CONDITIONS:
        per_bp = np.cumsum(diff[cond])[:n_bins * BIN_BP]
        cov[cond] = per_bp.reshape(n_bins, BIN_BP).mean(axis=1) / lib[cond] * 1e6
    cov.to_csv(cache, sep="\t", index=False)
    print(f"Wrote {cache}")
    return cov


def load_genes(chrom, window):
    ref = pd.read_csv(REFSEQ, sep="\t", header=None, names=["chr", "strand", "start", "end", "name"])
    ref = ref[(ref.chr == chrom) & (ref.start < window[1]) & (ref.end > window[0])]
    ref = ref.assign(length=ref.end - ref.start).sort_values("length", ascending=False)
    return ref.drop_duplicates("name").sort_values("start")


def load_exons(gene, chrom):
    """Exons of the Ensembl transcript matching `gene`'s RefSeq span (see alt script)."""
    cache = CACHE_DIR / f"{gene['name']}_exons.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    rows = []
    with open(GTF) as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.split("\t")
            if f[2] != "exon" or f[6] != gene.strand:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            if start >= gene.end or end <= gene.start:
                continue
            attrs = dict(kv.strip().split(" ", 1) for kv in f[8].strip().rstrip(";").split(";") if kv.strip())
            rows.append((attrs["gene_id"].strip('"'), attrs["transcript_id"].strip('"'), start, end))
    ex = pd.DataFrame(rows, columns=["gene_id", "transcript_id", "start", "end"])
    tss = gene.start if gene.strand == "+" else gene.end
    first = ex.groupby("gene_id").start.min() if gene.strand == "+" else ex.groupby("gene_id").end.max()
    gene_id = (first - tss).abs().idxmin()
    ex = ex[ex.gene_id == gene_id]
    tx = ex.transcript_id.value_counts().idxmax()
    ex = ex[ex.transcript_id == tx].sort_values("start").reset_index(drop=True)
    print(f"{gene['name']} -> {gene_id} / {tx} ({len(ex)} exons)")
    ex.to_csv(cache, sep="\t", index=False)
    return ex


# --- Drawing ----------------------------------------------------------------
def draw_gene(ax, g, exons, y, color, window, label=True):
    span = window[1] - window[0]
    x0, x1 = max(g.start, window[0]), min(g.end, window[1])
    ax.plot([x0, x1], [y, y], color=color, lw=1.0, solid_capstyle="butt", zorder=2)
    marker = ">" if g.strand == "+" else "<"
    xs = np.arange(x0 + span * 0.03, x1 - span * 0.01, span * 0.04)
    ax.plot(xs, np.full_like(xs, y, dtype=float), linestyle="none", marker=marker,
            ms=2.5, color=color, mew=0, zorder=2)
    if exons is not None:
        for e in exons.itertuples():
            if e.end > window[0] and e.start < window[1]:
                # keep exons visible at wide scales
                w = max(e.end - e.start, span * 0.0015)
                ax.add_patch(Rectangle((e.start, y - 0.22), w, 0.44, color=color, lw=0, zorder=3))
    if label:
        ax.text((x0 + x1) / 2, y - 0.35, g["name"], ha="center", va="top",
                fontproperties=IFONT, color=color)


def coord_axis(ax, window, y, chrom):
    ax.set_xlim(*window)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_position(("axes", y))
    ax.xaxis.set_ticks_position("bottom")
    ax.tick_params(axis="x", labeltop=False, labelbottom=True)
    ax.set_yticks([])
    set_kb_ticks(ax, window)
    style_axes(ax)
    for lbl in ax.get_xticklabels():
        lbl.set_fontproperties(FONT)
    ax.annotate(f"CHR {chrom[3:]}", xy=(0, y), xycoords="axes fraction", xytext=(0, 2),
                textcoords="offset points", ha="left", va="bottom", fontproperties=FONT, color=T.struct)


def set_kb_ticks(ax, window):
    span = window[1] - window[0]
    step = next(s for s in (1000, 2000, 5000, 10000, 20000, 50000, 100000) if span / s <= 6)
    ticks = np.arange(np.ceil(window[0] / step) * step, window[1], step)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / 1e3:,.0f} kb" for t in ticks])


def draw_link(ax, gene, tss_x, peak_mid, window):
    blank(ax)
    ax.set_ylim(0, 0.5)
    box_y, box_h = 0.14, 0.1
    box_w = (window[1] - window[0]) * 0.006
    ax.add_patch(Rectangle((tss_x - box_w / 2, box_y), box_w, box_h, color=T.gene, lw=0, zorder=3))
    top = box_y + box_h
    arc = MPath([(tss_x, top), (tss_x, top + 0.28), (peak_mid, top + 0.28), (peak_mid, top)],
                [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
    ax.add_patch(PathPatch(arc, fill=False, color=T.loop, lw=0.6, zorder=2))
    ax.text(tss_x, box_y - 0.07, gene, ha="center", va="top", fontproperties=IFONT, color=T.gene)
    return box_y, box_h


def peak_square(ax, peak_mid, box_y, box_h, window):
    """Square peak marker, sized after the layout is frozen."""
    px = ax.get_window_extent()
    sq_w = box_h * px.height / px.width * (window[1] - window[0])
    ax.add_patch(Rectangle((peak_mid - sq_w / 2, box_y), sq_w, box_h, color=T.peak, lw=0, zorder=3))


def peak_label(ax, chrom, peak, tfs, x0, x1, box_y, box_h, window):
    """Coordinates + TFs beside the peak band (figure-x edges x0/x1), on the roomier side.

    On the left (TSS) side the arc rises over the label, so it hangs down from the
    square's top edge instead of being centred on it, with all TFs on one line.
    """
    right = (sum(peak) / 2 - window[0]) / (window[1] - window[0]) < 0.6
    per_line = 4 if right else len(tfs)
    lines = [f"CHR {chrom[3:]}: {peak[0]:,}–{peak[1]:,}"]
    lines += [("TFs: " if i == 0 else "") + ", ".join(tfs[i:i + per_line])
              for i in range(0, len(tfs), per_line)]
    y = box_y + box_h / 2 if right else box_y + box_h
    ax.annotate("\n".join(lines), xy=(x1 if right else x0, y), xycoords=("figure fraction", "data"),
                xytext=(3 if right else -3, 0), textcoords="offset points",
                ha="left" if right else "right", va="center" if right else "top",
                ma="left" if right else "right", fontproperties=FONT, color=T.peak)


def draw_atac(ax, cov):
    top_val = max(cov[c].max() for c in CONDITIONS)
    ytick = float(f"{top_val:.1g}")
    edges = np.append(cov.bin_start.values, cov.bin_start.iloc[-1] + BIN_BP)
    # CTRL behind, CR on top
    for z, cond in enumerate(CONDITIONS):
        ax.stairs(cov[cond].values, edges, fill=True, color=T.atac[cond],
                  alpha=ATAC_ALPHA, lw=0, zorder=2 + z, label=cond)
    ax.set_ylim(0, top_val * 1.1)
    ax.set_yticks([0, ytick])
    ax.set_yticklabels(["0", f"{ytick:g}"])
    ax.set_ylabel("CPM", fontproperties=FONT, color=T.struct)
    style_axes(ax)
    for lbl in ax.get_yticklabels() + ax.get_xticklabels():
        lbl.set_fontproperties(FONT)
    leg = ax.legend(loc="upper left", frameon=False, prop=FONT, handlelength=1,
                    handleheight=0.8, borderaxespad=0)
    for txt in leg.get_texts():
        txt.set_color(T.struct)
    sns.despine(ax=ax)


def draw_zoom(ax, cov, peak, tss_x, window):
    """Inset of the CPM around the peak, joined to a dotted box on the main track."""
    zoom = (peak[0] - ZOOM_FLANK_BP, peak[1] + ZOOM_FLANK_BP)
    sub = cov[(cov.bin_start >= zoom[0]) & (cov.bin_start < zoom[1])]
    top_val = max(sub[c].max() for c in CONDITIONS) * 1.15

    # sit in the gap between TSS and peak, on the peak's side
    span = window[1] - window[0]
    frac = lambda x: (x - window[0]) / span  # noqa: E731
    w, h = ZOOM_W, ZOOM_H
    gap = 0.06
    if frac(tss_x) < frac(peak[0]):
        x = min(frac(zoom[0]) - gap - w, 1 - w)
    else:
        x = max(frac(zoom[1]) + gap, 0)
    axins = ax.inset_axes([x, ZOOM_Y, w, h])
    axins.patch.set_visible(False)

    edges = np.append(sub.bin_start.values, sub.bin_start.iloc[-1] + BIN_BP)
    axins.axvspan(*peak, color=T.band, alpha=T.band_alpha, lw=0, zorder=0)
    # condition with the higher max CPM behind, so the lower one stays visible
    for z, cond in enumerate(sorted(CONDITIONS, key=lambda c: -sub[c].max())):
        axins.stairs(sub[cond].values, edges, fill=True, color=T.atac[cond],
                     alpha=ATAC_ALPHA, lw=0, zorder=2 + z)
    axins.set_xlim(*zoom)
    axins.set_ylim(0, top_val)
    axins.set_xticks([])
    axins.set_yticks([])
    for spine in axins.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color(T.struct)

    ind = ax.indicate_inset_zoom(axins, edgecolor=T.struct, lw=0.5, alpha=1)
    ind.rectangle.set_linestyle((0, (1, 1.5)))
    # connectors: lower-left, upper-left, lower-right, upper-right. matplotlib picks the
    # visible pair from pre-layout positions, which is unstable; fix it by inset side.
    inset_left = frac(tss_x) < frac(peak[0])
    keep = (0, 3) if inset_left else (1, 2)
    for i, con in enumerate(ind.connectors):
        con.set_linestyle((0, (1, 1.5)))
        con.set_linewidth(0.5)
        con.set_visible(i in keep)
    return axins


def draw_genes(ax, genes, name, exons, window):
    blank(ax)
    ax.set_ylim(-1, 1)
    span = window[1] - window[0]
    for _, g in genes.iterrows():
        is_target = g["name"] == name
        # at wide scales, only label neighbours long enough to hold their name
        shown = min(g.end, window[1]) - max(g.start, window[0])
        draw_gene(ax, g, exons if is_target else None, 0.1, T.gene if is_target else T.other,
                  window, label=is_target or shown > span * 0.08)


def fig_x(fig, ax, x):
    return fig.transFigure.inverted().transform(ax.transData.transform((x, 0)))[0]


def plot_circuit(gene, c, window, cov, genes, exons):
    tss_x, peak = c["tss"], c["peak"]
    peak_mid = sum(peak) / 2
    fig, axes = plt.subplots(
        4, 1, sharex=True, dpi=600, layout="constrained",
        figsize=(FIG_W_MM * MM_TO_IN, FIG_W_MM * 1.5 / 4 * MM_TO_IN),
        gridspec_kw={"height_ratios": [1, 1, 2, 1]},
    )
    ax_coord, ax_loop, ax_atac, ax_gene = axes
    for ax in axes:
        ax.patch.set_visible(False)
    COORD_Y = 0.8
    coord_axis(ax_coord, window, COORD_Y, c["chrom"])
    box_y, box_h = draw_link(ax_loop, gene, tss_x, peak_mid, window)
    draw_atac(ax_atac, cov)
    if c.get("zoom", True):
        draw_zoom(ax_atac, cov, peak, tss_x, window)
    ax_atac.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_atac.spines["bottom"].set_visible(False)
    draw_genes(ax_gene, genes, gene, exons, window)

    fig.draw_without_rendering()
    fig.set_layout_engine("none")
    x0, x1 = fig_x(fig, ax_coord, peak[0]), fig_x(fig, ax_coord, peak[1])
    # the peak is sub-mm for long links; give the band a minimum visible width
    min_w = 0.6 / FIG_W_MM
    if x1 - x0 < min_w:
        mid = (x0 + x1) / 2
        x0, x1 = mid - min_w / 2, mid + min_w / 2
    c_pos, g_pos = ax_coord.get_position(), ax_gene.get_position()
    y_top = c_pos.y0 + COORD_Y * c_pos.height
    peak_square(ax_loop, peak_mid, box_y, box_h, window)
    peak_label(ax_loop, c["chrom"], peak, c["tfs"], x0, x1, box_y, box_h, window)
    fig.patches.append(Rectangle((x0, g_pos.y0), x1 - x0, y_top - g_pos.y0, transform=fig.transFigure,
                                 figure=fig, color=T.band, alpha=T.band_alpha, lw=0, zorder=-1))
    return fig


def main(selected):
    CACHE_DIR.mkdir(exist_ok=True)
    libdf = load_lib()
    print(libdf.to_string(index=False))
    lib = libdf.set_index("condition").lib.to_dict()

    PDF_DIR.mkdir(exist_ok=True)
    for gene in selected:
        c = CIRCUITS[gene]
        lo, hi = sorted([c["tss"], *c["peak"]])[0], sorted([c["tss"], *c["peak"]])[-1]
        pad = max(MIN_PAD, (hi - lo) * PAD_FRAC)
        window = (int(lo - pad) // 1000 * 1000, int(np.ceil((hi + pad) / 1000) * 1000))
        genes = load_genes(c["chrom"], window)
        target = genes.set_index("name").loc[gene].copy()
        target["name"] = gene
        ref_tss = target.start if target.strand == "+" else target.end
        assert ref_tss == c["tss"], f"{gene}: RefSeq TSS {ref_tss} != circuit TSS {c['tss']}"
        exons = load_exons(target, c["chrom"])

        cov = load_coverage(gene, c["chrom"], window, lib)
        print(f"{gene}: window {c['chrom']}:{window[0]:,}-{window[1]:,}")
        for theme in ("light", "dark"):
            set_theme(theme)
            fig = plot_circuit(gene, c, window, cov, genes, exons)
            save(fig, PDF_DIR / "circuit" / f"{gene.lower()}_circuit{T.suffix}.pdf")
            plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1:] or list(CIRCUITS))
