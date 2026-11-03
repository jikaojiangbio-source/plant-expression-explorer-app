"""Static report figures for the tomato jasmonate hairy-root case study.

Colour roles (validated with the dataviz palette checker, light surface):
categorical slots 1-3 for gene families (blue TPS, orange R2R3-MYB, aqua
JAZ), a blue/red diverging pair for down/up, and neutral greys for context.
Identity is never carried by colour alone: every highlighted family also has
a legend entry and the key genes carry direct text labels.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

import case_study_lib as lib  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
CONTEXT = "#c9c8c2"
TPS = "#2a78d6"
MYB = "#eb6834"
JAZ = "#1baf7a"
UP = "#e34948"
DOWN = "#2a78d6"
NEUTRAL_MID = "#f0efec"
MOCK = "#2a78d6"
JA = "#eb6834"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": INK_2,
        "axes.labelcolor": INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "svg.hashsalt": "case-study",
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=200, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def pca_figure(expression: pd.DataFrame, metadata: pd.DataFrame, path: Path, top_n: int = 500) -> None:
    matrix = expression.set_index("gene_id")
    variable = matrix.loc[matrix.var(axis=1).sort_values(ascending=False).index[:top_n]]
    centred = variable.sub(variable.mean(axis=1), axis=0).to_numpy().T
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    scores = u * s
    explained = s**2 / np.sum(s**2) * 100
    # Fix the sign so the jasmonate samples sit on the positive PC1 side.
    sign = 1 if scores[metadata["condition"].eq(lib.TEST_CONDITION).to_numpy(), 0].mean() > 0 else -1
    scores[:, 0] *= sign

    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for condition, colour in ((lib.REFERENCE_CONDITION, MOCK), (lib.TEST_CONDITION, JA)):
        mask = metadata["condition"].eq(condition).to_numpy()
        ax.scatter(scores[mask, 0], scores[mask, 1], s=64, color=colour, edgecolor=SURFACE,
                   linewidth=1.5, label=condition, zorder=3)
        for (x, y), name in zip(scores[mask, :2], metadata.loc[mask, "sample_id"], strict=True):
            ax.annotate(name, (x, y), xytext=(6, 4), textcoords="offset points", fontsize=7.5, color=INK_2)
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")
    ax.set_title(f"Samples separate by treatment (top {top_n} most variable genes)")
    ax.legend(loc="best", fontsize=8)
    _save(fig, path)


def volcano_figure(results: pd.DataFrame, families: pd.DataFrame, jaz_ids: set[str], path: Path) -> None:
    data = results.dropna(subset=["padj"]).copy()
    data["y"] = -np.log10(data["padj"].clip(lower=1e-300))
    cap = 80.0
    data["y_plot"] = data["y"].clip(upper=cap)
    tps_ids = set(families.loc[families["family"] == lib.FAMILY_TPS, "gene_id"])
    myb_ids = set(families.loc[families["family"] == lib.FAMILY_R2R3_MYB, "gene_id"])

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.scatter(data["log2FoldChange"], data["y_plot"], s=4, color=CONTEXT, linewidth=0,
               rasterized=True, label=f"all tested genes (n = {len(data):,})", zorder=1)
    layers = (
        (myb_ids, MYB, "R2R3-MYB (Pfam)"),
        (tps_ids, TPS, "terpene synthase (Pfam)"),
        (jaz_ids, JAZ, "JAZ (Chini et al. 2017)"),
    )
    for ids, colour, label in layers:
        subset = data.loc[data.index.isin(ids)]
        responsive = subset["regulation_lfc1"].isin([lib.REGULATION_UP, lib.REGULATION_DOWN])
        ax.scatter(subset.loc[responsive, "log2FoldChange"], subset.loc[responsive, "y_plot"], s=46,
                   color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=3, label=f"{label}, |log2FC| > 1 test significant")
        ax.scatter(subset.loc[~responsive, "log2FoldChange"], subset.loc[~responsive, "y_plot"], s=40,
                   facecolor=SURFACE, edgecolor=colour, linewidth=1.4, zorder=3, label=f"{label}, not significant")
    for x in (-1, 1):
        ax.axvline(x, color=INK_2, linewidth=0.8, linestyle=(0, (3, 3)), zorder=2)
    named = data.loc[data.index.isin(tps_ids) & data["atlas_gene_name"].notna() & data["regulation_lfc1"].eq(lib.REGULATION_UP)]
    for gene_id, row in named.iterrows():
        ax.annotate(str(row["atlas_gene_name"]), (row["log2FoldChange"], row["y_plot"]), xytext=(6, 2),
                    textcoords="offset points", fontsize=7.5, color=INK)
    ax.set_xlabel("log2 fold change, jasmonate vs mock (DESeq2 MLE)")
    ax.set_ylabel(f"-log10 adjusted p (standard Wald test; capped at {cap:.0f})")
    ax.set_title("Jasmonate response in tomato hairy roots, E-MTAB-11876")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7, ncol=2,
              handletextpad=0.3, columnspacing=1.2)
    _save(fig, path)


def family_bar_figure(rows: list[tuple[str, int, int, int]], path: Path) -> None:
    """``rows`` = (label, up, down, tested), drawn top to bottom in the given order."""

    rows = rows[::-1]

    fig, ax = plt.subplots(figsize=(6.0, 2.8))
    y = np.arange(len(rows))
    up_pct = np.array([r[1] / r[3] * 100 for r in rows])
    down_pct = np.array([r[2] / r[3] * 100 for r in rows])
    ax.barh(y, up_pct, height=0.55, color=UP, label="induced (up)")
    ax.barh(y, -down_pct, height=0.55, color=DOWN, label="repressed (down)")
    for yi, (name, up, down, tested), u, d in zip(y, rows, up_pct, down_pct, strict=True):
        ax.text(u + 1.2, yi, f"{up}/{tested}", va="center", fontsize=7.5, color=INK)
        ax.text(-d - 1.2, yi, f"{down}/{tested}", va="center", ha="right", fontsize=7.5, color=INK)
    ax.axvline(0, color=INK_2, linewidth=0.8)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlim(-30, 112)
    ax.set_xticks([-20, 0, 20, 40, 60, 80, 100], ["20%", "0", "20%", "40%", "60%", "80%", "100%"])
    ax.set_xlabel("share of expressed, tested family members (DESeq2 test of |log2FC| > 1, FDR < 0.05)")
    ax.set_title("TPS, R2R3-MYB and JAZ genes are over-represented among JA-induced genes")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=7.5)
    _save(fig, path)


def tps_heatmap_figure(results: pd.DataFrame, expression: pd.DataFrame, families: pd.DataFrame,
                       metadata: pd.DataFrame, path: Path) -> None:
    tps_ids = families.loc[families["family"] == lib.FAMILY_TPS, "gene_id"]
    matrix = expression.set_index("gene_id").reindex([g for g in tps_ids if g in set(expression["gene_id"])])
    matrix = matrix[list(metadata["sample_id"])]
    order = results.loc[matrix.index, "log2FoldChange"].sort_values(ascending=False).index
    matrix = matrix.loc[order]
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1).replace(0, np.nan), axis=0)

    cmap = LinearSegmentedColormap.from_list("blue_grey_red", [DOWN, NEUTRAL_MID, UP])
    fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(z) + 1.4))
    image = ax.imshow(z.to_numpy(), cmap=cmap, vmin=-2, vmax=2, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(z.shape[1]), z.columns, rotation=30, ha="right", fontsize=7.5)
    labels = []
    for gene_id in z.index:
        name = results.loc[gene_id, "atlas_gene_name"]
        status = results.loc[gene_id, "regulation_lfc1"]
        tag = " *" if status == lib.REGULATION_UP else ""
        labels.append(f"{gene_id}{' ' + str(name) if isinstance(name, str) else ''}{tag}")
    ax.set_yticks(range(len(z)), labels, fontsize=7.5)
    for i, gene_id in enumerate(z.index):
        ax.text(z.shape[1] - 0.4, i, f"{results.loc[gene_id, 'log2FoldChange']:+.1f}", va="center",
                fontsize=7.5, color=INK)
    ax.text(z.shape[1] - 0.4, -0.75, "log2FC", fontsize=7.5, color=INK_2, va="bottom")
    ax.set_xlim(-0.5, z.shape[1] + 0.6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colourbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    colourbar.set_label("row z-score of log2 normalised counts", fontsize=7.5)
    colourbar.ax.tick_params(labelsize=7)
    ax.set_title(f"Expressed terpene synthase genes ({len(z)} of {len(tps_ids)} in the proteome)\n"
                 "* = significant in the |log2FC| > 1 test", pad=18)
    _save(fig, path)


def atlas_agreement_figure(results: pd.DataFrame, atlas: pd.DataFrame, path: Path) -> None:
    lcol = next(c for c in atlas.columns if c.endswith(".log2foldchange"))
    theirs = atlas.set_index("Gene ID")[lcol]
    shared = results.index.intersection(theirs.index)
    x = results.loc[shared, "log2FoldChange"]
    y = theirs.loc[shared]
    mask = x.notna() & y.notna()
    r = float(np.corrcoef(x[mask], y[mask])[0, 1])

    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.scatter(x[mask], y[mask], s=3, color=TPS, alpha=0.35, linewidth=0, rasterized=True)
    lim = float(np.nanmax(np.abs(np.concatenate([x[mask], y[mask]])))) + 0.5
    ax.plot([-lim, lim], [-lim, lim], color=INK_2, linewidth=0.8, linestyle=(0, (3, 3)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("log2FC, this analysis (pyDESeq2)")
    ax.set_ylabel("log2FC, Expression Atlas (DESeq2 1.22, rounded to 0.1)")
    ax.set_title(f"Agreement with Expression Atlas\nPearson r = {r:.4f}, n = {int(mask.sum()):,} genes")
    _save(fig, path)


def make_all(results: pd.DataFrame, expression: pd.DataFrame, metadata: pd.DataFrame,
             families: pd.DataFrame, atlas: pd.DataFrame, out_dir: Path,
             jaz_ids: set[str] | None = None,
             family_rows: list[tuple[str, int, int, int]] | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pca_figure(expression, metadata, out_dir / "fig1_pca.png")
    volcano_figure(results, families, jaz_ids or set(), out_dir / "fig2_volcano.png")
    if family_rows:
        family_bar_figure(family_rows, out_dir / "fig3_family_response.png")
    tps_heatmap_figure(results, expression, families, metadata, out_dir / "fig4_tps_heatmap.png")
    atlas_agreement_figure(results, atlas, out_dir / "fig5_atlas_agreement.png")
