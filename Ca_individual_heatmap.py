from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages

# input settings
INPUT_CSV = "centered_matrix.csv"
OUTPUT_PNG = "calcium_individual_heatmap.png"
OUTPUT_PDF = "calcium_individual_heatmap.pdf"
COLOR_LIMIT = 2.0  # Lower = stronger red/blue colors; higher = softer colors.
HEATMAP_COLORS = ["#2166AC", "#F7F7F7", "#B2182B"]
# ==================================

SAMPLE_ORDER = [
    "N1", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
    "O1", "O2", "O3", "O5", "O6", "O7", "O8", "O9", "O10",
    "D1", "D2", "D3", "D5", "D6", "D7", "D8S",
    "HF1", "HF2", "HF3", "HF4", "HF5", "HF6", "HF7", "HF8", "HF9", "HF10",
]

COHORTS = [
    ("Normal", 0, 8),
    ("Obese", 9, 17),
    ("Diabetic", 18, 24),
    ("HFrEF", 25, 34),
]


def main():
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(INPUT_CSV)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Cannot find {input_path}. Put the CSV beside this script or pass a file path."
        )

    df = pd.read_csv(input_path)

    for col in ["gene_symbol", "major_class"]:
        if col not in df.columns:
            raise ValueError(f"Missing required annotation column: {col}")

    missing_samples = [sample for sample in SAMPLE_ORDER if sample not in df.columns]
    if missing_samples:
        raise ValueError(
            "Missing required sample columns: " + ", ".join(missing_samples)
        )

    # ----------------------------------------------------------
    # KEEP ONLY MAGNESIUM-ANNOTATED GENES
    # Includes multi-class genes such as:
    # "Magnesium channels; Calcium channels"
    # ----------------------------------------------------------
    df = df[
        df["major_class"].astype(str).str.contains(
            "Calcium channels",
            case=False,
            na=False
        )
    ].copy()

    # Optional readable ordering.
    sort_columns = [
        col for col in ["subcategory", "gene_symbol"]
        if col in df.columns
    ]

    if sort_columns:
        df = df.sort_values(sort_columns, kind="stable").reset_index(drop=True)

    # Create numeric matrix from the 35 individual samples.
    matrix = df[SAMPLE_ORDER].apply(
        pd.to_numeric,
        errors="coerce"
    ).to_numpy(dtype=float)

    if not np.isfinite(matrix).all():
        raise ValueError(
            "All 35 sample columns must contain finite numeric heatmap values."
        )

    # Keep the supplied row order after magnesium filtering.
    row_labels = [
        f"{symbol}  [{str(category).replace(' channels', '').replace('; ', ' + ')}]"
        for symbol, category in zip(df["gene_symbol"], df["major_class"])
    ]

    n_genes = len(df)

    fig_height = max(10, 0.45 * n_genes + 6)
    fig = plt.figure(figsize=(20, fig_height))

    gs = fig.add_gridspec(
        3, 3,
        height_ratios=[0.55, 0.18, 14],
        width_ratios=[0.28, 1.0, 0.035],
        hspace=0.06,
        wspace=0.025,
    )

    fig.subplots_adjust(
        left=0.13,
        right=0.96,
        top=0.975,
        bottom=0.05
    )

    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")

    title_ax.text(
        0.5,
        0.70,
        "Person-by-Person Calcium Gene Expression Heatmap",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold"
    )

    title_ax.text(
        0.5,
        0.20,
        f"All retained calcium-annotated genes (n = {n_genes}) • "
        "each value is centered to the 9-sample Normal-group mean",
        ha="center",
        va="center",
        fontsize=9
    )

    labels_ax = fig.add_subplot(gs[2, 0])
    labels_ax.set_xlim(0, 1)
    labels_ax.set_ylim(n_genes - 0.5, -0.5)
    labels_ax.axis("off")

    font_size = max(6.0, min(10.0, 190 / max(n_genes, 1)))

    for row_index, label in enumerate(row_labels):
        labels_ax.text(
            0.99,
            row_index,
            label,
            ha="right",
            va="center",
            fontsize=font_size
        )

    strip_ax = fig.add_subplot(gs[1, 1])

    cohort_codes = np.array([[0] * 9 + [1] * 9 + [2] * 7 + [3] * 10])

    strip_ax.imshow(
        cohort_codes,
        aspect="auto",
        cmap=plt.get_cmap("tab10", 4),
        vmin=0,
        vmax=3
    )

    strip_ax.set_xlim(-0.5, len(SAMPLE_ORDER) - 0.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])

    for boundary in [9, 18, 25]:
        strip_ax.axvline(
            boundary - 0.5,
            color="black",
            linewidth=1.0
        )

    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    for cohort, start, end in COHORTS:
        strip_ax.text(
            (start + end) / 2,
            1.55,
            cohort,
            transform=strip_ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )

    heat_ax = fig.add_subplot(gs[2, 1])

    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        HEATMAP_COLORS,
        N=256
    )

    norm = TwoSlopeNorm(
        vmin=-COLOR_LIMIT,
        vcenter=0,
        vmax=COLOR_LIMIT
    )

    image = heat_ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm
    )

    heat_ax.set_xticks(np.arange(len(SAMPLE_ORDER)))
    heat_ax.set_xticklabels(
        SAMPLE_ORDER,
        rotation=90,
        fontsize=8
    )

    heat_ax.tick_params(axis="x", length=0, pad=3)
    heat_ax.set_yticks([])
    heat_ax.set_xlim(-0.5, len(SAMPLE_ORDER) - 0.5)
    heat_ax.set_ylim(n_genes - 0.5, -0.5)

    for boundary in [9, 18, 25]:
        heat_ax.axvline(
            boundary - 0.5,
            color="black",
            linewidth=1.0
        )

    cbar_ax = fig.add_subplot(gs[2, 2])

    cbar = fig.colorbar(
        image,
        cax=cbar_ax
    )

    cbar.set_label(
        "Centered expression\n(log2-CPM deviation from the Normal mean)",
        fontsize=9
    )

    cbar.ax.tick_params(labelsize=8)

    fig.text(
        0.61,
        0.012,
        "Blue = lower than the Normal-group mean; white = near the Normal-group mean; red = higher. "
        f"Color range: −{COLOR_LIMIT:g} to +{COLOR_LIMIT:g}; values beyond this range are color-saturated.",
        ha="center",
        va="bottom",
        fontsize=8
    )

    fig.savefig(
        OUTPUT_PNG,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.15
    )

    with PdfPages(OUTPUT_PDF) as pdf:
        pdf.savefig(
            fig,
            bbox_inches="tight",
            pad_inches=0.15
        )

    plt.close(fig)

if __name__ == "__main__":
    main()