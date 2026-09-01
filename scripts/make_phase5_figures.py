from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

from make_deep_benchmark_figures import (
    BLACK,
    BLUE,
    GRAY,
    LIGHT_GRAY,
    ORANGE,
    PROGRAM_LABELS,
    PURPLE,
    RED,
    VERY_LIGHT,
    clean_axes,
    panel_label,
    save,
)


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

PROGRAM_LABEL = {
    "RAM2019_ENDO_2": "Endo2",
    "RAM2019_ENDO_6_SAENDO1": "SAEndo1",
}
PROGRAM_COLOR = {"Endo2": BLUE, "SAEndo1": PURPLE}
METHOD_LABEL = {"singscore": "Rank", "standardized_mean": "z-mean"}
METHOD_MARKER = {"singscore": "o", "standardized_mean": "s"}
DATASET_LABEL = {
    "GSE202379": "Clinical NASH cirrhosis",
    "GSE244832": "MASH F2–F4",
    "GSE256398_human": "MASH cirrhosis",
}


def p5(repo: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(repo / "results" / "phase5" / name)


def figure_9(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.2, 10.0), constrained_layout=True)
    grid = GridSpec(3, 2, figure=fig, height_ratios=[1.25, 1.0, 1.0])

    # A: exact abundance-versus-intensity decomposition.
    ax = fig.add_subplot(grid[0, 0])
    dec = p5(repo, "state_abundance_intensity_decomposition.csv")
    dec = dec[dec["program_id"].isin(PROGRAM_LABEL)].copy()
    dec["program"] = dec["program_id"].map(PROGRAM_LABEL)
    dec["method"] = dec["score_method"].map(METHOD_LABEL)
    order = []
    for program in ("Endo2", "SAEndo1"):
        for dataset in DATASET_LABEL:
            for method in ("Rank", "z-mean"):
                order.append((program, dataset, method))
    dec["order"] = dec.apply(lambda row: order.index((row["program"], row["dataset_id"], row["method"])), axis=1)
    dec = dec.sort_values("order")
    y = np.arange(len(dec))
    for index, row in enumerate(dec.itertuples()):
        positive = 0.0
        negative = 0.0
        for value, color in ((row.abundance_component, ORANGE), (row.intensity_component, BLUE)):
            left = positive if value >= 0 else negative
            ax.barh(index, value, left=left, height=0.68, color=color, edgecolor="white", linewidth=0.35)
            if value >= 0:
                positive += value
            else:
                negative += value
    labels = [
        f"{row.program} | {DATASET_LABEL[row.dataset_id]} | {row.method}"
        for row in dec.itertuples()
    ]
    ax.set_yticks(y, labels, fontsize=6.0)
    ax.invert_yaxis()
    ax.axvline(0, color=BLACK, lw=0.7)
    ax.set_xlabel("Exact Kitagawa component (score units)")
    ax.legend(
        handles=[
            Line2D([0], [0], color=ORANGE, lw=6, label="State abundance"),
            Line2D([0], [0], color=BLUE, lw=6, label="Within-state intensity"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        borderaxespad=0,
        fontsize=6.7,
    )
    clean_axes(ax)
    panel_label(ax, "A", x=-0.32)
    dec.to_csv(source / "figure_9a_state_decomposition.csv", index=False)

    # B: state-specific effects and frozen marker-matching boundary.
    ax = fig.add_subplot(grid[0, 1])
    intensity = p5(repo, "state_intensity_effects.csv")
    intensity = intensity[
        intensity["analysis_tier"].eq("primary")
        & intensity["program_id"].isin(PROGRAM_LABEL)
    ].copy()
    random = p5(repo, "primary_state_random_benchmark.csv")
    keys = ["dataset_id", "source_state", "program_id", "score_method", "cell_gate"]
    intensity = intensity.merge(
        random[keys + ["above_random_95th_percentile"]], on=keys, how="left"
    )
    rows = intensity[["dataset_id", "source_state"]].drop_duplicates().sort_values(["dataset_id", "source_state"])
    row_keys = list(rows.itertuples(index=False, name=None))
    columns = [
        ("RAM2019_ENDO_2", "singscore"),
        ("RAM2019_ENDO_2", "standardized_mean"),
        ("RAM2019_ENDO_6_SAENDO1", "singscore"),
        ("RAM2019_ENDO_6_SAENDO1", "standardized_mean"),
    ]
    matrix = np.full((len(row_keys), len(columns)), np.nan)
    annotations = np.full(matrix.shape, "", dtype=object)
    for i, (dataset, state) in enumerate(row_keys):
        for j, (program, method) in enumerate(columns):
            q = intensity[
                intensity["dataset_id"].eq(dataset)
                & intensity["source_state"].eq(state)
                & intensity["program_id"].eq(program)
                & intensity["score_method"].eq(method)
            ]
            if q.empty:
                continue
            row = q.iloc[0]
            matrix[i, j] = row["hedges_g"]
            suffix = ("*" if row["robust_ci95_low"] > 0 else "") + (
                "†" if bool(row["above_random_95th_percentile"]) else ""
            )
            annotations[i, j] = f"{row['hedges_g']:.2f}{suffix}"
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-2.2,
        vmax=2.2,
        annot=annotations,
        fmt="",
        linewidths=0.35,
        linecolor="white",
        mask=np.isnan(matrix),
        cbar_kws={"label": "State-specific Hedges g", "shrink": 0.72},
        annot_kws={"fontsize": 6.1},
    )
    ax.set_xticklabels(["Endo2\nRank", "Endo2\nz-mean", "SAEndo1\nRank", "SAEndo1\nz-mean"], rotation=0, fontsize=6.6)
    ax.set_yticklabels([f"{DATASET_LABEL[d]} | {s}" for d, s in row_keys], rotation=0, fontsize=5.8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.text(
        0,
        -0.16,
        "* HC3 lower 95% CI > 0; † above matched-random 95th percentile.\nNo positive state pair met the frozen cross-cohort marker-match rule.",
        transform=ax.transAxes,
        fontsize=6.4,
        va="top",
    )
    panel_label(ax, "B", x=-0.30)
    intensity.to_csv(source / "figure_9b_state_intensity.csv", index=False)

    # C: alcohol hepatitis versus cirrhosis context.
    ax = fig.add_subplot(grid[1, 0])
    alcohol = p5(repo, "alcohol_context_program_effects.csv")
    benchmark = p5(repo, "alcohol_context_random_benchmark.csv")
    alcohol = alcohol[
        alcohol["analysis_tier"].eq("primary")
        & alcohol["program_id"].isin(PROGRAM_LABEL)
    ].merge(
        benchmark[
            ["contrast", "program_id", "score_method", "cell_gate", "above_random_95th_percentile"]
        ],
        on=["contrast", "program_id", "score_method", "cell_gate"],
        how="left",
    )
    alcohol["program"] = alcohol["program_id"].map(PROGRAM_LABEL)
    contrast_order = ["alcohol_hepatitis_vs_healthy", "alcohol_cirrhosis_vs_alcohol_hepatitis"]
    alcohol["order"] = alcohol.apply(
        lambda row: contrast_order.index(row["contrast"]) * 4
        + (0 if row["program"] == "Endo2" else 2)
        + (0 if row["score_method"] == "singscore" else 1),
        axis=1,
    )
    alcohol = alcohol.sort_values("order")
    y = np.arange(len(alcohol))
    for index, row in enumerate(alcohol.itertuples()):
        color = PROGRAM_COLOR[row.program]
        marker = METHOD_MARKER[row.score_method]
        face = color if row.above_random_95th_percentile else "white"
        ax.errorbar(
            row.hedges_g,
            index,
            xerr=[[row.hedges_g - row.robust_ci95_low], [row.robust_ci95_high - row.hedges_g]],
            fmt=marker,
            ms=5.2,
            mfc=face,
            mec=color,
            ecolor=color,
            elinewidth=1.0,
            capsize=2,
        )
    ax.axvline(0, color=GRAY, lw=0.8, ls="--")
    ax.set_yticks(
        y,
        [
            f"{'AH vs healthy' if row.contrast == contrast_order[0] else 'Alcohol cirrhosis vs AH'} | {row.program} | {METHOD_LABEL[row.score_method]}"
            for row in alcohol.itertuples()
        ],
        fontsize=6.1,
    )
    ax.invert_yaxis()
    ax.set_xlabel("Hedges g (HC3 95% CI)")
    clean_axes(ax)
    panel_label(ax, "C", x=-0.34)
    alcohol.to_csv(source / "figure_9c_alcohol_context.csv", index=False)

    # D: 2022 independent donor-by-region spatial differences.
    ax = fig.add_subplot(grid[1, 1])
    spatial_scores = p5(repo, "spatial_2022_donor_region_scores.csv")
    spatial_scores = spatial_scores[spatial_scores["gene_set_id"].isin(["Endo2_primary", "SAEndo1_primary"])].copy()
    spatial_scores["label"] = spatial_scores.apply(
        lambda row: f"{row['gene_set_id'].replace('_primary', '')}\n{METHOD_LABEL[row['score_method']]}", axis=1
    )
    labels = ["Endo2\nRank", "Endo2\nz-mean", "SAEndo1\nRank", "SAEndo1\nz-mean"]
    for xpos, label in enumerate(labels):
        data = spatial_scores[spatial_scores["label"].eq(label)]["scar_minus_parenchyma"].to_numpy(float)
        jitter = np.linspace(-0.11, 0.11, len(data))
        color = PROGRAM_COLOR[label.split("\n")[0]]
        ax.scatter(np.full(len(data), xpos) + jitter, data, color=color, s=22, alpha=0.8, edgecolor="white", linewidth=0.35)
        ax.plot([xpos - 0.20, xpos + 0.20], [np.median(data), np.median(data)], color=BLACK, lw=1.6)
    ax.axhline(0, color=GRAY, lw=0.8, ls="--")
    ax.set_xticks(np.arange(len(labels)), labels, fontsize=6.7)
    ax.set_ylabel("Scar − parenchyma score")
    ax.set_xlabel("Eight independent cirrhotic explants")
    clean_axes(ax)
    panel_label(ax, "D")
    spatial_scores.to_csv(source / "figure_9d_spatial_2022_donor_scores.csv", index=False)

    # E: 2022 matched-random calibration.
    ax = fig.add_subplot(grid[2, 0])
    spatial_random = p5(repo, "spatial_2022_random_benchmark.csv")
    spatial_random = spatial_random[spatial_random["gene_set_id"].isin(["Endo2_primary", "SAEndo1_primary"])].copy()
    spatial_random["label"] = spatial_random.apply(
        lambda row: f"{row['gene_set_id'].replace('_primary', '')}\n{METHOD_LABEL[row['score_method']]}", axis=1
    )
    spatial_random = spatial_random.set_index("label").loc[labels].reset_index()
    x = np.arange(len(labels))
    for xpos, row in zip(x, spatial_random.itertuples()):
        color = PROGRAM_COLOR[row.label.split("\n")[0]]
        ax.plot([xpos, xpos], [row.random_95th_percentile, row.observed_median_scar_minus_parenchyma], color=LIGHT_GRAY, lw=1.4)
        ax.scatter(xpos, row.random_95th_percentile, marker="D", s=30, color=GRAY, zorder=3)
        ax.scatter(xpos, row.observed_median_scar_minus_parenchyma, marker="o", s=42, color=color, edgecolor="white", linewidth=0.45, zorder=4)
    ax.set_xticks(x, labels, fontsize=6.7)
    ax.set_ylabel("Median donor scar effect")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white", label="Observed", markersize=6),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=GRAY, label="Random 95th percentile", markersize=5),
        ],
        frameon=False,
        fontsize=6.7,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.21),
        ncol=2,
        borderaxespad=0,
    )
    clean_axes(ax)
    panel_label(ax, "E")
    spatial_random.to_csv(source / "figure_9e_spatial_2022_random.csv", index=False)

    # F: 2025 seven-sample scar-cluster recurrence.
    ax = fig.add_subplot(grid[2, 1])
    recurrence = p5(repo, "spatial_2025_scar_top100_random_benchmark.csv")
    recurrence_order = ["Endo2_primary", "SAEndo1_primary", "Endo2_secondary", "SAEndo1_secondary"]
    recurrence = recurrence.set_index("gene_set_id").loc[recurrence_order].reset_index()
    y = np.arange(len(recurrence))
    for ypos, row in zip(y, recurrence.itertuples()):
        color = BLUE if row.gene_set_id.startswith("Endo2") else PURPLE
        ax.plot([row.random_rank_weight_95th, row.observed_rank_weight], [ypos, ypos], color=LIGHT_GRAY, lw=1.6)
        ax.scatter(row.random_rank_weight_95th, ypos, marker="D", s=30, color=GRAY)
        ax.scatter(row.observed_rank_weight, ypos, s=46, color=color, edgecolor="white", linewidth=0.45)
        ax.text(row.observed_rank_weight + 0.22, ypos, f"{int(row.observed_samples_with_hit)}/7 samples", va="center", fontsize=6.4)
    ax.set_yticks(y, ["Endo2 primary", "SAEndo1 primary", "Endo2 secondary", "SAEndo1 secondary"], fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("Scar top-100 recurrence rank weight")
    ax.set_xlim(left=0)
    clean_axes(ax)
    panel_label(ax, "F", x=-0.25)
    recurrence.to_csv(source / "figure_9f_spatial_2025_recurrence.csv", index=False)

    save(fig, output, "figure_9_state_context_spatial_boundary")


def supplementary_figure_10(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.8, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)

    # A: all endothelial cross-cohort state matches.
    ax = fig.add_subplot(grid[0, 0])
    matches = p5(repo, "cross_cohort_state_matches.csv")
    matches = matches[matches["lineage"].eq("endothelial")].copy()
    ax.scatter(
        matches["marker_effect_spearman_rho"],
        matches["top50_marker_jaccard"],
        s=18,
        color=np.where(matches["supported_match"], BLUE, LIGHT_GRAY),
        alpha=0.85,
        edgecolor="white",
        linewidth=0.25,
    )
    ax.axvline(0, color=GRAY, lw=0.7, ls="--")
    ax.axhline(0.10, color=GRAY, lw=0.7, ls="--")
    ax.set_xlabel("Non-program marker-effect Spearman ρ")
    ax.set_ylabel("Top-50 marker Jaccard")
    ax.text(0.98, 0.94, f"{int(matches['supported_match'].sum())} supported pairs", transform=ax.transAxes, ha="right", va="top", fontsize=6.8)
    clean_axes(ax)
    panel_label(ax, "A")
    matches.to_csv(source / "supplementary_figure_10a_state_matches.csv", index=False)

    # B: all-program alcohol context effects.
    ax = fig.add_subplot(grid[0, 1])
    alcohol = p5(repo, "alcohol_context_program_effects.csv")
    alcohol = alcohol[
        alcohol["analysis_tier"].eq("primary")
        & alcohol["contrast"].eq("alcohol_cirrhosis_vs_alcohol_hepatitis")
    ].copy()
    pivot = alcohol.pivot(index="program_id", columns="score_method", values="hedges_g")
    pivot = pivot.sort_values("standardized_mean", ascending=False)
    sns.heatmap(
        pivot[["singscore", "standardized_mean"]],
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "Hedges g", "shrink": 0.7},
    )
    ax.set_yticklabels([PROGRAM_LABELS.get(value, value) for value in pivot.index], rotation=0, fontsize=5.5)
    ax.set_xticklabels(["Rank", "z-mean"], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    panel_label(ax, "B", x=-0.24)
    alcohol.to_csv(source / "supplementary_figure_10b_alcohol_all_programs.csv", index=False)

    # C: author cell2location scar-associated endothelial fractions.
    ax = fig.add_subplot(grid[1, 0])
    coloc = p5(repo, "spatial_2025_endothelial_colocation.csv")
    coloc = coloc[coloc["metric"].eq("EC4_6_fraction_of_EC1_6")].copy()
    for index, row in enumerate(coloc.itertuples()):
        ax.plot([0, 1], [row.other_median, row.scar_median], color=PROGRAM_COLOR["Endo2"], alpha=0.75, lw=1.2)
        ax.scatter([0, 1], [row.other_median, row.scar_median], color=PROGRAM_COLOR["Endo2"], s=34, edgecolor="white", linewidth=0.4)
        ax.text(1.04, row.scar_median, f"Sample {row.donor_id}", va="center", fontsize=6.7)
    ax.set_xticks([0, 1], ["Other spots", "Author scar cluster"])
    ax.set_ylabel("EC4–6 fraction among EC1–6")
    ax.set_xlim(-0.2, 1.35)
    ax.set_ylim(0, 1.0)
    clean_axes(ax)
    panel_label(ax, "C")
    coloc.to_csv(source / "supplementary_figure_10c_spatial_ec_colocation.csv", index=False)

    # D: primary versus secondary 2025 scar recurrence by sample.
    ax = fig.add_subplot(grid[1, 1])
    hits = p5(repo, "spatial_2025_scar_top100_hits.csv")
    matrix = hits.pivot(index="gene_set_id", columns="sample_id", values="top100_hit_count")
    row_order = ["Endo2_primary", "SAEndo1_primary", "Endo2_secondary", "SAEndo1_secondary"]
    matrix = matrix.loc[row_order, list("ABCDEFG")]
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="Blues",
        vmin=0,
        vmax=max(4, int(matrix.to_numpy().max())),
        annot=True,
        fmt="g",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Program genes in scar top 100", "shrink": 0.72},
        annot_kws={"fontsize": 6.8},
    )
    ax.set_yticklabels(["Endo2 primary", "SAEndo1 primary", "Endo2 secondary", "SAEndo1 secondary"], rotation=0, fontsize=6.6)
    ax.set_xticklabels(list("ABCDEFG"), rotation=0)
    ax.set_xlabel("Independent spatial samples")
    ax.set_ylabel("")
    panel_label(ax, "D", x=-0.24)
    hits.to_csv(source / "supplementary_figure_10d_spatial_scar_hits.csv", index=False)

    save(fig, output, "supplementary_figure_10_phase5_audit")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "figures"
    source = repo / "results" / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    figure_9(repo, output, source)
    supplementary_figure_10(repo, output, source)


if __name__ == "__main__":
    main()
