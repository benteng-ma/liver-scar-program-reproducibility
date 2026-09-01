from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import seaborn as sns

from make_deep_benchmark_figures import (
    BLACK,
    BLUE,
    GREEN,
    GRAY,
    LIGHT_GRAY,
    LINEAGE_COLORS,
    LINEAGE_LABELS,
    METHOD_COLORS,
    ORANGE,
    PROGRAM_LABELS,
    PURPLE,
    RED,
    VERY_LIGHT,
    clean_axes,
    panel,
    panel_label,
    save,
)


plt.rcParams.update({"font.family": "Arial", "font.size": 8.0, "pdf.fonttype": 42, "ps.fonttype": 42})

CONTRAST_LABELS = {
    "mash_cirrhosis_vs_healthy": "MASH cirrhosis\nvs healthy",
    "alcohol_cirrhosis_vs_healthy": "Alcohol cirrhosis\nvs healthy",
    "mash_fibrosis_vs_masld_f0": "MASH fibrosis\nvs MASLD F0",
    "mash_cirrhosis_vs_masld_f0": "MASH cirrhosis\nvs MASLD F0",
    "mash_vs_alcohol_cirrhosis_etiology": "MASH vs alcohol\ncirrhosis",
}


def program_order(repo: Path) -> list[str]:
    inventory = pd.read_csv(repo / "literature" / "program_inventory.csv")
    order = []
    for lineage in ("endothelial", "macrophage_monocyte", "mesenchymal_hsc_myofibroblast"):
        order.extend(sorted(inventory.loc[inventory["cell_lineage"].eq(lineage), "program_id"].unique()))
    return order


def figure_7(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.2, 9.0), constrained_layout=True)
    grid = GridSpec(3, 2, figure=fig, height_ratios=[0.88, 1.35, 1.0])

    # A: cohort geometry.
    ax = fig.add_subplot(grid[0, 0])
    donors = pd.read_csv(repo / "metadata" / "gse256398_donor_manifest.csv")
    group_order = ["healthy", "masld_f0", "mash_fibrosis", "mash_cirrhosis", "alcohol_cirrhosis", "alcohol_hepatitis"]
    labels = ["Healthy", "MASLD F0", "MASH F2-3/F3", "MASH cirrhosis", "Alcohol cirrhosis", "Alcohol hepatitis"]
    colors = [GRAY, "#7FB3D5", "#F5B041", RED, PURPLE, "#B39DDB"]
    counts = donors["disease_group"].value_counts().reindex(group_order).fillna(0)
    y = np.arange(len(group_order))
    ax.barh(y, counts, color=colors, edgecolor="white", height=0.72)
    for yy, value in zip(y, counts):
        ax.text(value + 0.08, yy, f"n={int(value)}", va="center", fontsize=7.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 7)
    ax.set_xlabel("Independent human donors")
    ax.text(0.99, 0.04, "159,784 QC-retained nuclei\n28,895 target-lineage nuclei\n26 donors × 3 pseudobulks", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.4, color=BLACK, bbox=dict(boxstyle="round,pad=0.3", facecolor=VERY_LIGHT, edgecolor=LIGHT_GRAY))
    panel_label(ax, "A")
    clean_axes(ax)
    donors[["donor_id", "disease_group", "etiology", "fibrosis_stage", "metabolic_order"]].to_csv(source / "figure_7a_gse256398_donors.csv", index=False)

    # B: sixth-cohort effect landscape.
    ax = fig.add_subplot(grid[1, 0])
    effects = pd.read_csv(repo / "results" / "phase3" / "gse256398_program_effects.csv")
    random = pd.read_csv(repo / "results" / "random_controls" / "gse256398_random_module_benchmark.csv")
    keys = ["contrast", "program_id", "score_method"]
    merged = effects.merge(random[keys + ["above_random_95th_percentile"]], on=keys)
    program_context = merged.groupby(["program_id", "contrast"], as_index=False).agg(
        mean_g=("hedges_g", "mean"),
        dual_random=("above_random_95th_percentile", "sum"),
        dual_positive_ci=("robust_ci95_low", lambda x: int((x > 0).sum())),
        dual_positive=("hedges_g", lambda x: int((x > 0).sum())),
    )
    order = program_order(repo)
    contrast_order = list(CONTRAST_LABELS)
    heat = program_context.pivot(index="program_id", columns="contrast", values="mean_g").reindex(index=order, columns=contrast_order)
    sns.heatmap(heat, ax=ax, cmap="vlag", center=0, vmin=-3, vmax=3, linewidths=0.25, linecolor="white", cbar_kws={"label": "Mean Hedges g\n(two scores)", "shrink": 0.7})
    ax.set_xticklabels([CONTRAST_LABELS[value] for value in contrast_order], rotation=25, ha="right")
    ax.set_yticklabels([PROGRAM_LABELS[value] for value in order], rotation=0, fontsize=6.6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    lookup = program_context.set_index(["program_id", "contrast"])
    for row, program in enumerate(order):
        for column, contrast in enumerate(contrast_order):
            if (program, contrast) not in lookup.index:
                continue
            value = lookup.loc[(program, contrast)]
            if value["dual_random"] == 2:
                ax.text(column + 0.5, row + 0.53, "●", ha="center", va="center", fontsize=5.2, color="white" if abs(heat.loc[program, contrast]) > 1.3 else BLACK)
            if value["dual_random"] == 2 and value["dual_positive_ci"] == 2 and value["dual_positive"] == 2:
                ax.add_patch(patches.Rectangle((column + 0.04, row + 0.04), 0.92, 0.92, fill=False, lw=1.2, edgecolor=BLACK))
    ax.text(0, -0.34, "● both scores > matched-random 95th; box adds positive HC3 intervals", transform=ax.transAxes, fontsize=7.0)
    panel_label(ax, "B", x=-0.14, y=1.04)
    program_context.to_csv(source / "figure_7b_gse256398_effect_landscape.csv", index=False)

    # C: selected cross-context concordance.
    ax = fig.add_subplot(grid[0, 1])
    concordance = pd.read_csv(repo / "results" / "phase3" / "phase3_selected_pair_concordance.csv")
    pair_order = [
        "mash_cross_cohort_stage_mismatched",
        "same_assay_cirrhosis_etiology_contrast",
        "same_assay_mash_stage_contrast",
        "mash_cirrhosis_endpoint_aligned_cross_cohort",
        "advanced_fibrosis_cross_assay_reference",
    ]
    pair_labels = ["MASH cross-cohort", "Same-assay etiology", "Same-assay MASH stage", "MASH cirrhosis cross-cohort", "Advanced-fibrosis reference"]
    column_order = [(lineage, method) for lineage in LINEAGE_LABELS for method in ("singscore", "standardized_mean")]
    matrix = pd.DataFrame(index=pair_order, columns=[f"{lineage}|{method}" for lineage, method in column_order], dtype=float)
    annotations = matrix.copy().astype(object)
    for row in concordance.itertuples(index=False):
        key = f"{row.lineage}|{row.score_method}"
        matrix.loc[row.pair_id, key] = row.spearman_rho
        annotations.loc[row.pair_id, key] = f"{row.spearman_rho:.2f}\n{row.sign_agreement:.0%}"
    sns.heatmap(matrix, ax=ax, cmap="vlag", center=0, vmin=-1, vmax=1, annot=annotations, fmt="", annot_kws={"fontsize": 5.8}, linewidths=0.4, linecolor="white", cbar_kws={"label": "Spearman rho", "shrink": 0.72})
    ax.set_yticklabels(pair_labels, rotation=0, fontsize=6.6)
    ax.set_xticklabels([f"{LINEAGE_LABELS[lineage].split('/')[0]}\n{'rank' if method == 'singscore' else 'z-mean'}" for lineage, method in column_order], rotation=0, fontsize=5.9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.text(0, -0.24, "Cells show rank rho / sign agreement", transform=ax.transAxes, fontsize=6.8)
    panel_label(ax, "C")
    concordance.to_csv(source / "figure_7c_selected_pair_concordance.csv", index=False)

    # D: two independent ordinal progression analyses.
    ax = fig.add_subplot(grid[1, 1])
    trends_244 = pd.read_csv(repo / "results" / "phase3" / "gse244832_metabolic_progression.csv")
    trends_244["cohort"] = "GSE244832: normal→MASL→MASH"
    trends_244["fdr"] = trends_244["fdr_within_lineage_method"]
    trends_256 = pd.read_csv(repo / "results" / "phase3" / "gse256398_metabolic_ordinal_trends.csv")
    trends_256["cohort"] = "GSE256398: F0→fibrosis→cirrhosis"
    trends_256["fdr"] = trends_256["permutation_fdr_within_lineage_method"]
    trends = pd.concat([trends_244, trends_256], ignore_index=True, sort=False)
    pivot = trends.pivot_table(index=["cohort", "lineage", "program_id"], columns="score_method", values=["spearman_rho", "fdr"], aggfunc="first").reset_index()
    markers = {"GSE244832: normal→MASL→MASH": "o", "GSE256398: F0→fibrosis→cirrhosis": "s"}
    for cohort, values in pivot.groupby("cohort"):
        for lineage, lineage_values in values.groupby("lineage"):
            ax.scatter(lineage_values[("spearman_rho", "singscore")], lineage_values[("spearman_rho", "standardized_mean")], s=38, marker=markers[cohort], color=LINEAGE_COLORS[lineage], alpha=0.78, edgecolor="white", linewidth=0.5)
            significant = lineage_values[((lineage_values[("fdr", "singscore")] < 0.05) | (lineage_values[("fdr", "standardized_mean")] < 0.05))]
            for row in significant.itertuples(index=False):
                x = getattr(row, "_3") if False else None
            for _, row in significant.iterrows():
                ax.annotate(PROGRAM_LABELS[row[("program_id", "")]], (row[("spearman_rho", "singscore")], row[("spearman_rho", "standardized_mean")]), xytext=(3, 3), textcoords="offset points", fontsize=5.8)
    ax.axhline(0, color=LIGHT_GRAY, lw=0.8)
    ax.axvline(0, color=LIGHT_GRAY, lw=0.8)
    ax.plot([-1, 1], [-1, 1], ls="--", lw=0.7, color=GRAY)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("Singscore ordinal rho")
    ax.set_ylabel("Standardized-mean ordinal rho")
    handles = [plt.Line2D([0], [0], marker=marker, color="none", markerfacecolor=GRAY, label=cohort, markersize=6) for cohort, marker in markers.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, borderaxespad=0, frameon=False, fontsize=6.6)
    panel_label(ax, "D", x=-0.14, y=1.04)
    clean_axes(ax)
    trends.to_csv(source / "figure_7d_metabolic_ordinal_trends.csv", index=False)

    # E: pathway transfer.
    ax = fig.add_subplot(grid[2, 0])
    pathway = pd.read_csv(repo / "results" / "phase3" / "reactome_pathway_transfer_pairwise.csv")
    definitions = [
        ("Same-cohort\ncirrhosis etiologies", "GSE256398_human::alcohol_cirrhosis_vs_healthy", "GSE256398_human::mash_cirrhosis_vs_healthy"),
        ("Cross-cohort\nMASH", "GSE244832::mash_f2f4_group_vs_normal_sensitivity", "GSE256398_human::mash_cirrhosis_vs_healthy"),
        ("Cross-cohort\ncirrhosis", "GSE202379::advanced_f3f4_vs_f0_non_end_stage", "GSE256398_human::mash_cirrhosis_vs_healthy"),
    ]
    plot_rows = []
    for label, left, right in definitions:
        match = pathway[((pathway["context_left"].eq(left) & pathway["context_right"].eq(right)) | (pathway["context_left"].eq(right) & pathway["context_right"].eq(left)))].copy()
        match["comparison"] = label
        plot_rows.append(match)
    pathway_plot = pd.concat(plot_rows, ignore_index=True)
    x = np.arange(len(definitions))
    width = 0.23
    for index, lineage in enumerate(LINEAGE_LABELS):
        values = pathway_plot[pathway_plot["lineage"].eq(lineage)].set_index("comparison")["spearman_rho_nes"].reindex([item[0] for item in definitions])
        ax.bar(x + (index - 1) * width, values, width=width, color=LINEAGE_COLORS[lineage], label=LINEAGE_LABELS[lineage])
    ax.axhline(0, color=BLACK, lw=0.7)
    ax.set_xticks(x, [item[0] for item in definitions])
    ax.set_ylim(-0.6, 0.9)
    ax.set_ylabel("Reactome NES rank rho")
    ax.legend(frameon=False, ncol=3, fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.20), borderaxespad=0)
    panel_label(ax, "E")
    clean_axes(ax)
    pathway_plot.to_csv(source / "figure_7e_reactome_transfer.csv", index=False)

    # F: weight robustness.
    ax = fig.add_subplot(grid[2, 1])
    weights = pd.read_csv(repo / "results" / "phase3" / "report_card_weight_sensitivity.csv").sort_values("top_five_probability", ascending=True).tail(9)
    y = np.arange(len(weights))
    colors = [LINEAGE_COLORS[value] for value in weights["lineage"]]
    ax.barh(y, weights["top_five_probability"], color=colors, alpha=0.9)
    ax.set_yticks(y, [PROGRAM_LABELS[value] for value in weights["program_id"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Probability of top-five rank\n(100,000 Dirichlet weight draws)")
    for yy, value in zip(y, weights["top_five_probability"]):
        ax.text(min(value + 0.02, 0.96), yy, f"{value:.0%}", va="center", fontsize=6.8)
    panel_label(ax, "F")
    clean_axes(ax)
    weights.to_csv(source / "figure_7f_report_card_weight_sensitivity.csv", index=False)

    save(fig, output, "figure_7_post_lock_external_enrichment")


def supplementary_1(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    cluster_audit = pd.read_csv(repo / "results" / "qc" / "gse256398_cluster_marker_audit.csv")
    cluster_mapping = pd.read_csv(repo / "metadata" / "gse256398_cluster_mapping.csv")
    clusters = cluster_audit.merge(
        cluster_mapping[["cluster", "harmonized_lineage", "target_included", "mapping_status"]],
        on="cluster",
        how="left",
    )
    ax = fig.add_subplot(grid[0, 0])
    target = clusters["target_included"].astype(str).str.lower().eq("yes")
    ax.scatter(clusters.loc[~target, "top_minus_second_score"], clusters.loc[~target, "winner_anchors_z_gt_0_5"], s=np.sqrt(clusters.loc[~target, "n_cells"]) * 3, color=LIGHT_GRAY, edgecolor=GRAY, linewidth=0.3)
    for lineage, values in clusters[target].groupby("harmonized_lineage"):
        ax.scatter(values["top_minus_second_score"], values["winner_anchors_z_gt_0_5"], s=np.sqrt(values["n_cells"]) * 3, color=LINEAGE_COLORS[lineage], label=LINEAGE_LABELS[lineage], edgecolor="white", linewidth=0.4)
    ax.axvline(0.5, ls="--", color=GRAY, lw=0.8)
    ax.axhline(2, ls="--", color=GRAY, lw=0.8)
    ax.set_xlabel("Top-minus-second lineage score")
    ax.set_ylabel("Winning marker anchors (z>0.5)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=3, borderaxespad=0)
    panel(ax, "A", "Frozen, outcome-blind cluster identity rule")
    clean_axes(ax)

    eligibility = pd.read_csv(repo / "metadata" / "gse256398_donor_lineage_eligibility.csv")
    ax = fig.add_subplot(grid[0, 1])
    medians = eligibility.groupby(["disease_group", "harmonized_lineage"])["n_cells"].median().unstack().reindex(["healthy", "masld_f0", "mash_fibrosis", "mash_cirrhosis", "alcohol_cirrhosis", "alcohol_hepatitis"])
    sns.heatmap(np.log10(medians + 1), annot=medians.astype(int), fmt="d", cmap="Blues", ax=ax, cbar_kws={"label": "log10(median cells+1)"})
    ax.set_yticklabels(["Healthy", "MASLD F0", "MASH fibrosis", "MASH cirrhosis", "Alcohol cirrhosis", "Alcohol hepatitis"], rotation=0)
    ax.set_xticklabels(["Endothelial", "Macrophage", "Mesenchymal"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    panel(ax, "B", "Median retained target cells per donor")

    ax = fig.add_subplot(grid[1, 0])
    coverage = pd.read_csv(repo / "results" / "qc" / "gse256398_program_coverage.csv").sort_values("coverage")
    ax.barh(np.arange(len(coverage)), coverage["coverage"], color=[LINEAGE_COLORS[value] for value in coverage["lineage"]])
    ax.set_yticks(np.arange(len(coverage)), [PROGRAM_LABELS[value] for value in coverage["program_id"]], fontsize=6.5)
    ax.axvline(0.8, ls="--", color=BLACK, lw=0.8)
    ax.set_xlim(0.75, 1.01)
    ax.set_xlabel("Measured program-gene fraction")
    panel(ax, "C", "All 19 programs pass primary feature coverage")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    gates = pd.read_csv(repo / "results" / "qc" / "gse256398_donor_gate_summary.csv")
    summary = gates.groupby(["contrast", "lineage"])[["donors_30"]].min().reset_index()
    heat = summary.pivot(index="contrast", columns="lineage", values="donors_30").reindex(list(CONTRAST_LABELS))
    sns.heatmap(heat, annot=True, fmt=".0f", cmap="Greens", vmin=0, vmax=6, ax=ax, cbar_kws={"label": "Minimum eligible donors/group"})
    ax.set_yticklabels([CONTRAST_LABELS[value].replace("\n", " ") for value in heat.index], rotation=0, fontsize=6.8)
    ax.set_xticklabels(["Endothelial", "Macrophage", "Mesenchymal"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    panel(ax, "D", "All fixed binary contrasts pass the 30-cell donor gate")
    save(fig, output, "supplementary_figure_1_gse256398_qc")
    clusters.to_csv(source / "supplementary_figure_1a_cluster_qc.csv", index=False)
    medians.to_csv(source / "supplementary_figure_1b_cell_counts.csv")
    coverage.to_csv(source / "supplementary_figure_1c_coverage.csv", index=False)
    summary.to_csv(source / "supplementary_figure_1d_gates.csv", index=False)


def supplementary_2(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    overlap = pd.read_csv(repo / "results" / "phase3" / "program_gene_overlap.csv")
    order = program_order(repo)
    matrix = pd.DataFrame(np.nan, index=order, columns=order)
    for row in overlap.itertuples(index=False):
        matrix.loc[row.program_left, row.program_right] = row.jaccard
        matrix.loc[row.program_right, row.program_left] = row.jaccard
    ax = fig.add_subplot(grid[:, 0])
    sns.heatmap(matrix, cmap="mako", vmin=0, vmax=0.25, mask=matrix.isna(), ax=ax, square=True, cbar_kws={"label": "Jaccard overlap", "shrink": 0.55})
    ax.set_xticklabels([PROGRAM_LABELS[value] for value in order], rotation=90, fontsize=5.8)
    ax.set_yticklabels([PROGRAM_LABELS[value] for value in order], rotation=0, fontsize=5.8)
    panel(ax, "A", "Published programs share few genes but form correlated scores")

    effective = pd.read_csv(repo / "results" / "phase3" / "program_effective_test_counts.csv")
    ax = fig.add_subplot(grid[0, 1])
    sns.boxplot(data=effective, x="lineage", y="effective_fraction", order=list(LINEAGE_LABELS), palette=LINEAGE_COLORS, ax=ax, width=0.55, fliersize=0)
    sns.stripplot(data=effective, x="lineage", y="effective_fraction", order=list(LINEAGE_LABELS), color=BLACK, size=3, alpha=0.55, ax=ax)
    ax.set_xticklabels(["Endothelial", "Macrophage", "Mesenchymal"], rotation=15)
    ax.set_ylim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("Effective tests / nominal programs")
    panel(ax, "B", "Donor-score correlation reduces effective dimensionality")
    clean_axes(ax)

    weights = pd.read_csv(repo / "results" / "phase3" / "report_card_weight_sensitivity.csv").sort_values("rank_median").head(10)
    ax = fig.add_subplot(grid[1, 1])
    y = np.arange(len(weights))
    ax.errorbar(weights["rank_median"], y, xerr=[weights["rank_median"] - weights["rank_ci95_low"], weights["rank_ci95_high"] - weights["rank_median"]], fmt="o", color=BLUE, ecolor=GRAY, capsize=2)
    ax.set_yticks(y, [PROGRAM_LABELS[value] for value in weights["program_id"]])
    ax.invert_yaxis()
    ax.set_xlabel("Rank under 100,000 weight draws")
    panel(ax, "C", "Report-card ordering is partly weight-sensitive")
    clean_axes(ax)
    save(fig, output, "supplementary_figure_2_redundancy_weights")
    overlap.to_csv(source / "supplementary_figure_2a_overlap.csv", index=False)
    effective.to_csv(source / "supplementary_figure_2b_effective_tests.csv", index=False)
    weights.to_csv(source / "supplementary_figure_2c_weight_ranks.csv", index=False)


def supplementary_3(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.1), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    univ = pd.read_csv(repo / "results" / "phase3" / "gse202379_histology_univariate.csv")
    univ = univ[univ["analysis_set"].eq("non_end_stage_primary_exploratory")]
    for index, method in enumerate(("singscore", "standardized_mean")):
        ax = fig.add_subplot(grid[0, index])
        values = univ[univ["score_method"].eq(method)]
        heat = values.pivot(index="program_id", columns="histology_axis", values="spearman_rho").reindex(index=program_order(repo), columns=["fibrosis", "steatosis", "ballooning", "inflammation"])
        sns.heatmap(heat, cmap="vlag", center=0, vmin=-0.7, vmax=0.7, ax=ax, linewidths=0.2, linecolor="white", cbar_kws={"label": "Spearman rho", "shrink": 0.7})
        ax.set_yticklabels([PROGRAM_LABELS[value] for value in heat.index], rotation=0, fontsize=5.8)
        ax.set_xticklabels(["Fibrosis", "Steatosis", "Ballooning", "Inflammation"], rotation=20, ha="right")
        ax.set_xlabel("")
        ax.set_ylabel("")
        panel(ax, "A" if index == 0 else "B", f"GSE202379 histology specificity: {'rank score' if method == 'singscore' else 'standardized mean'}")

    for index, (path, title) in enumerate((("gse244832_metabolic_progression.csv", "GSE244832 normal→MASL→MASH"), ("gse256398_metabolic_ordinal_trends.csv", "GSE256398 F0→fibrosis→cirrhosis"))):
        ax = fig.add_subplot(grid[1, index])
        trends = pd.read_csv(repo / "results" / "phase3" / path)
        fdr_col = "fdr_within_lineage_method" if "fdr_within_lineage_method" in trends else "permutation_fdr_within_lineage_method"
        pivot = trends.pivot_table(index=["lineage", "program_id"], columns="score_method", values=["spearman_rho", fdr_col], aggfunc="first").reset_index()
        for lineage, values in pivot.groupby("lineage"):
            ax.scatter(values[("spearman_rho", "singscore")], values[("spearman_rho", "standardized_mean")], color=LINEAGE_COLORS[lineage], s=42, alpha=0.78, edgecolor="white", linewidth=0.5, label=LINEAGE_LABELS[lineage])
            significant = values[(values[(fdr_col, "singscore")] < 0.05) | (values[(fdr_col, "standardized_mean")] < 0.05)]
            for _, row in significant.iterrows():
                ax.annotate(PROGRAM_LABELS[row[("program_id", "")]], (row[("spearman_rho", "singscore")], row[("spearman_rho", "standardized_mean")]), xytext=(3, 3), textcoords="offset points", fontsize=6)
        ax.axhline(0, color=LIGHT_GRAY, lw=0.7)
        ax.axvline(0, color=LIGHT_GRAY, lw=0.7)
        ax.plot([-1, 1], [-1, 1], ls="--", color=GRAY, lw=0.7)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xlabel("Singscore rho")
        ax.set_ylabel("Standardized-mean rho")
        if index == 1:
            ax.legend(frameon=False, fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, borderaxespad=0)
        panel(ax, "C" if index == 0 else "D", title)
        clean_axes(ax)
    save(fig, output, "supplementary_figure_3_histology_progression")
    univ.to_csv(source / "supplementary_figure_3ab_histology.csv", index=False)


def supplementary_4(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    transfer = pd.read_csv(repo / "results" / "phase3" / "reactome_pathway_transfer_pairwise.csv")
    for index, lineage in enumerate(LINEAGE_LABELS):
        ax = fig.add_subplot(grid[index // 2, index % 2])
        values = transfer[transfer["lineage"].eq(lineage)].copy()
        contexts = sorted(set(values["context_left"]) | set(values["context_right"]))
        matrix = pd.DataFrame(np.eye(len(contexts)), index=contexts, columns=contexts)
        for row in values.itertuples(index=False):
            matrix.loc[row.context_left, row.context_right] = row.spearman_rho_nes
            matrix.loc[row.context_right, row.context_left] = row.spearman_rho_nes
        short = [value.replace("GSE256398_human::", "256:").replace("GSE244832::", "244:").replace("GSE202379::", "202:").replace("mash_f2f4_group_vs_normal_sensitivity", "MASH").replace("advanced_f3f4_vs_f0_non_end_stage", "F3-4").replace("mash_cirrhosis_vs_healthy", "MASH-C").replace("alcohol_cirrhosis_vs_healthy", "ALD-C").replace("mash_fibrosis_vs_masld_f0", "MASH-F") for value in contexts]
        sns.heatmap(matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", annot_kws={"fontsize": 6}, ax=ax, cbar=index == 0, cbar_kws={"label": "NES rank rho", "shrink": 0.65})
        ax.set_xticklabels(short, rotation=35, ha="right", fontsize=6)
        ax.set_yticklabels(short, rotation=0, fontsize=6)
        panel(ax, chr(ord("A") + index), f"Reactome transfer: {LINEAGE_LABELS[lineage]}")

    ax = fig.add_subplot(grid[1, 1])
    recurrence = pd.read_csv(repo / "results" / "phase3" / "reactome_recurrent_pathways.csv")
    positive = recurrence[(recurrence["direction"].eq("positive")) & (recurrence["significant_contexts"].ge(2))].sort_values("minimum_absolute_nes", ascending=False).head(15).sort_values("minimum_absolute_nes")
    ax.barh(np.arange(len(positive)), positive["minimum_absolute_nes"], color=[LINEAGE_COLORS[value] for value in positive["lineage"]])
    ax.set_yticks(np.arange(len(positive)), [textwrap.fill(value, 33) for value in positive["pathway_name"]], fontsize=5.8)
    ax.set_xlabel("Minimum |NES| across recurrent contexts")
    panel(ax, "D", "Positive pathways recurring in at least two contexts")
    clean_axes(ax)
    save(fig, output, "supplementary_figure_4_reactome_transfer")
    transfer.to_csv(source / "supplementary_figure_4abc_pathway_transfer.csv", index=False)
    positive.to_csv(source / "supplementary_figure_4d_recurrent_positive_pathways.csv", index=False)


def supplementary_5(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    grid_data = pd.read_csv(repo / "results" / "phase3" / "threshold_sensitivity_grid.csv")
    for index, gate in enumerate((20, 30)):
        ax = fig.add_subplot(grid[0, index])
        values = grid_data[grid_data["cell_gate"].eq(gate)].pivot(index="coverage_threshold", columns="random_percentile_threshold", values="dual_score_positive_ci_random_contexts")
        sns.heatmap(values, annot=True, fmt="d", cmap="YlOrRd", ax=ax, cbar=False)
        ax.set_xlabel("Matched-random percentile")
        ax.set_ylabel("Program coverage")
        panel(ax, "A" if index == 0 else "B", f"Positive-CI random-specific contexts: {gate}-cell gate")

    ax = fig.add_subplot(grid[1, 0])
    curves = pd.read_csv(repo / "results" / "phase3" / "balanced_sample_size_power_curves.csv")
    for effect, values in curves.groupby("standardized_effect"):
        ax.plot(values["balanced_n_per_group"], values["power"], label=f"d={effect:g}")
    ax.axhline(0.8, ls="--", color=BLACK, lw=0.8)
    ax.set_xlabel("Balanced donors per group")
    ax.set_ylabel("Two-sided power")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, ncol=5, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.18), borderaxespad=0)
    panel(ax, "C", "Prospective power depends steeply on effect magnitude")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    precision = pd.read_csv(repo / "results" / "phase3" / "contrast_precision_and_mde.csv")
    order = ["GSE202379", "GSE244832", "GSE256398_human", "GSE290642_human", "GSE210077_Watson6"]
    sns.boxplot(data=precision, x="dataset_id", y="minimum_detectable_cohens_d_80_power", order=order, color="#A9CCE3", fliersize=0, ax=ax)
    ax.set_xticklabels(["202379", "244832", "256398", "290642", "Watson6"], rotation=20)
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Minimum detectable standardized effect")
    panel(ax, "D", "Current cohorts are powered only for large effects")
    clean_axes(ax)
    save(fig, output, "supplementary_figure_5_threshold_precision")
    grid_data.to_csv(source / "supplementary_figure_5ab_threshold_grid.csv", index=False)
    curves.to_csv(source / "supplementary_figure_5c_power_curves.csv", index=False)
    precision.to_csv(source / "supplementary_figure_5d_precision.csv", index=False)


def supplementary_6(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, width_ratios=[1.45, 1.0])
    boot = pd.read_csv(repo / "results" / "phase3" / "donor_bootstrap_program_rank_stability.csv")
    context_labels = {
        ("GSE202379", "clinical_cirrhosis_vs_healthy"): "202:Cirrhosis",
        ("GSE202379", "advanced_f3f4_vs_f0_non_end_stage"): "202:F3-4",
        ("GSE244832", "mash_f2f4_group_vs_normal_sensitivity"): "244:MASH",
        ("GSE256398_human", "mash_cirrhosis_vs_healthy"): "256:MASH-C",
        ("GSE256398_human", "alcohol_cirrhosis_vs_healthy"): "256:ALD-C",
        ("GSE256398_human", "mash_fibrosis_vs_masld_f0"): "256:MASH-F",
        ("GSE256398_human", "mash_cirrhosis_vs_masld_f0"): "256:MASH-C/F0",
        ("GSE256398_human", "mash_vs_alcohol_cirrhosis_etiology"): "256:Etiology",
    }
    boot["context"] = [context_labels.get((dataset, contrast), contrast) for dataset, contrast in zip(boot["dataset_id"], boot["contrast"])]
    dual = boot.groupby(["program_id", "context"], as_index=False).agg(
        min_positive_probability=("positive_effect_probability", "min"),
        min_top_five_probability=("top_five_probability", "min"),
        maximum_rank_interval_width=("rank_ci95_high", "max"),
    )
    order = program_order(repo)
    context_order = list(dict.fromkeys(context_labels.values()))
    ax = fig.add_subplot(grid[:, 0])
    heat = dual.pivot(index="program_id", columns="context", values="min_positive_probability").reindex(index=order, columns=context_order)
    sns.heatmap(heat, cmap="YlGnBu", vmin=0, vmax=1, ax=ax, linewidths=0.2, linecolor="white", cbar_kws={"label": "Minimum probability across scores", "shrink": 0.6})
    ax.set_xticklabels(context_order, rotation=35, ha="right", fontsize=6.2)
    ax.set_yticklabels([PROGRAM_LABELS[value] for value in order], rotation=0, fontsize=6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    panel(ax, "A", "Donor bootstrap shows context-specific direction stability")

    ax = fig.add_subplot(grid[0, 1])
    widths = boot.assign(rank_width=boot["rank_ci95_high"] - boot["rank_ci95_low"])
    order_dataset = ["GSE202379", "GSE244832", "GSE256398_human"]
    sns.boxplot(data=widths, x="dataset_id", y="rank_width", order=order_dataset, color="#D6EAF8", fliersize=0, ax=ax)
    ax.set_xticks(np.arange(3), ["202379", "244832", "256398"])
    ax.set_xlabel("Dataset")
    ax.set_ylabel("95% bootstrap rank-interval width")
    panel(ax, "B", "Small donor groups leave broad rank uncertainty")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    selected = (
        dual.groupby("program_id", as_index=False)
        .agg(
            median_min_top_five_probability=("min_top_five_probability", "median"),
            contexts_top_five_probability_ge_080=("min_top_five_probability", lambda x: int((x >= 0.80).sum())),
        )
        .sort_values("median_min_top_five_probability", ascending=False)
        .head(10)
        .sort_values("median_min_top_five_probability")
    )
    lineage_lookup = pd.read_csv(repo / "literature" / "program_inventory.csv").groupby("program_id")["cell_lineage"].first().to_dict()
    ax.barh(
        np.arange(len(selected)),
        selected["median_min_top_five_probability"],
        color=[LINEAGE_COLORS[lineage_lookup[value]] for value in selected["program_id"]],
    )
    ax.set_yticks(np.arange(len(selected)), [PROGRAM_LABELS[value] for value in selected["program_id"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Median across contexts of minimum\ntop-five probability across scores")
    panel(ax, "C", "Rank stability remains program- and context-dependent")
    clean_axes(ax)
    save(fig, output, "supplementary_figure_6_bootstrap_stability")
    boot.to_csv(source / "supplementary_figure_6_bootstrap_stability.csv", index=False)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "figures"
    source = repo / "results" / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    figure_7(repo, output, source)
    supplementary_1(repo, output, source)
    supplementary_2(repo, output, source)
    supplementary_3(repo, output, source)
    supplementary_4(repo, output, source)
    supplementary_5(repo, output, source)
    supplementary_6(repo, output, source)
    print("Phase 3 Figure 7 and Supplementary Figures S1-S6 created")


if __name__ == "__main__":
    main()
