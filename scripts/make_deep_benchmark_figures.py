from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import seaborn as sns


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.2,
        "axes.titlesize": 9.2,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.7,
    }
)
sns.set_theme(style="white", context="paper")

BLUE = "#2B6CB0"
RED = "#C94C4C"
ORANGE = "#D98C20"
GREEN = "#2F855A"
PURPLE = "#7251A3"
GRAY = "#777777"
LIGHT_GRAY = "#E6E6E6"
VERY_LIGHT = "#F5F5F5"
BLACK = "#222222"

METHOD_COLORS = {"singscore": BLUE, "standardized_mean": RED}
LINEAGE_COLORS = {
    "endothelial": BLUE,
    "macrophage_monocyte": ORANGE,
    "mesenchymal_hsc_myofibroblast": GREEN,
}
LINEAGE_LABELS = {
    "endothelial": "Endothelial",
    "macrophage_monocyte": "Macrophage/monocyte",
    "mesenchymal_hsc_myofibroblast": "Mesenchymal/HSC",
}
PROGRAM_LABELS = {
    "RAM2019_ENDO_1": "Endo1",
    "RAM2019_ENDO_2": "Endo2",
    "RAM2019_ENDO_3": "Endo3",
    "RAM2019_ENDO_4": "Endo4",
    "RAM2019_ENDO_5": "Endo5",
    "RAM2019_ENDO_6_SAENDO1": "SAEndo1",
    "RAM2019_ENDO_7_SAENDO2": "SAEndo2",
    "RAM2019_MAC_SIG_A_SAM": "SAM-A",
    "RAM2019_MAC_SIG_B_SAM": "SAM-B",
    "RAM2019_MAC_SIG_C_KC": "KC-C",
    "RAM2019_MAC_SIG_D_TMO": "TMo-D",
    "RAM2019_MAC_SIG_E_TMO": "TMo-E",
    "RAM2019_MAC_SIG_F_CDC1": "cDC1-F",
    "RAM2019_MES_HSC": "HSC",
    "RAM2019_MES_MESOTHELIAL": "Mesothelial",
    "RAM2019_MES_SAMES": "SAMes",
    "RAM2019_MES_VSMC": "VSMC",
    "RAM2019_SAMES_A": "SAMes-A",
    "RAM2019_SAMES_B": "SAMes-B",
}

REPRESENTATIVE_ENDPOINTS = {
    "GSE202379": "advanced_f3f4_vs_f0_non_end_stage",
    "GSE290642_human": "f4_vs_f0_reconstructed_label_sensitivity",
    "GSE244832": "mash_f2f4_group_vs_normal_sensitivity",
    "GSE210077_Watson6": "mixed_f2f4_fibrosis_vs_healthy_sensitivity",
    "GSE181483_human": "cirrhosis_vs_healthy_directional",
}


def panel_label(
    ax: plt.Axes,
    letter: str,
    *,
    x: float = -0.11,
    y: float = 1.02,
) -> None:
    """Add a journal-style panel letter without an in-panel narrative heading."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def panel(ax: plt.Axes, letter: str, title: str) -> None:
    """Add a panel letter and title for supplementary figures that need local headings."""
    ax.text(-0.11, 1.07, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    ax.set_title(title, loc="left", pad=6, fontweight="bold")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=3)


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.png", dpi=360, bbox_inches="tight", pad_inches=0.12, facecolor="white")
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.12, facecolor="white")
    plt.close(fig)


def read_random(repo: Path) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted((repo / "results" / "random_controls").glob("*_random_module_benchmark.csv"))]
    return pd.concat(frames, ignore_index=True, sort=False)


def representative_effects(repo: Path) -> pd.DataFrame:
    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    keep = pd.Series(False, index=effects.index)
    for dataset, contrast in REPRESENTATIVE_ENDPOINTS.items():
        keep |= effects["dataset_id"].eq(dataset) & effects["contrast"].eq(contrast)
    return effects[keep].copy()


def figure_1(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 7.7), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)

    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    ax.axis("off")
    steps = [
        ("Lineage identity\npositive controls", GREEN),
        ("Canonical disease\nresponse controls", ORANGE),
        ("19 frozen scar\nprograms", BLUE),
        ("Held-out Core5 +\nreport card", PURPLE),
    ]
    xs = np.linspace(0.02, 0.78, 4)
    for i, ((label, color), x) in enumerate(zip(steps, xs)):
        rect = patches.FancyBboxPatch(
            (x, 0.44), 0.19, 0.28, boxstyle="round,pad=0.01", facecolor="white", edgecolor=color,
            linewidth=1.5, transform=ax.transAxes
        )
        ax.add_patch(rect)
        ax.text(x + 0.095, 0.58, label, ha="center", va="center", transform=ax.transAxes)
        if i < 3:
            ax.annotate("", xy=(x + 0.24, 0.58), xytext=(x + 0.195, 0.58), xycoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
    ax.text(0.02, 0.18, "Discovery: GSE202379 + GSE244832   |   Held out: GSE290642 + Watson6 + GSE181483",
            transform=ax.transAxes, color=GRAY, fontsize=7.6)
    pd.DataFrame(
        {"step": np.arange(1, 5), "analysis": [item[0].replace("\n", " ") for item in steps],
         "purpose": ["technical sensitivity", "cohort disease-axis sensitivity", "external transportability", "rescue and prioritization"]}
    ).to_csv(source / "figure_1a_workflow.csv", index=False)

    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    cohort = pd.DataFrame(
        [
            ["GSE202379", 47, "snRNA", "Author", "Primary"],
            ["GSE290642", 24, "scRNA", "Reconstructed", "Sensitivity"],
            ["GSE244832", 18, "snRNA", "Mapped clusters", "Sensitivity"],
            ["Watson6", 6, "snRNA", "Author", "Sensitivity"],
            ["GSE181483", 4, "scRNA", "Reconstructed", "Directional"],
        ], columns=["dataset", "donors", "assay", "annotation", "evidence"]
    )
    cohort.to_csv(source / "figure_1b_cohort_summary.csv", index=False)
    colors = cohort["evidence"].map({"Primary": RED, "Sensitivity": ORANGE, "Directional": BLUE})
    ax.barh(np.arange(len(cohort)), cohort["donors"], color=colors, height=0.6)
    for y, row in cohort.iterrows():
        ax.text(row.donors + 0.7, y, f"{row.assay} | {row.annotation}", va="center", fontsize=7.2)
    ax.set_yticks(np.arange(len(cohort)), cohort["dataset"])
    ax.invert_yaxis()
    ax.set_xlim(0, 62)
    ax.set_xlabel("Independent human donors")
    clean_axes(ax)

    performance = pd.read_csv(repo / "results" / "deep_benchmark" / "identity_control_performance.csv")
    performance.to_csv(source / "figure_1c_donor_gate_counts.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    for annotation, marker in [("author", "o"), ("reconstructed", "s")]:
        subset = performance[performance["annotation"].eq(annotation)]
        for method, color in METHOD_COLORS.items():
            values = subset[subset["score_method"].eq(method)]
            ax.scatter(values["top_score_lineage_accuracy"], values["macro_one_vs_rest_auc"],
                       s=58, marker=marker, color=color, edgecolor="white", linewidth=0.5, alpha=0.9)
            for row in values.itertuples(index=False):
                if method != "standardized_mean" or row.dataset_id not in {"GSE202379", "GSE210077_Watson6"}:
                    continue
                label = str(row.dataset_id).replace("GSE210077_Watson6", "Watson6")
                offset = (-0.105, 0.018) if row.dataset_id == "GSE202379" else (0.014, -0.018)
                ax.text(
                    row.top_score_lineage_accuracy + offset[0],
                    row.macro_one_vs_rest_auc + offset[1],
                    label,
                    fontsize=6.8,
                )
    ax.axvline(0.80, color=GRAY, linestyle="--", lw=0.8)
    ax.axhline(0.90, color=GRAY, linestyle="--", lw=0.8)
    ax.set_xlim(0.60, 1.03)
    ax.set_ylim(0.60, 1.03)
    ax.set_xlabel("Top-score lineage accuracy")
    ax.set_ylabel("Macro one-vs-rest AUROC")
    clean_axes(ax)
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=BLACK, label="Author labels"),
        plt.Line2D([0], [0], marker="s", linestyle="", color=BLACK, label="Reconstructed labels"),
        plt.Line2D([0], [0], color=BLUE, lw=2, label="singscore"),
        plt.Line2D([0], [0], color=RED, lw=2, label="standardized mean"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2, borderaxespad=0)

    disease = pd.read_csv(repo / "results" / "deep_benchmark" / "disease_control_effects.csv")
    disease["cohort"] = disease["dataset_id"].replace({"GSE210077_Watson6": "Watson6", "GSE290642_human": "GSE290642", "GSE181483_human": "GSE181483"})
    disease["module"] = disease["control_id"].replace({
        "CANONICAL_SCAR_MACROPHAGE": "Scar macrophage",
        "CANONICAL_LSEC_CAPILLARIZATION": "LSEC capillarization",
        "CANONICAL_HSC_ACTIVATION": "HSC activation",
    })
    disease["column"] = disease["cohort"] + "\n" + disease["score_method"].map({"singscore": "Rank", "standardized_mean": "Z-mean"})
    disease.to_csv(source / "figure_1d_coverage_counts.csv", index=False)
    pivot = disease.pivot_table(index="module", columns="column", values="hedges_g", aggfunc="first")
    order = [
        "GSE202379\nRank", "GSE202379\nZ-mean", "GSE244832\nRank", "GSE244832\nZ-mean",
        "GSE290642\nRank", "GSE290642\nZ-mean", "Watson6\nRank", "Watson6\nZ-mean",
        "GSE181483\nRank", "GSE181483\nZ-mean",
    ]
    pivot = pivot.reindex(columns=order)
    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    sns.heatmap(pivot.clip(-3, 3), cmap="vlag", center=0, vmin=-3, vmax=3, mask=pivot.isna(),
                linewidths=0.4, linecolor="white", cbar_kws={"label": "Hedges g", "shrink": 0.75}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=42)
    save(fig, output, "figure_1_study_design_and_evidence")


def figure_2(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 7.8), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)

    eligibility = pd.read_csv(repo / "metadata" / "gse202379_donor_lineage_eligibility.csv")
    eligible = eligibility[eligibility["eligible_primary_30"].astype(bool)].copy()
    eligible["stage"] = "F" + eligible["Fibrosis.score..F0.4."].astype(int).astype(str)
    eligible.loc[eligible["Disease.status"].eq("end stage"), "stage"] = "End-stage"
    counts = eligible.groupby(["stage", "harmonized_lineage"], as_index=False)["canonical_donor_id"].nunique().rename(columns={"canonical_donor_id": "donors"})
    counts.to_csv(source / "figure_2a_stage_donor_counts.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    order = ["F0", "F1", "F2", "F3", "F4", "End-stage"]
    pivot = counts.pivot(index="stage", columns="harmonized_lineage", values="donors").fillna(0).reindex(order)
    x = np.arange(len(order))
    width = 0.25
    for i, lineage in enumerate(LINEAGE_LABELS):
        ax.bar(x + (i - 1) * width, pivot.get(lineage, 0), width, color=LINEAGE_COLORS[lineage], label=LINEAGE_LABELS[lineage])
    ax.set_xticks(x, order, rotation=20, ha="right")
    ax.set_ylabel("Donors passing 30-cell gate")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.19), borderaxespad=0)
    clean_axes(ax)

    primary = pd.read_csv(repo / "results" / "primary" / "gse202379_primary_effects.csv")
    primary.to_csv(source / "figure_2b_primary_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    display = primary.sort_values("hedges_g").reset_index(drop=True)
    y = np.arange(len(display))
    colors = display["lineage"].map(LINEAGE_COLORS)
    ax.hlines(y, display["robust_ci95_low"], display["robust_ci95_high"], color=colors, alpha=0.45, lw=1.1)
    ax.scatter(display["hedges_g"], y, c=colors, s=18, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color=BLACK, lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Hedges g with HC3 95% interval (34 rows)")
    clean_axes(ax)

    trends = pd.read_csv(repo / "results" / "exploratory" / "gse202379_stage_trends.csv")
    trends = trends[trends["analysis_set"].eq("non_end_stage_primary_exploratory")]
    trends.to_csv(source / "figure_2c_primary_specificity.csv", index=False)
    pair = trends.pivot_table(index=["program_id", "lineage"], columns="score_method", values="spearman_rho", aggfunc="first").dropna().reset_index()
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    ax.scatter(pair["singscore"], pair["standardized_mean"], c=pair["lineage"].map(LINEAGE_COLORS), s=38,
               edgecolor="white", linewidth=0.4)
    ax.axhline(0, color=GRAY, lw=0.7)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.plot([-1, 1], [-1, 1], linestyle="--", color=LIGHT_GRAY, lw=0.8)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("Stage rho: singscore")
    ax.set_ylabel("Stage rho: standardized mean")
    clean_axes(ax)

    random = pd.read_csv(repo / "results" / "random_controls" / "gse202379_random_module_benchmark.csv")
    merged = primary.merge(
        random[["contrast", "lineage", "program_id", "score_method", "real_effect_percentile", "above_random_95th_percentile"]],
        on=["contrast", "lineage", "program_id", "score_method"], how="left"
    )
    merged.to_csv(source / "figure_2d_stage_trends.csv", index=False)
    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    for method, marker in [("singscore", "o"), ("standardized_mean", "s")]:
        values = merged[merged["score_method"].eq(method)]
        ax.scatter(values["hedges_g"], values["real_effect_percentile"], marker=marker,
                   c=values["lineage"].map(LINEAGE_COLORS), s=34, edgecolor="white", linewidth=0.4)
    ax.axhline(0.95, color=BLACK, linestyle="--", lw=0.8)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.set_xlabel("Hedges g")
    ax.set_ylabel("Percentile among matched random modules")
    ax.set_ylim(0, 1.03)
    clean_axes(ax)
    save(fig, output, "figure_2_primary_and_stage")


def figure_3(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 7.8), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    within = pd.read_csv(repo / "results" / "exploratory" / "within_cohort_score_method_concordance.csv")
    across = pd.read_csv(repo / "results" / "exploratory" / "cross_cohort_program_concordance.csv")
    within.to_csv(source / "figure_3a_within_cohort_concordance.csv", index=False)
    across.to_csv(source / "figure_3b_cross_cohort_concordance.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    plot = pd.DataFrame({
        "rho": pd.concat([within["spearman_rho"], across["spearman_rho"]], ignore_index=True),
        "comparison": ["Within cohort\nscore methods"] * len(within) + ["Across cohorts\nsame method"] * len(across),
    })
    sns.boxplot(data=plot, x="comparison", y="rho", color=LIGHT_GRAY, width=0.55, fliersize=0, ax=ax)
    sns.stripplot(data=plot, x="comparison", y="rho", hue="comparison", palette=[GREEN, ORANGE], size=4, jitter=0.18, legend=False, ax=ax)
    ax.axhline(0, color=GRAY, lw=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Spearman rho")
    clean_axes(ax)

    strata = pd.read_csv(repo / "results" / "deep_benchmark" / "transfer_failure_stratified_summary.csv")
    strata.to_csv(source / "figure_3c_sign_agreement.csv", index=False)
    labels = {
        "all_pairs": "All pairs",
        "exclude_watson": "Exclude Watson",
        "minimum_coverage_ge_0_80_exclude_watson": "≥80% coverage\nexclude Watson",
        "comparable_advanced_endpoint": "Matched advanced\nfibrosis endpoint",
    }
    shown = strata[strata["stratum"].isin(labels)].copy()
    shown["label"] = shown["stratum"].map(labels)
    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    x = np.arange(len(shown))
    ax.bar(x - 0.18, shown["median_pairwise_program_spearman"], width=0.36, color=BLUE, label="Median program rho")
    ax.bar(x + 0.18, shown["program_pair_sign_agreement"], width=0.36, color=GREEN, label="Sign agreement")
    ax.set_xticks(x, shown["label"], rotation=20, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Transportability metric")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=2, borderaxespad=0)
    clean_axes(ax)

    predictions = pd.read_csv(repo / "results" / "deep_benchmark" / "leave_one_cohort_out_predictions.csv")
    predictions.to_csv(source / "figure_3d_representative_effects.csv", index=False)
    summary = pd.read_csv(repo / "results" / "deep_benchmark" / "leave_one_cohort_out_prediction_summary.csv")
    overall = summary[summary["held_out_dataset"].eq("ALL")].iloc[0]
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    colors = predictions["held_out_dataset"].eq("GSE210077_Watson6").map({True: RED, False: BLUE})
    ax.scatter(predictions["predicted_effect_z_from_other_cohorts"], predictions["observed_effect_z"],
               c=colors, s=24, alpha=0.7, edgecolor="white", linewidth=0.3)
    lim = 2.7
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", color=GRAY, lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Predicted standardized effect from other cohorts")
    ax.set_ylabel("Observed held-out standardized effect")
    ax.text(0.03, 0.95, f"rho={overall.spearman_rho:.2f}; predictive $R^2$={overall.predictive_r_squared:.2f}",
            transform=ax.transAxes, va="top")
    clean_axes(ax)

    decomposition = pd.read_csv(repo / "results" / "deep_benchmark" / "transfer_failure_decomposition_summary.csv")
    show = decomposition[
        decomposition["analysis_type"].eq("categorical_median_contrast")
        & decomposition["outcome"].eq("absolute_effect_difference")
    ].copy()
    show["label"] = show["descriptor"].map({
        "same_assay": "Same assay",
        "both_author_labelled": "Both author-labelled",
        "comparable_advanced_endpoint": "Matched advanced endpoint",
        "includes_watson": "Includes Watson",
    })
    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    y = np.arange(len(show))
    ax.hlines(y, show["bootstrap_ci95_low"], show["bootstrap_ci95_high"], color=GRAY, lw=2)
    ax.scatter(show["estimate_true_minus_false"], y, color=[BLUE, ORANGE, GREEN, RED], s=44, zorder=3)
    ax.axvline(0, color=BLACK, lw=0.8)
    ax.set_yticks(y, show["label"])
    ax.set_xlabel("Median absolute-effect difference\n(true minus false; clustered bootstrap 95% CI)")
    clean_axes(ax)
    save(fig, output, "figure_3_reproducibility_vs_transportability")


def figure_4(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 7.8), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    target = ["RAM2019_MAC_SIG_A_SAM", "RAM2019_MAC_SIG_B_SAM", "RAM2019_MAC_SIG_E_TMO"]
    scores = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz")
    scores = scores[
        scores["program_id"].isin(target) & scores["score_method"].eq("standardized_mean")
        & scores["disease_group"].isin(["normal", "MASH"]) & scores["eligible_30"].astype(bool)
    ].copy()
    scores["program"] = scores["program_id"].map(PROGRAM_LABELS)
    scores["group"] = scores["disease_group"].map({"normal": "Normal", "MASH": "MASH F2-F4"})
    scores.to_csv(source / "figure_4a_mash_macrophage_donor_scores.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    sns.boxplot(data=scores, x="program", y="score", hue="group", palette={"Normal": LIGHT_GRAY, "MASH F2-F4": ORANGE},
                width=0.65, fliersize=0, ax=ax)
    sns.stripplot(data=scores, x="program", y="score", hue="group", dodge=True,
                  palette={"Normal": GRAY, "MASH F2-F4": ORANGE}, size=4, alpha=0.85, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, borderaxespad=0)
    ax.set_xlabel("")
    ax.set_ylabel("Standardized-mean score")
    clean_axes(ax)

    effects = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_sensitivity_effects.csv")
    effects = effects[effects["program_id"].isin(target)].copy()
    effects["program"] = effects["program_id"].map(PROGRAM_LABELS)
    effects.to_csv(source / "figure_4b_mash_macrophage_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    for i, program in enumerate(["SAM-A", "SAM-B", "TMo-E"]):
        values = effects[effects["program"].eq(program)]
        for j, method in enumerate(["singscore", "standardized_mean"]):
            row = values[values["score_method"].eq(method)].iloc[0]
            y = i + (-0.12 if j == 0 else 0.12)
            ax.hlines(y, row.robust_ci95_low, row.robust_ci95_high, color=METHOD_COLORS[method], lw=1.5)
            ax.scatter(row.hedges_g, y, color=METHOD_COLORS[method], s=42, zorder=3)
    ax.axvline(0, color=BLACK, lw=0.8)
    ax.set_yticks(np.arange(3), ["SAM-A", "SAM-B", "TMo-E"])
    ax.set_xlabel("Hedges g with HC3 95% interval")
    clean_axes(ax)

    random = pd.read_csv(repo / "results" / "random_controls" / "gse244832_random_module_benchmark.csv")
    random = random[random["program_id"].isin(target)].copy()
    random["program"] = random["program_id"].map(PROGRAM_LABELS)
    random.to_csv(source / "figure_4c_mash_random_percentiles.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    for method, marker in [("singscore", "o"), ("standardized_mean", "s")]:
        values = random[random["score_method"].eq(method)]
        ax.scatter(values["program"], values["real_effect_percentile"], marker=marker, color=METHOD_COLORS[method], s=55,
                   edgecolor="white", linewidth=0.4, label=method)
    ax.axhline(0.95, color=BLACK, linestyle="--", lw=0.8)
    ax.set_ylim(0.75, 1.01)
    ax.set_ylabel("Observed-effect percentile")
    ax.set_xlabel("")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, borderaxespad=0)
    clean_axes(ax)

    card = pd.read_csv(repo / "results" / "deep_benchmark" / "program_transportability_report_card.csv").head(8).copy()
    card["program"] = card["program_id"].map(PROGRAM_LABELS)
    card.to_csv(source / "figure_4d_mash_evidence_gates.csv", index=False)
    domains = [
        ("measurement_domain_0_20", "Measurement", BLUE),
        ("score_method_domain_0_20", "Method", GREEN),
        ("directional_transfer_domain_0_20", "Direction", ORANGE),
        ("matched_random_specificity_domain_0_20", "Random specificity", PURPLE),
        ("endpoint_evidence_domain_0_20", "Endpoint", RED),
    ]
    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    y = np.arange(len(card))
    left = np.zeros(len(card))
    for column, label, color in domains:
        ax.barh(y, card[column], left=left, color=color, label=label, height=0.62)
        left += card[column].to_numpy()
    ax.set_yticks(y, card["program"])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Exploratory readiness total (0-100)")
    # Keep the five-domain key outside the plotting region so it cannot cover
    # the leading stacked bar. A single row also preserves the full vertical
    # area for the eight ranked programs.
    ax.legend(
        frameon=False,
        ncol=5,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        borderaxespad=0,
        columnspacing=1.0,
        handletextpad=0.35,
        fontsize=8,
    )
    clean_axes(ax)
    save(fig, output, "figure_4_mash_macrophage_specificity")


def figure_5(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 7.9), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    watson = pd.read_csv(repo / "results" / "deep_benchmark" / "watson_identity_control_donor_audit.csv")
    watson.to_csv(source / "figure_5a_gse290642_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    values = watson[watson["score_method"].eq("standardized_mean")]
    sns.stripplot(data=values, x="disease_group", y="lineage_top_score_accuracy", hue="disease_group",
                  palette={"healthy": GREEN, "fibrosis": RED}, size=8, jitter=0.08, legend=False, ax=ax)
    ax.plot([0, 1], [values[values.disease_group.eq("healthy")]["lineage_top_score_accuracy"].mean(),
                     values[values.disease_group.eq("fibrosis")]["lineage_top_score_accuracy"].mean()], color=GRAY, lw=1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("")
    ax.set_ylabel("Three-lineage top-score accuracy")
    clean_axes(ax)

    contrasts = pd.read_csv(repo / "results" / "deep_benchmark" / "watson_identity_control_contrasts.csv")
    contrasts.to_csv(source / "figure_5b_watson_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    standardized = watson[watson["score_method"].eq("standardized_mean")]
    x = np.arange(len(standardized))
    ax.bar(x, standardized["mean_matched_minus_best_off_margin"],
           color=standardized["disease_group"].map({"healthy": GREEN, "fibrosis": RED}))
    ax.axhline(0, color=BLACK, lw=0.8)
    ax.set_xticks(x, standardized["donor_id"], rotation=35, ha="right")
    ax.set_ylabel("Matched identity minus best off-lineage")
    clean_axes(ax)
    ax.text(0.02, 0.94, "Exact 3+3 P=0.10; context is confounded", transform=ax.transAxes,
            fontsize=7.1, color=GRAY, bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 1.5})

    discovery = pd.read_csv(repo / "results" / "deep_benchmark" / "core_discovery_gene_effects.csv")
    cores = pd.read_csv(repo / "results" / "deep_benchmark" / "minimal_core_membership.csv")
    discovery.to_csv(source / "figure_5c_all_random_percentiles.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    eligible = discovery[discovery["passes_measurement"].astype(bool)]
    ax.scatter(eligible["gse202379_stage_spearman_rho"], eligible["gse244832_mash_vs_normal_hedges_g"],
               c=eligible["lineage"].map(LINEAGE_COLORS), s=24, alpha=0.6, edgecolor="white", linewidth=0.3)
    ax.axvline(0.10, color=GRAY, linestyle="--", lw=0.8)
    ax.axhline(0.20, color=GRAY, linestyle="--", lw=0.8)
    label_offsets = {
        "TFF3": (3, 7), "PLPP1": (3, -10), "FTL": (3, 7), "CPE": (3, -10), "FTH1": (3, 7),
        "MDK": (3, 7), "CST3": (3, -10), "TM4SF1": (3, 7), "SERPINF1": (3, -10), "TMSB10": (3, 7),
    }
    for row in cores.itertuples(index=False):
        ax.annotate(
            row.gene_symbol,
            (row.gse202379_stage_spearman_rho, row.gse244832_mash_vs_normal_hedges_g),
            xytext=label_offsets.get(row.gene_symbol, (3, 5)),
            textcoords="offset points",
            fontsize=6.5,
            fontweight="bold",
        )
    ax.set_xlabel("GSE202379 gene-stage Spearman rho")
    ax.set_ylabel("GSE244832 gene-level Hedges g")
    clean_axes(ax)

    validation = pd.read_csv(repo / "results" / "deep_benchmark" / "minimal_core_random_benchmark.csv")
    validation["dataset"] = validation["dataset_id"].replace({
        "GSE210077_Watson6": "Watson6", "GSE290642_human": "GSE290642", "GSE181483_human": "GSE181483"
    })
    validation["core"] = validation["lineage"].map({"endothelial": "Endothelial Core5", "mesenchymal_hsc_myofibroblast": "Mesenchymal Core5"})
    validation.to_csv(source / "figure_5d_comparability_matrix.csv", index=False)
    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    validation["x"] = validation["dataset"] + "\n" + validation["core"]
    order = list(dict.fromkeys(validation["x"]))
    for method, marker in [("singscore", "o"), ("standardized_mean", "s")]:
        values = validation[validation["score_method"].eq(method)].copy()
        positions = [order.index(value) for value in values["x"]]
        ax.scatter(positions, values["real_effect_percentile"], marker=marker, color=METHOD_COLORS[method], s=52,
                   edgecolor="white", linewidth=0.4, label=method)
    ax.axhline(0.95, color=BLACK, linestyle="--", lw=0.8)
    ax.set_xticks(np.arange(len(order)), order, rotation=35, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Held-out percentile among 1,000 random Core5 modules")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2, borderaxespad=0)
    clean_axes(ax)
    save(fig, output, "figure_5_cohort_assay_stress_tests")


def figure_6(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.2), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[1.2, 0.8])
    card = pd.read_csv(repo / "results" / "deep_benchmark" / "program_transportability_report_card.csv").copy()
    card["program"] = card["program_id"].map(PROGRAM_LABELS)
    card.to_csv(source / "figure_6a_advanced_meta.csv", index=False)
    domains = [
        ("measurement_domain_0_20", "Measurement", BLUE),
        ("score_method_domain_0_20", "Method", GREEN),
        ("directional_transfer_domain_0_20", "Direction", ORANGE),
        ("matched_random_specificity_domain_0_20", "Random specificity", PURPLE),
        ("endpoint_evidence_domain_0_20", "Endpoint", RED),
    ]
    ax = fig.add_subplot(grid[0, 0])
    panel_label(ax, "A")
    shown = card.sort_values("transportability_readiness_total_0_100", ascending=True)
    y = np.arange(len(shown))
    left = np.zeros(len(shown))
    for column, label, color in domains:
        ax.barh(y, shown[column], left=left, color=color, label=label, height=0.72)
        left += shown[column].to_numpy()
    ax.set_yticks(y, shown["program"])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Exploratory readiness total (0-100)")
    ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.10), borderaxespad=0, fontsize=7.2)
    clean_axes(ax)

    matrix = card.set_index("program")[[column for column, _, _ in domains]].copy()
    matrix.columns = [label for _, label, _ in domains]
    matrix = matrix.loc[card.sort_values("transportability_readiness_total_0_100", ascending=False)["program"]]
    matrix.to_csv(source / "figure_6b_evidence_matrix.csv")
    ax = fig.add_subplot(grid[0, 1])
    panel_label(ax, "B")
    sns.heatmap(matrix, cmap="YlGnBu", vmin=0, vmax=20, linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Domain score (0-20)", "shrink": 0.75}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=35)

    evidence = pd.read_csv(repo / "results" / "exploratory" / "program_evidence_attrition_matrix.csv")
    gates = [
        ("Evaluable", "evaluable_independent"),
        ("Positive somewhere", "positive_both_scores_any_cohort"),
        ("Above random somewhere", "random_specific_both_scores_any_cohort"),
        ("Positive primary CI", "positive_primary_interval_both_scores"),
        ("Replicated", "within_cell_state_replicated"),
        ("Pan-cirrhotic", "pan_cirrhotic_transportable"),
        ("Assay robust", "assay_robust"),
    ]
    gate_data = pd.DataFrame({"gate": [x[0] for x in gates], "programs": [int(evidence[x[1]].astype(bool).sum()) for x in gates]})
    gate_data.to_csv(source / "figure_6c_evidence_counts.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel_label(ax, "C")
    colors = [GREEN, GREEN, ORANGE, RED, RED, RED, RED]
    ax.bar(np.arange(len(gate_data)), gate_data["programs"], color=colors)
    ax.set_xticks(np.arange(len(gate_data)), gate_data["gate"], rotation=35, ha="right")
    ax.set_ylim(0, 20)
    ax.set_ylabel("Programs (n=19)")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel_label(ax, "D")
    ax.axis("off")
    conclusions = [
        ("Technical sensitivity", "Identity controls pass in GSE202379", GREEN),
        ("Conditional positive", "Matched advanced endpoint: rho 0.69; 91.7% signs", BLUE),
        ("Remaining failure", "Held-out effect magnitude R²≈0; Core5 not random-specific", ORANGE),
        ("Translation boundary", "Not a staging biomarker, mechanism, or treatment target", RED),
    ]
    for i, (head, body, color) in enumerate(conclusions):
        y = 0.83 - i * 0.22
        ax.add_patch(patches.FancyBboxPatch((0.03, y - 0.08), 0.94, 0.15, boxstyle="round,pad=0.01",
                                           facecolor="white", edgecolor=color, linewidth=1.3, transform=ax.transAxes))
        ax.text(0.06, y, head, color=color, fontweight="bold", fontsize=7.7,
                transform=ax.transAxes, va="center")
        ax.text(0.40, y, textwrap.fill(body, 39), color=BLACK, fontsize=7.5,
                transform=ax.transAxes, va="center")
    pd.DataFrame(conclusions, columns=["level", "conclusion", "color"]).drop(columns="color").to_csv(
        source / "figure_6d_conclusion_boundary.csv", index=False
    )
    save(fig, output, "figure_6_evidence_synthesis")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "figures"
    source = repo / "results" / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    figure_1(repo, output, source)
    figure_2(repo, output, source)
    figure_3(repo, output, source)
    figure_4(repo, output, source)
    figure_5(repo, output, source)
    figure_6(repo, output, source)
    print("Deep benchmark Figures 1-6 written to", output)


if __name__ == "__main__":
    main()
