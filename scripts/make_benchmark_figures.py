from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd
import seaborn as sns


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.2,
        "axes.titlesize": 9.2,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.2,
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
GRAY = "#7A7A7A"
LIGHT_GRAY = "#E8E8E8"
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
PROGRAM_ORDER = list(PROGRAM_LABELS)

REPRESENTATIVE_ENDPOINTS = [
    ("GSE202379", "advanced_f3f4_vs_f0_non_end_stage", "GSE202379\nF3-F4"),
    ("GSE290642_human", "f4_vs_f0_reconstructed_label_sensitivity", "GSE290642\nF4 sens."),
    ("GSE244832", "mash_f2f4_group_vs_normal_sensitivity", "GSE244832\nMASH sens."),
    ("GSE210077_Watson6", "mixed_f2f4_fibrosis_vs_healthy_sensitivity", "Watson6\nmixed sens."),
    ("GSE181483_human", "cirrhosis_vs_healthy_directional", "GSE181483\n2+2 dir."),
]


def save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=360, bbox_inches="tight", pad_inches=0.12, facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.12, facecolor="white")
    plt.close(fig)


def panel(ax: plt.Axes, letter: str, title: str) -> None:
    ax.text(-0.11, 1.07, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")
    ax.set_title(title, loc="left", pad=6, fontweight="bold")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=3)


def read_random(repo: Path) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted((repo / "results" / "random_controls").glob("*_random_module_benchmark.csv"))]
    data = pd.concat(frames, ignore_index=True, sort=False)
    data["above_random_95th_percentile"] = data["above_random_95th_percentile"].astype(bool)
    return data


def representative_effects(effects: pd.DataFrame, score_method: str = "singscore") -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    for dataset, contrast, label in REPRESENTATIVE_ENDPOINTS:
        subset = effects[
            effects["dataset_id"].eq(dataset)
            & effects["contrast"].eq(contrast)
            & effects["score_method"].eq(score_method)
        ].copy()
        subset["cohort_endpoint"] = label
        selected.append(subset)
    return pd.concat(selected, ignore_index=True)


def figure_1(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.6, 7.6), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[0.82, 1.18])

    ax = fig.add_subplot(grid[0, 0])
    panel(ax, "A", "Prespecified donor-level transportability benchmark")
    ax.axis("off")
    labels = [
        ("19 published\nprograms frozen", BLUE),
        ("Donor x lineage\npseudobulks", GREEN),
        ("Cohort-specific\neffects + controls", ORANGE),
        ("Replication and\ntransportability gates", RED),
    ]
    xs = np.linspace(0.03, 0.78, len(labels))
    for index, ((text, color), x) in enumerate(zip(labels, xs)):
        rect = patches.FancyBboxPatch(
            (x, 0.43), 0.19, 0.28, boxstyle="round,pad=0.012,rounding_size=0.012",
            facecolor="white", edgecolor=color, linewidth=1.5, transform=ax.transAxes,
        )
        ax.add_patch(rect)
        ax.text(x + 0.095, 0.57, text, ha="center", va="center", transform=ax.transAxes, fontsize=8.2)
        if index < len(labels) - 1:
            ax.annotate("", xy=(x + 0.235, 0.57), xytext=(x + 0.195, 0.57), xycoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.1))
    ax.text(0.02, 0.18, "No cell-level disease test  |  no cross-study expression integration  |  no missing-gene zero filling",
            transform=ax.transAxes, fontsize=7.6, color=GRAY)
    workflow = pd.DataFrame(
        {"step": [1, 2, 3, 4], "label": [x[0].replace("\n", " ") for x in labels], "boundary": [
            "program membership fixed before validation", "independent donor is inferential unit",
            "two scores, robust intervals, permutations, matched random modules", "endpoints and evidence tiers preserved",
        ]}
    )
    workflow.to_csv(source_dir / "figure_1a_workflow.csv", index=False)

    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "Five independent resources contribute unequal evidence")
    cohort = pd.DataFrame(
        [
            ["GSE202379", 47, "snRNA", "Author", "Formal primary"],
            ["GSE290642", 24, "scRNA", "Reconstructed", "Sensitivity"],
            ["GSE244832", 18, "snRNA", "Author clusters", "Sensitivity"],
            ["Watson6", 6, "snRNA", "Author", "Sensitivity"],
            ["GSE181483", 4, "scRNA", "Reconstructed", "Directional"],
        ], columns=["dataset", "donors", "assay", "annotation", "highest_evidence_role"]
    )
    cohort.to_csv(source_dir / "figure_1b_cohort_summary.csv", index=False)
    y = np.arange(len(cohort))
    colors = cohort["highest_evidence_role"].map({"Formal primary": RED, "Sensitivity": ORANGE, "Directional": BLUE})
    ax.barh(y, cohort["donors"], color=colors, height=0.58)
    for yi, row in cohort.iterrows():
        ax.text(row["donors"] + 0.8, yi, f"{row['assay']} | {row['annotation']}", va="center", fontsize=7.3)
    ax.set_yticks(y, cohort["dataset"])
    ax.invert_yaxis()
    ax.set_xlabel("Independent human donors")
    ax.set_xlim(0, 62)
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 0])
    panel(ax, "C", "Donor-lineage eligibility collapses after the 30-cell gate")
    eligibility_specs = [
        ("GSE202379", "metadata/gse202379_donor_lineage_eligibility.csv", "eligible_primary_30"),
        ("GSE290642", "metadata/gse290642_donor_lineage_eligibility.csv", "eligible_30"),
        ("GSE244832", "metadata/gse244832_donor_lineage_eligibility.csv", "eligible_30"),
        ("Watson6", "metadata/gse210077_watson6_donor_lineage_eligibility.csv", "eligible_30"),
        ("GSE181483", "metadata/gse181483_donor_lineage_eligibility.csv", "eligible_30"),
    ]
    gate_rows = []
    for dataset, relative, gate in eligibility_specs:
        frame = pd.read_csv(repo / relative)
        for lineage, values in frame.groupby("harmonized_lineage"):
            gate_rows.append({"dataset": dataset, "lineage": lineage, "eligible_donors": int(values[gate].astype(bool).sum())})
    gates = pd.DataFrame(gate_rows)
    gates.to_csv(source_dir / "figure_1c_donor_gate_counts.csv", index=False)
    pivot = gates.pivot(index="dataset", columns="lineage", values="eligible_donors").reindex(cohort["dataset"])
    pivot = pivot.reindex(columns=list(LINEAGE_LABELS))
    bottom = np.zeros(len(pivot))
    for lineage in pivot.columns:
        values = pivot[lineage].fillna(0).to_numpy()
        ax.bar(np.arange(len(pivot)), values, bottom=bottom, label=LINEAGE_LABELS[lineage], color=LINEAGE_COLORS[lineage])
        bottom += values
    ax.set_xticks(np.arange(len(pivot)), pivot.index, rotation=25, ha="right")
    ax.set_ylabel("Eligible donor-lineage units")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel(ax, "D", "Feature coverage is high for RNA assays but not for the spatial panel")
    coverage_files = [
        ("GSE202379", "gse202379_program_coverage.csv"),
        ("GSE290642", "gse290642_program_coverage.csv"),
        ("GSE244832", "gse244832_program_coverage.csv"),
        ("Watson6", "gse210077_watson6_program_coverage.csv"),
        ("GSE181483", "gse181483_program_coverage.csv"),
    ]
    coverage_rows = []
    for dataset, filename in coverage_files:
        frame = pd.read_csv(repo / "results" / "qc" / filename)
        for tier, count in frame["coverage_tier"].value_counts().items():
            coverage_rows.append({"dataset": dataset, "coverage_tier": tier, "programs": int(count)})
    coverage_rows.extend([
        {"dataset": "MERFISH panel", "coverage_tier": "primary", "programs": 0},
        {"dataset": "MERFISH panel", "coverage_tier": "not_evaluated", "programs": 19},
    ])
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(source_dir / "figure_1d_coverage_counts.csv", index=False)
    order = ["GSE202379", "GSE290642", "GSE244832", "Watson6", "GSE181483", "MERFISH panel"]
    tier_order = ["primary", "flagged", "sensitivity", "not_evaluated"]
    tier_colors = {"primary": GREEN, "flagged": ORANGE, "sensitivity": BLUE, "not_evaluated": LIGHT_GRAY}
    cov_pivot = coverage.pivot_table(index="dataset", columns="coverage_tier", values="programs", aggfunc="sum", fill_value=0).reindex(order)
    bottom = np.zeros(len(cov_pivot))
    for tier in tier_order:
        values = cov_pivot[tier].to_numpy() if tier in cov_pivot else np.zeros(len(cov_pivot))
        ax.bar(np.arange(len(cov_pivot)), values, bottom=bottom, color=tier_colors[tier], label=tier.replace("_", " ").title())
        bottom += values
    ax.set_xticks(np.arange(len(cov_pivot)), cov_pivot.index, rotation=25, ha="right")
    ax.set_ylabel("Frozen programs (n=19)")
    ax.set_ylim(0, 20)
    ax.legend(frameon=False, ncol=2, loc="upper left")
    clean_axes(ax)
    save(fig, output_dir, "figure_1_study_design_and_evidence")


def figure_2(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.6, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[0.85, 1.15])
    scores = pd.read_csv(repo / "results" / "primary" / "gse202379_donor_program_scores.csv.gz")
    base = pd.read_csv(repo / "metadata" / "gse202379_donor_lineage_eligibility.csv")
    base = base[base["eligible_primary_30"].astype(bool)].copy()
    base["stage_label"] = "F" + base["Fibrosis.score..F0.4."].astype(int).astype(str)
    base.loc[base["Disease.status"].eq("end stage"), "stage_label"] = "End-stage F4"
    counts = base.groupby(["stage_label", "harmonized_lineage"], as_index=False)["canonical_donor_id"].nunique().rename(columns={"canonical_donor_id": "donors"})
    counts.to_csv(source_dir / "figure_2a_stage_donor_counts.csv", index=False)

    ax = fig.add_subplot(grid[0, 0])
    panel(ax, "A", "Eligible donors span F0-F4, but lineage depth is uneven")
    stage_order = ["F0", "F1", "F2", "F3", "F4", "End-stage F4"]
    count_pivot = counts.pivot(index="stage_label", columns="harmonized_lineage", values="donors").fillna(0).reindex(stage_order)
    x = np.arange(len(stage_order))
    width = 0.25
    for index, lineage in enumerate(LINEAGE_LABELS):
        ax.bar(x + (index - 1) * width, count_pivot.get(lineage, 0), width=width, color=LINEAGE_COLORS[lineage], label=LINEAGE_LABELS[lineage])
    ax.set_xticks(x, stage_order, rotation=20, ha="right")
    ax.set_ylabel("Donors passing 30-cell gate")
    ax.legend(frameon=False, ncol=1, loc="upper right")
    clean_axes(ax)

    primary = pd.read_csv(repo / "results" / "primary" / "gse202379_primary_effects.csv")
    primary["display"] = np.where(
        primary["contrast"].eq("clinical_cirrhosis_vs_healthy"),
        "Clinical\nEndothelial",
        "F3-F4\n" + primary["lineage"].map(LINEAGE_LABELS),
    )
    primary.to_csv(source_dir / "figure_2b_primary_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "All 34 formal-primary robust intervals include zero")
    primary["row"] = primary["program_id"].map(PROGRAM_LABELS)
    columns = []
    for label in ["Clinical\nEndothelial", "F3-F4\nEndothelial", "F3-F4\nMesenchymal/HSC"]:
        for method in ["singscore", "standardized_mean"]:
            columns.append((label, method, f"{label}\n{'Rank' if method == 'singscore' else 'Z-mean'}"))
    primary_programs = list(dict.fromkeys(primary["program_id"]))
    row_order = [p for p in PROGRAM_ORDER if p in primary_programs]
    matrix = pd.DataFrame(index=[PROGRAM_LABELS[p] for p in row_order], columns=[c[2] for c in columns], dtype=float)
    for label, method, display in columns:
        subset = primary[primary["display"].eq(label) & primary["score_method"].eq(method)].set_index("row")
        matrix.loc[subset.index, display] = subset["hedges_g"]
    sns.heatmap(matrix.clip(-2, 2), cmap="vlag", center=0, vmin=-2, vmax=2, mask=matrix.isna(), linewidths=0.3,
                linecolor="white", cbar_kws={"label": "Hedges g", "shrink": 0.75}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=35)

    random = pd.read_csv(repo / "results" / "random_controls" / "gse202379_random_module_benchmark.csv")
    direction = pd.read_csv(repo / "results" / "random_controls" / "gse202379_direction_randomization_benchmark.csv")
    merged = primary.merge(
        random[["dataset_id", "contrast", "lineage", "program_id", "score_method", "real_effect_percentile", "above_random_95th_percentile"]],
        on=["dataset_id", "contrast", "lineage", "program_id", "score_method"], how="left",
    ).merge(
        direction[["dataset_id", "contrast", "lineage", "program_id", "score_method", "above_direction_randomized_95th"]],
        on=["dataset_id", "contrast", "lineage", "program_id", "score_method"], how="left",
    )
    merged.to_csv(source_dir / "figure_2c_primary_specificity.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel(ax, "C", "Program specificity does not overcome donor-level imprecision")
    for method, marker in [("singscore", "o"), ("standardized_mean", "s")]:
        subset = merged[merged["score_method"].eq(method)]
        ax.scatter(subset["hedges_g"], subset["real_effect_percentile"], marker=marker, s=32,
                   c=subset["lineage"].map(LINEAGE_COLORS), edgecolor="white", linewidth=0.4, alpha=0.9, label=method)
    ax.axhline(0.95, color=BLACK, linestyle="--", lw=0.9)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.set_xlabel("Hedges g")
    ax.set_ylabel("Percentile among 1,000 matched modules")
    ax.set_ylim(-0.02, 1.03)
    clean_axes(ax)
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=BLACK, label="singscore"),
        plt.Line2D([0], [0], marker="s", linestyle="", color=BLACK, label="standardized mean"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    ax.text(0.02, 0.03, "No point also had a positive HC3 interval", transform=ax.transAxes, color=GRAY, fontsize=7.5)

    trends = pd.read_csv(repo / "results" / "exploratory" / "gse202379_stage_trends.csv")
    trends = trends[trends["analysis_set"].eq("non_end_stage_primary_exploratory")].copy()
    trends.to_csv(source_dir / "figure_2d_stage_trends.csv", index=False)
    ax = fig.add_subplot(grid[1, 1])
    panel(ax, "D", "No program tracks F0-F4 stage with both scores after FDR control")
    for method, marker in [("singscore", "o"), ("standardized_mean", "s")]:
        subset = trends[trends["score_method"].eq(method)]
        ax.scatter(subset["spearman_rho"], -np.log10(subset["fdr_bh"].clip(lower=1e-6)), marker=marker, s=34,
                   c=subset["lineage"].map(LINEAGE_COLORS), edgecolor="white", linewidth=0.4, alpha=0.9)
    ax.axhline(-np.log10(0.05), color=BLACK, linestyle="--", lw=0.9)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.set_xlabel("Spearman rho with fibrosis stage")
    ax.set_ylabel("-log10(BH FDR)")
    clean_axes(ax)
    top = trends.nsmallest(3, "fdr_bh")
    for _, row in top.iterrows():
        ax.annotate(PROGRAM_LABELS[row["program_id"]], (row["spearman_rho"], -np.log10(row["fdr_bh"])),
                    xytext=(3, 4), textcoords="offset points", fontsize=6.8)
    ax.text(0.98, 0.03, "0/18 evaluable programs dual-score FDR-positive", transform=ax.transAxes, ha="right", color=RED, fontweight="bold")
    save(fig, output_dir, "figure_2_primary_and_stage")


def figure_3(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.4), constrained_layout=True)
    grid = GridSpec(2, 3, figure=fig, height_ratios=[0.72, 1.28])
    within = pd.read_csv(repo / "results" / "exploratory" / "within_cohort_score_method_concordance.csv")
    across = pd.read_csv(repo / "results" / "exploratory" / "cross_cohort_program_concordance.csv")
    within.to_csv(source_dir / "figure_3a_within_cohort_concordance.csv", index=False)
    across.to_csv(source_dir / "figure_3b_cross_cohort_concordance.csv", index=False)

    ax = fig.add_subplot(grid[0, 0])
    panel(ax, "A", "Within-cohort score agreement")
    display = within.copy()
    display["label"] = display["dataset_id"].replace({"GSE210077_Watson6": "Watson6", "GSE290642_human": "GSE290642", "GSE181483_human": "GSE181483"}) + " | " + display["lineage"].map(LINEAGE_LABELS)
    display = display.sort_values("spearman_rho")
    y = np.arange(len(display))
    ax.scatter(display["spearman_rho"], y, c=display["lineage"].map(LINEAGE_COLORS), s=34)
    ax.set_yticks(y, display["label"])
    ax.set_xlim(-1.02, 1.02)
    ax.axvline(0, color=GRAY, lw=0.7)
    ax.axvline(display["spearman_rho"].median(), color=BLACK, linestyle="--", lw=0.9)
    ax.set_xlabel("Spearman rho: singscore vs standardized mean")
    clean_axes(ax)
    ax.text(0.02, 0.98, f"Median rho = {display['spearman_rho'].median():.2f}", transform=ax.transAxes, va="top", fontweight="bold")

    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "Cross-cohort program agreement")
    order = list(LINEAGE_LABELS)
    sns.boxplot(data=across, x="lineage", y="spearman_rho", order=order, color=VERY_LIGHT, width=0.58, fliersize=0, ax=ax)
    sns.stripplot(data=across, x="lineage", y="spearman_rho", order=order, hue="score_method", palette=METHOD_COLORS,
                  dodge=True, alpha=0.8, size=4, ax=ax)
    ax.axhline(0, color=GRAY, lw=0.7)
    ax.set_xticks(range(3), ["Endothelial", "Macrophage", "Mesenchymal"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Cross-cohort Spearman rho")
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    clean_axes(ax)
    ax.text(0.02, 0.98, f"Median rho = {across['spearman_rho'].median():.2f}", transform=ax.transAxes, va="top", fontweight="bold")

    ax = fig.add_subplot(grid[0, 2])
    panel(ax, "C", "Program-direction agreement")
    agreement = pd.concat([
        within.assign(comparison="Within cohort")[['comparison', 'sign_agreement']],
        across.assign(comparison="Across cohorts")[['comparison', 'sign_agreement']],
    ], ignore_index=True)
    agreement.to_csv(source_dir / "figure_3c_sign_agreement.csv", index=False)
    sns.boxplot(data=agreement, x="comparison", y="sign_agreement", order=["Within cohort", "Across cohorts"],
                palette=[GREEN, GRAY], width=0.52, fliersize=0, ax=ax)
    sns.stripplot(data=agreement, x="comparison", y="sign_agreement", order=["Within cohort", "Across cohorts"],
                  color=BLACK, alpha=0.5, size=3, ax=ax)
    ax.axhline(0.5, color=BLACK, linestyle="--", lw=0.9)
    ax.set_xlabel("")
    ax.set_ylabel("Program-direction agreement")
    ax.set_ylim(-0.03, 1.05)
    clean_axes(ax)
    medians = agreement.groupby("comparison")["sign_agreement"].median()
    ax.text(0.02, 0.98, f"Median: {medians['Within cohort']:.2f} vs {medians['Across cohorts']:.2f}",
            transform=ax.transAxes, va="top", fontweight="bold")

    effects = pd.read_csv(repo / "results" / "meta" / "cross_cohort_effect_matrix.csv")
    selected = representative_effects(effects, "singscore")
    selected.to_csv(source_dir / "figure_3d_representative_effects.csv", index=False)
    ax = fig.add_subplot(grid[1, :])
    panel(ax, "D", "Cohort context changes program magnitude and direction")
    columns = [item[2] for item in REPRESENTATIVE_ENDPOINTS]
    matrix = selected.pivot_table(index="program_id", columns="cohort_endpoint", values="hedges_g", aggfunc="first").reindex(index=PROGRAM_ORDER, columns=columns)
    annotations = matrix.map(lambda value: "" if pd.isna(value) else f"{value:.1f}")
    sns.heatmap(matrix.clip(-2, 2), cmap="vlag", center=0, vmin=-2, vmax=2, mask=matrix.isna(), annot=annotations,
                fmt="", linewidths=0.35, linecolor="white", cbar_kws={"label": "Hedges g (color clipped to +/-2)", "shrink": 0.72}, ax=ax)
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0)
    ax.axhline(7, color=BLACK, linewidth=1.1)
    ax.axhline(13, color=BLACK, linewidth=1.1)
    save(fig, output_dir, "figure_3_reproducibility_vs_transportability")


def figure_4(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.6, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[0.9, 1.1])
    target_programs = ["RAM2019_MAC_SIG_A_SAM", "RAM2019_MAC_SIG_B_SAM", "RAM2019_MAC_SIG_E_TMO"]
    donor = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_donor_program_scores.csv.gz")
    donor = donor[
        donor["program_id"].isin(target_programs)
        & donor["score_method"].eq("singscore")
        & donor["condition"].isin(["NORMAL", "NASH"])
        & donor["eligible_30"].astype(bool)
    ].copy()
    donor["program"] = donor["program_id"].map(PROGRAM_LABELS)
    donor["condition"] = donor["condition"].map({"NORMAL": "Normal", "NASH": "MASH"})
    donor.to_csv(source_dir / "figure_4a_mash_macrophage_donor_scores.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel(ax, "A", "Three macrophage programs rise in mixed-stage MASH")
    sns.boxplot(data=donor, x="program", y="score", hue="condition", order=["SAM-A", "SAM-B", "TMo-E"],
                hue_order=["Normal", "MASH"], palette={"Normal": LIGHT_GRAY, "MASH": ORANGE}, width=0.65, fliersize=0, ax=ax)
    sns.stripplot(data=donor, x="program", y="score", hue="condition", order=["SAM-A", "SAM-B", "TMo-E"],
                  hue_order=["Normal", "MASH"], dodge=True, palette={"Normal": GRAY, "MASH": "#8A4F00"}, size=3.5, alpha=0.8, ax=ax, legend=False)
    ax.set_xlabel("")
    ax.set_ylabel("Donor-level singscore")
    clean_axes(ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, loc="upper left")

    effects = pd.read_csv(repo / "results" / "sensitivity" / "gse244832_sensitivity_effects.csv")
    mac = effects[effects["lineage"].eq("macrophage_monocyte")].copy()
    mac.to_csv(source_dir / "figure_4b_mash_macrophage_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "Effect estimates are program-selective, not uniformly positive")
    programs = [p for p in PROGRAM_ORDER if p.startswith("RAM2019_MAC")]
    y = np.arange(len(programs))
    for offset, method in [(-0.09, "singscore"), (0.09, "standardized_mean")]:
        subset = mac[mac["score_method"].eq(method)].set_index("program_id").reindex(programs)
        ax.errorbar(subset["hedges_g"], y + offset,
                    xerr=[subset["hedges_g"] - subset["robust_ci95_low"], subset["robust_ci95_high"] - subset["hedges_g"]],
                    fmt="o" if method == "singscore" else "s", color=METHOD_COLORS[method], markersize=4, capsize=2,
                    label="singscore" if method == "singscore" else "standardized mean")
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_yticks(y, [PROGRAM_LABELS[p] for p in programs])
    ax.invert_yaxis()
    ax.set_xlabel("Hedges g (HC3 95% CI)")
    ax.legend(frameon=False, loc="lower right")
    clean_axes(ax)

    random = pd.read_csv(repo / "results" / "random_controls" / "gse244832_random_module_benchmark.csv")
    random.to_csv(source_dir / "figure_4c_mash_random_percentiles.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel(ax, "C", "Matched random modules isolate the selective signals")
    pct = random.pivot_table(index="program_id", columns="score_method", values="real_effect_percentile", aggfunc="first").reindex(PROGRAM_ORDER)
    annotations = pct.map(lambda value: "" if pd.isna(value) else f"{value:.2f}")
    sns.heatmap(pct, cmap=sns.light_palette(ORANGE, as_cmap=True), vmin=0.5, vmax=1.0, annot=annotations, fmt="",
                linewidths=0.35, linecolor="white", cbar_kws={"label": "Percentile", "shrink": 0.72}, ax=ax)
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xticklabels(["singscore", "standardized mean"], rotation=20, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    for row_index, program_id in enumerate(PROGRAM_ORDER):
        for col_index, method in enumerate(["singscore", "standardized_mean"]):
            if program_id in pct.index and method in pct.columns and pd.notna(pct.loc[program_id, method]) and pct.loc[program_id, method] > 0.95:
                ax.add_patch(patches.Rectangle((col_index, row_index), 1, 1, fill=False, edgecolor=BLACK, linewidth=1.1))

    ax = fig.add_subplot(grid[1, 1])
    panel(ax, "D", "Evidence narrows from broad positivity to three dual-score macrophage hypotheses")
    positive = int((effects["hedges_g"] > 0).sum())
    robust = int((effects["robust_ci95_low"] > 0).sum())
    above_random = int(random["above_random_95th_percentile"].astype(bool).sum())
    dual = 3
    gates = pd.DataFrame({
        "criterion": ["Positive effect", "Positive HC3 interval", "Above matched-random 95th", "Macrophage programs passing both scores"],
        "count": [positive, robust, above_random, dual],
        "denominator": [38, 38, 38, 6],
    })
    gates.to_csv(source_dir / "figure_4d_mash_evidence_gates.csv", index=False)
    bars = ax.barh(np.arange(4), gates["count"] / gates["denominator"], color=[LIGHT_GRAY, BLUE, ORANGE, RED])
    ax.set_yticks(np.arange(4), ["Positive", "Robust interval", "Random-specific", "Dual-score macrophage"])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of eligible rows or programs")
    for bar, count, denominator in zip(bars, gates["count"], gates["denominator"]):
        ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2, f"{count}/{denominator}", va="center", fontsize=8)
    clean_axes(ax)
    ax.text(0.02, -0.22, "Endpoint boundary: MASH donors span F2-F4 without donor-specific stage; this is not cirrhosis replication.",
            transform=ax.transAxes, color=RED, fontsize=7.4)
    save(fig, output_dir, "figure_4_mash_macrophage_specificity")


def figure_5(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.6, 8.0), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    gse290 = pd.read_csv(repo / "results" / "sensitivity" / "gse290642_sensitivity_effects.csv")
    gse290.to_csv(source_dir / "figure_5a_gse290642_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 0])
    panel(ax, "A", "Reconstructed-label GSE290642 does not yield random-specific transfer")
    gse290["column"] = gse290["contrast"].map({
        "f4_vs_f0_reconstructed_label_sensitivity": "F4 vs F0",
        "all_fibrosis_vs_f0_reconstructed_label_sensitivity": "F1-F4 vs F0",
    }) + "\n" + gse290["score_method"].map({"singscore": "Rank", "standardized_mean": "Z-mean"})
    matrix = gse290.pivot_table(index="program_id", columns="column", values="hedges_g", aggfunc="first").reindex([p for p in PROGRAM_ORDER if p.startswith("RAM2019_ENDO")])
    annotations = matrix.map(lambda value: "" if pd.isna(value) else f"{value:.1f}")
    sns.heatmap(matrix.clip(-2, 2), cmap="vlag", center=0, vmin=-2, vmax=2, annot=annotations, fmt="", linewidths=0.4,
                linecolor="white", cbar_kws={"label": "Hedges g", "shrink": 0.7}, ax=ax)
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in matrix.index], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=25)
    ax.text(0.98, -0.22, "0/28 rows exceeded matched-random 95th", transform=ax.transAxes, ha="right", color=RED, fontweight="bold")

    watson = pd.read_csv(repo / "results" / "sensitivity" / "gse210077_watson6_sensitivity_effects.csv")
    watson.to_csv(source_dir / "figure_5b_watson_effects.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "Watson6 shows a global negative shift across lineages")
    display = watson.copy()
    display["lineage_label"] = display["lineage"].map(LINEAGE_LABELS)
    sns.boxplot(data=display, x="lineage_label", y="hedges_g", hue="score_method", palette=METHOD_COLORS, showfliers=False, ax=ax)
    sns.stripplot(data=display, x="lineage_label", y="hedges_g", hue="score_method", palette=METHOD_COLORS, dodge=True,
                  alpha=0.55, size=3, ax=ax, legend=False)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Hedges g")
    ax.tick_params(axis="x", labelrotation=20)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ["singscore", "standardized mean"], frameon=False, loc="lower left")
    clean_axes(ax)
    ax.text(0.98, 0.89, "38/38 effects negative\n0/38 random-specific", transform=ax.transAxes, ha="right", va="top", color=RED, fontweight="bold")

    random = read_random(repo)
    random.to_csv(source_dir / "figure_5c_all_random_percentiles.csv", index=False)
    ax = fig.add_subplot(grid[1, 0])
    panel(ax, "C", "Random-module calibration distinguishes selective from global shifts")
    order = ["GSE202379", "GSE290642_human", "GSE244832", "GSE210077_Watson6", "GSE181483_human"]
    labels = ["GSE202379", "GSE290642", "GSE244832", "Watson6", "GSE181483"]
    sns.boxplot(data=random, x="dataset_id", y="real_effect_percentile", hue="score_method", order=order,
                palette=METHOD_COLORS, showfliers=False, ax=ax)
    sns.stripplot(data=random, x="dataset_id", y="real_effect_percentile", hue="score_method", order=order,
                  palette=METHOD_COLORS, dodge=True, alpha=0.22, size=2, ax=ax, legend=False)
    ax.axhline(0.95, color=BLACK, linestyle="--", lw=0.9)
    ax.set_xticks(range(len(order)), labels, rotation=25, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Real-effect percentile")
    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], legend_labels[:2], frameon=False, loc="lower left")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel(ax, "D", "No public resource satisfies all requirements for decisive external validation")
    comparability = pd.DataFrame(
        [
            ["GSE202379", 1, 1, 1, 1, 0],
            ["GSE290642", 0, 1, 1, 1, 0],
            ["GSE244832", 1, 0, 1, 1, 0],
            ["Watson6", 1, 1, 1, 1, 0],
            ["GSE181483", 0, 0, 0, 1, 0],
        ], columns=["dataset", "author_labels", "individual_stage", "at_least_3_per_group", "random_calibrated", "comparable_second_primary"]
    )
    comparability.to_csv(source_dir / "figure_5d_comparability_matrix.csv", index=False)
    matrix = comparability.set_index("dataset")
    matrix.columns = ["Author labels", "Individual stage", ">=3/group", "Random calibration", "Second formal primary"]
    sns.heatmap(matrix, cmap=sns.color_palette([LIGHT_GRAY, GREEN], as_cmap=True), vmin=0, vmax=1, cbar=False,
                annot=matrix.map(lambda x: "Yes" if x else "No"), fmt="", linewidths=0.7, linecolor="white", ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=30)
    ax.tick_params(axis="y", labelrotation=0)
    save(fig, output_dir, "figure_5_cohort_assay_stress_tests")


def figure_6(repo: Path, output_dir: Path, source_dir: Path) -> None:
    fig = plt.figure(figsize=(11.7, 8.2), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, width_ratios=[0.95, 1.35], height_ratios=[1.05, 0.95])
    meta = pd.read_csv(repo / "results" / "meta" / "advanced_endothelial_sensitivity_meta.csv")
    meta.to_csv(source_dir / "figure_6a_advanced_meta.csv", index=False)
    subgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[0, 0], wspace=0.12)
    programs = list(dict.fromkeys(meta["program_id"]))
    for index, (method, title) in enumerate([("singscore", "Singscore"), ("standardized_mean", "Standardized mean")]):
        ax = fig.add_subplot(subgrid[0, index])
        if index == 0:
            panel(ax, "A", "Singscore")
        else:
            ax.set_title(title, loc="left", pad=6, fontweight="bold")
        subset = meta[meta["score_method"].eq(method)].set_index("program_id").reindex(programs)
        y = np.arange(len(programs))
        ax.errorbar(subset["fixed_hedges_g"], y - 0.08,
                    xerr=[subset["fixed_hedges_g"] - subset["fixed_ci95_low"], subset["fixed_ci95_high"] - subset["fixed_hedges_g"]],
                    fmt="o", color=BLUE, markersize=3.6, capsize=2, label="Fixed")
        ax.errorbar(subset["random_reml_hedges_g"], y + 0.08,
                    xerr=[subset["random_reml_hedges_g"] - subset["random_reml_ci95_low"], subset["random_reml_ci95_high"] - subset["random_reml_hedges_g"]],
                    fmt="s", color=RED, markersize=3.4, capsize=2, label="REML")
        ax.axvline(0, color=GRAY, lw=0.7)
        ax.set_yticks(y, [PROGRAM_LABELS[p] for p in programs] if index == 0 else [])
        ax.invert_yaxis()
        ax.set_xlabel("Pooled Hedges g")
        clean_axes(ax)
        if index == 0:
            ax.legend(frameon=False, loc="lower right")
        if index == 1:
            ax.text(0.98, 0.02, "k=2; sensitivity only", transform=ax.transAxes, ha="right", color=RED, fontsize=7)

    evidence = pd.read_csv(repo / "results" / "exploratory" / "program_evidence_attrition_matrix.csv")
    evidence.to_csv(source_dir / "figure_6b_evidence_matrix.csv", index=False)
    ax = fig.add_subplot(grid[0, 1])
    panel(ax, "B", "Evidence is lost before formal replication")
    columns = [
        "evaluable_independent", "positive_both_scores_any_cohort", "random_specific_both_scores_any_cohort",
        "positive_primary_interval_both_scores", "advanced_sensitivity_meta_available", "within_cell_state_replicated",
        "pan_cirrhotic_transportable", "assay_robust",
    ]
    labels = ["Evaluable", "Both-score +", "Random-specific", "Primary CI +", "Comparable meta", "Replicated", "Pan-cirrhotic", "Assay robust"]
    matrix = evidence.set_index("program_id")[columns].astype(int).reindex(PROGRAM_ORDER)
    sns.heatmap(matrix, cmap=sns.color_palette([LIGHT_GRAY, GREEN], as_cmap=True), vmin=0, vmax=1, cbar=False,
                linewidths=0.45, linecolor="white", ax=ax)
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.axhline(7, color=BLACK, linewidth=1.0)
    ax.axhline(13, color=BLACK, linewidth=1.0)

    ax = fig.add_subplot(grid[1, 0])
    panel(ax, "C", "Evidence cascade")
    counts = matrix.sum(axis=0).rename_axis("criterion").reset_index(name="programs")
    counts["label"] = labels
    counts.to_csv(source_dir / "figure_6c_evidence_counts.csv", index=False)
    bars = ax.barh(np.arange(len(counts)), counts["programs"], color=[GREEN, GREEN, ORANGE, RED, BLUE, RED, RED, RED])
    ax.set_yticks(np.arange(len(counts)), counts["label"])
    ax.invert_yaxis()
    ax.set_xlim(0, 20)
    ax.set_xlabel("Programs (n=19)")
    for bar, value in zip(bars, counts["programs"]):
        ax.text(value + 0.35, bar.get_y() + bar.get_height() / 2, str(int(value)), va="center", fontsize=8)
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    panel(ax, "D", "What the negative validation result means")
    ax.axis("off")
    supported = [
        "Within-cohort program ranking is reproducible across scoring methods",
        "Selected macrophage programs form a focused MASH follow-up hypothesis",
        "Cohort context and data geometry limit cross-study transfer",
        "Donor-level external validation standards can prevent false universality",
    ]
    unsupported = [
        "Universal liver scar-cell state",
        "Fibrosis-stage biomarker",
        "Pan-etiologic diagnostic signature",
        "Causal mechanism or therapeutic target",
    ]
    ax.add_patch(patches.Rectangle((0.02, 0.12), 0.46, 0.72, facecolor="#EEF7F1", edgecolor=GREEN, linewidth=1.1, transform=ax.transAxes))
    ax.add_patch(patches.Rectangle((0.52, 0.12), 0.46, 0.72, facecolor="#FBEEEE", edgecolor=RED, linewidth=1.1, transform=ax.transAxes))
    ax.text(0.05, 0.78, "Supported", transform=ax.transAxes, color=GREEN, fontweight="bold", fontsize=9)
    ax.text(0.55, 0.78, "Not supported", transform=ax.transAxes, color=RED, fontweight="bold", fontsize=9)
    for index, text in enumerate(supported):
        wrapped = "+ " + textwrap.fill(text, width=37, subsequent_indent="  ")
        ax.text(0.05, 0.68 - index * 0.14, wrapped, transform=ax.transAxes, va="top", fontsize=6.8, linespacing=1.15)
    for index, text in enumerate(unsupported):
        wrapped = "- " + textwrap.fill(text, width=32, subsequent_indent="  ")
        ax.text(0.55, 0.68 - index * 0.14, wrapped, transform=ax.transAxes, va="top", fontsize=6.8, linespacing=1.15)
    boundary = pd.DataFrame({"supported": supported, "not_supported": unsupported})
    boundary.to_csv(source_dir / "figure_6d_conclusion_boundary.csv", index=False)
    ax.text(0.02, 0.02, "Clinical implication: these programs should not be used for cross-cohort staging or biomarker claims without new, matched donor validation.",
            transform=ax.transAxes, color=BLACK, fontsize=7.6, fontweight="bold")
    save(fig, output_dir, "figure_6_evidence_synthesis")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output_dir = repo / "results" / "figures"
    source_dir = repo / "results" / "source_data"
    output_dir.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)
    for old in output_dir.glob("figure_*.*"):
        old.unlink()
    for old in source_dir.glob("figure_*.csv"):
        old.unlink()
    figure_1(repo, output_dir, source_dir)
    figure_2(repo, output_dir, source_dir)
    figure_3(repo, output_dir, source_dir)
    figure_4(repo, output_dir, source_dir)
    figure_5(repo, output_dir, source_dir)
    figure_6(repo, output_dir, source_dir)
    print(f"Wrote six multi-panel figures (PNG/PDF) and source tables to {output_dir}")


if __name__ == "__main__":
    main()
