from __future__ import annotations

from pathlib import Path

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


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

SHORT_COHORT = {
    "GSE202379": "Clinical NASH\ncirrhosis",
    "GSE244832": "MASH F2–F4",
    "GSE256398": "MASH\ncirrhosis",
}

SHORT_CONTEXT = {
    "mash_f2f4_group_vs_normal_sensitivity": "GSE244832 MASH F2–F4",
    "mash_cirrhosis_vs_healthy": "GSE256398 MASH cirrhosis",
    "alcohol_cirrhosis_vs_healthy": "GSE256398 alcohol cirrhosis",
}

COMPONENT_COLORS = {
    "Program": BLUE,
    "Context": ORANGE,
    "Program × context / residual": "#B8B8B8",
}


def _p4(repo: Path, filename: str) -> pd.DataFrame:
    return pd.read_csv(repo / "results" / "phase4" / filename)


def _program_lineage(repo: Path) -> dict[str, str]:
    inventory = pd.read_csv(repo / "literature" / "program_inventory.csv")
    return inventory.drop_duplicates("program_id").set_index("program_id")["cell_lineage"].to_dict()


def _method_label(value: str) -> str:
    return "Rank score" if value == "singscore" else "z-mean"


def _stable_coupling(repo: Path) -> pd.DataFrame:
    stability = _p4(repo, "cross_lineage_coupling_stability.csv")
    meta = _p4(repo, "cross_lineage_coupling_meta.csv")
    stable = stability.loc[stability["stable_cross_lineage_coupling"]].copy()
    meta_wide = meta.pivot_table(
        index=["program_left", "lineage_left", "program_right", "lineage_right"],
        columns="score_method",
        values=["meta_spearman_rho", "meta_fdr_within_method"],
        aggfunc="first",
    ).reset_index()
    meta_wide.columns = [
        "__".join([str(part) for part in col if str(part)]) if isinstance(col, tuple) else col
        for col in meta_wide.columns
    ]
    keys = ["program_left", "lineage_left", "program_right", "lineage_right"]
    return stable.merge(meta_wide, on=keys, how="left")


def _force_layout(nodes: list[str], edges: list[tuple[str, str, float]]) -> dict[str, np.ndarray]:
    """Small deterministic force-directed layout without an external graph dependency."""
    rng = np.random.default_rng(20260831)
    position = rng.normal(0, 0.5, size=(len(nodes), 2))
    lookup = {node: idx for idx, node in enumerate(nodes)}
    ideal = 0.75
    for iteration in range(700):
        force = np.zeros_like(position)
        for left in range(len(nodes)):
            delta = position[left] - position
            distance2 = np.sum(delta * delta, axis=1) + 1e-3
            distance2[left] = np.inf
            force[left] += np.sum(delta / distance2[:, None], axis=0) * 0.018
        for left, right, weight in edges:
            i, j = lookup[left], lookup[right]
            delta = position[j] - position[i]
            distance = max(float(np.linalg.norm(delta)), 1e-3)
            attraction = (distance - ideal) * (0.025 + 0.03 * weight) * delta / distance
            force[i] += attraction
            force[j] -= attraction
        cooling = 0.12 * (1 - iteration / 700) + 0.01
        position += np.clip(force, -cooling, cooling)
        position -= position.mean(axis=0)
    scale = np.max(np.abs(position), axis=0)
    position = position / np.where(scale == 0, 1, scale)
    return {node: position[idx] for idx, node in enumerate(nodes)}


def _coupling_network(ax: plt.Axes, stable: pd.DataFrame) -> None:
    lineage_by_node: dict[str, str] = {}
    edges: list[tuple[str, str, float]] = []
    degrees: dict[str, int] = {}
    for row in stable.itertuples(index=False):
        lineage_by_node[row.program_left] = row.lineage_left
        lineage_by_node[row.program_right] = row.lineage_right
        edges.append((row.program_left, row.program_right, row.minimum_absolute_rho))
        degrees[row.program_left] = degrees.get(row.program_left, 0) + 1
        degrees[row.program_right] = degrees.get(row.program_right, 0) + 1
    nodes = sorted(lineage_by_node)
    ordered_nodes = {
        "endothelial": [
            "RAM2019_ENDO_6_SAENDO1",
            "RAM2019_ENDO_7_SAENDO2",
            "RAM2019_ENDO_5",
            "RAM2019_ENDO_1",
        ],
        "macrophage_monocyte": [
            "RAM2019_MAC_SIG_A_SAM",
            "RAM2019_MAC_SIG_B_SAM",
            "RAM2019_MAC_SIG_E_TMO",
            "RAM2019_MAC_SIG_F_CDC1",
        ],
        "mesenchymal_hsc_myofibroblast": [
            "RAM2019_MES_MESOTHELIAL",
            "RAM2019_MES_SAMES",
            "RAM2019_SAMES_B",
            "RAM2019_MES_VSMC",
            "RAM2019_SAMES_A",
        ],
    }
    x_by_lineage = {"endothelial": -1.0, "macrophage_monocyte": 0.0, "mesenchymal_hsc_myofibroblast": 1.0}
    position: dict[str, np.ndarray] = {}
    for lineage, ordered in ordered_nodes.items():
        present = [node for node in ordered if node in nodes]
        for node, yvalue in zip(present, np.linspace(0.72, -0.72, len(present))):
            position[node] = np.array([x_by_lineage[lineage], yvalue])
    for left, right, weight in edges:
        xy_left, xy_right = position[left], position[right]
        ax.plot(
            [xy_left[0], xy_right[0]],
            [xy_left[1], xy_right[1]],
            color="#8C8C8C",
            lw=0.6 + 5.2 * weight,
            alpha=0.62,
            zorder=1,
        )
    for node in nodes:
        xy = position[node]
        ax.scatter(
            xy[0],
            xy[1],
            s=205 + 70 * degrees[node],
            color=LINEAGE_COLORS[lineage_by_node[node]],
            edgecolor="white",
            linewidth=1.0,
            zorder=2,
        )
        ax.text(xy[0], xy[1], PROGRAM_LABELS[node], ha="center", va="center", fontsize=6.1, fontweight="bold", zorder=3)
    for lineage, xpos in x_by_lineage.items():
        ax.text(xpos, 0.96, LINEAGE_LABELS[lineage], ha="center", va="bottom", fontsize=6.7, fontweight="bold", color=LINEAGE_COLORS[lineage])
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-1.02, 1.08)
    ax.axis("off")
    ax.text(
        0.99,
        0.02,
        "16 edges; width = minimum |ρ|\nacross 2 cohorts × 2 scores",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.8,
        color=BLACK,
    )


def figure_8(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.2, 9.2), constrained_layout=True)
    grid = GridSpec(3, 2, figure=fig, height_ratios=[1.35, 1.05, 1.0])

    coherence = _p4(repo, "endothelial_member_gene_coherence.csv")
    summary = _p4(repo, "endothelial_member_gene_summary.csv")
    coherent = coherence.loc[coherence["coherent_positive_member"]].copy()
    coherent["program_label"] = coherent["program_id"].map(PROGRAM_LABELS)
    coherent = coherent.sort_values(["program_id", "minimum_g"], ascending=[True, False])
    cohort_cols = ["g__GSE202379", "g__GSE244832", "g__GSE256398"]

    # A: positive member-gene leading edge.
    ax = fig.add_subplot(grid[0, 0])
    matrix = coherent.set_index("gene_symbol")[cohort_cols]
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        linewidths=0.22,
        linecolor="white",
        cbar_kws={"label": "Member-gene Hedges g", "shrink": 0.72},
    )
    ax.set_xticklabels([SHORT_COHORT[col.replace("g__", "")] for col in cohort_cols], rotation=0, fontsize=6.6)
    labels = [f"{gene}{' *' if supported else ''}" for gene, supported in zip(coherent["gene_symbol"], coherent["meta_supported_member"])]
    ax.set_yticks(np.arange(len(labels)) + 0.5, labels, rotation=0, fontsize=5.7)
    ax.set_xlabel("")
    ax.set_ylabel("")
    first_count = int(coherent["program_id"].eq("RAM2019_ENDO_2").sum())
    ax.axhline(first_count, color=BLACK, lw=1.1)
    ax.text(-0.42, 1 - first_count / (2 * len(coherent)), "Endo2", transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=7.0, fontweight="bold", color=BLUE)
    ax.text(-0.42, (len(coherent) - first_count) / (2 * len(coherent)), "SAEndo1", transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=7.0, fontweight="bold", color=BLUE)
    ax.text(0, -0.16, "* fixed-effect meta-analysis FDR < 0.05 within program", transform=ax.transAxes, fontsize=6.7)
    panel_label(ax, "A")
    coherent.to_csv(source / "figure_8a_member_gene_coherence.csv", index=False)

    # B: detection-matched random-module calibration.
    ax = fig.add_subplot(grid[0, 1])
    random = _p4(repo, "endothelial_member_gene_random_controls.csv.gz")
    plot_summary = []
    for idx, row in summary.reset_index(drop=True).iterrows():
        program = row["program_id"]
        values = random.loc[random["program_id"].eq(program), "coherent_positive_genes"]
        bins = np.arange(values.min() - 0.5, values.max() + 1.5, 1)
        ax.hist(values, bins=bins, density=True, histtype="stepfilled", alpha=0.30, color=[BLUE, PURPLE][idx], label=f"{PROGRAM_LABELS[program]} random")
        ax.axvline(row["coherent_positive_genes"], color=[BLUE, PURPLE][idx], lw=2.1)
        ax.text(
            row["coherent_positive_genes"] + 0.25,
            0.49 - idx * 0.11,
            f"{PROGRAM_LABELS[program]}: {int(row['coherent_positive_genes'])}/{int(row['eligible_shared_genes'])}\nempirical P={row['empirical_p_coherent_count']:.4f}",
            color=[BLUE, PURPLE][idx],
            fontsize=7.1,
            va="top",
        )
        plot_summary.append(row)
    ax.set_xlim(-0.5, max(summary["coherent_positive_genes"]) + 5)
    ax.set_xlabel("Genes positive in all three cohorts")
    ax.set_ylabel("Density across 10,000 matched modules")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2, borderaxespad=0, fontsize=6.7)
    panel_label(ax, "B")
    clean_axes(ax)
    random.to_csv(source / "figure_8b_random_module_controls.csv", index=False)
    summary.to_csv(source / "figure_8b_observed_summary.csv", index=False)

    # C: stable cross-lineage coupling network.
    ax = fig.add_subplot(grid[1, 0])
    stable = _stable_coupling(repo)
    _coupling_network(ax, stable)
    panel_label(ax, "C")
    stable.to_csv(source / "figure_8c_stable_cross_lineage_network.csv", index=False)

    # D: shared cirrhosis signal versus etiology divergence.
    ax = fig.add_subplot(grid[1, 1])
    geometry = _p4(repo, "etiology_program_classification.csv")
    geometry["absolute_shared_component"] = geometry["median_shared_component"].abs()
    for lineage, values in geometry.groupby("lineage"):
        ax.scatter(
            values["absolute_shared_component"],
            values["median_absolute_divergence"],
            s=np.where(values["shared_random_specific_backbone"], 78, 38),
            color=LINEAGE_COLORS[lineage],
            alpha=np.where(values["shared_directional_backbone"], 0.90, 0.42),
            edgecolor=np.where(values["shared_random_specific_backbone"], BLACK, "white"),
            linewidth=np.where(values["shared_random_specific_backbone"], 1.0, 0.5),
            label=LINEAGE_LABELS[lineage],
        )
    max_value = float(max(geometry["absolute_shared_component"].max(), geometry["median_absolute_divergence"].max())) + 0.35
    ax.plot([0, max_value], [0, max_value], ls="--", color=GRAY, lw=0.8)
    for program in ["RAM2019_ENDO_2", "RAM2019_MAC_SIG_A_SAM", "RAM2019_MES_MESOTHELIAL", "RAM2019_ENDO_6_SAENDO1"]:
        row = geometry.loc[geometry["program_id"].eq(program)].iloc[0]
        ax.annotate(PROGRAM_LABELS[program], (row["absolute_shared_component"], row["median_absolute_divergence"]), xytext=(4, 3), textcoords="offset points", fontsize=6.5)
    ax.set_xlim(0, max_value)
    ax.set_ylim(0, max_value)
    ax.set_xlabel("|Shared cirrhosis component| (median across scores)")
    ax.set_ylabel("Etiology-divergence magnitude")
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, borderaxespad=0)
    ax.text(0.03, 0.96, "Below diagonal: shared component dominates", transform=ax.transAxes, va="top", fontsize=6.7)
    panel_label(ax, "D")
    clean_axes(ax)
    geometry.to_csv(source / "figure_8d_etiology_geometry.csv", index=False)

    # E: composition adjustment.
    ax = fig.add_subplot(grid[2, 0])
    composition = _p4(repo, "composition_adjusted_program_effects.csv")
    stability = _p4(repo, "composition_stability_summary.csv")
    comp = composition.merge(
        stability[["dataset_id", "contrast", "program_id", "composition_stable", "both_adjusted_intervals_positive"]],
        on=["dataset_id", "contrast", "program_id"],
        how="left",
    )
    for lineage, values in comp.groupby("lineage"):
        ax.scatter(
            values["unadjusted_standardized_beta"],
            values["adjusted_standardized_beta"],
            s=24,
            color=LINEAGE_COLORS[lineage],
            alpha=0.48,
            edgecolor="white",
            linewidth=0.35,
            label=LINEAGE_LABELS[lineage],
        )
    positive = comp.loc[comp["both_adjusted_intervals_positive"]]
    ax.scatter(
        positive["unadjusted_standardized_beta"],
        positive["adjusted_standardized_beta"],
        s=62,
        facecolors="none",
        edgecolors=BLACK,
        linewidths=1.05,
        label="Dual-score positive HC3 intervals",
    )
    bound = float(np.ceil(max(comp["unadjusted_standardized_beta"].abs().max(), comp["adjusted_standardized_beta"].abs().max()) * 1.05))
    ax.plot([-bound, bound], [-bound, bound], ls="--", lw=0.8, color=GRAY)
    ax.axhline(0, color=LIGHT_GRAY, lw=0.7)
    ax.axvline(0, color=LIGHT_GRAY, lw=0.7)
    ax.set_xlim(-bound, bound)
    ax.set_ylim(-bound, bound)
    ax.set_xlabel("Unadjusted standardized coefficient")
    ax.set_ylabel("Adjusted standardized coefficient")
    ax.legend(frameon=False, fontsize=5.9, loc="upper center", bbox_to_anchor=(0.5, -0.21), ncol=2, borderaxespad=0)
    ax.text(0.98, 0.03, "41/76 context–program pairs retained\nsign and ≥70% magnitude; 5 remained\npositive by both adjusted HC3 intervals", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5)
    panel_label(ax, "E", x=-0.17, y=1.04)
    clean_axes(ax)
    comp.to_csv(source / "figure_8e_composition_adjustment.csv", index=False)

    # F: variance partition.
    ax = fig.add_subplot(grid[2, 1])
    variance = _p4(repo, "program_context_variance_components.csv")
    variance["label"] = variance.apply(lambda row: f"{LINEAGE_LABELS[row['lineage']].replace('Macrophage/monocyte', 'Macrophage').replace('Mesenchymal/HSC', 'Mesenchymal')}\n{_method_label(row['score_method'])}", axis=1)
    x = np.arange(len(variance))
    bottom = np.zeros(len(variance))
    components = [
        ("Program", "program_fraction"),
        ("Context", "context_fraction"),
        ("Program × context / residual", "interaction_residual_fraction"),
    ]
    for label, column in components:
        ax.bar(x, variance[column], bottom=bottom, color=COMPONENT_COLORS[label], width=0.72, label=label, edgecolor="white", linewidth=0.5)
        bottom += variance[column].to_numpy()
    ax.set_xticks(x, variance["label"], rotation=24, ha="right", fontsize=6.4)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of effect-matrix variance")
    ax.legend(frameon=False, fontsize=6.2, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.26), borderaxespad=0)
    ax.text(0.02, 0.96, "Median context = 55.2%\nMedian program = 4.7%", transform=ax.transAxes, va="top", fontsize=6.9, bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=LIGHT_GRAY))
    panel_label(ax, "F")
    clean_axes(ax)
    variance.to_csv(source / "figure_8f_variance_partition.csv", index=False)

    save(fig, output, "figure_8_biological_structure")


def supplementary_figure_7(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.0, 9.1), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, width_ratios=[1.12, 0.88], height_ratios=[1.4, 0.85])
    coherence = _p4(repo, "endothelial_member_gene_coherence.csv")
    cohort_cols = ["g__GSE202379", "g__GSE244832", "g__GSE256398"]
    coherence = coherence.sort_values(["program_id", "minimum_g"], ascending=[True, False])

    ax = fig.add_subplot(grid[0, 0])
    matrix = coherence.set_index("gene_symbol")[cohort_cols]
    sns.heatmap(matrix, ax=ax, cmap="vlag", center=0, vmin=-3, vmax=3, linewidths=0.12, linecolor="white", cbar_kws={"label": "Member-gene Hedges g", "shrink": 0.65})
    labels = [
        f"{row.gene_symbol}{' †' if row.coherent_positive_member else ''}{' *' if row.meta_supported_member else ''}"
        for row in coherence.itertuples(index=False)
    ]
    ax.set_yticks(np.arange(len(labels)) + 0.5, labels, rotation=0, fontsize=4.8)
    ax.set_xticklabels([SHORT_COHORT[col.replace("g__", "")] for col in cohort_cols], rotation=0, fontsize=6.5)
    split = int(coherence["program_id"].eq("RAM2019_ENDO_2").sum())
    ax.axhline(split, color=BLACK, lw=1.0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.text(0, -0.12, "† positive in all cohorts; * fixed-effect FDR < 0.05", transform=ax.transAxes, fontsize=6.6)
    panel(ax, "A", "Complete member-gene audit for Endo2 and SAEndo1")

    ax = fig.add_subplot(grid[0, 1])
    for program, values in coherence.groupby("program_id"):
        ax.scatter(values["fixed_effect"], -np.log10(values["fixed_fdr_within_program"].clip(lower=1e-10)), color=BLUE if program == "RAM2019_ENDO_2" else PURPLE, s=np.where(values["meta_supported_member"], 45, 22), alpha=np.where(values["meta_supported_member"], 0.90, 0.42), edgecolor=np.where(values["meta_supported_member"], BLACK, "white"), linewidth=0.5, label=PROGRAM_LABELS[program])
    ax.axhline(-np.log10(0.05), ls="--", color=GRAY, lw=0.8)
    ax.axvline(0, color=LIGHT_GRAY, lw=0.7)
    label_genes = {"TFF3", "EFEMP1", "LGALS3", "GSN", "RBP7", "PLPP1", "PLVAP", "VWA1"}
    offsets = {
        "TFF3": (4, -11), "EFEMP1": (4, 2), "LGALS3": (4, 11),
        "GSN": (4, -9), "RBP7": (-31, 7), "PLPP1": (4, 2), "PLVAP": (4, 2), "VWA1": (4, 2),
    }
    for row in coherence.loc[coherence["meta_supported_member"] & coherence["gene_symbol"].isin(label_genes)].itertuples(index=False):
        ax.annotate(row.gene_symbol, (row.fixed_effect, -np.log10(max(row.fixed_fdr_within_program, 1e-10))), xytext=offsets[row.gene_symbol], textcoords="offset points", fontsize=5.8)
    ax.set_xlabel("Fixed-effect meta-analytic Hedges g")
    ax.set_ylabel("−log10 within-program FDR")
    ax.legend(frameon=False, fontsize=6.7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, borderaxespad=0)
    panel(ax, "B", "Meta-supported genes define a compact positive core")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 0])
    loco = _p4(repo, "endothelial_member_gene_loco.csv")
    cohorts = list(SHORT_COHORT)
    x = np.arange(len(cohorts))
    width = 0.34
    for idx, program in enumerate(["RAM2019_ENDO_2", "RAM2019_ENDO_6_SAENDO1"]):
        values = loco.loc[loco["program_id"].eq(program)].set_index("held_out_dataset").reindex(cohorts)
        ax.bar(x + (idx - 0.5) * width, values["held_out_sign_retention"], width=width, color=[BLUE, PURPLE][idx], label=PROGRAM_LABELS[program])
        for xx, value in zip(x + (idx - 0.5) * width, values["held_out_sign_retention"]):
            ax.text(xx, value + 0.025, f"{value:.0%}", ha="center", fontsize=6.2)
    ax.set_xticks(x, [SHORT_COHORT[value].replace("\n", " ") for value in cohorts], rotation=15, ha="right", fontsize=6.4)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Held-out positive-sign retention")
    ax.legend(frameon=False, fontsize=6.7, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, borderaxespad=0)
    panel(ax, "C", "Leave-one-cohort-out discovery retains direction in the held-out cohort")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    summary = _p4(repo, "endothelial_member_gene_summary.csv")
    xpos = np.arange(len(summary))
    ax.bar(xpos, summary["coherent_positive_genes"], color=[BLUE, PURPLE], width=0.58, label="Observed")
    ax.errorbar(xpos, summary["random_coherent_median"], yerr=[summary["random_coherent_median"] - summary["random_coherent_median"], summary["random_coherent_95th"] - summary["random_coherent_median"]], fmt="o", color=BLACK, capsize=4, label="Random median to 95th")
    for xx, row in zip(xpos, summary.itertuples(index=False)):
        ax.text(xx, row.coherent_positive_genes + 0.4, f"P={row.empirical_p_coherent_count:.4f}", ha="center", fontsize=6.5)
    ax.set_xticks(xpos, [PROGRAM_LABELS[value] for value in summary["program_id"]])
    ax.set_ylabel("Coherent-positive genes")
    ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, borderaxespad=0)
    panel(ax, "D", "Observed coherence clears the matched-random 95th percentile")
    clean_axes(ax)

    coherence.to_csv(source / "supplementary_figure_7ab_member_gene_audit.csv", index=False)
    loco.to_csv(source / "supplementary_figure_7c_loco.csv", index=False)
    summary.to_csv(source / "supplementary_figure_7d_random_summary.csv", index=False)
    save(fig, output, "supplementary_figure_7_member_gene_audit")


def _meta_matrix(meta: pd.DataFrame, method: str, order: list[str]) -> pd.DataFrame:
    matrix = pd.DataFrame(np.nan, index=order, columns=order)
    for row in meta.loc[meta["score_method"].eq(method)].itertuples(index=False):
        matrix.loc[row.program_left, row.program_right] = row.meta_spearman_rho
        matrix.loc[row.program_right, row.program_left] = row.meta_spearman_rho
    return matrix


def supplementary_figure_8(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.2, 9.2), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 0.85])
    meta = _p4(repo, "cross_lineage_coupling_meta.csv")
    stability = _p4(repo, "cross_lineage_coupling_stability.csv")
    inventory = pd.read_csv(repo / "literature" / "program_inventory.csv")
    order = inventory.drop_duplicates("program_id").sort_values(["cell_lineage", "program_id"])["program_id"].tolist()
    stable_keys = set(zip(stability.loc[stability["stable_cross_lineage_coupling"], "program_left"], stability.loc[stability["stable_cross_lineage_coupling"], "program_right"]))

    for letter, method, location in [("A", "singscore", grid[0, 0]), ("B", "standardized_mean", grid[0, 1])]:
        ax = fig.add_subplot(location)
        matrix = _meta_matrix(meta, method, order)
        sns.heatmap(matrix, ax=ax, cmap="vlag", center=0, vmin=-0.75, vmax=0.75, square=True, linewidths=0.18, linecolor="white", mask=matrix.isna(), cbar_kws={"label": "Meta Spearman ρ", "shrink": 0.68})
        labels = [PROGRAM_LABELS[value] for value in order]
        ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=5.3)
        ax.set_yticklabels(labels, rotation=0, fontsize=5.3)
        for i, left in enumerate(order):
            for j, right in enumerate(order):
                if (left, right) in stable_keys or (right, left) in stable_keys:
                    ax.add_patch(patches.Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, fill=False, edgecolor=BLACK, lw=0.85))
        panel(ax, letter, f"Cross-lineage coupling meta-map: {_method_label(method)}")

    ax = fig.add_subplot(grid[1, 0])
    wide = meta.pivot_table(index=["program_left", "lineage_left", "program_right", "lineage_right"], columns="score_method", values="meta_spearman_rho", aggfunc="first").reset_index()
    stable_lookup = stability.set_index(["program_left", "lineage_left", "program_right", "lineage_right"])["stable_cross_lineage_coupling"]
    wide["stable"] = [bool(stable_lookup.loc[(row.program_left, row.lineage_left, row.program_right, row.lineage_right)]) for row in wide.itertuples(index=False)]
    ax.scatter(wide["singscore"], wide["standardized_mean"], s=np.where(wide["stable"], 48, 20), color=np.where(wide["stable"], PURPLE, LIGHT_GRAY), alpha=np.where(wide["stable"], 0.9, 0.5), edgecolor=np.where(wide["stable"], BLACK, "white"), linewidth=0.55)
    ax.plot([-0.8, 0.8], [-0.8, 0.8], ls="--", color=GRAY, lw=0.8)
    key_pairs = {
        ("RAM2019_ENDO_7_SAENDO2", "RAM2019_MAC_SIG_B_SAM"): (12, 16),
        ("RAM2019_ENDO_6_SAENDO1", "RAM2019_MAC_SIG_A_SAM"): (-95, -15),
        ("RAM2019_MAC_SIG_E_TMO", "RAM2019_MES_MESOTHELIAL"): (12, -17),
    }
    for row in wide.loc[wide["stable"]].itertuples(index=False):
        key = (row.program_left, row.program_right)
        if key not in key_pairs:
            continue
        ax.annotate(
            f"{PROGRAM_LABELS[row.program_left]}–{PROGRAM_LABELS[row.program_right]}",
            (row.singscore, row.standardized_mean),
            xytext=key_pairs[key],
            textcoords="offset points",
            fontsize=5.8,
            arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.55),
        )
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    ax.set_xlabel("Rank-score meta ρ")
    ax.set_ylabel("z-mean meta ρ")
    ax.text(0.03, 0.95, "Purple: stable across both cohorts and scores", transform=ax.transAxes, va="top", fontsize=6.7)
    panel(ax, "C", "Stable edges are method-concordant")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    stable = stability.loc[stability["stable_cross_lineage_coupling"]]
    degree = pd.concat([stable[["program_left", "lineage_left"]].rename(columns={"program_left": "program_id", "lineage_left": "lineage"}), stable[["program_right", "lineage_right"]].rename(columns={"program_right": "program_id", "lineage_right": "lineage"})]).value_counts(["program_id", "lineage"]).rename("degree").reset_index().sort_values(["degree", "program_id"], ascending=[False, True])
    ax.barh(np.arange(len(degree)), degree["degree"], color=[LINEAGE_COLORS[value] for value in degree["lineage"]])
    ax.set_yticks(np.arange(len(degree)), [PROGRAM_LABELS[value] for value in degree["program_id"]], fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("Stable-network degree")
    ax.set_xticks(range(0, int(degree["degree"].max()) + 1))
    panel(ax, "D", "SAEndo2, SAM-B, TMo-E, and SAMes-B are coupling hubs")
    clean_axes(ax)

    meta.to_csv(source / "supplementary_figure_8ab_coupling_meta.csv", index=False)
    wide.to_csv(source / "supplementary_figure_8c_method_concordance.csv", index=False)
    degree.to_csv(source / "supplementary_figure_8d_hub_degree.csv", index=False)
    save(fig, output, "supplementary_figure_8_cross_lineage_coupling")


def supplementary_figure_9(repo: Path, output: Path, source: Path) -> None:
    fig = plt.figure(figsize=(12.0, 8.6), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig)
    components = _p4(repo, "etiology_shared_divergent_components.csv")

    ax = fig.add_subplot(grid[0, 0])
    for lineage, values in components.groupby("lineage"):
        ax.scatter(values["g_mash_cirrhosis_vs_healthy"], values["g_alcohol_cirrhosis_vs_healthy"], color=LINEAGE_COLORS[lineage], marker="o", s=32, alpha=0.62, edgecolor="white", linewidth=0.4, label=LINEAGE_LABELS[lineage])
    bound = float(np.ceil(max(components["g_mash_cirrhosis_vs_healthy"].abs().max(), components["g_alcohol_cirrhosis_vs_healthy"].abs().max())))
    ax.plot([-bound, bound], [-bound, bound], ls="--", color=GRAY, lw=0.8)
    both_random = components.loc[components["mash_above_random_95th"] & components["alcohol_above_random_95th"]]
    ax.scatter(both_random["g_mash_cirrhosis_vs_healthy"], both_random["g_alcohol_cirrhosis_vs_healthy"], s=70, facecolors="none", edgecolors=BLACK, linewidths=1.0)
    for row in both_random.itertuples(index=False):
        ax.annotate(PROGRAM_LABELS[row.program_id], (row.g_mash_cirrhosis_vs_healthy, row.g_alcohol_cirrhosis_vs_healthy), xytext=(3, 3), textcoords="offset points", fontsize=6.2)
    ax.set_xlim(-bound, bound)
    ax.set_ylim(-bound, bound)
    ax.set_xlabel("MASH cirrhosis vs healthy Hedges g")
    ax.set_ylabel("Alcohol cirrhosis vs healthy Hedges g")
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, borderaxespad=0)
    panel(ax, "A", "Same-assay etiologic contrasts expose shared and divergent programs")
    clean_axes(ax)

    ax = fig.add_subplot(grid[0, 1])
    selected = components.loc[components["program_id"].isin(["RAM2019_ENDO_2", "RAM2019_MAC_SIG_A_SAM", "RAM2019_ENDO_6_SAENDO1", "RAM2019_MES_MESOTHELIAL"])].copy()
    selected["label"] = selected.apply(lambda row: f"{PROGRAM_LABELS[row['program_id']]} · {_method_label(row['score_method'])}", axis=1)
    selected = selected.sort_values(["program_id", "score_method"])
    y = np.arange(len(selected))
    ax.errorbar(selected["shared_cirrhosis_component"], y - 0.11, xerr=[selected["shared_cirrhosis_component"] - selected["shared_bootstrap_low"], selected["shared_bootstrap_high"] - selected["shared_cirrhosis_component"]], fmt="o", color=BLUE, ecolor=BLUE, capsize=2.5, markersize=4.2, label="Shared component")
    ax.errorbar(selected["etiology_divergence_component"], y + 0.11, xerr=[selected["etiology_divergence_component"] - selected["divergence_bootstrap_low"], selected["divergence_bootstrap_high"] - selected["etiology_divergence_component"]], fmt="s", color=RED, ecolor=RED, capsize=2.5, markersize=3.8, label="Etiology divergence")
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_yticks(y, selected["label"], fontsize=6.3)
    ax.invert_yaxis()
    ax.set_xlabel("Component estimate with 95% donor bootstrap interval")
    ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2, borderaxespad=0)
    panel(ax, "B", "Bootstrap intervals calibrate the shared-versus-divergent decomposition")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 0])
    stability = _p4(repo, "composition_stability_summary.csv")
    stability["context_label"] = stability.apply(lambda row: SHORT_CONTEXT.get(row["contrast"], row["contrast"]), axis=1)
    stability["positive"] = stability["both_adjusted_intervals_positive"]
    for lineage, values in stability.groupby("lineage"):
        jitter = np.linspace(-0.035, 0.035, len(values))
        ax.scatter(values["minimum_magnitude_retention"].clip(upper=2.5), values["composition_stable"].astype(int) + jitter, s=np.where(values["positive"], 62, 24), color=LINEAGE_COLORS[lineage], alpha=np.where(values["positive"], 0.95, 0.48), edgecolor=np.where(values["positive"], BLACK, "white"), linewidth=0.6, label=LINEAGE_LABELS[lineage])
    ax.axvline(0.70, ls="--", color=GRAY, lw=0.8)
    ax.set_yticks([0, 1], ["Not stable", "Sign + ≥70% retained"])
    ax.set_xlabel("Minimum adjusted/unadjusted magnitude ratio across scores (clipped at 2.5)")
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, borderaxespad=0)
    ax.text(0.98, 0.05, "Large outlined points: both adjusted HC3 intervals > 0", transform=ax.transAxes, ha="right", fontsize=6.5)
    panel(ax, "C", "Composition sensitivity separates persistence from interval-level support")
    clean_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    descriptor = _p4(repo, "context_descriptor_permutation_regression.csv")
    predictors = ["same_endpoint", "same_etiology", "same_assay", "same_annotation"]
    labels = ["Same endpoint", "Same etiology", "Same assay", "Same annotation"]
    offsets = {("endothelial", "singscore"): -0.21, ("endothelial", "standardized_mean"): -0.07, ("macrophage_monocyte", "singscore"): 0.07, ("macrophage_monocyte", "standardized_mean"): 0.21, ("mesenchymal_hsc_myofibroblast", "singscore"): 0.35, ("mesenchymal_hsc_myofibroblast", "standardized_mean"): 0.49}
    base = np.arange(len(predictors))
    for (lineage, method), values in descriptor.groupby(["lineage", "score_method"]):
        values = values.set_index("predictor").reindex(predictors)
        y = base + offsets[(lineage, method)] - 0.14
        ax.scatter(values["coefficient_on_spearman_rho"], y, s=28, color=LINEAGE_COLORS[lineage], marker="o" if method == "singscore" else "s", alpha=0.78)
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_yticks(base + 0.04, labels)
    ax.set_xlabel("Coefficient on cross-context Spearman ρ")
    ax.text(0.98, 0.04, "0/24 descriptors significant after FDR correction", transform=ax.transAxes, ha="right", fontsize=6.7, bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=LIGHT_GRAY))
    panel(ax, "D", "No single recorded descriptor explains transport topology after correction")
    clean_axes(ax)

    components.to_csv(source / "supplementary_figure_9ab_etiology_components.csv", index=False)
    stability.to_csv(source / "supplementary_figure_9c_composition_stability.csv", index=False)
    descriptor.to_csv(source / "supplementary_figure_9d_descriptor_regression.csv", index=False)
    save(fig, output, "supplementary_figure_9_context_structure")


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    output = repo / "results" / "figures"
    source = repo / "results" / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    figure_8(repo, output, source)
    supplementary_figure_7(repo, output, source)
    supplementary_figure_8(repo, output, source)
    supplementary_figure_9(repo, output, source)
    print("Phase 4 figures and source-data files written.")


if __name__ == "__main__":
    main()
