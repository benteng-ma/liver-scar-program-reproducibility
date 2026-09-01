"""Build Phase 0 metadata and the frozen, traceable program inventory.

This script performs metadata-only curation.  It never estimates a disease
effect and never combines expression matrices across studies.
"""

from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "phase0_metadata"
META = ROOT / "metadata"
LIT = ROOT / "literature"
META.mkdir(exist_ok=True)
LIT.mkdir(exist_ok=True)


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def geo_samples(accession: str) -> list[dict[str, object]]:
    root = ET.parse(RAW / f"{accession}_family.xml").getroot()
    rows: list[dict[str, object]] = []
    for sample in root.iter():
        if local(sample.tag) != "Sample":
            continue
        title = next(
            ((node.text or "").strip() for node in sample.iter() if local(node.tag) == "Title"),
            "",
        )
        chars = {
            node.attrib.get("tag", "").strip().lower(): (node.text or "").strip()
            for node in sample.iter()
            if local(node.tag) == "Characteristics"
        }
        rows.append({"sample_id": sample.attrib.get("iid", ""), "title": title, "chars": chars})
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


dataset_fields = [
    "dataset_id", "publication", "pmid", "doi", "repository", "accession", "assay",
    "tissue", "cell_enrichment", "disease", "etiology", "fibrosis_stage",
    "compensation_status", "n_reported_samples", "n_independent_donors",
    "donor_id_available", "processed_counts_available", "original_annotation_available",
    "spatial_data_available", "program_discovery_source", "primary_validation_candidate",
    "verified", "limitations",
]


datasets = [
    dict(dataset_id="GSE136103", publication="Ramachandran et al., Nature 2019",
         pmid="31597160", doi="10.1038/s41586-019-1631-3", repository="GEO", accession="GSE136103",
         assay="scRNA-seq", tissue="human liver non-parenchymal cells", cell_enrichment="CD45+ and CD45- FACS fractions",
         disease="healthy and cirrhosis", etiology="NAFLD; alcohol; PBC", fibrosis_stage="cirrhosis",
         compensation_status="not reported", n_reported_samples=20, n_independent_donors=10,
         donor_id_available="yes", processed_counts_available="yes", original_annotation_available="yes",
         spatial_data_available="no", program_discovery_source="yes", primary_validation_candidate="internal recovery only",
         verified="yes", limitations="5 healthy/5 cirrhotic; macrophages 5+5, endothelial and mesenchymal 4+3; technical fractions must be collapsed; discovery leakage prohibited"),
    dict(dataset_id="GSE181483_human", publication="Zhang et al., Science Translational Medicine 2021",
         pmid="34613814", doi="10.1126/scitranslmed.abd1206", repository="GEO", accession="GSE181483",
         assay="scRNA-seq", tissue="human liver non-parenchymal cells", cell_enrichment="NPC",
         disease="control and cirrhosis", etiology="not recoverable from GEO metadata", fibrosis_stage="cirrhosis",
         compensation_status="not reported", n_reported_samples=4, n_independent_donors=4,
         donor_id_available="sample-level only", processed_counts_available="yes", original_annotation_available="partial",
         spatial_data_available="no", program_discovery_source="no", primary_validation_candidate="directional support only",
         verified="yes", limitations="only 2 controls and 2 cirrhosis donors; remaining GEO samples are pig; not suitable for standalone inference"),
    dict(dataset_id="E-MTAB-10553_unique", publication="Wang et al., Scientific Reports 2021",
         pmid="34588551", doi="10.1038/s41598-021-98806-y", repository="ArrayExpress/BioStudies", accession="E-MTAB-10553",
         assay="scRNA-seq", tissue="human liver non-parenchymal cells", cell_enrichment="hepatic NPC",
         disease="normal only", etiology="none", fibrosis_stage="none",
         compensation_status="not applicable", n_reported_samples=6, n_independent_donors=6,
         donor_id_available="yes", processed_counts_available="yes", original_annotation_available="yes",
         spatial_data_available="no", program_discovery_source="no", primary_validation_candidate="healthy assay baseline only",
         verified="yes", limitations="all six newly profiled donors are normal; processed RDS also embeds public GSE136103 cells and must not be treated as an independent disease cohort"),
    dict(dataset_id="GSE202379", publication="Gribben et al., Nature 2024",
         pmid="38778114", doi="10.1038/s41586-024-07465-2", repository="GEO", accession="GSE202379",
         assay="snRNA-seq", tissue="whole human liver", cell_enrichment="nuclei; no lineage enrichment",
         disease="healthy; NAFLD; NASH; NASH cirrhosis; end-stage liver disease", etiology="metabolic; end-stage etiology not fully encoded in GEO",
         fibrosis_stage="SAF F0-F4 plus end-stage", compensation_status="not reported", n_reported_samples=59,
         n_independent_donors=47, donor_id_available="yes", processed_counts_available="yes",
         original_annotation_available="yes", spatial_data_available="no", program_discovery_source="no",
         primary_validation_candidate="yes", verified="yes",
         limitations="59 regions/samples are 47 donors; repeated lobes from PHL1 and end-stage donors must be collapsed; spectrum includes non-cirrhotic fibrosis"),
    dict(dataset_id="GSE210077_Watson6", publication="Watson et al., Nature Communications 2025",
         pmid="39747812", doi="10.1038/s41467-024-55325-4", repository="GEO/HCA/Dryad", accession="GSE210077; HCA 64809a52-f703-4aec-b3a5-eca808a971d0",
         assay="snRNA-seq + MERFISH", tissue="whole human liver", cell_enrichment="nuclei and intact tissue imaging",
         disease="3 healthy; 3 fibrotic", etiology="not reported", fibrosis_stage="healthy; F2; F3; F4",
         compensation_status="not reported", n_reported_samples=6, n_independent_donors=6,
         donor_id_available="yes", processed_counts_available="yes", original_annotation_available="yes",
         spatial_data_available="yes", program_discovery_source="no", primary_validation_candidate="spatial and cross-assay support",
         verified="yes", limitations="fixed six-donor Watson subset only; spatial donors are matched to snRNA donors; disease group contains one F4 donor"),
    dict(dataset_id="GSE210077_full_series", publication="evolving GEO series; multiple linked studies",
         pmid="39747812", doi="10.1038/s41467-024-55325-4", repository="GEO", accession="GSE210077",
         assay="snRNA-seq", tissue="whole human liver", cell_enrichment="nuclei",
         disease="mostly normal plus F2/F3/F4 fibrosis", etiology="not reported", fibrosis_stage="normal; F2; F3; F4",
         compensation_status="not reported", n_reported_samples=40, n_independent_donors=40,
         donor_id_available="sample-as-donor for GEO entries", processed_counts_available="yes", original_annotation_available="yes for linked processed objects",
         spatial_data_available="only Watson six-donor subset", program_discovery_source="partly", primary_validation_candidate="no as one undifferentiated cohort",
         verified="yes", limitations="GEO has 40 records named Liver-1 through Liver-41 with missing numbers; later additions belong to other studies; roles must be subset-specific"),
    dict(dataset_id="GSE244832", publication="Kim et al., Journal of Hepatology 2025",
         pmid="39522884", doi="10.1016/j.jhep.2024.10.044", repository="GEO", accession="GSE244832",
         assay="snRNA-seq + snATAC-seq", tissue="whole human liver", cell_enrichment="nuclei",
         disease="normal; MASL; MASH; MetALD", etiology="metabolic and metabolic-alcohol overlap", fibrosis_stage="heterogeneous; donor-level mapping in processed metadata",
         compensation_status="not reported", n_reported_samples=18, n_independent_donors=18,
         donor_id_available="yes", processed_counts_available="yes", original_annotation_available="yes",
         spatial_data_available="no", program_discovery_source="yes for hA1 HSC program", primary_validation_candidate="candidate after condition mapping",
         verified="yes", limitations="GEO sample pages omit condition; processed 694 MB archive contains cell metadata; hA1 program cannot be externally validated in this same dataset"),
    dict(dataset_id="GSE290642_human", publication="Hu et al., Cell 2026",
         pmid="41794026", doi="10.1016/j.cell.2026.02.001", repository="GEO", accession="GSE290642",
         assay="scRNA-seq", tissue="human liver non-parenchymal cells", cell_enrichment="NPC",
         disease="F0 controls and F1-F4 fibrosis", etiology="not encoded in GEO", fibrosis_stage="F0-F4",
         compensation_status="not reported", n_reported_samples=24, n_independent_donors=24,
         donor_id_available="sample-level only", processed_counts_available="yes", original_annotation_available="not separately deposited",
         spatial_data_available="no", program_discovery_source="yes for ROCK2 angiocrine states", primary_validation_candidate="fibrosis-spectrum sensitivity",
         verified="yes", limitations="21 fibrosis donors but only three F4; most are F1/F1-2; etiology unresolved; 12 pig samples excluded"),
    dict(dataset_id="PRJNA833766_new4", publication="Bai et al., JHEP Reports 2023",
         pmid="37867598", doi="10.1016/j.jhepr.2023.100883", repository="SRA", accession="PRJNA833766",
         assay="scRNA-seq", tissue="human liver", cell_enrichment="CD45+ FACS",
         disease="2 healthy; 2 HBV cirrhosis newly generated", etiology="HBV", fibrosis_stage="cirrhosis",
         compensation_status="compensated/decompensated labels require run-table recovery", n_reported_samples=4,
         n_independent_donors=4, donor_id_available="BioSample", processed_counts_available="no; raw only",
         original_annotation_available="publication-integrated object not confirmed", spatial_data_available="no",
         program_discovery_source="no", primary_validation_candidate="macrophage directional only",
         verified="yes", limitations="raw data only; published 6 healthy/5 HBV total reuses HRA001730 and HRA000069; only four donors are new; immune-only"),
    dict(dataset_id="GSE185477", publication="Andrews et al., Hepatology Communications 2022",
         pmid="34792289", doi="10.1002/hep4.1854", repository="GEO", accession="GSE185477",
         assay="paired scRNA-seq and snRNA-seq; Visium", tissue="healthy human liver", cell_enrichment="whole tissue/NPC protocol comparison",
         disease="healthy only", etiology="none", fibrosis_stage="none", compensation_status="not applicable",
         n_reported_samples=5, n_independent_donors=5, donor_id_available="yes", processed_counts_available="yes",
         original_annotation_available="yes", spatial_data_available="yes", program_discovery_source="no",
         primary_validation_candidate="assay baseline only", verified="yes", limitations="four paired sc/sn donors plus a separate spatial donor; no disease contrast"),
    dict(dataset_id="ALC_HBV_2023_access_pending", publication="Zhang et al., Frontiers in Endocrinology 2023",
         pmid="36817578", doi="10.3389/fendo.2023.1132085", repository="publication supplement; accession unresolved",
         accession="none recovered", assay="scRNA-seq", tissue="human liver NPC", cell_enrichment="NPC; hepatocyte depletion",
         disease="3 healthy; 2 alcohol cirrhosis; 3 HBV cirrhosis", etiology="alcohol; HBV", fibrosis_stage="cirrhosis",
         compensation_status="not reported", n_reported_samples=8, n_independent_donors=8,
         donor_id_available="publication-level", processed_counts_available="not publicly recovered",
         original_annotation_available="figures/supplement only", spatial_data_available="no", program_discovery_source="no",
         primary_validation_candidate="access pending", verified="yes", limitations="paper examined 15 tissue specimens but scRNA used only eight; public expression object/accession not recovered in Phase 0"),
]
write_csv(META / "dataset_manifest.csv", datasets, dataset_fields)


donor_fields = [
    "dataset_id", "donor_id", "sample_id", "disease_group", "etiology", "fibrosis_stage",
    "compensated_or_decompensated", "age", "sex", "tissue_source", "cell_enrichment", "assay",
    "technical_fraction", "possible_overlap", "included", "exclusion_reason", "notes",
]
donors: list[dict[str, object]] = []

# GSE136103: retain sample rows but explicitly collapse fractions to ten donors.
for row in geo_samples("GSE136103"):
    title, chars = str(row["title"]), row["chars"]
    if title.startswith("Blood") or title.startswith("Mouse"):
        continue
    m = re.match(r"(Healthy|Cirrhotic)(\d+)", title, re.I)
    donor_id = f"{m.group(1).capitalize()}{m.group(2)}" if m else "UNRESOLVED"
    donors.append(dict(dataset_id="GSE136103", donor_id=donor_id, sample_id=row["sample_id"],
        disease_group="healthy" if donor_id.startswith("Healthy") else "cirrhosis",
        etiology=chars.get("cause of liver disease", ""), fibrosis_stage="none" if donor_id.startswith("Healthy") else "cirrhosis",
        compensated_or_decompensated="not reported", age="not reported", sex=chars.get("sex", ""),
        tissue_source="liver", cell_enrichment=chars.get("population", ""), assay="scRNA-seq",
        technical_fraction=chars.get("cell subtype", "") + ("; technical replicate" if title.endswith("A") or title.endswith("B") else ""),
        possible_overlap="none within project; discovery source", included="yes",
        exclusion_reason="", notes="Multiple CD45 fractions/replicates remain one biological donor."))

for row in geo_samples("GSE181483"):
    title = str(row["title"])
    if not title.startswith("human_"):
        continue
    disease = "healthy" if "ctrl" in title else "cirrhosis"
    donors.append(dict(dataset_id="GSE181483_human", donor_id=title, sample_id=row["sample_id"],
        disease_group=disease, etiology="unreported", fibrosis_stage="none" if disease == "healthy" else "cirrhosis",
        compensated_or_decompensated="not reported", age="not reported", sex="not reported", tissue_source="liver",
        cell_enrichment="NPC", assay="scRNA-seq", technical_fraction="none", possible_overlap="none identified",
        included="directional only", exclusion_reason="n<3 per group for standalone inferential model",
        notes="Pig controls/treatments excluded."))

# ArrayExpress unique donor rows (SDRF has one row per sequencing file/lane).
sdrf = pd.read_csv(RAW / "E-MTAB-10553.sdrf.txt", sep="\t", dtype=str)
for _, row in sdrf.drop_duplicates("Characteristics[individual]").iterrows():
    donors.append(dict(dataset_id="E-MTAB-10553_unique", donor_id=row["Characteristics[individual]"],
        sample_id=row["Source Name"], disease_group=row["Characteristics[disease]"], etiology="none",
        fibrosis_stage="none", compensated_or_decompensated="not applicable", age=row["Characteristics[age]"],
        sex=row["Characteristics[sex]"], tissue_source=row["Characteristics[organism part]"],
        cell_enrichment=row["Characteristics[cell type]"], assay="scRNA-seq", technical_fraction="sequencing lanes collapsed",
        possible_overlap="processed RDS contains additional reused GSE136103 cells", included="healthy baseline only",
        exclusion_reason="no disease donors", notes="Only the six donor-labelled cells are unique; never count embedded GSE136103 cells again."))

# GSE202379: preserve every sample/region, with donor as biological unit.
for row in geo_samples("GSE202379"):
    chars = row["chars"]
    status = chars.get("disease status", "")
    saf = chars.get("saf score", "")
    fmatch = re.search(r"F(\d)", saf)
    stage = f"F{fmatch.group(1)}" if fmatch else ("end-stage" if "end stage" in status.lower() else "unresolved")
    donors.append(dict(dataset_id="GSE202379", donor_id=chars.get("patient id", "UNRESOLVED"), sample_id=row["sample_id"],
        disease_group=status, etiology="metabolic" if ("NASH" in status or "NAFLD" in status) else "unresolved",
        fibrosis_stage=stage, compensated_or_decompensated="not reported", age=chars.get("age", ""),
        sex=chars.get("gender", ""), tissue_source="liver; " + chars.get("liver lobe", "region not reported"),
        cell_enrichment="nuclei", assay="snRNA-seq", technical_fraction="lobe/region",
        possible_overlap="same donor has multiple lobes where donor_id repeats", included="yes; regions must be aggregated/nested",
        exclusion_reason="", notes="Donor, not lobe, is the biological unit."))

# Fixed Watson six-donor HCA subset.
watson = [
    ("AM042", "GSM6416567/Liver-13", "healthy", "none", "51-60", "F"),
    ("AM061", "GSM6416578/Liver-18", "healthy", "none", "81-90", "F"),
    ("AM048", "GSM6416569/Liver-14", "healthy", "none", "61-70", "M"),
    ("AM031", "GSM8493744/Liver-32", "fibrosis", "F4", "51-60", "M"),
    ("AM062", "GSM8493745/Liver-33", "fibrosis", "F2", "41-50", "F"),
    ("AM072", "GSM8493746/Liver-34", "fibrosis", "F3", "61-70", "M"),
]
for donor, sample, disease, stage, age, sex in watson:
    donors.append(dict(dataset_id="GSE210077_Watson6", donor_id=donor, sample_id=sample,
        disease_group=disease, etiology="unreported", fibrosis_stage=stage,
        compensated_or_decompensated="not reported", age=age, sex=sex, tissue_source="liver",
        cell_enrichment="matched snRNA-seq and MERFISH", assay="snRNA-seq + MERFISH", technical_fraction="matched assays",
        possible_overlap="same donor contributes both assays by design", included="yes",
        exclusion_reason="", notes="HCA metadata verified mapping; assay observations must not be counted as separate donors."))

# GSE290642 human entries only; pig samples explicitly excluded from the manifest scope.
for row in geo_samples("GSE290642"):
    title, chars = str(row["title"]), row["chars"]
    if not (title.startswith("Fibrosis_") or title.startswith("Control_")):
        continue
    stage = chars.get("fibrosis stage", "")
    donors.append(dict(dataset_id="GSE290642_human", donor_id=title, sample_id=row["sample_id"],
        disease_group="control" if title.startswith("Control") else "fibrosis", etiology="unresolved",
        fibrosis_stage=stage, compensated_or_decompensated="not reported", age=chars.get("age", ""),
        sex=chars.get("gender", ""), tissue_source="liver", cell_enrichment="NPC", assay="scRNA-seq",
        technical_fraction="none", possible_overlap="none identified", included="yes; F1-F3 sensitivity, F4 cirrhosis subset",
        exclusion_reason="", notes="GEO sample is treated as donor; etiologic transfer cannot use this cohort until etiology is recovered."))

# GSE244832 donor IDs are real but condition mapping is not exposed on sample pages.
for row in geo_samples("GSE244832"):
    title = str(row["title"])
    if "snRNAseq" not in title:
        continue
    donor_id = title.split(",", 1)[0]
    donors.append(dict(dataset_id="GSE244832", donor_id=donor_id, sample_id=row["sample_id"],
        disease_group="pending processed cell metadata", etiology="pending", fibrosis_stage="pending",
        compensated_or_decompensated="not reported", age="pending", sex="pending", tissue_source="liver",
        cell_enrichment="nuclei", assay="snRNA-seq", technical_fraction="paired snATAC exists; not independent",
        possible_overlap="same donor has paired snATAC", included="pending condition-map recovery",
        exclusion_reason="not eligible for formal effect until processed metadata is ingested",
        notes="Publication reports 18 donors across normal/MASL/MASH/MetALD; no condition was guessed."))

write_csv(META / "donor_manifest.csv", donors, donor_fields)


# Dataset/lineage evidence is a Phase 0 feasibility statement, not a cell-count result.
lineages = ["macrophage_monocyte", "endothelial", "mesenchymal_hsc_myofibroblast"]
coverage_rows: list[dict[str, object]] = []
coverage_spec = {
    "GSE136103": {"macrophage_monocyte": (5, 5, "author counts/annotations"), "endothelial": (4, 3, "author Fig.4"), "mesenchymal_hsc_myofibroblast": (4, 3, "author Fig.5")},
    "GSE181483_human": {x: (2, 2, "NPC atlas; donor-state counts pending") for x in lineages},
    "GSE202379": {x: (4, 9, "whole-liver snRNA; ≥9 advanced/cirrhosis donors; per-state ≥30 audit pending") for x in lineages},
    "GSE210077_Watson6": {x: (3, 3, "matched snRNA/MERFISH; lineage labels in object") for x in lineages},
    "GSE244832": {x: (5, 13, "publication-level normal/disease totals; donor-state counts pending") for x in lineages},
    "GSE290642_human": {x: (3, 21, "NPC scRNA; three F4; per-state ≥30 audit pending") for x in lineages},
    "PRJNA833766_new4": {"macrophage_monocyte": (2, 2, "CD45+ only")},
}
for dataset_id, values in coverage_spec.items():
    for lineage in lineages:
        if lineage in values:
            nh, nd, evidence = values[lineage]
            available = "yes"
        else:
            nh, nd, evidence, available = 0, 0, "not biologically covered", "no"
        coverage_rows.append(dict(dataset_id=dataset_id, cell_lineage=lineage, lineage_available=available,
            approximate_healthy_donors=nh, approximate_disease_donors=nd,
            minimum_30_cells_verified="yes" if dataset_id == "GSE136103" else "pending matrix-level audit",
            evidence_basis=evidence, phase0_role="feasibility only"))
write_csv(META / "cell_lineage_coverage.csv", coverage_rows)

# One row per donor/sample/lineage so technical repeats are visible.
sample_lineage: list[dict[str, object]] = []
for row in donors:
    dset = str(row["dataset_id"])
    for lineage in lineages:
        c = next((x for x in coverage_rows if x["dataset_id"] == dset and x["cell_lineage"] == lineage), None)
        if c is None:
            continue
        sample_lineage.append(dict(dataset_id=dset, donor_id=row["donor_id"], sample_id=row["sample_id"],
            assay=row["assay"], cell_lineage=lineage, expected_from_design=c["lineage_available"],
            donor_level_cell_count="not computed in Phase 0", eligible_at_30_cells="pending",
            evidence_basis=c["evidence_basis"]))
write_csv(META / "sample_lineage_matrix.csv", sample_lineage)

assay_rows = []
for d in datasets:
    assay_rows.append(dict(dataset_id=d["dataset_id"], assay=d["assay"], scRNA="yes" if "scRNA" in d["assay"] else "no",
        snRNA="yes" if "snRNA" in d["assay"] else "no", spatial="yes" if d["spatial_data_available"] == "yes" else "no",
        processed_object=d["processed_counts_available"], independent_of_Ramachandran_programs="no" if d["dataset_id"] == "GSE136103" else ("partially" if d["dataset_id"] == "E-MTAB-10553_unique" else "yes"),
        phase0_use=d["primary_validation_candidate"]))
write_csv(META / "assay_coverage.csv", assay_rows)

spatial_rows = [
    dict(dataset_id="GSE210077_Watson6", platform="MERFISH", tissue="human liver", donors=6, healthy_donors=3,
         fibrotic_donors=3, fibrosis_stages="F2;F3;F4", gene_panel_size="317 genes (332 codebook entries including 15 blanks)",
         coordinates_available="yes (x,y)", cell_labels_available="yes", scar_or_fibrosis_annotation="condition and fibrosis stage; septal proximity requires derived geometry",
         matched_single_nucleus="yes, same six donors", processed_object="HCA h5ad 1.322 GB; cell-properties CSV 66.99 MB; gene list 2.3 KB",
         memory_32gb="feasible with backed/chunked h5ad or CSV chunks", smoke_test="schema and bounded-row read passed; full h5ad download deferred",
         limitations="targeted panel; one F4 donor; spatial pixels/cells are not biological replicates"),
    dict(dataset_id="GSE185477", platform="10x Visium", tissue="healthy human liver", donors=1, healthy_donors=1,
         fibrotic_donors=0, fibrosis_stages="none", gene_panel_size="transcriptome-wide",
         coordinates_available="yes", cell_labels_available="spot deconvolution/reference mapping", scar_or_fibrosis_annotation="no disease",
         matched_single_nucleus="different donor", processed_object="GEO processed files", memory_32gb="feasible",
         smoke_test="not selected", limitations="healthy baseline only; cannot support fibrotic localization"),
]
write_csv(META / "spatial_manifest.csv", spatial_rows)


# Frozen programs from exact Ramachandran supplementary tables.
program_rows: list[dict[str, object]] = []
common = dict(paper_id="Ramachandran2019_Nature", discovery_dataset="GSE136103",
              discovery_assay="scRNA-seq", disease_etiology="mixed nonviral cirrhosis",
              source_page="supplementary file", verified="TRUE", data_lineage_class="REUSED")

table9 = pd.read_excel(RAW / "PMC6876711_supplementary" / "EMS84316-supplement-Suppl_Table_9.xlsx",
                       sheet_name="Gene per Signature")
mac_names = {
    "Signature A": ("RAM2019_MAC_SIG_A_SAM", "scar-associated macrophage metagene A", "SCAR_ASSOCIATED_MACROPHAGE"),
    "Signature B": ("RAM2019_MAC_SIG_B_SAM", "scar-associated macrophage metagene B", "SCAR_ASSOCIATED_MACROPHAGE"),
    "Signature C": ("RAM2019_MAC_SIG_C_KC", "Kupffer-cell metagene C", "HOMEOSTATIC_MACROPHAGE"),
    "Signature D": ("RAM2019_MAC_SIG_D_TMO", "tissue-monocyte metagene D", "TISSUE_MONOCYTE"),
    "Signature E": ("RAM2019_MAC_SIG_E_TMO", "tissue-monocyte metagene E", "TISSUE_MONOCYTE"),
    "Signature F": ("RAM2019_MAC_SIG_F_CDC1", "cDC1 metagene F", "DENDRITIC_CELL"),
}
for column, (pid, name, ptype) in mac_names.items():
    for gene in table9[column].dropna().astype(str):
        gene = gene.strip()
        if not gene:
            continue
        program_rows.append(dict(program_id=pid, program_name=name, cell_lineage="macrophage_monocyte",
            gene_symbol=gene, direction="UP_IN_STATE", source_table="Supplementary Table 9",
            notes=f"{ptype}; exact author metagene list; external only outside GSE136103", **common))

for table_number, lineage, mapping in [
    (13, "endothelial", {
        "Endothelia (1)": ("RAM2019_ENDO_1", "endothelial cluster 1 / healthy LSEC marker program", "HEALTHY_LSEC"),
        "Endothelia (2)": ("RAM2019_ENDO_2", "endothelial cluster 2 marker program", "ENDOTHELIAL_STATE"),
        "Endothelia (3)": ("RAM2019_ENDO_3", "endothelial cluster 3 marker program", "ENDOTHELIAL_STATE"),
        "Endothelia (4)": ("RAM2019_ENDO_4", "endothelial cluster 4 marker program", "ENDOTHELIAL_STATE"),
        "Endothelia (5)": ("RAM2019_ENDO_5", "endothelial cluster 5 marker program", "ENDOTHELIAL_STATE"),
        "Endothelia (6)": ("RAM2019_ENDO_6_SAENDO1", "scar-associated endothelial 1 marker program", "SCAR_ENDOTHELIAL_PLVAP"),
        "Endothelia (7)": ("RAM2019_ENDO_7_SAENDO2", "scar-associated endothelial 2 marker program", "SCAR_ENDOTHELIAL_ACKR1"),
    }),
    (16, "mesenchymal_hsc_myofibroblast", {
        "Mesenchyme (VSMC)": ("RAM2019_MES_VSMC", "vascular smooth-muscle marker program", "VSMC"),
        "Mesenchyme (HSC)": ("RAM2019_MES_HSC", "hepatic stellate-cell marker program", "QUIESCENT_HSC"),
        "(Myo) Fibroblast": ("RAM2019_MES_SAMES", "scar-associated myofibroblast marker program", "ACTIVATED_HSC_MYOFIBROBLAST"),
        "Mesothelia": ("RAM2019_MES_MESOTHELIAL", "mesothelial marker program", "MESOTHELIAL"),
    }),
    (17, "mesenchymal_hsc_myofibroblast", {
        "SAMes (A)": ("RAM2019_SAMES_A", "scar-associated mesenchymal metagene A", "ACTIVATED_HSC_MYOFIBROBLAST"),
        "SAMes (B)": ("RAM2019_SAMES_B", "scar-associated mesenchymal metagene B", "ACTIVATED_HSC_MYOFIBROBLAST"),
    }),
]:
    data = pd.read_csv(RAW / "PMC6876711_supplementary" / f"EMS84316-supplement-Suppl_Table_{table_number}.csv")
    for cluster, (pid, name, ptype) in mapping.items():
        genes = data.loc[data["cluster"] == cluster, "gene"].dropna().astype(str)
        for gene in genes:
            program_rows.append(dict(program_id=pid, program_name=name, cell_lineage=lineage,
                gene_symbol=gene.strip(), direction="UP_IN_STATE", source_table=f"Supplementary Table {table_number}",
                notes=f"{ptype}; positive cluster marker (avg_logFC>0); external only outside GSE136103", **common))

program_fields = ["program_id", "paper_id", "program_name", "cell_lineage", "gene_symbol", "direction",
                  "discovery_dataset", "discovery_assay", "disease_etiology", "source_table", "source_page",
                  "verified", "data_lineage_class", "notes"]
write_csv(LIT / "program_inventory.csv", program_rows, program_fields)

program_summary = []
for pid, group in pd.DataFrame(program_rows).groupby("program_id", sort=False):
    program_summary.append(dict(program_id=pid, program_name=group.iloc[0]["program_name"],
        cell_lineage=group.iloc[0]["cell_lineage"], n_genes=group["gene_symbol"].nunique(),
        discovery_dataset="GSE136103", gse136103_lineage="REUSED/internal only",
        gse181483_lineage="INDEPENDENT; small", gse202379_lineage="INDEPENDENT snRNA candidate",
        gse210077_watson_lineage="INDEPENDENT snRNA/spatial candidate", gse244832_lineage="INDEPENDENT except own hA1 program",
        gse290642_lineage="INDEPENDENT fibrosis-spectrum candidate"))
write_csv(LIT / "program_lineage_matrix.csv", program_summary)

print(f"datasets={len(datasets)} donor/sample_rows={len(donors)} programs={len(program_summary)} program_gene_rows={len(program_rows)}")
