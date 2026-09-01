from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    manifest_path = repo / "metadata" / "donor_manifest.csv"
    source_path = (
        repo
        / "data"
        / "raw"
        / "GSE244832"
        / "processed_files"
        / "hLIVER_metadata.csv"
    )
    manifest = pd.read_csv(manifest_path, keep_default_na=False)
    source = pd.read_csv(source_path, index_col=0)
    donor_condition = (
        source[["orig.ident", "condition"]]
        .drop_duplicates()
        .sort_values("orig.ident")
    )
    if donor_condition["orig.ident"].duplicated().any() or len(donor_condition) != 18:
        raise RuntimeError("expected one condition for each of 18 GSE244832 donors")
    condition_map = dict(zip(donor_condition["orig.ident"], donor_condition["condition"]))
    disease_map = {"NORMAL": "normal", "NAFL": "MASL", "NASH": "MASH"}

    target = manifest["dataset_id"].eq("GSE244832")
    if target.sum() != 18:
        raise RuntimeError("donor manifest does not contain 18 GSE244832 rows")
    for index in manifest.index[target]:
        author_donor = manifest.at[index, "donor_id"].replace("JB_", "JB")
        condition = condition_map[author_donor]
        disease_group = disease_map[condition]
        manifest.at[index, "disease_group"] = disease_group
        manifest.at[index, "etiology"] = "none" if condition == "NORMAL" else "metabolic"
        manifest.at[index, "fibrosis_stage"] = (
            "group-level F2-F4 only" if condition == "NASH" else "not reported"
        )
        manifest.at[index, "age"] = "not reported"
        manifest.at[index, "sex"] = "not reported"
        manifest.at[index, "included"] = "sensitivity only"
        manifest.at[index, "exclusion_reason"] = (
            "no donor-level histologic stage; MASH group spans F2-F4"
        )
        manifest.at[index, "notes"] = (
            "Condition recovered directly from author hLIVER_metadata.csv; "
            "NORMAL/MASL/MASH effects must not be relabelled cirrhosis or F3-F4."
        )

    manifest.to_csv(manifest_path, index=False)
    observed = (
        manifest[target]
        .groupby("disease_group")["donor_id"]
        .nunique()
        .to_dict()
    )
    expected = {"MASH": 9, "MASL": 4, "normal": 5}
    if observed != expected:
        raise RuntimeError(f"condition count mismatch: {observed} != {expected}")
    print(observed)


if __name__ == "__main__":
    main()
