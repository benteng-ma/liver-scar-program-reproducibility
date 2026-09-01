# Data access and redistribution boundary

The analyses use public human liver datasets. Accession, publication, assay, donor-count, endpoint, annotation, and limitation fields are recorded in `metadata/dataset_manifest.csv`.

Raw matrices and third-party processed objects are not included in this repository. Retrieve them from GEO, ArrayExpress/BioStudies, Dryad, or the linked study repository using the source accession and citation. Preserve downloaded files as read-only inputs under `data/raw/` or `data/external/`; generated intermediates belong under `data/interim/` and `data/processed/`. These directories are ignored by Git.

Do not infer biological zero from missing values in derived tables. Missing fields indicate not applicable, unavailable, or not estimated according to the table schema. Public donor identifiers are retained only as source-study codes; no direct identifiers or re-identification keys are included.
