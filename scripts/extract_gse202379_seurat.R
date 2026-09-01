#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: extract_gse202379_seurat.R <input.rds.gz> <output_dir> <r_library>",
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
r_library <- normalizePath(args[[3L]], mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(r_library, .libPaths()))

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratObject))

message("Reading author Seurat object: ", input_path)
connection <- gzfile(input_path, open = "rb")
on.exit(close(connection), add = TRUE)
object <- readRDS(connection)
close(connection)
on.exit(NULL, add = FALSE)

metadata <- object[[]]
metadata$cell_id <- rownames(metadata)
metadata <- metadata[, c("cell_id", setdiff(names(metadata), "cell_id")), drop = FALSE]

metadata_path <- file.path(output_dir, "cell_metadata.csv.gz")
metadata_connection <- gzfile(metadata_path, open = "wt")
write.csv(metadata, metadata_connection, row.names = FALSE, na = "")
close(metadata_connection)

required_columns <- c(
  "Patient.ID",
  "cell.annotation",
  "Disease.status",
  "Fibrosis.score..F0.4."
)
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "author object is missing required metadata columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

lineage_map <- c(
  "Macrophages" = "macrophage_monocyte",
  "Endothelial" = "endothelial",
  "Stellate" = "mesenchymal_hsc_myofibroblast"
)
metadata$harmonized_lineage <- unname(lineage_map[as.character(metadata$cell.annotation)])
target_cells <- !is.na(metadata$harmonized_lineage)
if (!any(target_cells)) {
  stop("no cells map to the three frozen target lineages", call. = FALSE)
}

counts <- tryCatch(
  LayerData(object, assay = "RNA", layer = "counts"),
  error = function(error) {
    GetAssayData(object, assay = "RNA", slot = "counts")
  }
)
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}
if (ncol(counts) != nrow(metadata)) {
  stop("RNA count columns do not match metadata rows", call. = FALSE)
}
if (!identical(colnames(counts), metadata$cell_id)) {
  stop("RNA count column order does not match metadata cell order", call. = FALSE)
}

group_key <- paste(
  metadata$Patient.ID[target_cells],
  metadata$harmonized_lineage[target_cells],
  sep = "|||"
)
group_levels <- sort(unique(group_key))
group_index <- match(group_key, group_levels)
design <- sparseMatrix(
  i = which(target_cells),
  j = group_index,
  x = 1,
  dims = c(nrow(metadata), length(group_levels))
)
pseudobulk <- counts %*% design
group_ids <- sprintf("G%04d", seq_along(group_levels))
colnames(pseudobulk) <- group_ids

writeMM(pseudobulk, file.path(output_dir, "donor_lineage_raw_counts.mtx"))
write.csv(
  data.frame(gene = rownames(pseudobulk), stringsAsFactors = FALSE),
  file.path(output_dir, "genes.csv"),
  row.names = FALSE,
  na = ""
)

audit_columns <- intersect(
  c(
    "orig.ident",
    "Patient.ID",
    "manuscript.expt",
    "Disease.status",
    "Fibrosis.score..F0.4.",
    "Steatosis",
    "Ballooning",
    "Inflammation",
    "Age",
    "Sex",
    "Gender",
    "Etiology",
    "etiology"
  ),
  names(metadata)
)
manifest_rows <- lapply(seq_along(group_levels), function(index) {
  selected <- target_cells
  selected[target_cells] <- group_index == index
  parts <- strsplit(group_levels[[index]], "\\|\\|\\|")[[1L]]
  row <- list(
    group_id = group_ids[[index]],
    donor_id = parts[[1L]],
    harmonized_lineage = parts[[2L]],
    n_cells = sum(selected)
  )
  for (column in audit_columns) {
    values <- sort(unique(as.character(metadata[selected, column])))
    values <- values[!is.na(values) & nzchar(values)]
    row[[column]] <- paste(values, collapse = ";")
  }
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
})
group_manifest <- do.call(rbind, manifest_rows)
write.csv(
  group_manifest,
  file.path(output_dir, "donor_lineage_manifest.csv"),
  row.names = FALSE,
  na = ""
)

author_annotation_counts <- as.data.frame(
  table(
    donor_id = metadata$Patient.ID,
    author_annotation = metadata$cell.annotation,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
author_annotation_counts <- author_annotation_counts[author_annotation_counts$Freq > 0L, ]
write.csv(
  author_annotation_counts,
  file.path(output_dir, "donor_author_annotation_cell_counts.csv"),
  row.names = FALSE,
  na = ""
)

inventory <- c(
  paste("r_version", R.version.string, sep = "\t"),
  paste("SeuratObject_version", as.character(packageVersion("SeuratObject")), sep = "\t"),
  paste("Matrix_version", as.character(packageVersion("Matrix")), sep = "\t"),
  paste("object_class", paste(class(object), collapse = ";"), sep = "\t"),
  paste("object_version", as.character(object@version), sep = "\t"),
  paste("cells", nrow(metadata), sep = "\t"),
  paste("genes", nrow(counts), sep = "\t"),
  paste("rna_count_nonzero", length(counts@x), sep = "\t"),
  paste("target_cells", sum(target_cells), sep = "\t"),
  paste("donors", length(unique(metadata$Patient.ID)), sep = "\t"),
  paste("target_donor_lineage_groups", ncol(pseudobulk), sep = "\t"),
  paste("metadata_columns", paste(names(metadata), collapse = ";"), sep = "\t"),
  paste("author_annotations", paste(sort(unique(metadata$cell.annotation)), collapse = ";"), sep = "\t"),
  paste("disease_statuses", paste(sort(unique(metadata$Disease.status)), collapse = ";"), sep = "\t"),
  paste("fibrosis_scores", paste(sort(unique(metadata$Fibrosis.score..F0.4.)), collapse = ";"), sep = "\t")
)
writeLines(inventory, file.path(output_dir, "object_inventory.tsv"), useBytes = TRUE)

message(
  "Exported ", nrow(metadata), " cells, ", nrow(counts), " genes, and ",
  ncol(pseudobulk), " donor-lineage pseudobulks."
)
