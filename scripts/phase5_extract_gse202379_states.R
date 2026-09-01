#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: phase5_extract_gse202379_states.R <input.rds> <output_dir> <r_library>",
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
r_library <- normalizePath(args[[3L]], mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(r_library, .libPaths()))

suppressPackageStartupMessages(library(Matrix))

object <- readRDS(input_path)
object_attributes <- attributes(object)
metadata <- object_attributes[["meta.data"]]
if (!is.data.frame(metadata)) {
  stop("serialized Seurat metadata slot is not a data.frame", call. = FALSE)
}
metadata$cell_id <- rownames(metadata)

required <- c(
  "Patient.ID", "cell.annotation", "Disease.status",
  "Fibrosis.score..F0.4.", "SCT_snn_res.0.8"
)
missing <- setdiff(required, names(metadata))
if (length(missing) > 0L) {
  stop("missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

lineage_map <- c(
  "Macrophages" = "macrophage_monocyte",
  "Endothelial" = "endothelial",
  "Stellate" = "mesenchymal_hsc_myofibroblast"
)
metadata$harmonized_lineage <- unname(lineage_map[as.character(metadata$cell.annotation)])
target <- !is.na(metadata$harmonized_lineage)
metadata$source_state <- paste0("cluster_", as.character(metadata$SCT_snn_res.0.8))

assays <- object_attributes[["assays"]]
if (is.null(assays[["RNA"]])) {
  stop("serialized Seurat object has no RNA assay", call. = FALSE)
}
counts <- attributes(assays[["RNA"]])[["counts"]]
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
if (!identical(colnames(counts), metadata$cell_id)) {
  stop("RNA count columns do not match metadata cell order", call. = FALSE)
}

key <- paste(
  metadata$Patient.ID[target],
  metadata$harmonized_lineage[target],
  metadata$source_state[target],
  sep = "|||"
)
levels <- sort(unique(key))
index <- match(key, levels)
design <- sparseMatrix(
  i = which(target), j = index, x = 1,
  dims = c(nrow(metadata), length(levels))
)
aggregated <- counts %*% design
group_ids <- sprintf("S%05d", seq_along(levels))
colnames(aggregated) <- group_ids

writeMM(aggregated, file.path(output_dir, "phase5_donor_state_raw_counts.mtx"))
write.csv(
  data.frame(gene = rownames(aggregated), stringsAsFactors = FALSE),
  file.path(output_dir, "phase5_state_genes.csv"),
  row.names = FALSE,
  na = ""
)

audit_columns <- intersect(
  c(
    "orig.ident", "Patient.ID", "manuscript.expt", "Disease.status",
    "Fibrosis.score..F0.4.", "Steatosis", "Ballooning", "Inflammation",
    "Age", "Sex", "Gender", "Etiology", "etiology"
  ),
  names(metadata)
)

rows <- lapply(seq_along(levels), function(i) {
  selected <- target
  selected[target] <- index == i
  parts <- strsplit(levels[[i]], "\\|\\|\\|")[[1L]]
  row <- list(
    state_group_id = group_ids[[i]],
    donor_id = parts[[1L]],
    harmonized_lineage = parts[[2L]],
    source_state = parts[[3L]],
    n_cells = sum(selected),
    library_size = sum(aggregated[, i])
  )
  for (column in audit_columns) {
    values <- sort(unique(as.character(metadata[selected, column])))
    values <- values[!is.na(values) & nzchar(values)]
    row[[column]] <- paste(values, collapse = ";")
  }
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
})
manifest <- do.call(rbind, rows)
write.csv(
  manifest,
  file.path(output_dir, "phase5_donor_state_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Exported ", ncol(aggregated), " donor-state pseudobulks from ",
  sum(target), " target cells and ", nrow(aggregated), " genes."
)
