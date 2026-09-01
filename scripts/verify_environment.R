required <- c("renv", "yaml", "jsonlite", "testthat")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
lines <- c(
  paste("verified_at", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste("r_version", R.version.string),
  vapply(required, function(pkg) paste(pkg, as.character(utils::packageVersion(pkg))), character(1)),
  "phase0_disease_effects FALSE"
)
writeLines(lines, "results/logs/environment_verification.txt")
cat(paste(lines, collapse = "\n"), "\n")

