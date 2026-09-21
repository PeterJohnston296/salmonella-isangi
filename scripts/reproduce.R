args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("salmonella-isangi.Rproj")) stop("Run from the repository root.", call. = FALSE)
source("scripts/lib/project.R")
stages <- c(cohorts = "scripts/analysis/01_cohorts.R",
  snps = "scripts/analysis/02_snp_distances.R", epidemiology = "scripts/figures/01_epidemiology.R",
  phenotype = "scripts/figures/02_susceptibility.R", phylogenies = "scripts/figures/03_phylogenies.R",
  burden = "scripts/analysis/03_amr_burden.R", matrices = "scripts/figures/04_snp_matrices.R",
  artwork = "scripts/figures/05_retained_artwork.R", composite = "scripts/figures/06_composite.R")
selected <- if (length(args)) strsplit(args[1], ",", fixed = TRUE)[[1]] else names(stages)
if (!all(selected %in% names(stages))) stop("Stages: ", paste(names(stages), collapse = ", "))
status <- list()
for (stage in selected) {
  message(stage)
  result <- tryCatch({sys.source(stages[[stage]], envir = globalenv()); "complete"},
    error = function(e) conditionMessage(e))
  status[[stage]] <- data.frame(stage = stage, result = result)
}
status <- do.call(rbind, status)
write.csv(status, file.path(OUT, "audit", "stage_results.csv"), row.names = FALSE)
source("environment/capture.R")
writeLines(OUT, file.path(ROOT, "results", "LATEST_RUN.txt"))
print(status, row.names = FALSE)
if (any(status$result != "complete")) quit(status = 1L)
