paths <- c(list.files("scripts", pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("environment", pattern = "[.]R$", recursive = TRUE, full.names = TRUE))
for (path in paths) invisible(parse(file = path))
source("scripts/lib/project.R")
stopifnot(all(d$checks), nrow(d$meta) == 345, length(d$cc_tree$ids) == 327,
  length(d$st_tree$ids) == 224, nrow(d$phenotype) == 42)
cat("R parsing and frozen-data checks passed.\n")
