packages <- c("ape", "phangorn", "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "here", "ragg", "cowplot", "magick", "jsonlite", "png", "digest", "igraph", "knitr", "rmarkdown", "scales", "pheatmap")
if ("--map" %in% commandArgs(trailingOnly = TRUE)) packages <- c(packages, "sf", "ggrepel")
needed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(needed)) install.packages(needed, repos = "https://cloud.r-project.org")
