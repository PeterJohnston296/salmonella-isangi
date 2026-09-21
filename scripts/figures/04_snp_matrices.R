cfg$matrix$max_cc25 <- max(d$full)
cfg$matrix$max_st335 <- max(d$st)
save_native(native$build_snp(d, cfg, "CC25"), "Figure_S20_CC25_tree_aligned_327")
save_native(native$build_snp(d, cfg, "ST335"), "Figure_S21_ST335_tree_aligned_224")
save_native(native$build_snp(d, cfg, "CC25", complete = TRUE), "CC25_complete_331_SNP_matrix")
for (pid in unique(d$panels$panel_id)) {
  pop <- if (startsWith(pid, "CC25")) "CC25" else "ST335"
  save_native(native$build_snp(d, cfg, pop, pid), paste0(pid, "_SNP_detail"))
}
local_matrix <- function(path, stem, size_mm, numbers = FALSE, clustered = FALSE) {
  M <- native$read_matrix(path)
  if (clustered) {
    h <- hclust(as.dist(M), method = "average")
    M <- M[h$order, h$order, drop = FALSE]
  }
  ids <- rownames(M)
  annotations <- d$meta[ids, c("country", "year", "source_group", "lineage"), drop = FALSE]
  names(annotations) <- c("Country", "Year", "Source", "Lineage")
  annotations[] <- lapply(annotations, function(x) {x <- as.character(x); x[is.na(x) | x == ""] <- "Other"; x})
  annotation_colours <- list(Country = cfg$palette$country, Source = cfg$palette$source, Lineage = cfg$palette$lineage)
  annotation_colours$Year <- setNames(grDevices::hcl.colors(length(unique(annotations$Year)), "Blues 3"), sort(unique(annotations$Year)))
  ph <- pheatmap::pheatmap(M, cluster_rows = FALSE, cluster_cols = FALSE,
    color = grDevices::colorRampPalette(c("#FFFFFF", "#D7E8F3", "#7FAACB", "#2F6997", "#173C5A"))(256),
    breaks = seq(0, max(M), length.out = 257), border_color = "#E7EBEF",
    annotation_col = if (clustered) NULL else annotations,
    annotation_colors = if (clustered) NULL else annotation_colours,
    display_numbers = numbers, number_format = "%.0f", fontsize_number = 5,
    fontsize = 8, fontsize_row = 6, fontsize_col = 6, angle_col = "90",
    legend_breaks = unique(round(seq(0, max(M), length.out = 6))),
    main = NA, silent = TRUE)
  save_drawing(function(w,h) grid::grid.draw(ph$gtable), stem, size_mm, size_mm, "supplementary")
  native$export_matrix(M, stem, file.path(OUT, "matrices"))
  write.csv(data.frame(position = seq_along(ids), isolate_id = ids, run_accession = d$meta[ids, "run_accession"]),
    file.path(OUT, "source_data", paste0(stem, "_index.csv")), row.names = FALSE)
}
local_matrix(file.path(DATA, "matrices", "Malawi_ST335_74.csv"), "Figure_S23_Malawi_ST335_74", 290)
local_matrix(file.path(DATA, "matrices", "South_Africa_outbreak_ST335_37.csv"), "Figure_S24_South_Africa_ST335_37", 245, TRUE, TRUE)
save_native(native$build_phenotype(d, cfg, "Comparator"), "Figure_S22_comparator_phenotypes")
