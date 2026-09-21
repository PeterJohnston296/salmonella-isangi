native$export_analysis(d, cfg, OUT)
m <- d$meta
for (nm in c("cc25", "st335")) {
  tr <- if (nm == "cc25") d$cc_tree else d$st_tree
  write.csv(data.frame(position = seq_along(tr$ids), m[tr$ids, c("isolate_id", "run_accession", "sample_accession", "sequence_type", "lineage")]),
    file.path(OUT, "source_data", paste0(nm, "_tree_accession_index.csv")), row.names = FALSE)
}
write.csv(m[m$phenotyped_Kenneth == 1, c("isolate_id", "run_accession", "sample_accession", "sequence_type", "CC25_panel", "ST335_panel")],
  file.path(OUT, "source_data", "phenotype_tree_linkage.csv"), row.names = FALSE)
