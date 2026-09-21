matrices <- list(CC25 = d$cc, ST335 = d$st, Malawi = d$malawi,
  South_Africa = native$read_matrix(file.path(DATA, "matrices", "South_Africa_outbreak_ST335_37.csv")))
rows <- list()
for (nm in names(matrices)) {
  M <- matrices[[nm]]
  native$export_matrix(M, nm, file.path(OUT, "matrices"), long = FALSE)
  ij <- which(upper.tri(M), arr.ind = TRUE)
  v <- M[ij]
  groups <- d$meta[rownames(M), "source_group"]
  pairs <- data.frame(isolate_a = rownames(M)[ij[, 1]], isolate_b = colnames(M)[ij[, 2]],
    source_a = groups[ij[, 1]], source_b = groups[ij[, 2]], SNP_distance = as.integer(v))
  write.csv(pairs, file.path(OUT, "source_data", paste0(nm, "_unique_pairs.csv")), row.names = FALSE)
  diag(M) <- Inf
  nearest <- apply(M, 1, min)
  ties <- vapply(seq_len(nrow(M)), function(i) paste(colnames(M)[M[i, ] == nearest[i]], collapse = ";"), character(1))
  write.csv(data.frame(isolate_id = rownames(M), nearest_SNP_distance = nearest, nearest_isolates = ties),
    file.path(OUT, "source_data", paste0(nm, "_nearest_neighbours.csv")), row.names = FALSE)
  q <- quantile(v, c(.25, .5, .75), names = FALSE)
  rows[[nm]] <- data.frame(population = nm, isolates = nrow(M), pairs = length(v),
    minimum = min(v), q1 = q[1], median = q[2], q3 = q[3], maximum = max(v),
    pairs_0_SNP = sum(v == 0), pairs_le_2 = sum(v <= 2), pairs_le_5 = sum(v <= 5), pairs_le_10 = sum(v <= 10))
}
write.csv(do.call(rbind, rows), file.path(OUT, "source_data", "SNP_population_summary.csv"), row.names = FALSE)
