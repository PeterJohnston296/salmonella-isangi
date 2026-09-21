agents <- native$ISANGI_AGENTS
labels <- c("AMP", "AZM", "CIP", "FOX", "MEM", "CHL", "PEF", "FEP", "SXT", "TGC", "IPM", "AMK", "CPD", "ATM", "GEN", "TZP", "CTX")
ids <- d$phenotype$isolate_id[d$phenotype$group == "Isangi"]
ph <- do.call(rbind, lapply(agents, function(a) data.frame(isolate_id = ids,
  run_accession = d$meta[ids, "run_accession"], agent = a,
  category = d$phenotype[ids, a], zone_mm = d$zones[ids, a])))
ph$row_label <- paste(ph$isolate_id, ph$run_accession, sep = " | ")
row_labels <- paste(ids, d$meta[ids, "run_accession"], sep = " | ")
ph$row_label <- factor(ph$row_label, levels = rev(row_labels))
ph$agent <- factor(ph$agent, levels = agents)
p <- ggplot2::ggplot(ph, ggplot2::aes(agent, row_label, fill = category)) +
  ggplot2::geom_tile(width = .97, height = .96, colour = "white", linewidth = .48) +
  ggplot2::geom_text(ggplot2::aes(label = zone_mm), size = 2.55, colour = "#20252A") +
  ggplot2::scale_fill_manual(values = c(S = "#DCE9D8", I = "#E8D59C", R = "#B94E5D"),
    breaks = c("S", "I", "R"), labels = c(S = "Susceptible", I = "Increased exposure", R = "Resistant"), name = NULL) +
  ggplot2::scale_x_discrete(labels = setNames(labels, agents), expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(x = "Antimicrobial agent", y = NULL) +
  ggplot2::theme_minimal(base_size = 8.5) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 7.3, colour = "#30353B"),
    axis.text.y = ggplot2::element_text(size = 6.75, colour = "#30353B"),
    legend.position = "top", legend.justification = "left", legend.text = ggplot2::element_text(size = 7.4),
    plot.margin = ggplot2::margin(3, 4, 4, 4))
save_plot(p, "Figure_2_phenotypic_susceptibility", 190, 125)
write.csv(ph, file.path(OUT, "source_data", "Figure_2_source_data.csv"), row.names = FALSE)
