needed <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "here", "ragg")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install R dependencies with Rscript environment/install.R.")

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(here)

STYLE <- list(
  width_mm = 125,
  height_mm = 110,
  font_family = "sans",
  base_size = 9,
  violin_bandwidth = 0.65,
  violin_width = 0.74,
  violin_alpha = 0.22,
  point_size = 0.85,
  point_alpha = 0.48,
  point_jitter = 0.15,
  seed = 335,
  show_study_points = TRUE,
  show_source_labels = TRUE
)

group_levels <- c("NCBI non-Isangi", "Isangi ST335", "Other Isangi")
group_colours <- c(
  "NCBI non-Isangi" = "#75818A",
  "Isangi ST335" = "#B74F60",
  "Other Isangi" = "#4B81A6"
)

axis_labels <- expression(
  atop("Non-Isangi", italic(Salmonella)),
  atop(italic(S.) ~ "Isangi", "ST335"),
  atop("Other", italic(S.) ~ "Isangi")
)

outdir <- file.path(OUT, "main")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

knitr::opts_chunk$set(
  dev = "ragg_png",
  fig.width = STYLE$width_mm / 25.4,
  fig.height = STYLE$height_mm / 25.4,
  dpi = 180
)

input_paths <- c(
  ncbi = here("data", "data_raw", "2026.01.13.ncbi_pathogen_salmonella_amr_serovar.tsv.tsv.gz"),
  amr = here("data", "data_raw", "all_isangi_assemblies_amrfinder2_20062024.csv"),
  lookup = here("data", "data_raw", "isangi_st_lookup.csv")
)

if (any(!file.exists(input_paths))) {
  stop("Missing input file(s):\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"))
}

ncbi <- readr::read_tsv(input_paths["ncbi"], show_col_types = FALSE, progress = FALSE)
isangi_amr <- readr::read_csv(input_paths["amr"], show_col_types = FALSE, progress = FALSE)
isangi_lookup <- readr::read_csv(
  input_paths["lookup"], show_col_types = FALSE, progress = FALSE,
  col_types = cols(.default = col_character())
)

require_columns <- function(x, required, label) {
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "))
}
require_columns(ncbi, c("#AMR genotypes", "Computed types", "Run"), "NCBI table")
require_columns(isangi_amr, c("Isolate_ID", "Element type", "Element subtype", "Gene symbol"), "AMRFinder table")
require_columns(isangi_lookup, c("isolate_id", "run_accession", "ST"), "ST lookup")

stopifnot(
  nrow(isangi_lookup) == 345L,
  !anyNA(isangi_lookup$isolate_id),
  !anyDuplicated(isangi_lookup$isolate_id),
  !anyNA(isangi_lookup$ST),
  sum(isangi_lookup$ST == "335") == 224L
)

ncbi_non_isangi <- ncbi %>%
  filter(!str_detect(coalesce(`Computed types`, ""), regex("serotype=Isangi", ignore_case = TRUE))) %>%
  mutate(
    amr_string = coalesce(`#AMR genotypes`, ""),
    amr_gene_count =
      str_count(amr_string, fixed("=COMPLETE")) -
      as.integer(str_detect(amr_string, "(^|,)mdsA=COMPLETE(,|$)")) -
      as.integer(str_detect(amr_string, "(^|,)mdsB=COMPLETE(,|$)")),
    amr_gene_count = pmax(amr_gene_count, 0L)
  ) %>%
  transmute(isolate_id = Run, group = "NCBI non-Isangi", amr_gene_count)

isangi_gene_counts <- isangi_amr %>%
  mutate(isolate_id = str_remove(Isolate_ID, "_assembly_contigs$")) %>%
  filter(
    `Element type` == "AMR", `Element subtype` == "AMR",
    !is.na(`Gene symbol`), `Gene symbol` != ""
  ) %>%
  distinct(isolate_id, `Gene symbol`) %>%
  count(isolate_id, name = "amr_gene_count")

amr_ids <- str_remove(unique(isangi_amr$Isolate_ID), "_assembly_contigs$")
unmatched_ids <- setdiff(amr_ids, isangi_lookup$isolate_id)
if (length(unmatched_ids)) {
  stop("AMRFinder IDs absent from the lookup: ", paste(unmatched_ids, collapse = ", "))
}

isangi_counts <- isangi_lookup %>%
  left_join(isangi_gene_counts, by = "isolate_id") %>%
  mutate(
    count_before_zero_fill = amr_gene_count,
    amr_gene_count = replace_na(amr_gene_count, 0L),
    present_in_feature_file = isolate_id %in% amr_ids,
    ST_assigned = grepl("^[0-9]+$", ST),
    group = if_else(ST == "335", "Isangi ST335", "Other Isangi")
  )

plot_data <- bind_rows(
  ncbi_non_isangi,
  isangi_counts %>% select(isolate_id, group, amr_gene_count)
) %>%
  mutate(group = factor(group, levels = group_levels))

stopifnot(
  sum(isangi_counts$group == "Other Isangi") == 121L,
  all(is.finite(plot_data$amr_gene_count)),
  all(plot_data$amr_gene_count >= 0),
  all(plot_data$amr_gene_count == round(plot_data$amr_gene_count)),
  !anyNA(plot_data$group)
)

study_points <- plot_data %>% filter(group != "NCBI non-Isangi")

group_summary <- plot_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_amr_genes = median(amr_gene_count),
    q1 = unname(quantile(amr_gene_count, 0.25, type = 7)),
    q3 = unname(quantile(amr_gene_count, 0.75, type = 7)),
    IQR = IQR(amr_gene_count),
    min = min(amr_gene_count),
    max = max(amr_gene_count),
    zero_count = sum(amr_gene_count == 0),
    .groups = "drop"
  ) %>%
  mutate(source = if_else(group == "NCBI non-Isangi", "NCBI comparator", "Study collection"))

group_summary

isangi_counts %>% count(group, ST, name = "n")

valid_runs <- ncbi_non_isangi$isolate_id[
  !is.na(ncbi_non_isangi$isolate_id) & ncbi_non_isangi$isolate_id != ""
]
input_audit <- tibble(
  item = c(
    "NCBI comparator records",
    "NCBI records lacking a run identifier",
    "Repeated non-missing NCBI run identifiers",
    "Study isolates",
    "Study isolates with unassigned ST",
    "Study lookup isolates absent from the feature-only AMRFinder file"
  ),
  value = c(
    nrow(ncbi_non_isangi),
    sum(is.na(ncbi_non_isangi$isolate_id) | ncbi_non_isangi$isolate_id == ""),
    sum(duplicated(valid_runs)),
    nrow(isangi_counts),
    sum(!isangi_counts$ST_assigned),
    sum(!isangi_counts$present_in_feature_file)
  )
)
input_audit

axis_max <- max(5, ceiling(max(plot_data$amr_gene_count) / 5) * 5)
header_y <- axis_max + 3.8
source_y <- axis_max + 1.7

violin_arguments <- list(
  mapping = aes(fill = group, colour = group),
  trim = TRUE,
  bw = STYLE$violin_bandwidth,
  adjust = 1,
  scale = "width",
  width = STYLE$violin_width,
  alpha = STYLE$violin_alpha,
  linewidth = 0.35,
  show.legend = FALSE
)
if (utils::packageVersion("ggplot2") >= "3.5.0") {
  violin_arguments$bounds <- c(0, Inf)
}

figure_3b <- ggplot(plot_data, aes(x = group, y = amr_gene_count)) +
  do.call(ggplot2::geom_violin, violin_arguments)

if (STYLE$show_study_points) {
  figure_3b <- figure_3b +
    geom_point(
      data = study_points,
      aes(colour = group),
      position = position_jitter(width = STYLE$point_jitter, height = 0, seed = STYLE$seed),
      size = STYLE$point_size,
      alpha = STYLE$point_alpha,
      stroke = 0,
      show.legend = FALSE
    )
}

figure_3b <- figure_3b +

  geom_linerange(
    data = group_summary,
    aes(x = group, ymin = q1, ymax = q3),
    inherit.aes = FALSE, linewidth = 1.65, colour = "white"
  ) +
  geom_linerange(
    data = group_summary,
    aes(x = group, ymin = q1, ymax = q3),
    inherit.aes = FALSE, linewidth = 0.95, colour = "#28323B"
  ) +
  geom_point(
    data = group_summary,
    aes(x = group, y = median_amr_genes),
    inherit.aes = FALSE, shape = 21, size = 2.25,
    stroke = 0.45, fill = "white", colour = "#28323B"
  ) +
  geom_text(
    data = group_summary,
    aes(x = group, label = paste0("n = ", scales::comma(n))),
    y = header_y,
    inherit.aes = FALSE, size = 3.15,
    family = STYLE$font_family, colour = "#28323B"
  ) +
  scale_fill_manual(values = group_colours, guide = "none") +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.55), drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, axis_max, by = 5),
    minor_breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(ylim = c(-0.8, axis_max + 5.2), clip = "on") +
  labs(x = NULL, y = "Number of AMR genes") +
  theme_classic(base_size = STYLE$base_size, base_family = STYLE$font_family) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "#EBEDF0", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(colour = "#616B75", linewidth = 0.3),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "#616B75", linewidth = 0.3),
    axis.ticks.length = grid::unit(1.2, "mm"),
    axis.text.x = element_text(size = 9, colour = "#28323B", margin = margin(t = 7)),
    axis.text.y = element_text(size = 8.2, colour = "#404A54"),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 8)),
    legend.position = "none",
    plot.margin = margin(t = 4, r = 5, b = 5, l = 4)
  )

if (STYLE$show_source_labels) {
  figure_3b <- figure_3b +
    geom_text(
      data = group_summary,
      aes(x = group, label = source),
      y = source_y,
      inherit.aes = FALSE, size = 2.55,
      family = STYLE$font_family, colour = "#707982"
    )
}

figure_3b

tiff_path <- file.path(outdir, "Figure_3B_AMR_gene_burden_600dpi.tiff")
png_path <- file.path(outdir, "Figure_3B_AMR_gene_burden_preview.png")

save_tiff <- function(p, path) {
  ragg::agg_tiff(
    filename = path, width = STYLE$width_mm, height = STYLE$height_mm,
    units = "mm", res = 600, compression = "lzw", background = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
}
save_png <- function(p, path) {
  ragg::agg_png(
    filename = path, width = STYLE$width_mm, height = STYLE$height_mm,
    units = "mm", res = 220, background = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
}

save_tiff(figure_3b, tiff_path)
save_png(figure_3b, png_path)

saveRDS(figure_3b, file.path(outdir, "Figure_3B_AMR_gene_burden_plot.rds"))

write_csv(group_summary, file.path(outdir, "Figure_3B_AMR_gene_burden_summary.csv"))
write_csv(isangi_counts, file.path(outdir, "Figure_3B_study_isolate_counts.csv"))
write_csv(input_audit, file.path(outdir, "Figure_3B_input_audit.csv"))
write_csv(
  plot_data %>% count(group, amr_gene_count, name = "n") %>%
    group_by(group) %>% mutate(proportion = n / sum(n)) %>% ungroup(),
  file.path(outdir, "Figure_3B_count_frequencies.csv")
)

writeLines(
  c(
    "Figure 3B. AMR gene-count distributions using the existing counting workflow.",
    paste(
      "Violin shapes show kernel-density estimates, using a common bandwidth of",
      STYLE$violin_bandwidth,
      "genes and equal maximum widths, trimmed to each group's observed range."
    ),
    "White points indicate medians and black vertical lines indicate interquartile ranges.",
    "Small coloured points show all study isolates; horizontal jitter separates overlapping counts.",
    "The NCBI comparator is not shown as individual points. n denotes records included in each group.",
    paste(
      "Other Isangi comprises",
      sum(isangi_counts$group == "Other Isangi" & isangi_counts$ST_assigned),
      "isolates assigned a non-ST335 sequence type and",
      sum(isangi_counts$group == "Other Isangi" & !isangi_counts$ST_assigned),
      "with an unassigned sequence type."
    ),
    "AMR counting rules and upstream database/version differences are unchanged from the preceding QMD."
  ),
  file.path(outdir, "Figure_3B_display_notes.txt")
)
capture.output(sessionInfo(), file = file.path(outdir, "Figure_3B_sessionInfo.txt"))

stopifnot(file.exists(tiff_path), file.info(tiff_path)$size > 0)
cat("TIFF saved to:\n", tiff_path, "\n\nPreview saved to:\n", png_path, "\n", sep = "")
