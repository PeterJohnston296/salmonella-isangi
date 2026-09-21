library(here)
library(tidyverse)
library(ape)
library(phangorn)
library(ggtree)
library(treeio)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)

tree_file        <- here::here("data", "data_raw",   "RAxML_best_trees", "RAxML_bestTree.ebg25_quality_assured_tree")
metadata_file    <- here::here("data", "data_clean", "EBG_25_metadata.csv")
cluster_file     <- here::here("data", "data_clean", "ebg25.snp-sites.fastbaps_lineages.csv")
mlst_file        <- here::here("data", "data_clean", "all_isolates_mlst_24052024.csv")
amr_classes_file <- here::here("data", "data_clean", "all_isangi_assemblies_amrfinder2_20062024.number_of_classes.csv")
snp_matrix_file  <- here::here("data", "data_clean", "ebg_25_snpmatrix.tsv")

out_fig_dir <- here::here("outputs", "figures")
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(tree_file)
tree <- midpoint(tree)

metadata_raw <- readr::read_csv(metadata_file, show_col_types = FALSE)

metadata <- metadata_raw |>
  dplyr::rename(Collection.Year = `Collection Year`) |>
  mutate(Isolate_ID = gsub("#", "_", Isolate_ID)) |>
  filter(Isolate_ID %in% tree$tip.label)

cluster_data <- readr::read_csv(cluster_file, show_col_types = FALSE) |>
  dplyr::rename(
    Isolate_ID = Isolates,
    Level.1    = `Level 1`,
    Level.2    = `Level 2`
  ) |>
  mutate(Isolate_ID = gsub("#", "_", Isolate_ID)) |>
  filter(Isolate_ID %in% tree$tip.label)

mlst <- readr::read_csv(mlst_file, col_names = FALSE, show_col_types = FALSE) |>
  mutate(
    V1 = gsub("_assembly_contigs\\.fasta$", "", X1),
    V1 = gsub("#", "_", V1),
    V3 = X3
  )

metadata <- metadata |>
  left_join(mlst, by = c("Isolate_ID" = "V1"))

amr_classes <- readr::read_csv(amr_classes_file, show_col_types = FALSE) |>
  mutate(Isolate_ID = gsub("#", "_", Isolate_ID))

metadata <- metadata |>
  left_join(amr_classes, by = "Isolate_ID") |>
  mutate(Unique_Class_Count = tidyr::replace_na(Unique_Class_Count, 0))

country_order <- c(
  "Malawi", "South Africa", "Nigeria", "United Kingdom",
  "Brazil", "United States", "Mexico", "Taiwan", "Other"
)

country_cols <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "grey70"
)
names(country_cols) <- country_order

country <- metadata |>
  transmute(Isolate_ID,
            Country = if_else(is.na(Country) | Country == "", "Other", Country))

top_countries <- names(sort(table(country$Country), decreasing = TRUE))[1:9]
country$Country[!country$Country %in% top_countries] <- "Other"

country <- as.data.frame(country)
rownames(country) <- country$Isolate_ID
country$Isolate_ID <- NULL

source_cols <- c(
  "Malawi water sample"   = "#FDB462",
  "QECH clinical isolate" = "#80B1D3",
  "Chatinkha nursery"     = "#FB8072",
  "Eastern Cape Outbreak" = "#B3DE69",
  "Other"                 = "white"
)

comment <- metadata |>
  transmute(Isolate_ID,
            Comment = if_else(is.na(Comment) | Comment == "", "Other", Comment))

comment <- as.data.frame(comment)
rownames(comment) <- comment$Isolate_ID
comment$Isolate_ID <- NULL

year <- metadata |>
  transmute(
    Isolate_ID,
    Collection.Year = as.numeric(Collection.Year)
  ) |>
  mutate(
    YearGroup = cut(
      Collection.Year,
      breaks = seq(2000, 2024, by = 4),
      labels = c("2000-2004", "2005-2008", "2009-2012",
                 "2013-2016", "2017-2020", "2021-2024"),
      right = TRUE
    ),
    YearGroup = as.character(YearGroup),
    YearGroup = if_else(is.na(YearGroup), "Unknown", YearGroup)
  ) |>
  select(Isolate_ID, YearGroup) |>
  as.data.frame()

rownames(year) <- year$Isolate_ID
year$Isolate_ID <- NULL

year_cols <- c(
  "2000-2004" = "#f0f9e8",
  "2005-2008" = "#ccebc5",
  "2009-2012" = "#a8ddb5",
  "2013-2016" = "#7bccc4",
  "2017-2020" = "#4eb3d3",
  "2021-2024" = "#2b8cbe",
  "Unknown"   = "#08589e"
)

amr_subclasses <- metadata |>
  transmute(Isolate_ID, Unique_Class_Count) |>
  mutate(
    Unique_Class_Count = case_when(
      Unique_Class_Count == 0 ~ "0",
      Unique_Class_Count %in% 1:3   ~ "1-3",
      Unique_Class_Count %in% 4:6   ~ "4-6",
      Unique_Class_Count %in% 7:10  ~ "7-10",
      Unique_Class_Count %in% 11:13 ~ "11-13",
      TRUE ~ NA_character_
    )
  ) |>
  as.data.frame()

rownames(amr_subclasses) <- amr_subclasses$Isolate_ID
amr_subclasses$Isolate_ID <- NULL

amr_cols <- c(
  "0"      = "#f2f0f7",
  "1-3"    = "#bcbddc",
  "4-6"    = "#9e9ac8",
  "7-10"   = "#6a51a3",
  "11-13"  = "#4a1486"
)

st_levels <- c("216", "335", "5028")
st_pal <- c(
  "216"  = "#FFF3A3",
  "335"  = "#5AB4AC",
  "5028" = "#d95f0e"
)

lineage_levels <- c("1", "2", "3", "4", "5")
lineage_pal <- c(
  "1" = "#F4E39C",
  "2" = "#D4E9E2",
  "3" = "#E7C1C8",
  "4" = "#C9D3E8",
  "5" = "#E1D4C6"
)

meta_all <- metadata |>
  left_join(cluster_data, by = "Isolate_ID") |>
  select(Isolate_ID, Country, Comment,
         Collection.Year, Unique_Class_Count, V3, Level.1)

circ_base <- ggtree(tree, layout = "circular") %<+% meta_all

circ_hl <- circ_base

for (lev in lineage_levels) {
  tips <- meta_all |> filter(Level.1 == lev) |> pull(Isolate_ID)
  tip_nums <- match(tips, tree$tip.label)
  tip_nums <- tip_nums[!is.na(tip_nums)]
  if (length(tip_nums) >= 2) {
    node <- ape::getMRCA(tree, tip_nums)
    if (!is.na(node)) {
      circ_hl <- circ_hl +
        geom_hilight(node = node, fill = lineage_pal[lev], alpha = 0.3)
    }
  }
}

tip_df <- circ_base$data |> filter(isTip)

circ <- circ_hl +
  geom_point(
    data = tip_df,
    aes(x = x, y = y, fill = factor(V3)),
    shape  = 21,
    size   = 2,
    colour = "black",
    stroke = 0.3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = st_pal,
    name   = "Sequence Type",
    breaks = st_levels,
    labels = st_levels
  ) +
  theme(legend.position = "left")

lineage_legend_df <- tibble(
  x = Inf,
  y = Inf,
  Lineage = factor(lineage_levels, levels = lineage_levels)
)

circ <- circ +
  geom_point(
    data = lineage_legend_df,
    aes(x = x, y = y, colour = Lineage),
    shape = 22,
    size  = 4,

    fill  = lineage_pal[as.character(lineage_legend_df$Lineage)],
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_colour_manual(
    values = lineage_pal,
    name   = "Lineage",
    guide  = guide_legend(
      override.aes = list(
        shape = 22,
        size  = 5,

        fill  = lineage_pal
      )
    )
  )

add_heatmap_layer <- function(p, df, offset, width, fill_values, name,
                              breaks = NULL, labels = NULL,
                              legend_order = NULL,
                              legend_position = "left",
                              border_color = "white") {

  p2 <- p + ggnewscale::new_scale_fill()

  gheatmap(
    p2, df,
    offset   = offset,
    width    = width,
    colnames = FALSE,
    color    = border_color
  ) +
    scale_fill_manual(
      values = fill_values,
      name   = name,
      breaks = breaks,
      labels = labels
    ) +
    guides(fill = guide_legend(order = legend_order)) +
    theme(legend.position = legend_position)
}

ring_specs <- list(
  list(df = country,        offset = 0.0001, width = 0.10,
       fill_values = country_cols, name = "Country",
       breaks = country_order, labels = country_order,
       legend_order = 3, border = "white"),

  list(df = comment,        offset = 0.015, width = 0.10,
       fill_values = source_cols, name = "Source",
       breaks = names(source_cols), labels = names(source_cols),
       legend_order = 4, border = "grey80"),

  list(df = year,           offset = 0.030, width = 0.10,
       fill_values = year_cols, name = "Year",
       breaks = names(year_cols), labels = names(year_cols),
       legend_order = 5, border = "white"),

  list(df = amr_subclasses, offset = 0.045, width = 0.10,
       fill_values = amr_cols, name = "No. AMR sub-classes",
       breaks = names(amr_cols), labels = names(amr_cols),
       legend_order = 6, border = "white")
)

circ_final <- purrr::reduce(
  ring_specs,
  .init = circ,
  .f = function(p, s)
    add_heatmap_layer(
      p, s$df, s$offset, s$width,
      s$fill_values, s$name,
      s$breaks, s$labels,
      s$legend_order, "left",
      border_color = s$border
    )
)

snp_matrix_raw <- readr::read_tsv(snp_matrix_file, show_col_types = FALSE)

max_snp_dist <- snp_matrix_raw |>
  dplyr::select(where(is.numeric)) |>
  unlist() |>
  max(na.rm = TRUE)

max_tree_dist <- max(cophenetic(tree))
tree_units_per_snp <- max_tree_dist / max_snp_dist

scale_snp <- 50L
scale_len <- scale_snp * tree_units_per_snp

p_scale <- ggplot() +
  geom_segment(
    aes(x = 0, xend = scale_len, y = 0, yend = 0),
    linewidth = 0.6
  ) +
  geom_text(
    aes(x = scale_len / 2, y = -0.4 * scale_len),
    label = paste0(scale_snp, " SNPs"),
    size  = 3
  ) +
  xlim(-2 * scale_len, 8 * scale_len) +
  ylim(-scale_len,  scale_len) +
  theme_void()

p_final <- circ_final / p_scale + plot_layout(heights = c(10, 1))

ggsave(
  filename = file.path(out_fig_dir, "ebg25_circular_metadata_lineages_scalebar.tiff"),
  plot     = p_final,
  width    = 22,
  height   = 14,
  dpi      = 600
)

library(ggtree)
library(ggnewscale)
library(patchwork)
library(dplyr)
library(tidyr)

meta_all <- metadata |>
  left_join(cluster_data,  by = "Isolate_ID") |>
  left_join(mlst,          by = c("Isolate_ID" = "V1")) |>
  left_join(amr_classes,   by = "Isolate_ID") |>
  mutate(
    V3 = as.character(V3),
    Unique_Class_Count = tidyr::replace_na(Unique_Class_Count, 0L)
  )

st_levels <- c("216", "335", "5028")
st_pal <- c(
  "216"  = "#FFF3A3",
  "335"  = "#5AB4AC",
  "5028" = "#C6B7D8"
)

lineage_levels <- c("1", "2", "3", "4", "5")
lineage_pal <- c(
  "1" = "#F4E39C",
  "2" = "#D4E9E2",
  "3" = "#E7C1C8",
  "4" = "#C9D3E8",
  "5" = "#E1D4C6"
)

year_cols <- c(
  "2000-2004" = "#f0f9e8",
  "2005-2008" = "#ccebc5",
  "2009-2012" = "#a8ddb5",
  "2013-2016" = "#7bccc4",
  "2017-2020" = "#4eb3d3",
  "2021-2024" = "#2b8cbe",
  "Unknown"   = "#08589e"
)

amr_cols <- c(
  "0"      = "#f2f0f7",
  "1-3"    = "#bcbddc",
  "4-6"    = "#9e9ac8",
  "7-10"   = "#6a51a3",
  "11-13"  = "#4a1486"
)

heat_meta <- meta_all |>
  transmute(
    Isolate_ID,
    Collection.Year = as.numeric(Collection.Year),
    YearGroup = cut(
      Collection.Year,
      breaks = seq(2000, 2024, by = 4),
      labels = c("2000-2004", "2005-2008", "2009-2012",
                 "2013-2016", "2017-2020", "2021-2024"),
      right = TRUE
    ),
    YearGroup = as.character(YearGroup),
    YearGroup = if_else(is.na(YearGroup), "Unknown", YearGroup),
    AMR_bin = case_when(
      Unique_Class_Count == 0       ~ "0",
      Unique_Class_Count %in% 1:3   ~ "1-3",
      Unique_Class_Count %in% 4:6   ~ "4-6",
      Unique_Class_Count %in% 7:10  ~ "7-10",
      Unique_Class_Count %in% 11:13 ~ "11-13",
      TRUE                          ~ NA_character_
    )
  )

heat_df <- heat_meta |>
  select(Isolate_ID, YearGroup, AMR_bin) |>
  as.data.frame()

rownames(heat_df) <- heat_df$Isolate_ID
heat_df$Isolate_ID <- NULL

p_tree <- ggtree(tree, layout = "rectangular") %<+% meta_all

p_tree <- p_tree +
  geom_tree(aes(colour = factor(Level.1, levels = lineage_levels)),
            linewidth = 0.4) +
  scale_colour_manual(
    values = lineage_pal,
    name   = "Lineage"
  )

p_tree <- p_tree +
  geom_tippoint(
    aes(fill = factor(V3, levels = st_levels)),
    shape  = 21,
    size   = 1.8,
    colour = "black",
    stroke = 0.2
  ) +
  scale_fill_manual(
    values = st_pal,
    name   = "Sequence Type"
  )

p_tree <- p_tree +
  geom_treescale(
    x = 0,
    y = -5,
    width = max(tree$edge.length) / 4,
    fontsize = 3,
    linesize = 0.4
  ) +
  theme(
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 9),
    axis.text       = element_blank(),
    axis.ticks      = element_blank()
  )

p_tree_heat <- p_tree +
  ggnewscale::new_scale_fill()

p_tree_heat <- gheatmap(
  p_tree_heat,
  heat_df,
  offset   = 0.3,
  width    = 0.15,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 90,
  colnames_offset_y = 0.3,
  color    = NA
)

p_tree_heat <- p_tree_heat +
  scale_fill_manual(
    values = c(year_cols, amr_cols),
    name   = "Metadata",
    breaks = c(names(year_cols), names(amr_cols)),
    labels = c(
      paste0("Year ", names(year_cols)),
      paste0("AMR ", names(amr_cols))
    )
  ) +
  theme(
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(5, 20, 5, 5)
  )

ggsave(
  filename = here::here("outputs", "figures", "ebg25_rectangular_lineage_ST_heatmap.png"),
  plot     = p_tree_heat,
  width    = 16,
  height   = 9,
  dpi      = 600
)
