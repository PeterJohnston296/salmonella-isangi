library(tidyverse)
library(lubridate)
library(patchwork)
library(here)

metadata <- readr::read_csv(
  here::here("data", "data_clean", "malawi_metadata_for_plotting.csv"),
  show_col_types = FALSE
)

metadata_cols <- 10
snp_matrix <- metadata[, -(1:metadata_cols)] %>% as.data.frame()
rownames(snp_matrix) <- metadata$Isolate
snp_matrix <- as.matrix(snp_matrix)
storage.mode(snp_matrix) <- "numeric"

ref_isolate <- metadata %>% filter(Class == "Blood culture") %>% pull(Isolate) %>% .[1]
ref_index <- which(metadata$Isolate == ref_isolate)

metadata <- metadata %>%
  mutate(
    Class            = if_else(Class == "Chatinkha nursery", "Neonatal unit", Class),
    Distance_to_ref  = as.vector(snp_matrix[, ref_index]),
    `Collection date` = dmy(`Collection date`)
  )

class_colors <- c(
  "River water"   = "#1f78b4",
  "Blood culture" = "#e31a1c",
  "Neonatal unit" = "#bdbdbd",
  "CSF"           = "#6a3d9a"
)

plot_full <- ggplot() +
  geom_point(
    data = metadata %>% filter(Class == "Neonatal unit"),
    aes(x = `Collection date`, y = Distance_to_ref, color = Class),
    size = 4, alpha = 0.7,
    position = position_jitter(width = 20)
  ) +
  geom_point(
    data = metadata %>% filter(Class != "Neonatal unit"),
    aes(x = `Collection date`, y = Distance_to_ref, color = Class),
    size = 4, alpha = 0.7,
    position = position_jitter(width = 20)
  ) +
  scale_color_manual(values = class_colors) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(
    title = "2018–2023",
    x = "Collection Date",
    y = "SNP Distance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot_2020 <- ggplot() +
  geom_point(
    data = metadata %>% filter(Class == "Neonatal unit"),
    aes(x = `Collection date`, y = Distance_to_ref, color = Class),
    size = 4, alpha = 0.7,
    position = position_jitter(width = 20)
  ) +
  geom_point(
    data = metadata %>% filter(Class != "Neonatal unit"),
    aes(x = `Collection date`, y = Distance_to_ref, color = Class),
    size = 4, alpha = 0.7,
    position = position_jitter(width = 20)
  ) +
  scale_color_manual(values = class_colors) +
  scale_x_date(
    limits = c(as.Date("2020-01-01"), as.Date("2020-12-31")),
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    title = "2020",
    x = "Collection Date",
    y = "SNP Distance",
    color = "Sample Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

combined_plot <- plot_full / plot_2020 +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "i")

print(combined_plot)

ggsave(
  filename = here::here("outputs", "figures", "final_snp_combined_with_tags.tiff"),
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 600,
  units = "in",
  device = "tiff",
  bg = "white"
)
