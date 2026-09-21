# Load libraries
library(tidyverse)
library(ggridges)
library(scales)


# Read the data
ncbi_data <- read_delim("2026.01.13.ncbi_pathogen_salmonella_amr_serovar.tsv.tsv", delim = "\t")
enterobase_data <- read_delim("isangi_enterobase.txt", delim = "\t")

# Process the data
enterobase_accessions <- enterobase_data %>%
  mutate(
    accession = str_extract(
      `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Bases;Average Length;Status)`,
      "^[^;]+"
    )
  ) %>%
  select(accession, ST)

processed_data <- ncbi_data %>%
  mutate(
    amr_gene_count = str_count(`#AMR genotypes`, "COMPLETE"),
    is_isangi = str_detect(`Computed types`, regex("Isangi", ignore_case = TRUE))
  ) %>%
  left_join(enterobase_accessions, by = c("Run" = "accession")) %>%
  mutate(
    ST = as.character(ST),
    isangi_group = case_when(
      !is_isangi ~ "Salmonella\n(non-Isangi)",
      is_isangi & ST == "335" ~ "S. Isangi\n(ST335)",
      is_isangi ~ "S. Isangi\n(non-ST335)",
      TRUE ~ NA_character_
    )
  )

processed_data %>% group_by(is_isangi) %>% summarise(n = n())

# Violin plot with different smoothing
violin_counts <- processed_data %>%
  filter(!is.na(amr_gene_count), !is.na(isangi_group)) %>%
  mutate(
    isangi_group = factor(
      isangi_group,
      levels = c("Salmonella\n(non-Isangi)", "S. Isangi\n(ST335)", "S. Isangi\n(non-ST335)")
    )
  ) %>%
  count(isangi_group, name = "n")

violin_plot <- processed_data %>%
  filter(!is.na(amr_gene_count), !is.na(isangi_group)) %>%
  mutate(
    isangi_group = factor(
      isangi_group,
      levels = c("Salmonella\n(non-Isangi)", "S. Isangi\n(ST335)", "S. Isangi\n(non-ST335)")
    )
  ) %>%
  ggplot(aes(x = isangi_group, y = amr_gene_count, fill = isangi_group)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 6) +
  geom_text(
    data = violin_counts,
    aes(x = isangi_group, y = Inf, label = paste0("n = ", n)),
    vjust = 1.3,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  labs(x = NULL, y = "Number of AMR genes", fill = "Group") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave("amr_gene_violin_plot.png", violin_plot, width = 8, height = 5, dpi = 300)

## then, plot it as a normalised histogram with the absolute numbers annotated on.

# remove rows with missing values
filtered_data <- processed_data %>% 
  filter(!is.na(amr_gene_count), !is.na(isangi_group))

# counts and proportions
plot_data <- filtered_data %>% 
  count(isangi_group, amr_gene_count, name = "n") %>% 
  group_by(isangi_group) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup()

# normalised histogram with counts
ggplot(plot_data, aes(amr_gene_count, prop, fill = isangi_group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.15, size = 3) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Number of AMR genes",
       y = "Proportion within group",
       fill = "Group") +
  theme_minimal()
