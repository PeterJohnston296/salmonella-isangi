library(here)
library(tidyverse)

bsi <- read_csv(here("data/data_raw/BSI_2015_2025.csv"))
isangi <- read_csv(here("data/data_raw/isangi_bioinformatically_confirmed.csv"), col_names = FALSE)
potential_isangi <- read_csv(here("data/data_raw/list_of_potential_isangi_BCs.csv"), col_names = TRUE)

isangi$X1 <- substr(isangi$X1, 1, 6)

setdiff(isangi$X1, potential_isangi$IDS)

setdiff(potential_isangi$IDS, isangi$X1)

library(dplyr)
library(stringr)

bsi <- bsi %>%
  mutate(organism = str_squish(organism))

bsi <- bsi %>%
  filter(organism != "Salmonella typhi")

bsi <- bsi %>%
  mutate(XDR = as.integer(
    `Ampicillin 10` == "Resistant" &
    `Cefpodoxime 10 (Ceftriaxone)` == "Resistant" &
    `Chloramphenicol 30` == "Resistant" &
    `Cotrimoxazole 25` == "Resistant" &
    (
      `Ciprofloxacin 5` == "Resistant" | `Pefloxacin 5` == "Resistant"
    )
  ))

XDR <- bsi %>%
  filter(XDR == 1)

setdiff(isangi$X1, XDR$`Request Number`)

twenty_one <- semi_join(bsi, isangi, by = c("Request Number" = "X1"))

library(ComplexUpset)

ggplot(twenty_one, aes(x = DSP)) +
  geom_bar() +
  labs(title = "Confirmed Isangi Bloodstream Isolates Over Time")

library(ComplexUpset)

upset(twenty_one,
      intersect = c("Ampicillin 10", "Cefpodoxime 10 (Ceftriaxone)",
                    "Chloramphenicol 30", "Cotrimoxazole 25", "Ciprofloxacin 5", "Pefloxacin 5"),
      name = "Resistance Pattern")

library(dplyr)
library(tidyr)
library(ggplot2)

resistance_heatmap_data <- twenty_one %>%
  select(`Request Number`,
         `Ampicillin 10`,
         `Augmentin 30`,
         `Cefpodoxime 10 (Ceftriaxone)`,
         `Chloramphenicol 30`,
         `Cotrimoxazole 25`,
         `Ciprofloxacin 5`,
         `Pefloxacin 5`) %>%
  pivot_longer(
    cols = -`Request Number`,
    names_to = "Antibiotic",
    values_to = "Result"
  )

resistance_heatmap_data <- resistance_heatmap_data %>%
  mutate(
    Result = replace_na(Result, "Not tested"),
    Result = factor(Result, levels = c("Sensitive", "Intermediate", "Resistant", "Not tested")),
    Antibiotic = factor(Antibiotic, levels = c(
      "Ampicillin 10",
      "Augmentin 30",
      "Cefpodoxime 10 (Ceftriaxone)",
      "Chloramphenicol 30",
      "Cotrimoxazole 25",
      "Ciprofloxacin 5",
      "Pefloxacin 5"
    ))
  )

ggplot(resistance_heatmap_data, aes(x = Antibiotic, y = `Request Number`, fill = Result)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(
    "Sensitive" = "steelblue",
    "Intermediate" = "gold",
    "Resistant" = "firebrick",
    "Not tested" = "grey80"
  )) +
  labs(
    title = "Antibiotic Resistance Profiles of WGS-confirmed *Salmonella Isangi* Isolates",
    x = "Antibiotic",
    y = "Request Number"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

ggsave(
  filename = here::here("outputs", "figures", "Phenotypic_resistance_Isangis.tiff"),
  dpi = 300,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

bsi <- bsi %>%
  mutate(
    collection_date = dmy(DSP),
    year = year(collection_date)
  )

xdr_by_year <- bsi %>%
  filter(!is.na(year)) %>%
  group_by(year) %>%
  summarise(
    total_isolates = n(),
    xdr_isolates = sum(XDR == 1, na.rm = TRUE),
    prop_xdr = xdr_isolates / total_isolates
  ) %>%
  ungroup()

xmin <- 2018 + (as.numeric(ymd("2018-05-28") - ymd("2018-01-01")) / 365.25)
xmax <- 2023 + (as.numeric(ymd("2023-02-09") - ymd("2023-01-01")) / 365.25)

ggplot(xdr_by_year, aes(x = year, y = prop_xdr)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  annotate("text",
           x = (xmin + xmax) / 2,
           y = 0.17,
           label = "Isangi outbreak window",
           size = 4, fontface = "italic") +
  geom_line(size = 1.3, color = "firebrick") +
  geom_point(size = 3, color = "firebrick") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.2)) +
  scale_x_continuous(breaks = seq(min(xdr_by_year$year), max(xdr_by_year$year), 1)) +
  labs(
    title = "Proportion of XDR Isolates Over Time",
    subtitle = "Shaded region shows the Isangi outbreak period (28 May 2018 – 9 Feb 2023)",
    x = "Year",
    y = "Proportion XDR",
    caption = "Data: BSI dataset"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

ggsave(
  filename = here::here("outputs", "figures", "XDR_BC_over_time.tiff"),
  dpi = 300,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)

bsi <- bsi %>%
  mutate(
    collection_date = dmy(DSP),
    year = year(collection_date),
    MDR = as.integer(
      `Ampicillin 10` == "Resistant" &
      `Chloramphenicol 30` == "Resistant" &
      `Cotrimoxazole 25` == "Resistant"
    ),
    MDR_status = case_when(
      XDR == 1 ~ "XDR",
      MDR == 1 ~ "MDR",
      TRUE     ~ "Non-MDR"
    )
  )

bsi_summary <- bsi %>%
  filter(!is.na(year)) %>%
  group_by(year, MDR_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label = paste0(round(prop * 100), "%")
  ) %>%
  ungroup()

ggplot(bsi_summary, aes(x = year, y = n, fill = MDR_status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(prop > 0.04, label, "")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "Non-MDR" = "#a6bddb",
    "MDR"     = "#fd8d3c",
    "XDR"     = "#e31a1c"
  )) +
  scale_x_continuous(breaks = seq(min(bsi_summary$year), max(bsi_summary$year), 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Resistance Category of Bloodstream Isolates by Year",
    subtitle = "Based on phenotypic resistance to ampicillin, chloramphenicol, cotrimoxazole,\nfluoroquinolones and cefpodoxime (proxy for ceftriaxone)",
    x = "Collection Year",
    y = "Number of Isolates",
    fill = "Resistance Category"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = here::here("outputs", "figures", "XDR_MDR_stacked_bar.tiff"),
  dpi = 600,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)

bsi_prop <- bsi %>%
  mutate(
    collection_date = dmy(DSP),
    year = year(collection_date),
    MDR = as.integer(
      `Ampicillin 10` == "Resistant" &
      `Chloramphenicol 30` == "Resistant" &
      `Cotrimoxazole 25` == "Resistant"
    ),
    MDR_status = case_when(
      XDR == 1 ~ "XDR",
      MDR == 1 ~ "MDR",
      TRUE     ~ "Non-MDR"
    )
  ) %>%
  filter(!is.na(year)) %>%
  group_by(year, MDR_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label = paste0(round(prop * 100), "%")
  ) %>%
  ungroup()

ggplot(bsi_prop, aes(x = factor(year), y = prop, fill = MDR_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(
    aes(label = ifelse(prop > 0.04, label, "")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c(
    "Non-MDR" = "#a6bddb",
    "MDR"     = "#fd8d3c",
    "XDR"     = "#e31a1c"
  )) +
  labs(
    title = "Proportion of Resistance Categories by Year (100% Stacked)",
    subtitle = "Each bar represents the relative composition of Non-MDR, MDR, and XDR isolates per year",
    x = "Collection Year",
    y = "Proportion of Isolates",
    fill = "Resistance Category"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = here::here("outputs", "figures", "MDR_XDR_proportional_bar.tiff"),
  dpi = 600,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)

bsi_prop <- bsi %>%
  mutate(
    collection_date = dmy(DSP),
    year = year(collection_date),
    MDR = as.integer(
      `Ampicillin 10` == "Resistant" &
      `Chloramphenicol 30` == "Resistant" &
      `Cotrimoxazole 25` == "Resistant"
    ),
    MDR_status = case_when(
      XDR == 1 ~ "XDR",
      MDR == 1 ~ "MDR",
      TRUE     ~ "Non-MDR"
    )
  ) %>%
  filter(!is.na(year)) %>%
  group_by(year, MDR_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    prop = n / total
  ) %>%
  ungroup()

ggplot(bsi_prop, aes(x = factor(year), y = prop, fill = MDR_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = c(
    "Non-MDR" = "#a6bddb",
    "MDR"     = "#fd8d3c",
    "XDR"     = "#e31a1c"
  )) +
  labs(
    x = "Collection Year",
    y = "Percentage of Isolates",
    fill = "Resistance Category"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = here::here("outputs", "figures", "MDR_XDR_proportional_no_bar_labels.tiff"),
  dpi = 300,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(here)

date_end <- as.Date("2023-12-31")

bsi_limited <- bsi %>%
  mutate(
    collection_date = suppressWarnings(dmy(DSP)),
    year            = year(collection_date)
  ) %>%
  filter(!is.na(collection_date) & collection_date <= date_end) %>%
  mutate(

    XDR = as.integer(
      `Ampicillin 10` == "Resistant" &
      `Cefpodoxime 10 (Ceftriaxone)` == "Resistant" &
      `Chloramphenicol 30` == "Resistant" &
      `Cotrimoxazole 25` == "Resistant" &
      (`Ciprofloxacin 5` == "Resistant" | `Pefloxacin 5` == "Resistant")
    ),
    MDR = as.integer(
      `Ampicillin 10` == "Resistant" &
      `Chloramphenicol 30` == "Resistant" &
      `Cotrimoxazole 25` == "Resistant"
    ),
    MDR_status = dplyr::case_when(
      XDR == 1 ~ "XDR",
      MDR == 1 ~ "MDR",
      TRUE     ~ "Non-MDR"
    )
  )

bsi_prop <- bsi_limited %>%
  filter(!is.na(year)) %>%
  count(year, MDR_status, name = "n") %>%
  group_by(year) %>%
  mutate(
    total = sum(n),
    prop  = n / total
  ) %>%
  ungroup() %>%
  mutate(

    MDR_status = factor(MDR_status, levels = c("Non-MDR", "MDR", "XDR"))
  )

p_mdr_xdr_prop <- ggplot(bsi_prop, aes(x = factor(year), y = prop, fill = MDR_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c("Non-MDR" = "#a6bddb", "MDR" = "#fd8d3c", "XDR" = "#e31a1c"),
    breaks = c("XDR", "MDR", "Non-MDR")
  ) +
  labs(
    x = "Collection Year",
    y = "Percentage of Isolates",
    fill = "Resistance Category"
  ) +
  theme_minimal(base_size = 12)

print(p_mdr_xdr_prop)

ggsave(
  filename = here::here("outputs", "figures", "MDR_XDR_proportional_upto_2023.tiff"),
  plot     = p_mdr_xdr_prop,
  dpi      = 600,
  width    = 8,
  height   = 6,
  units    = "in",
  bg       = "white"
)

message("Years plotted: ", paste(sort(unique(bsi_prop$year)), collapse = ", "))
