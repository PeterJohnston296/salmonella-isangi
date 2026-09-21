library(tidyverse)
library(here)
library(lubridate)

ebg <- read_csv(here("data", "data_raw", "ebg.csv"))

genomes <- read_csv(here("data", "data_raw", "good_genomes.csv"))

metadata <-  read_csv(here("data", "data_raw", "metadata.csv"))

combined <- left_join(genomes, metadata, by = c("Accession" = "Isolate_ID"))

combined <- combined %>%
  dplyr::select(-c(V4, V5, V6, V7, V8, V9, V10, Comment.y))

combined <- combined %>%
  rename(sequence_type = V3
  )

table(combined$clonal_complex, useNA = "ifany")

table(combined$sequence_type, useNA = "ifany")

xdr_all <- combined %>%
  filter(`XDR (Traditional)` == 1)

xdr_335 <- combined %>%
  filter(`XDR (Traditional)` == 1 & sequence_type == 335)

xdr_216 <- combined %>%
  filter(`XDR (Traditional)` == 1 & sequence_type == 216)

combined %>%
  filter(Country == "South Africa") %>%
  tally()

combined %>%
  filter(Country == "South Africa", sequence_type == 335) %>%
  tally()

combined %>%
  filter(Country == "South Africa", sequence_type != 335) %>%
  distinct(sequence_type) %>%
  arrange(sequence_type)

combined %>%
  filter(Country == "South Africa", `XDR (Traditional)` == 1) %>%
  tally()

combined %>%
  filter(Country == "South Africa") %>%
  count(`XDR (Traditional)`, sequence_type, sort = TRUE)

combined %>%
  filter(Country == "South Africa") %>%
  group_by(`XDR (Traditional)`, sequence_type) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n))

combined %>%
  filter(sequence_type == 335) %>%
  summarise(
    total = n(),
    xdr = sum(`XDR (Traditional)` == 1),
    prop_xdr = xdr / total
  )

combined %>%
  mutate(ST_group = case_when(
    sequence_type == 335 ~ "ST335",
    sequence_type == 216 ~ "ST216",
    TRUE ~ "Other"
  )) %>%
  count(ST_group, sort = TRUE)

combined %>%
  filter(sequence_type == 335, Country %in% c("Malawi", "South Africa", "Mozambique")) %>%
  count()

   combined %>%
  filter(sequence_type == 216, !is.na(Country)) %>%
  summarise(n_countries = n_distinct(Country))

combined %>%
  filter(sequence_type == 216, !is.na(Continent)) %>%
  summarise(n_continents = n_distinct(Continent))

combined %>%
  filter(Country == "Malawi") %>%
  tally()

combined %>%
  filter(Country == "Malawi", sequence_type == 335) %>%
  tally()

filtered <- combined %>%
  filter(!(sequence_type == 216 & `XDR (Traditional)` == 1))

write_csv(
  combined,
  here("data", "data_clean", "all_good_isangi_genomes_metadata_amr_ebg.csv")
)
