# R script to parse CAT files and CoverM contig coverage outputs and summarise them at a high level.
# Author: Raphael Eisenhofer 18/12/2022

library(tidyverse)
library(janitor)

# Import arguments
# N.B. arg[1] = {group}, arg[2] = coverm file, arg[3] = official CAT output, arg[4] = unoffical CAT output
args = commandArgs(trailingOnly=TRUE)


# Import data
cat_official <- read_delim(args[3], delim = "\t") %>%
  clean_names() %>%
  select(-c("lineage_scores", "lineage", "classification")) %>%
  mutate(number_contig = str_replace_all(number_contig, "-.*", ""),
         across(cols = c(3, 4, 5, 6, 7, 8, 9), .fns = ~ str_replace_all(., ":.*", "")))

cat_unofficial <- read_delim(args[4], delim = "\t") %>%
  clean_names() %>%
  select(-c("lineage_scores", "lineage", "classification")) %>%
  filter(str_detect(full_lineage_names, "Bacteria", negate = TRUE)) %>%
  mutate(number_contig = str_replace_all(number_contig, "-.*", ""),
         kingdom = case_when(str_detect(full_lineage_names, "Fungi \\(kingdom\\)") ~ "fungi",
                             str_detect(full_lineage_names, "Metazoa \\(kingdom\\)") ~ "metazoa",
                             str_detect(full_lineage_names, "Viridiplantae \\(kingdom\\)") ~ "plants",
                             )
         )


coverm_official <- read_delim(args[2]) %>%
  clean_names() %>%
  select(!ends_with("mean") & !ends_with("count") & !ends_with("fraction")) %>%
  pivot_longer(., cols = !genome, names_to = "sample", values_to = "relative_abundance") %>%
  mutate(sample = str_replace_all(sample, "_m_relative_abundance_percent", "")) %>%
  inner_join(., cat_official, by = c("genome" = "number_contig")) %>%
  filter(superkingdom != "Bacteria" & superkingdom != "Archaea")

write_tsv(coverm_official, paste0("3_Outputs/non_bacterial/", args[1], "_CAT_full_table.tsv"))

coverm_unofficial <- read_delim(args[2]) %>%
  clean_names() %>%
  select(!ends_with("mean") & !ends_with("count") & !ends_with("fraction")) %>%
  pivot_longer(., cols = !genome, names_to = "sample", values_to = "relative_abundance") %>%
  mutate(sample = str_replace_all(sample, "_m_relative_abundance_percent", "")) %>%
  inner_join(., cat_unofficial, by = c("genome" = "number_contig"))


# Old code for where you have a coverm file per sample

# coverage_files <- list.files(path = "non_bacterial/coverm_outputs/", pattern = "*genome.tsv", full.names = TRUE)
# coverage_dfs_official <- coverage_files %>%
#   map(read_delim) %>%
#   map(clean_names) %>%
#   map_df(inner_join, cat_official, by = c("genome" = "number_contig")) %>%
#   select(!ends_with("mean") & !ends_with("count") & !ends_with("fraction")) %>%
#   pivot_longer(., cols = ends_with("percent"), names_to = "sample", values_to = "relative_abundance") %>%
#   filter(superkingdom != "Bacteria" & superkingdom != "Archaea") %>%
#   mutate(sample = str_replace_all(sample, "_.*", ""))

# coverage_dfs_unofficial <- coverage_files %>%
#   map(read_delim) %>%
#   map(clean_names) %>%
#   map_df(right_join, cat_unofficial, by = c("genome" = "number_contig")) %>%
#   select(!ends_with("mean") & !ends_with("count") & !ends_with("fraction")) %>%
#   pivot_longer(., cols = ends_with("percent"), names_to = "sample", values_to = "relative_abundance") %>%
#   mutate(sample = str_replace_all(sample, "_.*", ""))


# High level summaries
superkingdom_summary <- coverm_official %>%
  group_by(sample, superkingdom) %>%
  summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE)) %>%
  filter(relative_abundance > 0.2) %>%
  pivot_wider(., names_from = superkingdom, values_from = relative_abundance)

write_tsv(superkingdom_summary, paste0("3_Outputs/non_bacterial/", args[1], "_superkingdom_summary.tsv"))

kingdom_summary <- coverm_unofficial %>%
  group_by(sample, kingdom) %>%
  summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE)) %>%
  filter(relative_abundance > 0.2) %>%
  pivot_wider(., names_from = kingdom, values_from = relative_abundance)

write_tsv(kingdom_summary, paste0("3_Outputs/non_bacterial/", args[1], "_kingdom_summary.tsv"))

phylum_summary <- coverm_official %>%
  group_by(sample, superkingdom, phylum) %>%
  summarise(relative_abundance = sum(relative_abundance, na.rm = TRUE)) %>%
  filter(relative_abundance > 0.2) %>%
  pivot_wider(names_from = !c("sample", "relative_abundance"), values_from = relative_abundance)
  
write_tsv(phylum_summary, paste0("3_Outputs/non_bacterial/", args[1], "_phylum_summary.tsv"))



