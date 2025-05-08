### Setup env
library(tidyverse)

setwd("/mfd_distance_decay")

### Load metadata for 10 km grid reference samples
grid <- data.table::fread('output/2025-04-23_MFD_samples_grid_10km.tsv')

## Pull sample IDs
samples <- grid %>%
  pull(fieldsample_barcode)

## Load extended and aggregated metadata and filter to grid representatives
## select sample ID and seq ID 
metadata <- data.table::fread('data/2023-09-22_corrected_combined_metadata.csv') %>%
  filter(fieldsample_barcode %in% samples) %>%
  select(fieldsample_barcode, seq_id)

## Create list of file names based on seq IDs
files <- metadata %>%
  select(seq_id) %>%
  distinct() %>%
  mutate(file = str_c("arc_bac_", seq_id, "_forward.fq")) %>%
  select(file)

### Write list of files to output
data.table::fwrite(files, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_libs_classification.txt"))

## Load file to read ID key
key <- data.table::fread('data/2023-10-12_arcbac_MFD_read_key.csv')

## Filter key based on reduced metadata
sub.key <- key %>%
  filter(seq_id %in% metadata$seq_id)

rm(key)
gc()

## Create a list of patterns to search
patterns <- sub.key %>% 
  select(read_name)

### Write patterns to output
data.table::fwrite(patterns, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_forward_reads_classification.txt"))
