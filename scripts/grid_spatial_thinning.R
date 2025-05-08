### Setup env
library(tidyverse)
library(vegan)

### Import color palettes
source('scripts/MFD_colors.R')

setwd("/mfd_distance_decay")

### Import data
## Load genus-aggeregated observational table submitted to random 
## subsampling without replacement
data <- data.table::fread('data/2025-02-13_MFD_arcbac_genus_rarefaction_rel.csv')

## Load metadata file and create "complex" correponding to full MFDO1 string
## Remove samples with no GPS information or unreliable coordinates
data.grid <- readxl::read_excel('data/2025-04-14_mfd_db.xlsx') %>%
  filter(fieldsample_barcode %in% colnames(data),
         !coords_reliable == "No",
         !is.na(latitude),
         !is.na(mfd_hab1)) %>%
  # mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
  #        across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
  ungroup() %>%
  rename(long = longitude,
         lat = latitude)

### Create summary of MFDO1 categories with sufficient coverage 
## and/or number of representative samples
groups.summary <- data.grid %>%
  filter(!project_id == "P12_1") %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater")) %>%
  mutate(grid = "All samples") %>%
  arrange(desc(size))


### Select representative samples across DK based on community similarity
## Summary of number of samples per MFDO1 category when mapping to the 10 km grid
summary.10km <- data.grid %>%
  group_by(complex, cell.10km) %>%
  summarise(n = n())

## Groups with more than 1 representative sample per grid cell
group1.10km <- summary.10km %>%
  filter(n > 1)

## Create a list of MFDO1 categories with more than 1 representative sample per grid cell
list.10km <- data.grid %>%
  group_by(complex, cell.10km) %>%
  filter(n() > 1) %>%
  select(fieldsample_barcode) %>%
  group_split()

## Groups with just 1 representative sample per grid cell
group2.10km <- data.grid %>%
  group_by(complex, cell.10km) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(fieldsample_barcode)

## Create function to select samples, which shows the largest community similarity 
## to the samples of the same MFDO1 category within the same grid cell 
filter.grid <- function(data, list) {
  out <- lst()
  
  for (i in 1:length(list)) {
    filter <- list[[i]] %>%
      pull(fieldsample_barcode)
    
    tmp <- data %>%
      select(any_of(filter), Genus) %>%
      filter(rowSums(across(where(is.numeric)))!=0) %>%
      filter(!Genus == "Unclassified") %>%
      select(where(is.numeric), Genus) %>%
      column_to_rownames(var = "Genus") %>%
      t() %>%
      decostand(., method = "hellinger") %>%
      parallelDist::parDist(., method = "bray", threads = 100) %>%
      as.matrix() %>%
      data.frame() %>%
      rownames_to_column(var = "comp") %>%
      pivot_longer(!comp, names_to = "fieldsample_barcode", values_to = "dist") %>%
      filter(!comp == fieldsample_barcode) %>%
      group_by(fieldsample_barcode) %>%
      reframe(across(dist, ~mean(1-.))) %>%
      filter(dist == max(dist)) %>%
      slice_head(n = 1)
    
    out[[i]] <- tmp
  }
  return(out)
}

## Run function on the 10 km list and bind with grids with only one representative sample
data.10km <- filter.grid(data, list.10km) %>%
  bind_rows() %>%
  select(fieldsample_barcode) %>%
  rbind(group2.10km) %>%
  select(fieldsample_barcode)

## Subset the metadata to the selected representative 10 km grid samples
data.grid.10km <- data.grid %>%
  filter(fieldsample_barcode %in% data.10km$fieldsample_barcode) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
  ungroup() %>%
  mutate(group = 1)

## Create a new summary of the selected MFDO1 categories 
groups.summary.10km <- data.grid.10km %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater")) %>%
  mutate(grid = "10 km") %>%
  mutate(across(complex, ~factor(.))) %>%
  arrange(desc(size))

## Summary of number of samples per MFDO1 category when mapping to the 1 km grid
summary.1km <- data.grid %>%
  group_by(complex, cell.1km) %>%
  summarise(n = n())

## Groups with more than 1 representative sample per grid cell
group1.1km <- summary.1km %>%
  filter(n > 1)

## Groups with just 1 representative sample per grid cell
group2.1km <- data.grid %>%
  group_by(complex, cell.1km) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(fieldsample_barcode)

## Create a list of MFDO1 categories with more than 1 representative sample per grid cell
list.1km <- data.grid %>%
  group_by(complex, cell.1km) %>%
  filter(n() > 1) %>%
  select(fieldsample_barcode) %>%
  group_split()

## Run function on the 1 km list and bind with grids with only one representative sample
data.1km <- filter.grid(data, list.1km) %>%
  bind_rows() %>%
  select(fieldsample_barcode) %>%
  rbind(group2.1km) %>%
  select(fieldsample_barcode)

## Subset the metadata to the selected representative 10 km grid samples
data.grid.1km <- data.grid %>%
  filter(fieldsample_barcode %in% data.1km$fieldsample_barcode) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
  ungroup() %>%
  mutate(group = 1)

## Create a new summary of the selected MFDO1 categories
groups.summary.1km <- data.grid.1km %>%
  group_by(complex) %>%
  summarise(size = n()) %>%
  filter(!complex %in% c("Other, Urban, Landfill",
                         "Soil, Urban, Roadside",
                         "Water, Subterranean, Freshwater",
                         "Other, Urban, Drinking water",
                         "Water, Urban, Drinking water",
                         "Soil, Subterranean, Urban",
                         "Sediment, Subterranean, Saltwater")) %>%
  mutate(grid = "1 km") %>%
  mutate(across(complex, ~factor(.))) %>%
  arrange(desc(size))


### Write the metadata files to output
data.table::fwrite(data.grid.10km, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_samples_grid_10km.tsv"))
data.table::fwrite(data.grid.1km, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_samples_grid_1km.tsv"))


### Visualise result of spatial thinning
## Combine the 10 km and 1 km group summaries
plot.data <- groups.summary %>%
  rbind(groups.summary.10km, groups.summary.1km) %>%
  mutate(across(grid, ~factor(., levels = c("All samples", "1 km", "10 km"))),
         across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  group_by(grid) %>%
  mutate(mean_size = round(mean(size), 0))

## Plot the effect of spatial thinning on the sample size pf MFDO1 groups
plot.grid <- plot.data %>%
  ggplot(aes(x = complex, y = size, fill = complex, group = grid)) +
  geom_col(aes(color = grid),
           position = "dodge") +
  geom_hline(aes(yintercept = mean_size, color = grid)) +
  scale_fill_manual(values = mfdo1.palette) +
  scale_color_manual(values = c("grey", "black", "grey20")) +
  facet_grid(cols = vars(complex), scales = "free") +
  theme_minimal() +
  xlab('MFDO1') +
  ylab('Number of samples') +
  guides(fill = guide_legend(title = "MFDO1", override.aes = list(shape = 22, size = 5, alpha = 1)),
         color = guide_legend(title = "Grid", override.aes = list(shape = 22, size = 5, fill = NA))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_blank(),
        plot.margin = margin(0,0,0,0))

## Render plot
plot.grid

## Save plot
pdf(paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_samples_MFDO1_groups-grid.pdf"))
plot.grid
dev.off()


### Create individual maps to visualise the geographical spread of the samples
## Create empty 10 km data list
var.list.10km = list()

## Extract MFDO1 categories from the combined summaries
level.groups <- plot.data %>%
  pull(complex) %>%
  unique()

## Change data to list format based on the MFDO1 categories for 10 km data
for (i in 1:length(level.groups)) {
  
  var.list.10km[[i]] = data.grid.10km[data.grid.10km$complex == level.groups[i],]
  
}

## Create empty 10 km plot list
plot.10km = list()

## Create plot for each MFDO1 category based on 10 km data
for (i in 1:length(var.list.10km)) {
  p = mapDK::mapDK() +
    theme_bw() + 
    ggtitle(str_c(level.groups[i], "\n", var.list.10km[[i]] %>% nrow(), sep = ", samples = ")) +
    geom_point(data = var.list.10km[[i]], fill = "red", pch = 21) + 
    theme(legend.position = "right", 
          axis.text.x=element_text(size=20), 
          axis.text.y=element_text(size=20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size =20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20, face = "bold"), 
          plot.title = element_text(size=20))
  
  plot.10km[[i]] = p
}

## Create empty 1 km data list
var.list.1km = list()

## Change data to list format based on the MFDO1 categories for 1 km data
for (i in 1:length(level.groups)) {
  
  var.list.1km[[i]] = data.grid.1km[data.grid.1km$complex == level.groups[i],]
  
}

## Create empty 1 km plot list
plot.1km = list()

## Create plot for each MFDO1 category based on 1 km data
for (i in 1:length(var.list.1km)) {
  p = mapDK::mapDK() +
    theme_bw() + 
    ggtitle(str_c(level.groups[i], "\n", var.list.1km[[i]] %>% nrow(), sep = ", samples = ")) +
    geom_point(data = var.list.1km[[i]], fill = "red", pch = 21) + 
    theme(legend.position = "right", 
          axis.text.x=element_text(size=20), 
          axis.text.y=element_text(size=20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size =20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20, face = "bold"), 
          plot.title = element_text(size=20))
  
  plot.1km[[i]] = p
}


### Save as multi-page pdf, with each page representing af MFDO1 category

## 10 km
pdf(paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_map_MFDO1_groups-10km-grid.pdf"))
for (i in 1:length(var.list.10km)) {
  print(plot.10km[[i]])
}
dev.off()

## 1 km
pdf(paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_map_MFDO1_groups-1km-grid.pdf"))
for (i in 1:length(var.list.1km)) {
  print(plot.1km[[i]])
}
dev.off()

  