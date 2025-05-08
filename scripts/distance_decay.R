### Setup env
library(tidyverse)
library(patchwork)

# source('scripts/MFD_colors.R')

setwd("/mfd_abundance_tables")

###!! This script is very computational heavy to run and generates two massive output files !!###

### Import data files
## Read observational table
data <- data.table::fread('data/2025-02-13_MFD_arcbac_genus_rarefaction_rel.csv')

metadata <- readxl::read_excel('data/2025-04-14_mfd_db.xlsx') %>%
  filter(fieldsample_barcode %in% colnames(data)) %>%
  relocate(coords_reliable, .after = "longitude") %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "),
         across(complex, ~factor(.))) %>%
  mutate(across(mfd_sampletype, ~factor(., levels = c("Soil", "Sediment", "Water", "Other"))),
         across(mfd_areatype, ~factor(., levels = c("Natural", "Subterranean", "Agriculture", "Urban")))) %>%
  arrange(mfd_sampletype, mfd_areatype)

## Import metadata of grid representatives and subset based on MFDO1
###!! This section can be uncommented for a lighter dataset of 2k samples !!###
#metadata.grid <- data.table::fread('output/2024-05-10_samples-grid-10km.csv', na.strings = "") %>%
#  # mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
#  #        across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
#  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
#  group_by(complex) %>%
#  mutate(complex_size = n()) %>%
#  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
#  ungroup() %>%
#  filter(!complex %in% c("Water, Subterranean, Freshwater",
#                         "Water, Urban, Sandfilter",
#                         "Soil, Urban, Other",
#                         "Sediment, Urban, Other",
#                         "Soil, Urban, Roadside",
#                         "Soil, Subterranean, Urban")) %>%
#  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
#  filter(!is.na(complex),
#         !is.na(lat)) %>%
#  select(-project_id)

## Calculate Hellinger-transformed Bray-Curtis dissimilarity in parallel
bray.dist <- data %>%
  filter(!Genus == "Unclassified") %>%
  select(where(is.numeric), Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>% 
  vegan::decostand(., method = "hellinger") %>%
  parallelDist::parDist(., method = "bray", threads = 100) %>%
  as.matrix() %>%
  data.frame()

## Save full Bray-Curtis matrix
data.table::fwrite(bray.dist, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_full.csv"))

## Subset metadata based on coordinates and reliability
metadata.sub <- metadata %>%
  filter(!is.na(mfd_hab3),
         !is.na(longitude),
         !coords_reliable == "No")

## Extract sample IDs in the subset
samples.sub <- metadata.sub %>%
  pull(fieldsample_barcode)

###!! This section can be uncommented for a lighter dataset if 2k samples !!###
#samples.sub.grid <- metadata.grid %>%
#  pull(fieldsample_barcode)

rm(metadata, data)


### Calculate geographical distance
###!! This section can be uncommented for a lighter dataset if 2k samples !!###
## Geogrphical distance
geo.dist <- metadata.sub %>%
  # filter(fieldsample_barcode %in% samples.sub.grid) %>%
  filter(fieldsample_barcode %in% samples.sub) %>%
  column_to_rownames(var = "fieldsample_barcode") %>%
  select(latitude, longitude) %>% # Order very important
  codep::geodesics(method = "haversine") %>%
  as.matrix(.) %>%
  data.frame() %>%
  mutate(across(where(is.numeric), ~./1000))

## Transform to long format
geo.dist.long <- geo.dist %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  pivot_longer(!fieldsample_barcode, names_to = "comparison", values_to = "geo_dist")


### Evaluate individual projects with zero distances
## Remove samples which has zero geographical distance
geo.dist.long.sub <- geo.dist.long %>%
  filter(!fieldsample_barcode == comparison) %>%
  filter(!geo_dist == 0) %>%
  left_join(metadata.sub %>% select(fieldsample_barcode, project_id, coords_reliable, sitename, sampling_date, sampling_comment, mfd_sampletype:mfd_hab3))

## Find samples with zero distance
geo.dist.long.zero <- geo.dist.long %>%
  filter(!fieldsample_barcode == comparison) %>%
  filter(geo_dist == 0) %>%
  left_join(metadata.sub %>% select(fieldsample_barcode, project_id, coords_reliable, sitename, sampling_date, sampling_comment, mfd_sampletype:mfd_hab3)) %>%
  arrange(project_id)

rm(geo.dist.long)

## P04_3
geo.dist.long.zero.p04_3 <- geo.dist.long.zero %>%
  filter(project_id %in% "P04_3") %>%
  group_by(sitename) %>%
  slice(which.max(as.Date(sampling_date, '%Y-%m-%d'))) %>%
  ungroup()

## P11_2
geo.dist.long.zero.p11_2 <- geo.dist.long.zero %>%
  filter(project_id %in% "P11_2") %>%
  group_by(sitename) %>%
  slice(which.max(as.Date(sampling_date, '%Y-%m-%d'))) %>%
  ungroup()

## P12_5
geo.dist.long.zero.p12_5 <- geo.dist.long.zero %>%
  filter(project_id %in% "P12_5") %>%
  filter(str_detect(sampling_comment, "surface"))

## P13_3
geo.dist.long.zero.p13_3 <- geo.dist.long.zero %>%
  filter(project_id %in% "P13_3") %>%
  group_by(sitename) %>%
  slice(which.max(as.Date(sampling_date, '%Y-%m-%d'))) %>%
  ungroup()

## List of projects to remove in this analysis
remove <- c("P04_2", "P04_4", "P06_1", "P06_2", "P06_3", "P09_3", "P09_4", "P12_1", "P16_4", "P19_1")

## Check number of samples
# tmp <- geo.dist.long.zero %>%
#   filter(!project_id %in% c("P02_2", "P03_1", "P04_2", "P04_3", "P04_4", "P04_5",
#                             "P05_1", "P05_2", "P17_1", "P06_1", "P06_2", "P06_3",
#                             "P08_1", "P08_3", "P08_7", "P09_3", "P09_4", "P11_2",
#                             "P11_3", "P12_1", "P12_5", "P13_1", "P13_2", "P13_3",
#                             "P16_1", "P16_3", "P16_4", "P18_1", "P18_2", "P19_1",
#                             "P21_1", "P25_1"))

## Samples with zero distance that are OK
geo.dist.long.zero.ok <- geo.dist.long.zero %>%
  filter(!project_id %in% remove,
         !project_id %in% c("P04_3", "P11_2", "P12_5", "P13_3")) %>%
  rbind(geo.dist.long.zero.p04_3, geo.dist.long.zero.p11_2, geo.dist.long.zero.p12_5, geo.dist.long.zero.p13_3)

## Generate list of sample IDs after removal of zero distance cases
samples.filt <- metadata.sub %>% 
  filter(!project_id %in% remove,
         !project_id %in% c("P04_3", "P11_2", "P12_5", "P13_3"),
         longitude > 8 & longitude < 13,
         latitude > 54.5 & latitude < 58) %>%
  select(fieldsample_barcode) %>%
  rbind(geo.dist.long.zero.p04_3 %>% select(fieldsample_barcode), 
        geo.dist.long.zero.p11_2 %>% select(fieldsample_barcode), 
        geo.dist.long.zero.p12_5 %>% select(fieldsample_barcode), 
        geo.dist.long.zero.p13_3 %>% select(fieldsample_barcode)) %>%
  arrange(fieldsample_barcode) %>%
  pull(fieldsample_barcode) %>% 
  unique()

## Filter geographical distance matrix
geo.dist.filt <- geo.dist %>%
  filter(rownames(.) %in% samples.filt) %>%
  select(all_of(samples.filt))

bray.dist.filt <- bray.dist %>%
  filter(rownames(.) %in% samples.filt) %>%
  select(all_of(samples.filt))

## Write to disk
data.table::fwrite(geo.dist.filt, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_spatial_distance_filtered.csv"))

data.table::fwrite(bray.dist.filt, sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_BC_distance_filtered.csv"))


### Filter metadata
metadata.filt <- metadata.sub %>%
  # filter(fieldsample_barcode %in% samples.sub.grid) %>%
  filter(fieldsample_barcode %in% samples.filt) %>%
  select(fieldsample_barcode, mfd_sampletype:mfd_hab3) %>%
  arrange(fieldsample_barcode)

## Write to disk
data.table::fwrite(metadata.filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_metadata_dd_filtered_subset.tsv"))

## Clean
rm(geo.dist, geo.dist.long.sub, geo.dist.long.zero, geo.dist.long.zero.ok, geo.dist.long.zero.p04_3, 
   geo.dist.long.zero.p11_2, geo.dist.long.zero.p12_5, geo.dist.long.zero.p13_3, metadata.sub, bray.dist)


### Filter comparisons to self and duplicates + make long format
## Cange spatial to long format and ad pseudo-count 0f 0.1 m
geo.dist.long <- geo.dist.filt %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  pivot_longer(!fieldsample_barcode, names_to = "comp", values_to = "geo_dist") %>%
  filter(!fieldsample_barcode == comp) %>%
  mutate(tmp1 = str_remove(fieldsample_barcode, "MFD"),
         tmp2 = str_remove(comp, "MFD"),
         across(tmp1:tmp2, ~as.numeric(.)),
         tmp3 = case_when(tmp1 < tmp2 ~ tmp1,
                          tmp1 > tmp2 ~ tmp2),
         tmp4 = case_when(tmp1 > tmp2 ~ tmp1,
                          tmp1 < tmp2 ~ tmp2),
         across(fieldsample_barcode, ~str_c("MFD", str_pad(tmp3, width = 5, side = "left", pad = "0"))),
         across(comp, ~str_c("MFD", str_pad(tmp4, width = 5, side = "left", pad = "0")))) %>%
  select(-c(tmp1:tmp4)) %>%
  distinct() %>%
  mutate(across(geo_dist, ~.+0.0001)) %>%
  arrange(fieldsample_barcode, comp)


## Change BC to long format and change from BC dissimilarity to simmilarity
bray.dist.long <- bray.dist.filt %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  pivot_longer(!fieldsample_barcode, names_to = "comp", values_to = "bray_dist") %>%
  mutate(across(bray_dist, ~1-.)) %>%
  filter(!fieldsample_barcode == comp) %>%
  mutate(tmp1 = str_remove(fieldsample_barcode, "MFD"),
         tmp2 = str_remove(comp, "MFD"),
         across(tmp1:tmp2, ~as.numeric(.)),
         tmp3 = case_when(tmp1 < tmp2 ~ tmp1,
                          tmp1 > tmp2 ~ tmp2),
         tmp4 = case_when(tmp1 > tmp2 ~ tmp1,
                          tmp1 < tmp2 ~ tmp2),
         across(fieldsample_barcode, ~str_c("MFD", str_pad(tmp3, width = 5, side = "left", pad = "0"))),
         across(comp, ~str_c("MFD", str_pad(tmp4, width = 5, side = "left", pad = "0")))) %>%
  select(-c(tmp1:tmp4)) %>%
  distinct() %>%
  arrange(fieldsample_barcode, comp)

## Combine BC and spatial distances
df.long.comb <- bray.dist.long %>%
  left_join(geo.dist.long) %>%
  left_join(metadata.filt, by = "fieldsample_barcode") %>%
  rename(sample_sampletype = mfd_sampletype,
         sample_areatype = mfd_areatype,
         sample_hab1 = mfd_hab1,
         sample_hab2 = mfd_hab2,
         sample_hab3 = mfd_hab3) %>%
  left_join(metadata.filt, by = c("comp" = "fieldsample_barcode")) %>%
  rename(comp_sampletype = mfd_sampletype,
         comp_areatype = mfd_areatype,
         comp_hab1 = mfd_hab1,
         comp_hab2 = mfd_hab2,
         comp_hab3 = mfd_hab3) %>%
  mutate(comparison_sampletype = case_when(sample_sampletype == comp_sampletype ~ "Within",
                                           TRUE ~ "Between"),
         comparison_areatype = case_when(sample_areatype == comp_areatype & sample_sampletype == comp_sampletype ~ "Within",
                                         TRUE ~ "Between"),
         comparison_MFDO1 = case_when(sample_areatype == comp_areatype & sample_sampletype == comp_sampletype & sample_hab1 == comp_hab1 ~ "Within",
                                      TRUE ~ "Between"),
         comparison_MFDO2 = case_when(sample_areatype == comp_areatype & sample_sampletype == comp_sampletype & sample_hab1 == comp_hab1 & sample_hab2 == comp_hab2 ~ "Within",
                                      TRUE ~ "Between"),
         comparison_MFDO3 = case_when(sample_areatype == comp_areatype & sample_sampletype == comp_sampletype & sample_hab1 == comp_hab1 & sample_hab2 == comp_hab2 & sample_hab3 == comp_hab3 ~ "Within",
                                      TRUE ~ "Between")) %>%
  pivot_longer(!fieldsample_barcode:comp_hab3, names_to = "comparison_level", values_to = "comparison_type") %>%
  filter(!c(comparison_level == "comparison_sampletype" & is.na(sample_sampletype)),
         !c(comparison_level == "comparison_sampletype" & is.na(comp_sampletype)),
         !c(comparison_level == "comparison_areatype" & is.na(sample_areatype)),
         !c(comparison_level == "comparison_areatype" & is.na(comp_areatype)),
         !c(comparison_level == "comparison_MFDO1" & is.na(sample_hab1)),
         !c(comparison_level == "comparison_MFDO1" & is.na(comp_hab1)),
         !c(comparison_level == "comparison_MFDO2" & is.na(sample_hab2)),
         !c(comparison_level == "comparison_MFDO2" & is.na(comp_hab2)),
         !c(comparison_level == "comparison_MFDO3" & is.na(sample_hab3)),
         !c(comparison_level == "comparison_MFDO3" & is.na(comp_hab3)))

rm(geo.dist.long, bray.dist.long)

## Set levels for types of comparisons
levels <- c("comparison_sampletype", "comparison_areatype", "comparison_MFDO1", "comparison_MFDO2", "comparison_MFDO3")

## Filter the dataset
model.data <- df.long.comb %>%
  mutate(across(comparison_level, ~factor(., levels = levels)),
         across(comparison_type, ~factor(., levels = c("Within", "Between")))) %>%
  filter(geo_dist <= 300)

## Clean env
rm(df.long.comb)

gc()


### Distance-decay
## Write summary function to use with bins
midcut <- function(x, from, to, by){
  ## cut the data into bins...
  x=cut(x,seq(from,to,by),include.lowest=T)
  ## make a named vector of the midpoints, names=binnames
  vec=seq(from+by/2,to-by/2,by)
  names(vec)=levels(x)
  ## use the vector to map the names of the bins to the midpoint values
  unname(vec[x])
}

## Summarise data at 1 km binsize
data.sum.1km <- model.data %>%
  mutate(geo_bin = midcut(geo_dist, 0, 300, 1)) %>%
  group_by(geo_bin, comparison_level, comparison_type) %>%
  summarise(n = n(),
            mean = mean(bray_dist),
            sd = sd(bray_dist),
            median = median(bray_dist)) %>%
  ungroup()

## Summarise data at 5 km binsize
data.sum.5km <- model.data %>%
  mutate(geo_bin = midcut(geo_dist, 0, 300, 5)) %>%
  group_by(geo_bin, comparison_level, comparison_type) %>%
  summarise(n = n(),
            mean = mean(bray_dist),
            sd = sd(bray_dist),
            median = median(bray_dist)) %>%
  ungroup()

## Write summaries to disk
data.table::fwrite(data.sum.1km, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_dd_sum_binsize1km.tsv"))
data.table::fwrite(data.sum.5km, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_dd_sum_binsize5km.tsv"))

## Fit log decay on data nested by ontology levels and type of comparison
tidy.log <- model.data %>%
  nest(data = c(-comparison_level, -comparison_type)) %>%
  mutate(model = map(data, ~lm(bray_dist ~ log(geo_dist), data = .)), tidied = map(model, broom::tidy, conf.int = TRUE, conf.level = 0.95)) %>%
  unnest(tidied)

## Inspect models
tidy.log

## Create object without the model and data objects
tidy.log.df <- tidy.log %>%
  select(-c(model, data))

## Write model statisitcs to file
data.table::fwrite(tidy.log.df, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_dd_model_stats.tsv"))


## Write summary function to calculate R^2 for the models
r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

## Create empty data list
list <- lst()

## Loop over the models and predict based on the 5 km binned data
for (i in 1:nrow(tidy.log.df)) {
  
  new.data <- data.sum.5km %>%
    filter(comparison_level == tidy.log$comparison_level[i],
           comparison_type == tidy.log$comparison_type[i]) %>%
    rename(geo_dist = geo_bin,
           bray_dist = mean)
  
  pred.data <- predict(tidy.log$model[[i]], newdata = new.data)
  
  res <- cbind(new.data, pred.data) %>%
    rename(prediction = pred.data) %>%
    group_by(comparison_level, comparison_type) %>%
    summarise(r_squared = r_squared(bray_dist, prediction))
  
  list[[i]] <- res
}

## Bind results into one dataframe
result <- bind_rows(list)

## Combine with 5 km summary
data.sum.5km <- data.sum.5km %>%
  left_join(distinct(result))

## Make wide object with model terms
models <- tidy.log.df %>%
  select(1:4) %>%
  group_by(comparison_level, comparison_type) %>%
  pivot_wider(names_from = "term", values_from = "estimate") %>%
  rename(a = `log(geo_dist)`,
         b = `(Intercept)`) %>%
  ungroup() %>%
  mutate(x = 110,
         y = 0.95)

## Create dataframe of predicted values of BC similarity
data.pred <- models %>%
  group_by(comparison_level, comparison_type) %>%
  mutate(geo_dist = map2(0, 300, ~ seq(from = .x, to = .y, by = 1))) %>%
  unnest(cols = c(geo_dist)) %>%
  ungroup() %>%
  mutate(across(geo_dist, ~case_when(geo_dist == 0 ~.+0.0001,
                                     TRUE ~ geo_dist))) %>%
  mutate(bray_dist = a * log(geo_dist) + b)

## Write to disk
data.table::fwrite(data.pred, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_dd_model_predictions.tsv"))



### Visualisation
## Plot distance-decay summarised at 5 km bins
plot.decay.sum <- ggplot(data = data.pred, aes(x = geo_dist, y = bray_dist, group = comparison_type)) +
  geom_errorbar(data = data.sum.5km, aes(x = geo_bin, y = mean, ymin = mean-sd, ymax = mean+sd), alpha = 0.25, color = "black") +
  geom_point(data = data.sum.5km, aes(x = geo_bin, y = mean, fill = comparison_type), shape = 21, color = "black", alpha = 1) +
  scale_x_continuous(breaks = seq(0, 300, 30), name = "Geographial distance, km") +
  scale_y_continuous("Bray-Curtis similarity", limits = c(-0.01, 0.6), breaks = seq(0, 0.6, 0.1)) +
  facet_grid(cols = vars(comparison_level),
             labeller = labeller(comparison_level = c("comparison_sampletype" = "Sample type", 
                                                      "comparison_areatype" = "Area type", 
                                                      "comparison_MFDO1" = "MFDO1", 
                                                      "comparison_MFDO2" = "MFDO2", 
                                                      "comparison_MFDO3" = "MFDO3"))) +
  guides(color = guide_legend(title = "Comparison data"),
         fill = guide_legend(title = "Comparison data", override.aes = list(size = 4))) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank())

## Create filtered dataset for plotting with hexbin
data.all.filt <- model.data %>% 
  mutate(geo_bin = midcut(geo_dist, 0, 300, 1)) %>%
  # group_by(geo_bin, comparison_level, comparison_type) %>%
  # ungroup() %>%
  left_join(data.sum.1km, by = c("geo_bin", "comparison_level", "comparison_type"))

## Write to disk
data.table::fwrite(data.all.filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_dd_long_format.tsv"))

rm(model.data)

gc()

## Plot distance-decay with hexbin including fitted models
plot.decay.hex <- ggplot(data = data.pred, aes(x = geo_dist, y = bray_dist)) +
  geom_hex(data = data.all.filt, aes(x = geo_dist, y = bray_dist, weight = 1/n), bins = 300) +
  geom_line(data = data.pred, aes(x = geo_dist, y = bray_dist)) +
  geom_text(data = models, aes(x = x, y = y,
                               label = paste('italic(y)~`=`~', round(b, 3), '~', round(a, 3), '~log(italic(x))')), 
            parse = TRUE, size = 6, color = "black") +
  scale_x_continuous(breaks = seq(0, 300, 30), name = "Geographial distance, km") +
  scale_y_continuous("Bray-Curtis similarity", limits = c(-0.01, 1), breaks = seq(0, 1, 0.1)) +
  scale_fill_viridis_c(option = "magma", transform = "log10", breaks = c(0.00001, 0.0001, 0.001, 0.01)) +
  facet_grid(cols = vars(comparison_level),
             rows = vars(comparison_type)) +
  guides(fill = guide_colorbar(title = "Binned proportion", theme(legend.title.position = "top"))) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_blank(),
        panel.grid.minor.x = element_blank())

## Create density curves to visualise the data distributions for comparions within and between
## the different levels of the ontology
plot.decay.density <- ggplot() +
  geom_density(data = data.all.filt, aes(x = bray_dist, linetype = comparison_level), show.legend = FALSE) +
  stat_density(data = data.all.filt, aes(x = bray_dist, linetype = comparison_level),
               geom = "line", position = "identity", linewidth = 0) + 
  guides(linetype = guide_legend(title = "Comparison type", 
                                 labels = c("Sample type", "Area type", "MFDO1", "MFDO2", "MFDO3"),
                                 theme(legend.title.position = "top"), 
                                 override.aes = list(linewidth = 1))) +
  scale_x_continuous(limits = c(-0.01, 1), breaks = seq(0, 1, 0.1), name = "") +
  # scale_y_continuous("Bray-Curtis similarity", limits = c(-0.01, 1), breaks = seq(0, 1, 0.1)) +
  facet_grid(rows = vars(comparison_type)) +
  coord_flip() +
  theme_bw(base_size = 18) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank())

## Customise plot layout
design <- "
111114
222223
222223"

## Arrange plots
ggarranged <- plot.decay.sum + plot.decay.hex + plot.decay.density + guide_area() + 
  plot_layout(design = design, guides = "collect")

gc()

## Save plots to disk
png(file = 'output/distance_decay.png',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

pdf(file = 'output/distance_decay.pdf',
    width = 19,
    height = 12)
ggarranged
dev.off()

tiff(file = 'output/distance_decay.tiff',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

jpeg(file = 'output/distance_decay.jpeg',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

ggsave("output/distance_decay.svg", plot = ggarranged, width = 19, height = 12, units = "in", dpi = "retina")

# save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_MFD_distance_decay.RData"))
