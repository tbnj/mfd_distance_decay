### Setup env
library(tidyverse)
library(vegan)
library(parallelDist)
library(ggpubr)
library(rstatix)
library(ggpmisc)
library(patchwork)
library(hexbin)

source('scripts/MFD_colors.R')

setwd("/mfd_abundance_tables")

###!! This script is very computational heavy to run and generates two massive output files !!###

### Import data files
## Read observational table
table <- data.table::fread('data/2024-03-07_arcbac-rarefaction-rel.csv')

## Read sample metadata
metadata <- data.table::fread('data/2024-03-07_MFD_combined_metadata.csv', na.strings = "") %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub")))

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

## Subset metadata based on coordinates and reliability
metadata.sub <- metadata %>%
  filter(fieldsample_barcode %in% colnames(table),
         !is.na(mfd_hab1),
         !is.na(longitude),
         !coords_reliable == "No")

## Extract sample IDs in the subset
samples.sub <- metadata.sub %>%
  pull(fieldsample_barcode)

###!! This section can be uncommented for a lighter dataset if 2k samples !!###
#samples.sub.grid <- metadata.grid %>%
#  pull(fieldsample_barcode)

rm(metadata)


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

## Remove samples which has zero geographical distance
geo.dist.long.sub <- geo.dist.long %>%
  filter(!fieldsample_barcode == comparison) %>%
  filter(!geo_dist == 0) %>%
  left_join(metadata.sub %>% select(fieldsample_barcode, project_id, coords_reliable, sitename, sampling_date, mfd_sampletype:mfd_hab3))


### Evaluate individual projects with zero distances
## Find samples with zero distance
geo.dist.long.zero <- geo.dist.long %>%
  filter(!fieldsample_barcode == comparison) %>%
  filter(geo_dist == 0) %>%
  left_join(metadata.sub %>% select(fieldsample_barcode, project_id, coords_reliable, sitename, sampling_date, mfd_sampletype:mfd_hab3)) %>%
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
p12_5 <- readxl::read_excel("data/P12_5_minimal_metadata.xlsx") %>%
  filter(str_detect(sampling_comment, "surface")) %>%
  pull(fieldsample_barcode)

geo.dist.long.zero.p12_5 <- geo.dist.long.zero %>%
  filter(fieldsample_barcode %in% p12_5)

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
  pull(fieldsample_barcode) %>% 
  unique()

## Filter geographical distance matrix
geo.dist.filt <- geo.dist %>%
  filter(rownames(.) %in% samples.filt) %>%
  select(any_of(samples.filt))

## Write to disk
data.table::fwrite(geo.dist.filt, "output/2024-03-21_spatial-distance-filtered-subset.csv")

## Filter metadata
metadata.filt <- metadata.sub %>%
  # filter(fieldsample_barcode %in% samples.sub.grid) %>%
  filter(fieldsample_barcode %in% samples.filt) %>%
  select(fieldsample_barcode, mfd_sampletype:mfd_hab3)

## Write to disk
data.table::fwrite(metadata.filt, "output/2024-03-21_metadata-filtered-subset.csv")

## Clean
rm(geo.dist, geo.dist.long.sub, geo.dist.long.zero, geo.dist.long.zero.ok, geo.dist.long.zero.p04_3, 
   geo.dist.long.zero.p11_2, geo.dist.long.zero.p12_5, geo.dist.long.zero.p13_3, metadata.sub)


### Get total number of 16S reads in filtered dataset
table.counts <- data.table::fread('../format_data/output/2024-03-07_arcbac-rarefaction-count.csv')

## Calculate total counts
total.count <- table.counts %>%
  # select(any_of(samples.sub.grid), Genus) %>%
  select(any_of(samples.filt), Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  filter(!Genus == "Unclassified") %>%
  column_to_rownames(var = "Genus") %>%
  colSums() %>%
  sum()


### Bray-Curtis dissimilarity calculations
bray.dist <- table %>%
  # select(any_of(samples.sub.grid), Genus) %>%
  select(any_of(samples.filt), Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%
  filter(!Genus == "Unclassified") %>%
  column_to_rownames(var = "Genus") %>%
  t() %>%
  decostand(., method = "hellinger") %>%
  parallelDist::parDist(., method = "bray", threads = 100) %>%
  as.matrix() %>%
  data.frame()

data.table::fwrite(bray.dist, "output/2024-03-21_BC-distance-filtered-subset.csv")

rm(table)

## Transform to long format
geo.dist.long <- geo.dist.filt %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  pivot_longer(!fieldsample_barcode, names_to = "comp", values_to = "geo_dist") %>%
  filter(!fieldsample_barcode == comp) %>%
  mutate(across(geo_dist, ~.+0.0001))

## Change BC to long format and change from BC dissimilarity to simmilarity
bray.dist.long <- bray.dist %>%
  rownames_to_column(var = "fieldsample_barcode") %>%
  pivot_longer(!fieldsample_barcode, names_to = "comp", values_to = "bray_dist") %>%
  mutate(across(bray_dist, ~1-.)) %>%
  filter(!fieldsample_barcode == comp)

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
  filter(!c(comparison_level == "comparison_MFDO1" & is.na(sample_hab1)),
         !c(comparison_level == "comparison_MFDO1" & is.na(comp_hab1)),
         !c(comparison_level == "comparison_MFDO2" & is.na(sample_hab2)),
         !c(comparison_level == "comparison_MFDO2" & is.na(comp_hab2)),
         !c(comparison_level == "comparison_MFDO3" & is.na(sample_hab3)),
         !c(comparison_level == "comparison_MFDO3" & is.na(comp_hab3)))
         
rm(geo.dist.long, bray.dist.long)

## Set levels for types of comparisons
levels <- c("comparison_sampletype", "comparison_areatype", "comparison_MFDO1", "comparison_MFDO2", "comparison_MFDO3")

## Filter the dataset
data <- df.long.comb %>%
  mutate(across(comparison_level, ~factor(., levels = levels)),
         across(comparison_type, ~factor(., levels = c("Within", "Between")))) %>%
  filter(geo_dist <= 300)

## Evaluate the size of the dataset
data %>%
  pull(fieldsample_barcode) %>%
  unique() %>%
  length()

data %>%
  pull(comp) %>%
  unique() %>%
  length()

## Write to disk
data.table::fwrite(data, "output/2024-03-21_dd_long.csv")

## Clean env
rm(df.long.comb, metadata.grid)

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
data.sum.1km <- data %>%
  mutate(geo_bin = midcut(geo_dist, 0, 300, 1)) %>%
  group_by(geo_bin, comparison_level, comparison_type) %>%
  summarise(n = n(),
            mean = mean(bray_dist),
            sd = sd(bray_dist),
            median = median(bray_dist)) %>%
  ungroup()

## Summarise data at 5 km binsize
data.sum.5km <- data %>%
  mutate(geo_bin = midcut(geo_dist, 0, 300, 5)) %>%
  group_by(geo_bin, comparison_level, comparison_type) %>%
  summarise(n = n(),
            mean = mean(bray_dist),
            sd = sd(bray_dist),
            median = median(bray_dist)) %>%
  ungroup()

## Write summaries to disk
data.table::fwrite(data.sum.1km, "output/2024-03-21_dd-sum-binsize1.csv")
data.table::fwrite(data.sum.5km, "output/2024-03-21_dd-sum-binsize5.csv")

## Fit log decay on data nested by ontology levels and type of comparison
tidy.log <- data %>%
  nest(data = c(-comparison_level, -comparison_type)) %>%
  mutate(model = map(data, ~lm(bray_dist ~ log(geo_dist), data = .)), tidied = map(model, broom::tidy, conf.int = TRUE, conf.level = 0.95)) %>%
  unnest(tidied)

## Inspect models
tidy.log

## Create object without the model and data objects
tidy.log.df <- tidy.log %>%
  select(-c(model, data))

## Write model statisitcs to file
data.table::fwrite(tidy.log.df, "output/2024-03-21_tidy-stats.csv")


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
  mutate(x = 200,
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
data.table::fwrite(data.pred, "output/2024-03-21_dd_predictions.csv")


### Visualisation
## Plot distance-decay summarised at 5 km bins
plot.decay.sum <- ggplot(data = data.pred, aes(x = geo_dist, y = bray_dist, group = comparison_type)) +
  geom_errorbar(data = data.sum.5km, aes(x = geo_bin, y = mean, ymin = mean-sd, ymax = mean+sd), alpha = 0.25, color = "black") +
  geom_point(data = data.sum.5km, aes(x = geo_bin, y = mean, fill = comparison_type), shape = 21, color = "black", alpha = 1) +
  scale_x_continuous(breaks = seq(0, 300, 30), name = "Geographial distance, km") +
  scale_y_continuous("Bray-Curtis similarity", limits = c(-0.01, 0.6), breaks = seq(0, 0.6, 0.1)) +
  facet_grid(cols = vars(comparison_level)) +
  guides(color = guide_legend(title = "Comparison"),
         fill = guide_legend(title = "Comparison"),
         linetype = guide_legend(title = "Comparison")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank())

## Create filtered dataset for plotting with hexbin
data.all.filt <- data %>% 
  mutate(geo_bin = midcut(geo_dist, 0, 300, 1)) %>%
  group_by(geo_bin, comparison_level, comparison_type) %>%
  ungroup() %>%
  left_join(data.sum.1km, by = c("geo_bin", "comparison_level", "comparison_type"))

## Write to disk
data.table::fwrite(data.all.filt, "output/2024-03-21_dd_long_filt.csv")

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
  guides(fill = guide_colorbar(title = "Binned porportion")) +
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
  geom_density(data = data.all.filt, aes(x = bray_dist)) +
  scale_x_continuous(breaks = seq(-0.01, 1, 0.1), name = "") +
  # scale_y_continuous("Bray-Curtis similarity", limits = c(-0.01, 1), breaks = seq(0, 1, 0.1)) +
  facet_grid(rows = vars(comparison_type)) +
  coord_flip() +
  theme_bw(base_size = 18) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank())

## Customise plot layout
design <- "
11111#
222223
222223"

## Arrange plots
ggarranged <- plot.decay.sum + plot.decay.hex + plot.decay.density + 
  plot_layout(design = design, guides = "collect") & theme(legend.position = "bottom")

## Save plots to disk
png(file = 'output/dd-all.png',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

pdf(file = 'output/dd-all.pdf',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

tiff(file = 'output/dd-all.tiff',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

ggsave("output/dd-all.svg", plot = ggarranged, width = 19, height = 12, units = "in", dpi = "retina")
