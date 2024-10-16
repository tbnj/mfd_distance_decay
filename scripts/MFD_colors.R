library(wesanderson)

setwd('/mfd_distance_decay')

metadata <- readxl::read_excel('data/2024-02-13_mfd_db.xlsx') %>%
  select(mfd_sampletype:mfd_hab1)

sample.levels <- metadata %>%
  select(mfd_sampletype) %>%
  filter(!is.na(mfd_sampletype)) %>%
  distinct() %>%
  mutate(complex = mfd_sampletype) %>%
  arrange(mfd_sampletype) %>%
  pull(complex)

area.levels <- metadata %>%
  select(mfd_sampletype:mfd_areatype) %>%
  filter(!is.na(mfd_areatype)) %>%
  distinct() %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, sep = ", ")) %>%
  arrange(mfd_sampletype, mfd_areatype) %>%
  pull(complex)

mfdo1.levels <- metadata %>%
  select(mfd_sampletype:mfd_hab1) %>%
  filter(!is.na(mfd_hab1)) %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllus scrub", "Temperate heath and scrub"))) %>%
  distinct() %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  arrange(mfd_sampletype, mfd_areatype, mfd_hab1) %>%
  pull(complex)

rm(metadata)

## Ontology palettes
sediment.palette <- colorRampPalette(c(wes_palette("IsleofDogs2")[3], wes_palette("FantasticFox1")[1]))
plot(rep(1, 6), col = sediment.palette(6), pch = 19, cex = 3)
soil.palette <- colorRampPalette(c(wes_palette("AsteroidCity1")[4], wes_palette("AsteroidCity1")[1]))
plot(rep(1, 12), col = soil.palette(12), pch = 19, cex = 3)
water.palette <- colorRampPalette(c(wes_palette("Darjeeling2")[2], wes_palette("Zissou1")[2]))
plot(rep(1, 6), col = water.palette(6), pch = 19, cex = 3)

## Sampletype palette
sampletype.palette <- c(sediment.palette(1), soil.palette(1), water.palette(1))
names(sampletype.palette) <- sample.levels

sampletype.palette
plot(rep(1, 3), col = sampletype.palette, pch = 19, cex = 3)

## Areatype palette
areatype.palette <- c(sediment.palette(3), soil.palette(4), water.palette(3))
names(areatype.palette) <- area.levels

sampletype.palette
plot(rep(1, 3), col = sampletype.palette, pch = 19, cex = 3)

# MFDO1 palette
mfdo1.palette <- c(sediment.palette(6), soil.palette(12), water.palette(6))
names(mfdo1.palette) <- mfdo1.levels

mfdo1.palette
plot(rep(1, 24), col = mfdo1.palette, pch = 19, cex = 3)









