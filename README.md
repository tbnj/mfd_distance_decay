# mfd_distance_decay
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Metagenomic 16S data 
`scripts/maps.R` creates maps for each MFDO1 ontology level to visualise spatial separation.


`scripts/grid_spatial_thinning.R` performs spatial thinning of the dataset by selecting only one sample per cell of the 1 km and 10 km reference grid of Denmark. The sample most similar to other samples with the same MFDO1 classification within the same grid cell. 


`scripts/grid_reads_for_classification.R` combines sample ID from the 10 km spatially thinned set of samples with the read ID information to create a list of read IDs to be used in the evaluation of 16S rRNA gene reference databases. 


`scripts/distance_decay.R` performs the distance decay analysis across the five levels of the MFD Ontology. The script also renders the figures used in the manuscript. 

## Data
The scripts rely on data files available from the MFD [github](https://github.com/cmc-aau/mfd_metadata) and the MFD Zenodo [repo](https://zenodo.org/records/12605769). 
