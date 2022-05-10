# Nearctic Odonata Diversity

Code and data used to generate figures of "Diversity of Nearctic Dragonflies and Damselflies (Odonata)".
A 7zip containing a csv has all of our thresholded SDM outputs predicting the distribution of the 509 species in our study. Individual distributions could be examined by filtering this csv.

## A breif overview of what is found in each directory and subdirectory

## data
Contains data that was used to generate figures and SDM outputs.
- sdmOutputs.7z: Zip file containing sdm outputs of each species binded into a single csv.
- nObs_per_cell_aggregate15: Data used to generate Figure 8, where the number of observations were aggregated per cell. 
- rangeSize_byTraits: Data used to generate Figure 5, where the range size per species is stored. 
- regular_richnessCSVs: Subdirectory with richness csvs with richness values listed for each cell broken up by different traits or IUCN status.
- regions: Subdirectory with shapefiles used for mapping.

## Biodiverse Outputs
Overall richness and CWE values, including all species and all cells where Nearctic species are predicted to occur.

## scipts
Scripts used to generate figures

## sdmPipelineScripts
Scripts used to generate functions that were ultimately used in the sdm_pipeline.R script which was used to generate species-specific SDMs for each species of dragonfly found in the Nearctic. 
