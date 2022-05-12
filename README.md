# Nearctic Odonata Diversity

Code and data used to generate figures of "Diversity of Nearctic Dragonflies and Damselflies (Odonata)". In order to fully reproduce figures, you must first 7zip files found in the richEnd_Outputs subdirectory and the 7zip file titled regular_richnessCSVs.7z found in the /data subdirectory. 

This repository also hosts a 7zip containing a csv with all of our thresholded SDM outputs. This csv included the predicted distributions of the 509 species included in our study. Individual distributions can be examined by filtering this csv by species name. 

## A breif overview of what is found in each directory and subdirectory

## data
Contains data that was used to generate figures and SDM outputs
- sdmOutputs.7z: Zip file containing sdm outputs of each species binded into a single csv
- regular_richnessCSVs.7z: Zip with richness csvs of richness values listed for each cell broken up by different traits or IUCN status
- nObs_per_cell_aggregate15: Data used to generate Figure 8, where the number of observations were aggregated per cell
- rangeSize_byTraits: Data used to generate Figure 5, where the range size per species is stored
- regions: Subdirectory with shapefiles used for mapping

## Biodiverse Outputs
Overall richness and CWE values, including all species and all cells where Nearctic species are predicted to occur

## scipts
Scripts used to generate figures

## sdmPipelineScripts
Scripts hosting functions that were ultimately used in the sdm_pipeline.R script which was used to generate species-specific SDMs for each species of dragonfly found in the Nearctic

## Supplemental_Information
Excel datasheet hosting data associated with Supplementary Table 1 (Species Name, IUCN Status, Strictly Lotic, Strictly forest, Region)

## richEnd_Outputs
Contains four 7zip files that are needed to run /scripts/setup_forEndemismMapping.R script. 
