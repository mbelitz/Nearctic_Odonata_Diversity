#' function to make recularized csvs
library(raster)
library(dplyr)

regularize <- function(x){

  csv_df <- read.csv(x) %>%
    na.omit()

  r <- raster("Biodiverse_Ouputs/biod_ENDC_CWE.tif")
  rdf <- raster::as.data.frame(r, xy = T)%>%
    mutate(z = 1:nrow(.))
  rxyz <- rdf %>%
    dplyr::select(x,y,z) %>%
    rasterFromXYZ(.)

  spdf <- csv_df
  coordinates(spdf) <- ~ x + y

  e <- extract(rxyz, spdf)

  csv_df <- mutate(csv_df, z = e)

  lj <- left_join(csv_df, rdf, by = "z") %>%
    dplyr::rename(x = x.y, y = y.y) %>%
    mutate(propEnd = CWE / biod_ENDC_CWE) %>%
    dplyr::select(x, y, z, CWE, biod_ENDC_CWE, propEnd)

  return(lj)

}

## forest
f <- regularize('richEnd_outputs/Nearctic_forest_endemism')

#non forest
nf <- regularize("richEnd_Outputs/Nearctic_nonforest_endemism")

#lentic
len <- regularize("richEnd_Outputs/Nearctic_lentic_endemism")

#lotic
lot <- regularize("richEnd_Outputs/Nearctic_lotic_endemism")
