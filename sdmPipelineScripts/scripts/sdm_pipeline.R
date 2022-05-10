library(raster)
library(dplyr)
library(rnaturalearth)

setwd("/srv/mybook/mbelitz/Ode_SDMs/")

source("scripts/03_define_accessibleArea.R")
source("scripts/04_clip_modelLayers.R")
source("scripts/05_select_modelVariables.R")
source("scripts/07_save_SDMoutputs_LPT0.R")
source("scripts/08_save_SDMoutputs_TSS.R")

## read in data
nea <- data.table::fread("Occ_Data/nearctic_occurrenceRecords_12-15-21.csv")

nea <- nea %>%
  dplyr::filter(Note.for.Mike != "Remove") %>%
  dplyr::filter(Note.for.Mike != "remove") %>%
  dplyr::filter(!is.na(decimalLongitude)) %>%
  dplyr::select(scientific_name, decimalLongitude, decimalLatitude,
                validName, source)

palea <- data.table::fread("Occ_Data/AAADBMAPSPalaearctic20Oct21.txt")
palea <- palea %>%
  dplyr::rename(scientific_name = "Species name GEODE",
                decimalLongitude = "Longitude",
                decimalLatitude = "Latitude",
                source = "Source") %>%
  mutate(validName = stringr::str_replace(string = scientific_name, pattern = "_", replacement = " ")) %>%
  dplyr::select(scientific_name, decimalLongitude, decimalLatitude,
                validName, source)

occs <- rbind(nea, palea)

##  load pipeline
run_spp_pipeline <- function(binomial){

  cleanedOccs <- occs %>%
    dplyr::filter(validName == binomial) %>%
    na.omit()

  aa_shp <- define_accessibleArea(species_df = cleanedOccs, minBuff = 75000,
                                  saveImage = F,
                                  saveShapefile = FALSE)

  mod_vars <- clip_variableLayers(layerDir = "ClimateOnly/",
                                  accessibleArea = aa_shp)

  mod_vars <- aggregate(mod_vars, 5) #use this step to upaggregate resolution

  spp_df <- cleanedOccs

  coordinates(spp_df) <- ~ decimalLongitude + decimalLatitude

  ## spatial thinning
  area_sqkm <- raster::area(aa_shp)/1000000

  # my guess for now is 25km thinning for > 100,000km2, no thinning below that,
  # 50 km thinning for > 100,000km2 but < 1,000,000km2 and
  # 100 km thinning for > 1,000,000km2 but < 2,500,000km2 and
  # 200 km thinning for > 2,500,000km2
  if(area_sqkm < 100000){

    spp_df <- spp_df

  } else if(area_sqkm >= 100000 & area_sqkm < 250000){

    bio12 <- mod_vars[[2]]
    bio12 <- aggregate(bio12, 5)
    bio12mdf <- raster::as.data.frame(bio12, xy = T)
    bio12mdf_noNA <- na.omit(bio12mdf) %>%
      dplyr::rename(z = 3)
    bio12mdf_noNA <- mutate(bio12mdf_noNA, z = 1:nrow(bio12mdf_noNA))
    bio12_2 <- rasterFromXYZ(bio12mdf_noNA)

    e <- raster::extract(bio12_2, spp_df)

    spp_df$cell_id <- e
    spp_df_thinned <- st_as_sf(spp_df) %>%
      dplyr::group_by(cell_id) %>%
      slice(1)
    spp_df <- as_Spatial(spp_df_thinned)

  } else if(area_sqkm >= 250000 & area_sqkm < 1000000){

    bio12 <- mod_vars[[2]]
    bio12 <- aggregate(bio12, 10)
    bio12mdf <- raster::as.data.frame(bio12, xy = T)
    bio12mdf_noNA <- na.omit(bio12mdf) %>%
      dplyr::rename(z = 3)
    bio12mdf_noNA <- mutate(bio12mdf_noNA, z = 1:nrow(bio12mdf_noNA))
    bio12_2 <- rasterFromXYZ(bio12mdf_noNA)

    e <- raster::extract(bio12_2, spp_df)

    spp_df$cell_id <- e
    spp_df_thinned <- st_as_sf(spp_df) %>%
      dplyr::group_by(cell_id) %>%
      slice(1)
    spp_df <- as_Spatial(spp_df_thinned)

  } else if (area_sqkm >= 1000000 & area_sqkm < 2500000){

    bio12 <- mod_vars[[2]]
    bio12 <- aggregate(bio12, 20)
    bio12mdf <- raster::as.data.frame(bio12, xy = T)
    bio12mdf_noNA <- na.omit(bio12mdf) %>%
      dplyr::rename(z = 3)
    bio12mdf_noNA <- mutate(bio12mdf_noNA, z = 1:nrow(bio12mdf_noNA))
    bio12_2 <- rasterFromXYZ(bio12mdf_noNA)

    e <- raster::extract(bio12_2, spp_df)

    spp_df$cell_id <- e
    spp_df_thinned <- st_as_sf(spp_df) %>%
      dplyr::group_by(cell_id) %>%
      slice(1)
    spp_df <- as_Spatial(spp_df_thinned)

  } else{

    bio12 <- mod_vars[[2]]
    bio12 <- aggregate(bio12, 40)
    bio12mdf <- raster::as.data.frame(bio12, xy = T)
    bio12mdf_noNA <- na.omit(bio12mdf) %>%
      dplyr::rename(z = 3)
    bio12mdf_noNA <- mutate(bio12mdf_noNA, z = 1:nrow(bio12mdf_noNA))
    bio12_2 <- rasterFromXYZ(bio12mdf_noNA)

    e <- raster::extract(bio12_2, spp_df)

    spp_df$cell_id <- e
    spp_df_thinned <- st_as_sf(spp_df) %>%
      dplyr::group_by(cell_id) %>%
      slice(1)
    spp_df <- as_Spatial(spp_df_thinned)

  }


  max_model <- maxent(x = mod_vars, p = coordinates(spp_df),
                      progress = "text")

  predictors <- select_sdmVariables(pred_vars = mod_vars,
                                    maxent_mod = max_model,
                                    maxVIF = 5)

  eval1 <- ENMeval::ENMevaluate(occ = coordinates(spp_df),
                                env = predictors,
                                method = "block",
                                RMvalues = c(0.5, 1, 2, 3, 4),
                                fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                                parallel = TRUE, numCores = 5,
                                algorithm = 'maxent.jar')

  bw <- stringr::str_replace(string = binomial, pattern = " ", replacement = "_")
  fp <- file.path("Fixed_Sdms/Nearctic_sdmOutputs_R2/")
  fp2 <- file.path("Fixed_Sdms/Nearctic_sdmMaps_R2/")
  dir.create(paste(fp, bw, sep = ""))

  spp_df2 <- as.data.frame(spp_df) %>%
    dplyr::rename(decimalLongitude = coords.x1,
                  decimalLatitude = coords.x2)

  save_SDM_results(ENMeval_output = eval1,
                   AUCmin = 0.7,
                   resultDir = paste(fp, bw, "/", sep = ""),
                   spp = binomial,
                   occ_df = spp_df2)

  save(eval1,
       file = paste(fp, bw, "/",
                    bw, "_ENMeval", ".RData", sep = ""))

  p <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_sf(st_as_sf(aa_shp), mapping = aes(), fill = "orange", alpha = 0.5) +
    geom_point(cleanedOccs, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               shape = 1, size = 0.75) +
    coord_sf(xlim = c(min(cleanedOccs$decimalLongitude) - 5, max(cleanedOccs$decimalLongitude) + 5),
             ylim = c(min(cleanedOccs$decimalLatitude) - 5, max(cleanedOccs$decimalLatitude) + 5)) +
    ggtitle(paste(binomial, "n =", length(cleanedOccs$validName), sep = " ")) +
    theme_classic()

  r <- raster(paste(fp,bw,"/",bw,"_SDM.tif", sep = ""))
  p2 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = ClogLog)) +
    coord_sf(xlim = c(min(cleanedOccs$decimalLongitude) - 5, max(cleanedOccs$decimalLongitude) + 5),
             ylim = c(min(cleanedOccs$decimalLatitude) - 5, max(cleanedOccs$decimalLatitude) + 5)) +
    scale_fill_viridis_c() +
    theme_classic()

  r2 <- raster(paste(fp,bw,"/",bw,"_SDM_PA.tif", sep = ""))
  p3 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r2,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog))) +
    coord_sf(xlim = c(min(cleanedOccs$decimalLongitude) - 5, max(cleanedOccs$decimalLongitude) + 5),
             ylim = c(min(cleanedOccs$decimalLatitude) - 5, max(cleanedOccs$decimalLatitude) + 5)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    theme_classic()

  e <- egg::ggarrange(p2, p3, p, nrow = 2,
                      heights = c(1,0.65),
                      draw = F)

  ggsave(plot = e,
         filename = paste(fp2, bw, "_sdmMap.png", sep = ""),
         width = 8,
         height = 6)


  bp <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r2,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(cleanedOccs, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               size = 0.5, alpha = 0.75, shape = 1) +
    coord_sf(xlim = c(min(cleanedOccs$decimalLongitude) - 10, max(cleanedOccs$decimalLongitude) + 10),
             ylim = c(min(cleanedOccs$decimalLatitude) - 10, max(cleanedOccs$decimalLatitude) + 10)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle(paste(binomial, "n =", length(cleanedOccs$validName), sep = " ")) +
    theme_classic()

  ggsave(plot = bp, filename = paste(fp2, bw, "_sdmMapOccPoints.png", sep = ""))

  rm(p,r,p2,r2,p3,e,bp)
  gc()

  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()

}
