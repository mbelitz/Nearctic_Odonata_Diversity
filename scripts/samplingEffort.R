library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(sf)
library(raster)
library(rmapshaper)

#proj
proj_str <- "+proj=laea +lon_0=-96.32 +lat_0=8.27 +datum=WGS84 +units=m +no_defs"

## na shapefile
na <- ne_countries(scale = 10, country = c("Mexico", "United States of America",
                                           "Canada"), returnclass = 'sf')
na_proj <- st_transform(na, crs = proj_str)
world <- ne_countries(scale = 110, returnclass = 'sf')

americas <- ne_countries(scale = 10, returnclass = 'sf', 
                         continent = c("North America", "South America"))
americas_proj <- st_transform(americas, crs = proj_str)

# nearctic shapefile
nearctic <- shapefile("data/regions/tnc_terr_ecoregions.shp") %>% 
  st_as_sf() %>% 
  dplyr::filter(WWF_REALM2 == "Nearctic") 

nearctic <- nearctic %>%
  st_union() %>% 
  as_Spatial()

nearctic_proj <- spTransform(nearctic,CRSobj = proj_str)

neotropic <- shapefile("data/regions/tnc_terr_ecoregions.shp") %>% 
  st_as_sf() %>% 
  dplyr::filter(WWF_REALM2 == "Neotropic") 

neotropic <- neotropic %>%
  st_union() %>% 
  as_Spatial()

neotropic_proj <- spTransform(neotropic, CRSobj = proj_str)

#read in sampling effort data
se <- read.csv("data/nObs_per_cell_aggregate15.csv") %>% 
  na.omit()

#filter se to north america
se_xyz <- dplyr::select(se, x,y,nObs) %>% 
  rasterFromXYZ()
se_mask <- mask(se_xyz, na_proj)
crs(se_mask) <- crs(na_proj)
se_mask <- projectRaster(se_mask, crs = proj_str)

rdf <- raster::as.data.frame(se_mask, xy = T) %>% 
  filter(!is.na(nObs))

neotropic_proj_amer <- crop(neotropic_proj, na_proj)
neo_proj_sf <- st_as_sf(neotropic_proj_amer)
m <- st_intersection(neo_proj_sf, na_proj)


ggplot() +
  geom_tile(rdf, mapping = aes(x = x, y = y, fill = nObs)) +
  geom_sf(na_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.15) +
  geom_sf(m, mapping = aes(), 
          fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "mako", trans = "log",
                       breaks = c(1,10,200,3000),begin = 0.2) +
  coord_sf(c(-4500000, 3073990)) +
  labs(fill = "Occurrences") +
  theme_void()


ggsave(filename = "figures/nearctic/Figure4_nearctic.png", dpi = 450,
       width = 5, height = 5)


ggplot() +
  geom_tile(rdf, mapping = aes(x = x, y = y, fill = nObs)) +
  geom_sf(na_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.15) +
  geom_sf(m, mapping = aes(), 
          fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "inferno", trans = "log",
                       breaks = c(1,10,200,3000),begin = 0.2) +
  coord_sf(c(-4500000, 3073990)) +
  labs(fill = "Occurrences") +
  theme_void()


ggsave(filename = "figures/nearctic/Figure4_nearctic_infernoPalette.png", dpi = 450,
       width = 5, height = 5)
