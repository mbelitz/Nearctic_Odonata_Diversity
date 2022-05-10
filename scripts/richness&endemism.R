library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(sf)
library(raster)
library(rmapshaper)
library(ggspatial)

#proj
proj_str <- "+proj=laea +lon_0=-96.32 +lat_0=8.27 +datum=WGS84 +units=m +no_defs"

## na shapefile
na <- ne_countries(scale = 10, country = c("Mexico", "United States of America",
                                           "Canada"), returnclass = 'sf')
na_proj <- st_transform(na, crs = proj_str)
world <- ne_countries(scale = 110, returnclass = 'sf')

americas <- ne_countries(scale = 10, returnclass = 'sf',
                         continent = c("North America"))
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

lake <- shapefile('data/ne_10m_lakes.shp') %>%
  st_as_sf() %>%
  dplyr::filter(scalerank < 3)

nearctic_clip <- st_as_sf(nearctic) %>%
  ms_erase(., lake) %>%
  as_Spatial()

world_clip <- world %>%
  ms_erase(., lake) %>%
  as_Spatial()

americas_clip <- americas %>%
  ms_erase(., lake) %>%
  as_Spatial()

# read in richness tifs
rich <- raster("Biodiverse_Outputs/biod_Richness.tif")

# now americas
rich_mask <- mask(rich, americas_clip)
crs(rich_mask) <- crs(world)
memory.limit(size=56000)
rich_mask <- projectRaster(rich_mask, crs = proj_str)
rich_mask_df <- raster::as.data.frame(rich_mask, xy = T) %>%
  na.omit()

neotropic_proj_amer <- crop(neotropic_proj, americas_proj)
neo_proj_sf <- st_as_sf(neotropic_proj_amer)
m <- st_intersection(neo_proj_sf, americas_proj)

p1 <-ggplot() +
  geom_tile(rich_mask_df, mapping = aes(x = x , y = y,
                                        fill = biod_Richness,
                                        color = biod_Richness)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.125) +
 geom_sf(m, mapping = aes(),
         fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "plasma",
                       guide_legend(override.aes = list(size = 2),
                                    title = "Richness")) +
  scale_color_viridis_c(option = "plasma",
                        guide_legend(override.aes = list(size = 2),
                                                       title = "Richness")) +
  labs(fill = "Richness",
       color = "Richness") +
  theme_void() +
  theme(legend.position = c(0.75,0.5),
        legend.key.height = unit(0.5, "cm")) +
 coord_sf(xlim = c(-4500000, 7073990), ylim = c(-60e4, 8181452),
          crs = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ggsave(p1, filename = "figures/nearctic/totalRichness.png",
       dpi = 450, width = 6, height = 8)

# now we do the same thing but for endemism
# read in richness tifs
# read in endemism tifs
end <- raster("Biodiverse_Outputs/biod_ENDC_CWE.tif")

# now americas
end_mask <- mask(end, americas_clip)
crs(end_mask) <- crs(world)
memory.limit(size=56000)
end_mask <- projectRaster(end_mask, crs = proj_str)
end_mask_df <- raster::as.data.frame(end_mask, xy = T) %>%
  na.omit()

# plot the endemism! :D

p2 <- ggplot() +
  geom_tile(end_mask_df, mapping = aes(x = x , y = y,
                                        fill = biod_ENDC_CWE,
                                        color = biod_ENDC_CWE)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.125) +
  geom_sf(m, mapping = aes(),
          fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "turbo", trans = "log",
                       breaks = c(0.0000003, 0.000009, 0.00033),
                       labels = c("Low", "Mid", "High"),
                       guide_legend(override.aes = list(size = 2),
                                    title = "CWE"))  +
  scale_color_viridis_c(option = "turbo", trans = "log",
                       breaks = c(0.0000003, 0.000009, 0.00033),
                       labels = c("Low", "Mid", "High"),
                       guide_legend(override.aes = list(size = 2),
                                    title = "CWE"))  +
  theme_void() +
  theme(legend.position = c(0.75,0.5),
        legend.key.height = unit(0.5, "cm")) +
  coord_sf(xlim = c(-4500000, 7073990), ylim = c(-60e4, 8181452),
           crs = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


ggsave(p2, filename = "figures/nearctic/totalEndemism.png",
       dpi = 450, width = 6, height = 8)

cp <- cowplot::plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))

ggsave(cp, filename = "figures/nearctic/Figure1_cowplot.png", dpi = 450,
       width = 8, height =  3)

gp <- ggpubr::ggarrange(p1, p2, labels = c("A", "B"))

ggsave(gp, filename = "figures/nearctic/Figure1_ggarange.png", dpi = 450,
       width = 8, height =  3)
