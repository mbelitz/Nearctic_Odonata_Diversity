# Run lines 3 & 4 if an older version of sf is needed, newer versions may cause errors

#purl <- "http://cran.us.r-project.org/src/contrib/Archive/sf/sf_0.9-6.tar.gz"
#install.packages(purl, repos=NULL, type="source")

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

na_clip <- na %>%
  ms_erase(., lake) %>%
  as_Spatial()

## now on to IUCN listing groups
listed <- read.csv("data/regular_richnessCSVs/listed_richness.csv")
lc <- read.csv("data/regular_richnessCSVs/leastConcern_richness.csv")

## mapping options
# fist Americas extent
listed_rich <- dplyr::select(listed, x,y,richness) %>%
  rasterFromXYZ()
listed_rich_mask <- mask(listed_rich, americas_clip)
crs(listed_rich_mask) <- crs(world)
memory.limit(size=56000)
listed_rich_mask <- projectRaster(listed_rich_mask, crs = proj_str)

neotropic_proj_amer <- crop(neotropic_proj, americas_proj)
neo_proj_sf <- st_as_sf(neotropic_proj_amer)
m <- st_intersection(neo_proj_sf, americas_proj)


t_plot <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(listed_rich_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = richness, color = richness)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.125) +
  geom_sf(m, mapping = aes(),
          fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "plasma",
                       guide_legend(override.aes = list(size = 2),
                                    title = "Richness")) +
  scale_color_viridis_c(option = "plasma",
                        guide_legend(override.aes = list(size = 2),
                                     title = "Richness")) +
  coord_sf(xlim = c(-4500000, 7073990), ylim = c(20e4, 7581452)) +
  labs(fill = "Richness",
       color = "Richness") +
  ggtitle("Threatened") +
  theme_void() +
  theme(legend.position = c(0.75,0.5),
        legend.key.height = unit(0.5, "cm"))

ggsave(filename = "figures/richness_threatened_americasExtent.png",
       dpi = 450,  width = 5, height = 7)

# now do least concern
lc_rich <- dplyr::select(lc, x,y,richness) %>%
  rasterFromXYZ()
lc_rich_mask <- mask(lc_rich, americas_clip)
crs(lc_rich_mask) <- crs(world)
memory.limit(size=56000)
lc_rich_mask <- projectRaster(lc_rich_mask, crs = proj_str)

lc_plot <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(lc_rich_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = richness, color = richness)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75", size = 0.125) +
  geom_sf(m, mapping = aes(),
          fill = "grey55", alpha = 0.65, color = NA) +
  scale_fill_viridis_c(option = "plasma",
                       guide_legend(override.aes = list(size = 2),
                                    title = "Richness")) +
  scale_color_viridis_c(option = "plasma",
                        guide_legend(override.aes = list(size = 2),
                                     title = "Richness")) +
  coord_sf(xlim = c(-4500000, 7073990), ylim = c(20e4, 7581452)) +
  labs(fill = "Richness",
       color = "Richness") +
  ggtitle("Least concern") +
  theme_void() +
  theme(legend.position = c(0.75,0.5),
        legend.key.height = unit(0.5, "cm"))

ggsave(filename = "figures/richness_leastConcern_americasExtent.png",
       dpi = 450,  width = 5, height = 7)



## plot them all together
cp <- cowplot::plot_grid(lc_plot, t_plot,
                         nrow = 1, ncol = 2, labels = c("A", "B"))

ggsave(cp, filename = "figures/nearctic/Figure3_IUCN.png",
       dpi = 450, height = 3.5, width = 8)
