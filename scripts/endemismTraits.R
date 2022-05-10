# Run lines 3 & 4 if an older version of sf is needed, newer versions may cause errors

#purl <- "http://cran.us.r-project.org/src/contrib/Archive/sf/sf_0.9-6.tar.gz"
#install.packages(purl, repos=NULL, type="source")

library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(sf)
library(raster)
library(rmapshaper)
library(biscale)
library(cowplot)

#proj
proj_str <- "+proj=laea +lon_0=-96.32 +lat_0=8.27 +datum=WGS84 +units=m +no_defs"

## na shapefile
na <- ne_countries(scale = 10, country = c("Mexico", "Unitied States of America",
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

## read in forest only data
source('scripts/regularize_endemismCSVs.R')

## mapping options
# then western hemisphere
#forest
f_end <- dplyr::select(f, x,y, CWE) %>%
  rasterFromXYZ()
#non forest
nf_end <- dplyr::select(nf, x,y, CWE) %>%
  rasterFromXYZ()

# then nearctic only
#forest
f_end_mask <- mask(f_end, nearctic_clip)
crs(f_end_mask) <- crs(world)
memory.limit(size=56000)
f_end_mask <- projectRaster(f_end_mask, crs = proj_str)

fp <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(f_end_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = CWE, color = CWE)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA,
          color = "grey75", size = 0.15) +
  #scale_fill_gradient(low = "#E8E8E8",
  #                    high = "#2A5A5B",
  #                    limits = c(0,0.00005))  +
  #scale_color_gradient(low = "#E8E8E8",
  #                    high = "#2A5A5B",
  #                    limits = c(0,0.00005))  +
  scale_fill_viridis_c(option = "turbo",
                       limits = c(0,0.00005),
                       breaks = c(0.00001, 0.000025, 0.00004),
                       labels = c("Low", "Mid", "High"))  +
  scale_color_viridis_c(option = "turbo",
                        limits = c(0,0.00005),
                        breaks = c(0.00001, 0.000025, 0.00004),
                        labels = c("Low", "Mid", "High")) +
  labs(fill = "CWE",
       color = "CWE") +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  theme_void()

ggsave(fp, filename = "figures/nearctic/forestEndemism.png", dpi = 450)

# non forest
#nonforest
nf_end_mask <- mask(nf_end, nearctic_clip)
crs(nf_end_mask) <- crs(world)
memory.limit(size=56000)
nf_end_mask <- projectRaster(nf_end_mask, crs = proj_str)

nfp <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(nf_end_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = CWE, color = CWE)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA,
          color = "grey75", size = 0.15) +
  #scale_fill_gradient(low = "#E8E8E8",
  #                    high = "#2A5A5B",
  #                    limits = c(0,0.00006))  +
  #scale_color_gradient(low = "#E8E8E8",
  #                     high = "#2A5A5B",
  #                     limits = c(0,0.00006))  +
  scale_fill_viridis_c(option = "turbo",
                       limits = c(0,0.00005),
                       breaks = c(0.00001, 0.000025, 0.00004),
                       labels = c("Low", "Mid", "High"))  +
  scale_color_viridis_c(option = "turbo",
                        limits = c(0,0.00005),
                        breaks = c(0.00001, 0.000025, 0.00004),
                        labels = c("Low", "Mid", "High")) +
  labs(fill = "CWE",
       color = "CWE") +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  theme_void()

ggsave(nfp, filename = "figures/nearctic/nonforestEndemism.png", dpi = 450)

# then nearctic only
#lentic
len_end <- dplyr::select(len, x,y, CWE) %>%
  rasterFromXYZ()
len_end_mask <- mask(len_end, nearctic_clip)
crs(len_end_mask) <- crs(world)
len_end_mask <- projectRaster(len_end_mask, crs = proj_str)

len_p <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(len_end_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = layer, color = layer)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA,
          color = "grey75", size = 0.15) +
  #scale_fill_gradient(low = "#E8E8E8",
  #                    high = "#3B4994",
  #                    limits = c(0,0.00005))  +
  #scale_color_gradient(low = "#E8E8E8",
  #                     high = "#3B4994",
  #                     limits = c(0,0.00005))  +
  scale_fill_viridis_c(option = "turbo",
                       limits = c(0,0.00005),
                       breaks = c(0.00001, 0.000025, 0.00004),
                       labels = c("Low", "Mid", "High"))  +
  scale_color_viridis_c(option = "turbo",
                        limits = c(0,0.00005),
                        breaks = c(0.00001, 0.000025, 0.00004),
                        labels = c("Low", "Mid", "High")) +
  labs(fill = "CWE",
       color = "CWE") +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  theme_void()

ggsave(len_p, filename = "figures/nearctic/lenticEndemism.png", dpi = 450)

#lotic
lot_end <- dplyr::select(lot, x,y, CWE) %>%
  rasterFromXYZ()
lot_end_mask <- mask(lot_end, nearctic_clip)
crs(lot_end_mask) <- crs(world)
lot_end_mask <- projectRaster(lot_end_mask, crs = proj_str)

lot_p <- ggplot() +
  geom_tile(na.omit(raster::as.data.frame(lot_end_mask, xy = T)),
            mapping = aes(x = x , y = y, fill = CWE, color = CWE)) +
  geom_sf(americas_proj, mapping = aes(), fill = NA,
          color = "grey75", size = 0.15) +
  #scale_fill_gradient(low = "#E8E8E8",
  #                    high = "#3B4994",
  #                    limits = c(0,0.00006))  +
  #scale_color_gradient(low = "#E8E8E8",
  #                     high = "#3B4994",
  #                     limits = c(0,0.00006)) +
  scale_fill_viridis_c(option = "turbo",
                       limits = c(0,0.00005),
                       breaks = c(0.00001, 0.000025, 0.00004),
                       labels = c("Low", "Mid", "High"))  +
  scale_color_viridis_c(option = "turbo",
                        limits = c(0,0.00005),
                        breaks = c(0.00001, 0.000025, 0.00004),
                        labels = c("Low", "Mid", "High")) +
  labs(fill = "CWE",
       color = "CWE") +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  theme_void()

ggsave(lot_p, filename = "figures/nearctic/Endemism_loticSpp_nearcticExtent.png", dpi = 450)

## cow plot everything together
cp <- cowplot::plot_grid(fp, nfp,
                         len_p, lot_p, ncol = 2, nrow = 2,
                         labels = c("A", "B", "C", "D"),
                         label_size = 20)

ggsave(cp, filename = "figures/nearctic/FigureX_endemismTraits.png", dpi = 450,
       width = 14, height = 10)

ga <- ggpubr::ggarrange(fp, nfp,
                        len_p, lot_p,
                        ncol = 2, nrow = 2,
                        labels = c("Forest", "Nonforest", "Lentic", "Lotic"),
                        common.legend = TRUE, legend = "right"
                        )

ggsave(ga, filename = "figures/nearctic/FigureX_endemismTraits_GA.png", dpi = 450,
       width = 14, height = 10)
