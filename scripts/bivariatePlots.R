# Run lines 3 & 4 if an older version of sf is needed, newer versions may cause errors

#purl <- "http://cran.us.r-project.org/src/contrib/Archive/sf/sf_0.9-6.tar.gz"
#install.packages(purl, repos=NULL, type="source")

library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
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

na_clip <- na %>%
  ms_erase(., lake) %>%
  as_Spatial()

## read in lentic only data
len <- read.csv("data/regular_richnessCSVs/propRichness_lentic.csv")
lot <- read.csv("data/regular_richnessCSVs/propRichness_lotic.csv")

#generate bivariate data
bdf <- left_join(len, lot, by = c("z")) %>%
  dplyr::rename(for_richness = richness.x,
                nonfor_richness = richness.y,
                x = x.y, y = y.y) %>%
  dplyr::filter(!is.na(x), !is.na(y))

bdf_sf <- st_as_sf(bdf, coords = c("x","y"), crs = crs(na_clip))
near_clip_sf <- st_as_sf(nearctic_clip)
bdf_sf <- st_transform(bdf_sf, crs = st_crs(near_clip_sf))
bdf_sf_c <- st_join(bdf_sf, st_as_sf(near_clip_sf))
bdf_sf_c <- dplyr::filter(bdf_sf_c,!is.na(rmapshaperid))

lon <- st_coordinates(bdf_sf_c)[,1]
lat <- st_coordinates(bdf_sf_c)[,2]

bdf_df <- as.data.frame(st_drop_geometry(bdf_sf_c)) %>%
  mutate(x = lon, y = lat)

# create classes
data <- bi_class(bdf_df, x = nonfor_richness, y = for_richness, style = "quantile", dim = 3)

data_proj <- sf_project(from = "+proj=longlat +ellps=WGS84 +no_defs",
                        to = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                        pts = data[11:12])

t <- data %>%
  mutate(eq_lon = data_proj[,1], eq_lat = data_proj[,2])

ggplot(t, mapping = aes(x = eq_lon, y = eq_lat, fill = bi_class)) + geom_tile()

#non lentic
lot_rich <- dplyr::select(lot, x,y,propRich) %>%
  rasterFromXYZ()
lot_rich_mask <- mask(lot_rich, nearctic_clip)
crs(lot_rich_mask) <- crs(americas)
lot_rich_mask <- projectRaster(lot_rich_mask, crs = proj_str)

rdf <- raster::as.data.frame(lot_rich_mask, xy = T)%>%
  mutate(z = 1:nrow(.))
rxyz <- rdf %>%
  dplyr::select(x,y,z) %>%
  rasterFromXYZ(.)

spdf <- t
coordinates(spdf) <- ~ eq_lon + eq_lat

e <- extract(rxyz, spdf)

t <- dplyr::mutate(as.data.frame(t), z2 = e)

lj <- left_join(t, rdf, by = c("z2" = "z"))
lj2 <- lj %>%
  dplyr::rename(x = x.y, y = y.y) %>%
  dplyr::select(x, y, bi_class)

p <- ggplot() +
  geom_tile(lj2, mapping = aes(x = x, y = y, fill = bi_class, color = bi_class), show.legend = F) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75") +
  bi_scale_fill(pal = "DkBlue", dim = 3) +
  bi_scale_color(pal = "DkBlue", dim = 3) +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  labs(x = "", y ="") +
  bi_theme()

legend <- bi_legend(pal = "DkBlue",
                    dim = 3,
                    xlab = "Lotic rich",
                    ylab = "Lentic rich",
                    size = 12)

c <- ggdraw() +
  draw_plot(p, 0,0,1,1) +
  draw_plot(legend, x = 0.1, y =0.3,
            width = 0.3, height = 0.3)

ggsave(c, filename = "figures/nearctic/lenticLotic_richness.png", dpi = 450,
       width = 11, height = 7)


## read in endemism lentic/lotic only data
f <- read.csv("data/regular_endemismCSVs/propEndemism_lentic.csv")
nf <- read.csv("data/regular_endemismCSVs/propEndemism_lotic.csv")

#generate bivariate data
bdf <- left_join(f, nf, by = c("z")) %>%
  dplyr::rename(for_endemism = CWE.x,
                nonfor_endemism = CWE.y,
                x = x.y, y = y.y) %>%
  dplyr::filter(!is.na(x), !is.na(y))

bdf_sf <- st_as_sf(bdf, coords = c("x","y"), crs = crs(na_clip))
near_clip_sf <- st_as_sf(nearctic_clip)
bdf_sf <- st_transform(bdf_sf, crs = st_crs(near_clip_sf))
bdf_sf_c <- st_join(bdf_sf, st_as_sf(near_clip_sf))
bdf_sf_c <- dplyr::filter(bdf_sf_c,!is.na(rmapshaperid))

lon <- st_coordinates(bdf_sf_c)[,1]
lat <- st_coordinates(bdf_sf_c)[,2]

bdf_df <- as.data.frame(st_drop_geometry(bdf_sf_c)) %>%
  mutate(x = lon, y = lat)

# create classes
data <- bi_class(bdf_df, x = nonfor_endemism, y = for_endemism, style = "quantile", dim = 3)

data_proj <- sf_project(from = "+proj=longlat +ellps=WGS84 +no_defs",
                        to = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                        pts = data[11:12])

t <- data %>%
  mutate(eq_lon = data_proj[,1], eq_lat = data_proj[,2])

spdf <- t
coordinates(spdf) <- ~ eq_lon + eq_lat

e <- extract(rxyz, spdf)

t <- dplyr::mutate(as.data.frame(t), z2 = e)

lj <- left_join(t, rdf, by = c("z2" = "z"))
lj2 <- lj %>%
  dplyr::rename(x = x.y, y = y.y) %>%
  dplyr::select(x, y, bi_class)

p_LL_end <- ggplot() +
  geom_tile(lj2, mapping = aes(x = x, y = y, fill = bi_class, color = bi_class), show.legend = F) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75") +
  bi_scale_fill(pal = "DkBlue", dim = 3) +
  bi_scale_color(pal = "DkBlue", dim = 3) +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  labs(x = "", y ="") +
  bi_theme()

legend_LL_end <- bi_legend(pal = "DkBlue",
                    dim = 3,
                    xlab = "Lotic CWE",
                    ylab = "Lentic CWE",
                    size = 12)

c_LL_end <- ggdraw() +
  draw_plot(p_LL_end, 0,0,1,1) +
  draw_plot(legend_LL_end, x = 0.1, y =0.3,
            width = 0.3, height = 0.3)

ggsave(c_LL_end, filename = "figures/nearctic/lenticLotic_endemism.png", dpi = 450,
       width = 11, height = 7)

## now do forest/nonforest
## read in lentic only data
len <- read.csv("data/regular_richnessCSVs/propRichness_forest.csv")
lot <- read.csv("data/regular_richnessCSVs/propRichness_nonforest.csv")

#generate bivariate data
bdf <- left_join(len, lot, by = c("z")) %>%
  dplyr::rename(for_richness = richness.x,
                nonfor_richness = richness.y,
                x = x.y, y = y.y) %>%
  dplyr::filter(!is.na(x), !is.na(y))

bdf_sf <- st_as_sf(bdf, coords = c("x","y"), crs = crs(na_clip))
near_clip_sf <- st_as_sf(nearctic_clip)
bdf_sf <- st_transform(bdf_sf, crs = st_crs(near_clip_sf))
bdf_sf_c <- st_join(bdf_sf, st_as_sf(near_clip_sf))
bdf_sf_c <- dplyr::filter(bdf_sf_c,!is.na(rmapshaperid))

lon <- st_coordinates(bdf_sf_c)[,1]
lat <- st_coordinates(bdf_sf_c)[,2]

bdf_df <- as.data.frame(st_drop_geometry(bdf_sf_c)) %>%
  mutate(x = lon, y = lat)

# create classes
data <- bi_class(bdf_df, x = nonfor_richness, y = for_richness, style = "quantile", dim = 3)

data_proj <- sf_project(from = "+proj=longlat +ellps=WGS84 +no_defs",
                        to = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                        pts = data[11:12])

t <- data %>%
  mutate(eq_lon = data_proj[,1], eq_lat = data_proj[,2])

# regularlize
spdf <- t
coordinates(spdf) <- ~ eq_lon + eq_lat

e <- extract(rxyz, spdf)

t <- dplyr::mutate(as.data.frame(t), z2 = e)

lj <- left_join(t, rdf, by = c("z2" = "z"))
lj2 <- lj %>%
  dplyr::rename(x = x.y, y = y.y) %>%
  dplyr::select(x, y, bi_class)

p_for_richness <- ggplot() +
  geom_tile(lj2, mapping = aes(x = x, y = y, fill = bi_class, color = bi_class), show.legend = F) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75") +
  bi_scale_fill(pal = "DkCyan", dim = 3) +
  bi_scale_color(pal = "DkCyan", dim = 3) +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  labs(x = "", y ="") +
  bi_theme()

legend_for_richness <- bi_legend(pal = "DkCyan",
                    dim = 3,
                    xlab = "Nonforest rich",
                    ylab = "Forest rich",
                    size = 12)

c_for_richness <- ggdraw() +
  draw_plot(p_for_richness, 0,0,1,1) +
  draw_plot(legend_for_richness, x = 0.1, y =0.3,
            width = 0.3, height = 0.3)

ggsave(c_for_richness, filename = "figures/nearctic/forest-nonforest_richness.png", dpi = 450,
       width = 11, height = 7)


## read in endemism forest/nonforest data
f <- read.csv("data/regular_endemismCSVs/propEndemism_forest.csv")
nf <- read.csv("data/regular_endemismCSVs/propEndemism_nonforest.csv")

#generate bivariate data
bdf <- left_join(f, nf, by = c("z")) %>%
  dplyr::rename(for_endemism = CWE.x,
                nonfor_endemism = CWE.y,
                x = x.y, y = y.y) %>%
  dplyr::filter(!is.na(x), !is.na(y))

bdf_sf <- st_as_sf(bdf, coords = c("x","y"), crs = crs(na_clip))
near_clip_sf <- st_as_sf(nearctic_clip)
bdf_sf <- st_transform(bdf_sf, crs = st_crs(near_clip_sf))
bdf_sf_c <- st_join(bdf_sf, st_as_sf(near_clip_sf))
bdf_sf_c <- dplyr::filter(bdf_sf_c,!is.na(rmapshaperid))

lon <- st_coordinates(bdf_sf_c)[,1]
lat <- st_coordinates(bdf_sf_c)[,2]

bdf_df <- as.data.frame(st_drop_geometry(bdf_sf_c)) %>%
  mutate(x = lon, y = lat)

# create classes
data <- bi_class(bdf_df, x = nonfor_endemism, y = for_endemism, style = "quantile", dim = 3)

data_proj <- sf_project(from = "+proj=longlat +ellps=WGS84 +no_defs",
                        to = "+proj=laea +lat_0=8.27 +lon_0=-96.32 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                        pts = data[11:12])

t <- data %>%
  mutate(eq_lon = data_proj[,1], eq_lat = data_proj[,2])

spdf <- t
coordinates(spdf) <- ~ eq_lon + eq_lat

e <- extract(rxyz, spdf)

t <- dplyr::mutate(as.data.frame(t), z2 = e)

lj <- left_join(t, rdf, by = c("z2" = "z"))
lj2 <- lj %>%
  dplyr::rename(x = x.y, y = y.y) %>%
  dplyr::select(x, y, bi_class)

p_for_end <- ggplot() +
  geom_tile(lj2, mapping = aes(x = x, y = y, fill = bi_class, color = bi_class), show.legend = F) +
  geom_sf(americas_proj, mapping = aes(), fill = NA, color = "grey75") +
  bi_scale_fill(pal = "DkCyan", dim = 3) +
  bi_scale_color(pal = "DkCyan", dim = 3) +
  coord_sf(xlim = c(-4500000, 4000000), ylim = c(1000000,8000000)) +
  labs(x = "", y ="") +
  bi_theme()

legend_for_end <- bi_legend(pal = "DkCyan",
                           dim = 3,
                           xlab = "Nonforest CWE",
                           ylab = "Forest CWE",
                           size = 12)

c_for_end <- ggdraw() +
  draw_plot(p_for_end, 0,0,1,1) +
  draw_plot(legend_for_end, x = 0.1, y =0.3,
            width = 0.3, height = 0.3)

ggsave(c_for_end, filename = "figures/nearctic/forest-nonforest_endemism.png", dpi = 450,
       width = 11, height = 7)

cp <- cowplot::plot_grid(c_for_richness, c_for_end,
                         c, c_LL_end, ncol = 2, nrow = 2,
                          labels = c("A", "B", "C", "D"),
                         label_size = 20)

ggsave(cp, filename = "figures/nearctic/Figure2.png", dpi = 450,
       width = 14, height = 10)
