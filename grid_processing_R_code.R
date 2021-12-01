######################
#### Create Grids ####
######################

# this script can be used to create differing sizes of grids across an area
# to update the code to a different country then just find and replace the ISO3 code.

library(sf)
library(rgeos)
library(tidyverse)

# read in object the size and shape of the output you want
SEAS <- read_sf("O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/raw_data/Regional_Seas_simplified1km.shp")
MED <- SEAS %>% filter(Name == "Mediterranean")

# make grids
lv1_grid <- st_make_grid(MED, cellsize = c(0.5, 0.5), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('id' = 1:length(.)))
lv1_grid$geom <- "polygon"
lv1_grid$level_qdgc <- 1
lv1_grid$cellsizerees <- 0.5
centroids <- st_centroid(lv1_grid)
pts <- do.call(rbind, st_geometry(centroids))
lv1_grid$lon_center <- pts[,1]
lv1_grid$lat_center <- pts[,2]
qdgc_1 <- ifelse(lv1_grid$lon_center < 0,"E","W")
qdgc_2 <- str_pad(str_replace_all(round(abs(lv1_grid$lon_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
qdgc_3 <- ifelse(lv1_grid$lat_center < 0,"S","N")
qdgc_4 <- str_pad(str_replace_all(round(abs(lv1_grid$lat_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
lv1_grid$qdgc <- paste0(qdgc_1, qdgc_2, qdgc_3, qdgc_4)
lv1_grid$area_m2 <- as.numeric(st_area(lv1_grid))
# qdgc_5 <- rep(c(rep(c("A", "B"), each = 1, times = 93), rep(c("C", "D"), each = 1, times = 93)), each = 1, times = nrow(lv1_grid)/ (93*4))


lv2_grid <- st_make_grid(MED, cellsize = c(0.25, 0.25), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('id' = 1:length(.)))
lv2_grid$geom <- "polygon"
lv2_grid$level_qdgc <- 2
lv2_grid$cellsizerees <- 0.25
centroids <- st_centroid(lv2_grid)
pts <- do.call(rbind, st_geometry(centroids))
lv2_grid$lon_center <- pts[,1]
lv2_grid$lat_center <- pts[,2]
qdgc_1 <- ifelse(lv2_grid$lon_center < 0,"E","W")
qdgc_2 <- str_pad(str_replace_all(round(abs(lv2_grid$lon_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
qdgc_3 <- ifelse(lv2_grid$lat_center < 0,"S","N")
qdgc_4 <- str_pad(str_replace_all(round(abs(lv2_grid$lat_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
lv2_grid$qdgc <- paste0(qdgc_1, qdgc_2, qdgc_3, qdgc_4)
lv2_grid$area_m2 <- as.numeric(st_area(lv2_grid))


lv3_grid <- st_make_grid(MED, cellsize = c(0.125, 0.125), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('id' = 1:length(.)))
lv3_grid$geom <- "polygon"
lv3_grid$level_qdgc <- 3
lv3_grid$cellsizerees <- 0.125
centroids <- st_centroid(lv3_grid)
pts <- do.call(rbind, st_geometry(centroids))
lv3_grid$lon_center <- pts[,1]
lv3_grid$lat_center <- pts[,2]
qdgc_1 <- ifelse(lv3_grid$lon_center < 0,"E","W")
qdgc_2 <- str_pad(str_replace_all(round(abs(lv3_grid$lon_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
qdgc_3 <- ifelse(lv3_grid$lat_center < 0,"S","N")
qdgc_4 <- str_pad(str_replace_all(round(abs(lv3_grid$lat_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
lv3_grid$qdgc <- paste0(qdgc_1, qdgc_2, qdgc_3, qdgc_4)
lv3_grid$area_m2 <- as.numeric(st_area(lv3_grid))


lv4_grid <- st_make_grid(MED, cellsize = c(0.0625, 0.0625), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('id' = 1:length(.)))
lv4_grid$geom <- "polygon"
lv4_grid$level_qdgc <- 4
lv4_grid$cellsizerees <- 0.0625
centroids <- st_centroid(lv4_grid)
pts <- do.call(rbind, st_geometry(centroids))
lv4_grid$lon_center <- pts[,1]
lv4_grid$lat_center <- pts[,2]
qdgc_1 <- ifelse(lv4_grid$lon_center < 0,"E","W")
qdgc_2 <- str_pad(str_replace_all(round(abs(lv4_grid$lon_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
qdgc_3 <- ifelse(lv4_grid$lat_center < 0,"S","N")
qdgc_4 <- str_pad(str_replace_all(round(abs(lv4_grid$lat_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
lv4_grid$qdgc <- paste0(qdgc_1, qdgc_2, qdgc_3, qdgc_4)
lv4_grid$area_m2 <- as.numeric(st_area(lv4_grid))


lv5_grid <- st_make_grid(MED, cellsize = c(0.03125, 0.03125), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('id' = 1:length(.)))
lv5_grid$geom <- "polygon"
lv5_grid$level_qdgc <- 5
lv5_grid$cellsizerees <- 0.03125
centroids <- st_centroid(lv5_grid)
pts <- do.call(rbind, st_geometry(centroids))
lv5_grid$lon_center <- pts[,1]
lv5_grid$lat_center <- pts[,2]
qdgc_1 <- ifelse(lv5_grid$lon_center < 0,"E","W")
qdgc_2 <- str_pad(str_replace_all(round(abs(lv5_grid$lon_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
qdgc_3 <- ifelse(lv5_grid$lat_center < 0,"S","N")
qdgc_4 <- str_pad(str_replace_all(round(abs(lv5_grid$lat_center), 3), "[^[:alnum:]]", ""), width=5, side="left", pad="0")
lv5_grid$qdgc <- paste0(qdgc_1, qdgc_2, qdgc_3, qdgc_4)
lv5_grid$area_m2 <- as.numeric(st_area(lv5_grid))


# plot to check
plot(st_geometry(MED))
plot(st_geometry(lv1_grid), add = TRUE)

# set output path
out_path <- "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/scratch/MESA_fisheries/grids/"

# write grids to file
st_write(lv1_grid, paste(out_path, "MED_lv1_grid.shp", sep = ""))
st_write(lv2_grid, paste(out_path, "MED_lv2_grid.shp", sep = ""))
st_write(lv3_grid, paste(out_path, "MED_lv3_grid.shp", sep = ""))
st_write(lv4_grid, paste(out_path, "MED_lv4_grid.shp", sep = ""))
st_write(lv5_grid, paste(out_path, "MED_lv5_grid.shp", sep = ""))

# write to geopackage
st_write(lv1_grid, paste(out_path, "MED_Grids.gpkg", sep = ""), "lv1_grid")
st_write(lv2_grid, paste(out_path, "MED_Grids.gpkg", sep = ""), "lv2_grid", append = TRUE)
st_write(lv3_grid, paste(out_path, "MED_Grids.gpkg", sep = ""), "lv3_grid", append = TRUE)
st_write(lv4_grid, paste(out_path, "MED_Grids.gpkg", sep = ""), "lv4_grid", append = TRUE)
st_write(lv5_grid, paste(out_path, "MED_Grids.gpkg", sep = ""), "lv5_grid", append = TRUE)
