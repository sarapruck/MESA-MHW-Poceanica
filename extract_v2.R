#################################
#### p07352E ODYSSEA Overlay ####
#################################

library(rgdal)
library(raster)
library(sf)
library(tidyverse)
library(terra)

# set path to where data stored
path <- "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/"
workpath <- past(path, "work_in_progress")
resultpath <- paste(path, "outputs/Susceptibility")

# list of shapefiles
shapefile_list <- list.files(workpath, pattern = "\\.shp$", full.names = TRUE)
shapefile_list2 <- list.files(workpath, pattern = "\\.shp$")

stack <- raster(paste(resultpath, "/susceptibility_filled_round_v1711.tif"))

df_list <- list()

for (i in 1:length(shapefile_list)){
  
  sf <- read_sf(shapefile_list[i])
  sf_project <- st_transform(sf, crs(stack))
  # sf_project_cast <- st_cast(sf_project,"POLYGON")
  shapefile <- as_Spatial(sf_project)
  extract_max <- raster::extract(stack, shapefile, fun = max,  df = TRUE, na.rm=TRUE)
  df <- round(extract_max,0)
  df$shapefile <- shapefile_list2[i]
  extract_max$shapefile <- shapefile_list2[i]
  df[df== -Inf] <- 3
  sf_project_merge <- merge(sf_project, df, by = 0)
  
  max_values <- unique(df$susceptibility_filled_round_v1711)

  
  for (j in max_values){
  
    sf_max <- sf_project_merge %>% filter(susceptibility_filled_round_v1711 == j)
    write_sf(sf_max, paste0(path, "work_in_progress/v3/sus_",j,"_imp_",str_sub(shapefile_list2[i],1,4), "_", str_sub(shapefile_list2[i],-5,-1)))
    
    }
  
  
  df_list[i] <- extract_max

}



