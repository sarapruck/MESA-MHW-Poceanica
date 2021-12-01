# load libraries 
library(raster)
library(rgdal)
library(tidyverse)
library(ncdf4)
library(sf)
library(RColorBrewer)
library(ggplot2)
library(ggplot2)
library(rasterVis)

#putting variables at the beginning so they are not baked in:

#category 1: Max temperature
temp1 <- 24
temp2 <- 27
temp3 <- 29
temp4 <- 32
#temperature matrix in the following way. THis means that everything < temp1 will be assigned 1, between temp1 & temp2 = 2, etc.
matrix_temp <- matrix(c(-Inf, temp1, 1, temp1, temp2, 2, temp2, temp3, 3, temp3, temp4, 4, temp4, Inf, 5), ncol=3, byrow=TRUE)
#temp doesnt need a function, as max() is already sufficient for this analysis.

#category 2: Max HW length
days1 <- 4
days2 <- 10
days3 <- 20
days4 <- 30
#days risk matrix:
matrix_days <- matrix(c(-Inf, days1, 1, days1, days2, 2, days2, days3, 3, days3, days4, 4, days4, Inf, 5), ncol=3, byrow=TRUE)

#for this, we do need a function. we want to find out how many days in a row a temperature of above 25 degrees was reached; min set at 3 days
fn_maxHWlength <-function(x) {
  seq <- rle(x > temp1)
  df <- do.call(cbind, seq) %>% as.data.frame() %>% filter(values == 1)
  n = ifelse(max(df$lengths) > 0, max(df$lengths), 0)
  return(n)
}

#category 3: Time between heatwaves
time1 <- 60
time2 <- 40
time3 <- 20
time4 <- 5
matrix_time <- matrix(c(time1, Inf, 1, time2, time1, 2, time3, time2, 3, time4, time3, 4, 0, time4, 5), ncol=3, byrow=TRUE)

#we need a function that will look at the time between heatwaves, so where the temperature is below temp1
fn_timebetween <- function(x){
  seq <- rle(x < temp1)
  df <- do.call(cbind, seq) %>% as.data.frame() %>% filter(values == 1)
  n = ifelse(max(df$lengths) > 0, max(df$lengths), 0)
  return(n)
}

#no. of heatwaves in a year
annualhw1 <- 0
annualhw2 <- 1
annualhw3 <- 2
annualhw4 <- 3
matrix_annualhw <- matrix(c(0, annualhw1, 1, annualhw1, annualhw2, 2, annualhw2, annualhw3, 3, annualhw3, annualhw4, 4, annualhw4, Inf, 5), ncol=3, byrow=TRUE)
# mapping the count, not just the length! this will show how many times, the temperature was exceeded for that amount of days
fn_numHWyear <-function(x) {
  seq <- rle(x>temp1)
  df <- do.call(cbind, seq) %>% as.data.frame() %>% filter(values == 1)
  n = ifelse(length(df$values) > 0, length(df$values), 0)
  return(n)
}

#no of years with heatwaves
years1 <- 0
years2 <- 2
years3 <- 5
years4 <- 7
matrix_years <- matrix(c(0, years1, 1, years1, years2, 2, years2, years3, 3, years3, years4, 4, years4, Inf, 5), ncol=3, byrow=TRUE)


# alternative paths
path <- "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/raw_data/"
resultpath <- "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/outputs/"

# import the Regional Seas layer and extract the Mediterranean
RegSeas <- st_read(paste(path, "Regional_Seas_simplified1km.shp", sep = ""))
MedSea <- subset(RegSeas, Name == "Mediterranean")

# set extent of med
Med_extent <- raster::extent(MedSea)

# get folder list
folder_list <- list.dirs(path = paste(path, "dataset-2046-2055/RCP4_5", sep = ""))

#empty raster that has all the highest SST
overall_stack <- stack()

for (i in 2:length(folder_list)){
  
  # get list of files
  file_list <- list.files(folder_list[i], pattern = ".nc", full.names = TRUE)
  
  # create raster stack
  r <- raster::stack(file_list)
  
  # crop raster to med extent
  r_crop <- raster::crop(r,Med_extent)
  
  #year for the filename later
  year <- str_sub(folder_list[i],-4,-1)
  
  # calculate the maximum temperature for each pixel, and add it to a stack
  highestSST <- max(r_crop)

  #now calculate the max length of that year, and add to stack
  maxlength <- raster::calc(r_crop,fn_maxHWlength)
  maxlength[maxlength >= cellStats(maxlength, stat = max)] <- NA

  #same thing with time between heatwaves
  timebetween <- raster::calc(r_crop,fn_timebetween)

  #same thing with the number of annual heatwaves
  annualHW <- raster::calc(r_crop,fn_numHWyear)

  #reclassify
  riskmaxtemp <- reclassify(highestSST, matrix_temp)
  riskmaxlength <- reclassify(maxlength,matrix_days)
  riskmintime <- reclassify(timebetween,matrix_time)
  riskmaxannual <- reclassify(annualHW, matrix_annualhw)
  
  # calculate average susceptibility score for year
  mean <- mean(stack(riskmaxtemp, riskmaxlength, riskmintime, riskmaxannual))
  
  # create overall stack
  overall_stack <- stack(overall_stack, mean)
  
}

susceptibility <- calc(overall_stack, fun = max)

gap.fill = function(r, max.it = 1e4, tol = 1e-2, verbose=FALSE) {
  gaps = which(is.na(r)[])
  r.filled = r
  w = matrix(c(0,0.25,0,0.25,0,0.25,0,0.25,0), nc=3, nr=3)
  i = 0
  while(i < max.it) {
    i = i + 1
    new.vals = focal(r.filled, w=w, na.rm=TRUE)[gaps]
    max.residual = suppressWarnings(max(abs(r.filled[gaps] - new.vals), na.rm = TRUE))
    if (verbose) print(paste('Iteration', i, ': residual = ', max.residual))
    r.filled[gaps] = new.vals
    if (is.finite(max.residual) & max.residual <= tol) break
  }
  return(r.filled)
}

filled <- gap.fill(susceptibility)

filled_round <- round(filled, 0)
susceptibility_round <- round(susceptibility, 0)


writeRaster(susceptibility_round, "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/outputs/Susceptibility/susceptibility_round_v1711.tif", overwrite = TRUE)
writeRaster(filled_round, "O:/f01_projects_active/PanEuropean/p07352E_ODYSSEA/outputs/Susceptibility/susceptibility_filled_round_v1711.tif", overwrite = TRUE)

