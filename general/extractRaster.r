library(rgdal)
library(raster)
library(ggpubr)
library(ggplot2)
library(devtools)
library(data.table)
library(ggridges)
library(parallel)
library(sp)
library(sf)

###FUNCTION FOR EXTRACTING DATA FROM AREAS IN MAAKUNTA REGIONS BASED ON PERIOD, HARVEST, & CLIMATE SCENARIOS
    ##shapex is the raster of the area you would like to extract data from
    ##pathdata is the directory of files
    ##VARx is the variable of data to extract (Soil Carbon, Wood Harvested, etc.)
    ##periodx is the period of time for the scenario (2017-2025, 2026-2033, 2034-2050)
    ##harvscenx is the harvest scenario (Base, Low, MaxSust, NoHarv)
    ##climscenx is the climate scenario (CurrClim)

##Extraction Function
extractraster <- function(shapex, pathdata, VARx, periodx, harvscenx, climscenx){

    ##Set the path
    setwd(pathdata)
    
    ##Upload the country-wide raster
    con <- shapefile("/scratch/project_2000994/PREBASruns/ZonUnc/Testmaps/SuomenMaakuntajako_2021_10k.shp")

    ##Crop the country raster to just the desired area
    crop <- crop(con, shapex)

    ##Extract vector of region IDs present in the desired area
    regs <- c(crop$maakID)
    
    ##The number of regions within the area
    nreg <- length(regs)
    
    ##Create an empty list of the rasters to fill and merge
    xrast <- list()
    
    ##Loop to extract data from maakunta regions
    ##If more than 1 region, the rasters for the same time, variables, and scenarios are merged into one raster
    if(nreg > 1){
        for(j in 1:nreg) xrast[[j]] <- raster(paste0("rasters/maakunta",regs[j],"/",VARx,"_",periodx,"_",harvscenx,"_",climscenx,".tif"))
            rast1 <- do.call(merge, xrast)
            
    } else{
        rast1 <- raster(paste0("rasters/maakunta",regs[1],"/",VARx,"_",periodx,"_",harvscenx,"_",climscenx,".tif"))
        }
    
    ##The region raster is cropped to the desired area
    crop1 <- mask(crop(rast1, shapex),shapex)
    
    ##Save the new raster for that particular area
    writeRaster(crop1, paste0("outRast/", VARx,"_",periodx,"_",harvscenx,"_",climscenx,".tif"),overwrite=T)
}

###INPUTS FOR FUNCTION
##Path for file
pathdata <- "/scratch/project_2000994/PREBASruns/ZonUnc/Testmaps/"

##List of variables
VAR <- c("WenergyWood","Wharvested","Wtot","soilC")

##Periods of time
periods <- c("2017-2025", "2034-2050")

##Desired Area
shapex <- shapefile("/scratch/project_2000994/PREBASruns/ZonUnc/Testmaps/EVO/hyperlentoehdotus.shp")

###RUN FUNCTION
##Loop to input multiple variables and periods into the function
for(i in 1:length(VAR)){
    for(j in 1:length(periods)){
        extractraster(shapex, pathdata, VAR[i], periods[j], "Base", "CurrClim")
    }
}

