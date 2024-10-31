# define the settings
CSCrun <- T
vPREBAS = "master"   #### choose PREBAS version to run the model  "master" "v0.2.x" "newVersion"
setwd("/scratch/project_2000994/srinetri/regional/regRuns_currV/")
# load("/scratch/project_2000994/srinetri/calibrated parameters/pPRELES_calPcurrV.rdata")
# pPRELES <- pPRELES_calPcurrV
# load("/scratch/project_2000994/srinetri/calibrated parameters/pCROBAS_calPcurrV.rdata")
# pCROBAS <- pCROBAS_calPcurrV
r_no <- regions <- 18 ### forest center ID
# nSetRuns <- 10 #number of set runs
# harvScen <- "Base" ### c("Low","MaxSust","NoHarv","Base","Mitigation", "MitigationNoAdH","adapt","protect")
# harvInten <- "Base"
# regSets <- "maakunta" ### "forCent", "maakunta"
# minDharvX <- 15
# rcpfile <- rcps <- "CurrClim"  #"CurrClim" "CanESM2.rcp26.rdata" "CanESM2.rcp45.rdata" "CanESM2.rcp85.rdata")
# # climIDName <- "clim1800" ### climIDs for ISIMIP resolutions and clim files
# # outType <- "testRun"
# deadWoodCalc <- F
# HSIruns <- F

##load settings
devtools::source_url("https://raw.githubusercontent.com/rsrinet/IBCcarbon_runs/master/calRuns_currV/Rsrc/settings.r")
##load functions
source_url("https://raw.githubusercontent.com/rsrinet/IBCcarbon_runs/master/general/functions.r")

load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/data.all_maakunta_",r_no,".rdata"))
load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
data.all$segID <- data.all$maakuntaID

####procData
cord.ne = SpatialPoints(cbind(data.IDs$x,data.IDs$y), proj4string=CRS("+init=EPSG:3067"))
location<-as.data.frame(spTransform(cord.ne, CRS("+init=epsg:4326")))
rm(list=c("cord.ne")); gc()
lat <- location$coords.x2
rm(list=c("location")); gc()
data.IDs$latitude <- lat
rm(list=c("lat")); gc()

latitude <- aggregate(.~maakuntaID, data=data.IDs, mean)$latitude
data.all <- data.all[, latitude := latitude]
latitude.dt <- data.all[, c("segID", "latitude")]
save(latitude.dt, file=paste0("add_inputs/latitude_",r_no,".rdata"))
