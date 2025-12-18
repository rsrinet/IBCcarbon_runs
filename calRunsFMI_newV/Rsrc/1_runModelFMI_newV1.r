devtools::source_url("https://raw.githubusercontent.com/rsrinet/IBCcarbon_runs/master/calRunsFMI_newV/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/rsrinet/IBCcarbon_runs/master/general/functions_FMI.r")

# load data
if(TRUE){
  load(paste0("/scratch/project_2000994/PREBASruns/finRuns/input/maakunta/maakunta_",r_no,"_IDsTab.rdata"))
  data.IDs <- data.IDs[segID!=0]
  data.IDs$segID <- data.IDs$maakuntaID
  setkey(data.IDs,segID)
  setkey(data.all,segID)
  #setkey(data.IDs,maakuntaID)
  
  tabX <- merge(data.IDs,data.all)
  ntabX <- tabX[,.I[which.max(y)],by=segID]$V1
  data.all <- cbind(data.all, tabX[ntabX,c("x","y")])
  
  #set_thin_PROJ6_warnings(TRUE)
  xy <- data.all[,c("segID","x","y")]
  coordinates(xy) <- c("x","y")
  proj4string(xy) <- crsX
  #cord = SpatialPoints(xy, proj4string=CRS("+init=EPSG:3067"))
  location<-as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
  data.all$lat <- location$coords.x2#location$y
}

nSitesRun <- 10000
nSamples <- ceiling(dim(data.all)[1]/nSitesRun)
nSetRuns <- 10
sampleIDs <- split(1:nSamples,             # Applying split() function
                   cut(seq_along(1:nSamples),
                       nSetRuns,
                       labels = FALSE))[[setX]]
set.seed(1)
ops <- split(data.all, sample(1:nSamples, nrow(data.all), replace=T))

totArea <- sum(data.all$area)
# test
toMem <- ls()
###check and run missing sampleIDs 
# library('stringi')
# fileX <- list.files(path= "/scratch/project_2000994/PREBASruns/finRuns/outputDT/forCent12/", pattern = "age")
# sampleIDs <- which(!1:nSamples %in%  as.numeric(stri_extract_last(fileX, regex = "(\\d+)")))
# print(sampleIDs)
# sampleIDs <- c(66,342,395)
jx<-1

# Initialize soil & deadwood
if(harvScen == "Base" & harvInten == "Base"){
  out <- runModel(sampleID = jx, outType = "testRun", 
                  rcps = "CurrClim", climScen = 0, 
                  harvScen = "Base", harvInten = "Base", procDrPeat=T, 
                  thinFactX= thinFactX,
                  compHarvX = compHarvX,ageHarvPriorX = ageHarvPriorX,
                  forceSaveInitSoil=F, sampleX = ops[[jx]])
}

if(!CSCrun){
  out <- runModel(sampleID = jx, outType = "testRun", 
                  rcps = rcps, climScen = climScen, 
                  harvScen = harvScen, harvInten =harvInten, procDrPeat=T, 
                  thinFactX= thinFactX,
                  compHarvX = compHarvX,ageHarvPriorX = ageHarvPriorX,
                  forceSaveInitSoil=F, sampleX = ops[[jx]])
} else {
  mclapply(sampleIDs, function(jx) {
    runModel(sampleID = jx, outType = "testRun", 
             rcps = rcps, climScen = climScen, 
             harvScen = harvScen, harvInten =harvInten, procDrPeat=T, 
             thinFactX= thinFactX,
             compHarvX = compHarvX,ageHarvPriorX = ageHarvPriorX,
             forceSaveInitSoil=F, sampleX = ops[[jx]])}, 
    mc.cores = nCores,mc.silent=FALSE)      
}

break

mclapply(sampleIDs, function(jx) {
  runModel(
    jx,compHarvX = compHarvX,
    harvScen=harvScen,harvInten=harvInten,
    cons10run=cons10run,landClassUnman=landClassUnman,
    procDrPeat=procDrPeat
    # outModReStart = reStartMod, initSoilCreStart = reStartSoil,
    # funPreb = reStartRegionPrebas,reStartYear = reStartYearX
  )}, 
  mc.cores = nCores,mc.silent=FALSE)      

# models outputs to NAs, outputDT, initSoilC and plots
Sys.chmod(list.dirs("NAs"), "0777",use_umask=FALSE)
f <- list.files("NAs", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("outputDT"), "0777",use_umask=FALSE)
f <- list.files("outputDT", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("initSoilC"), "0777",use_umask=FALSE)
f <- list.files("initSoilC", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)

Sys.chmod(list.dirs("plots"), "0777",use_umask=FALSE)
f <- list.files("plots", all.files = TRUE, full.names = TRUE, recursive = TRUE)
Sys.chmod(f, (file.info(f)$mode | "0777"),use_umask=FALSE)
