# define the settings
CSCrun <- T
vPREBAS <- "newVersion"
setwd("/scratch/project_2000994/srinetri/regional/regional_newV")
load("input/pPRELES_newVcalP_CN.rdata")
load("input/pCROBAS_newVcalP_CN.rdata")
r_no <- regions <- 19 ### forest center ID
nSetRuns <- 10 #number of set runs
harvScen <- "Base" ### c("Low","MaxSust","NoHarv","Base","Mitigation", "MitigationNoAdH","adapt","protect")
harvInten <- "Base"
regSets <- "maakunta" ### "forCent", "maakunta"
minDharvX <- 15
rcpfile <- rcps <- "CurrClim"  #"CurrClim" "CanESM2.rcp26.rdata" "CanESM2.rcp45.rdata" "CanESM2.rcp85.rdata")
# climIDName <- "clim1800" ### climIDs for ISIMIP resolutions and clim files
# outType <- "testRun"
deadWoodCalc <- F
HSIruns <- F

##load general settings
source("Rsrc/settings_newV.r")
##load functions
source("Rsrc/functions_newV.r")

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
save(latitude.dt, file=paste0("/scratch/project_2000994/srinetri/regional/regional_newV/newV_inputs/latitude_",r_no,".rdata"))

load("/scratch/project_2000994/RCP/CurrClim.rdata")
smoothP0=0
maxYears <- 43
startingYear <- 1980
nSites <- nrow(data.all)
dat <- dat[id %in% data.all[, unique(id)]]
nClimID <- length(unique(dat$id))

dat[, pvm:= as.Date('1980-01-01') - 1 + rday ]
dat[, DOY:= as.numeric(format(pvm, "%j"))]
dat[, Year:= as.numeric(format(pvm, "%Y"))]
dat = dat[Year >= startingYear]
dat[DOY==366, DOY:=365]
  
dat[, Mon:= as.numeric(format(pvm, "%m"))]
Tmean = dat[, mean(TAir), by = Year]
Tsum = dat[, sum(ifelse(TAir>5, TAir-5, 0)), by=.(id, Year)][, mean(V1), by=Year]
PARmean = dat[, mean(PAR), by = Year]
VPDmean = dat[, mean(VPD), by = Year]
CO2mean = dat[, mean(CO2), by = Year]
Precip = dat[, sum(Precip), by = .(id, Year)][, mean(V1), by=Year]
Tampl = dat[, .(mean(TAir)), by = .(id, Year, Mon)][, (max(V1)-min(V1))/2, by=Year]
  
id = dat[,unique(id)]
PARtran = t( dcast(dat[, list(id, rday, PAR)], rday ~ id,
                     value.var="PAR")[, -1])
TAirtran = t( dcast(dat[, list(id, rday, TAir)], rday ~ id,
                      value.var="TAir")[, -1])
VPDtran = t( dcast(dat[, list(id, rday, VPD)], rday ~ id,
                     value.var="VPD")[, -1])
Preciptran = t( dcast(dat[, list(id, rday, Precip)], rday ~ id,
                        value.var="Precip")[, -1])
CO2tran = t( dcast(dat[, list(id, rday, CO2)], rday ~ id,
                     value.var="CO2")[, -1])
  
### NLimit
##process weather inputs for YASSO
weatherYasso <- array(0,dim=c(length(id),maxYears,3))
if(length(id)>1){
  weatherYasso[,,1] <- t(apply(TAirtran[,1:(maxYears*365)],1,aTmean,maxYears))
  weatherYasso[,,3] <- t(apply(TAirtran[,1:(maxYears*365)],1,aTampl,maxYears))
  weatherYasso[,,2] <- t(apply(Preciptran[,1:(maxYears*365)],1,aPrecip,maxYears))
  }else{
    weatherYasso[1,,1] <- aTmean(TAirtran[1,1:(maxYears*365)],maxYears) 
    weatherYasso[1,,3] <- aTampl(TAirtran[1,1:(maxYears*365)],maxYears) 
    weatherYasso[1,,2] <- aPrecip(Preciptran[1,1:(maxYears*365)],maxYears)
}

multiP0 <- array(NA,dim=c(length(id),maxYears,2))
for(climID in 1:length(id)){
  P0 <- PRELES(DOY=rep(1:365,maxYears),
              PAR=dat$PAR[id==id[climID]][1:(365*maxYears)],
              TAir=dat$TAir[id==id[climID]][1:(365*maxYears)],
              VPD=dat$VPD[id==id[climID]][1:(365*maxYears)],
              Precip=dat$Precip[id==id[climID]][1:(365*maxYears)],
              CO2=dat$CO2[id==id[climID]][1:(365*maxYears)],
              fAPAR=rep(1,(365*maxYears)),LOGFLAG=0,p=pPRELES_newVcalP_CN)$GPP
  P0 <- matrix(P0,365,maxYears)
  multiP0[climID,(1:maxYears),1] <- colSums(P0)
}
if(smoothP0==1 & maxYears > 1){
  multiP0[,1,2] <- multiP0[,1,1]
for(i in 2:maxYears) multiP0[,i,2] <- multiP0[,(i-1),2] + (multiP0[,i,1]-multiP0[,(i-1),2])/min(i,smoothYear)
# multiP0[,,2] <- matrix(rowMeans(multiP0[,,1]),nClimID,maxYears,byrow = F)
}else{multiP0[,,2] <- multiP0[,,1]
}
P0currClim <- apply(multiP0,1,mean)
P0.dt <- data.table(id=id, P0currClim=P0currClim)
save(P0.dt, file=paste0("/scratch/project_2000994/srinetri/regional/regional_newV/newV_inputs/P0currClim_",r_no,".rdata"))

fT <- fTfun(weatherYasso[,,1],weatherYasso[,,2],weatherYasso[,,3])
fT0AvgCurrClim <- apply(fT,1,mean)
fT0.dt <- data.table(id=id, fT0=fT0AvgCurrClim)
save(fT0.dt, file=paste0("/scratch/project_2000994/srinetri/regional/regional_newV/newV_inputs/fT0_",r_no,".rdata"))

