#Script to resample the trawling sampling each yeear between 1997-2018 to get the same sampling effort
#based on David Eme s 'script
# to do divide the grid by 2 or 3
# change the aggergation for taxonomy
# increase the time serie


### Setting working directory 
setwd("/Users/calbouy/Documents/Projets/Maestro/Romane/Unvertebrates_BoB/")

### Libraries loading
lib_vect <- c("sfheaders","Rcpp","rgdal","sf","vegan","betapart","parallel","maptools")
sapply(lib_vect,library,character.only=TRUE)

Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

### Data loading
load("data/BIOdatasetSEL.RData")
load("data/BioDataset_EVHOE.RData")

data_evhoe <- BioDataset[[5]]

#### Mapping preparation
ICES_rect <- st_read("data/ICES_rectangles/ShapeFiles/ICES.BayBiscay.Grid.CoastAccurate.shp")
ICES.rect.shp <- st_read("data/ICES_rectangles/ShapeFiles/ICES_Statistical_Rectangles_Eco.shp")

coast <- readShapePoly("data/ICES_rectangles/ShapeFiles/Bob_land.shp")
# multiple polygone of land, that need to be aggregated into one polygone.
land.id = rep("1", 8)
# merge polygone by polygone id.
coast <- unionSpatialPolygons(coast, land.id)
coast <- as(coast, "SpatialPolygonsDataFrame")
proj4string(coast) <- CRS(Proj)
coast <- st_as_sf(coast)

#ICES_rect <- st_crop(ICES.rect.shp,ICES_rect)
#ICES_rect <- st_difference(ICES_rect,coast)

################################################################################
###                     PART 1: division of the data                         ###
################################################################################
### Trawling stations per year
Trawls.info <- BioDataset$Traits

### Community matrix per station
mat_sp_stations <- BioDataset[[5]]

#change format xtab to real matrix
attr(mat_sp_stations, "class") <- NULL
attr(mat_sp_stations, "call") <- NULL

Worrmtax <- BioDataset$Taxo
mat_sp_station_number <- BioDataset$MATRIXstspdNombre

#change format xtab to real matrix
attr(mat_sp_stations, "class") <- NULL
attr(mat_sp_stations, "call") <- NULL

################################################################################
###                 PART 2: Modification of the data                         ###
################################################################################

#TRAWLS--------------
### Nb of trawl per year
table(Trawls.info$Annee) # 25 : TECHNICAL PROBLEM IN 2017 this is why number so low -- has to be excluded from the data 

#removal year 2017 with 25 trawls only 
Trawls.info.2 <- Trawls.info[which(Trawls.info$Annee != "2017"),]
Years <- sort(unique(Trawls.info.2$Annee))

################################################################################
###         PART 3: Finding in which ICES rectangle each trawl               ###
###           was made based on the coord intersections                      ###
################################################################################

### Plot the grid
plot(rnorm(1000),ylim=c(40,55),xlim=c(-12, 5),axes=T,xlab="",ylab="",lwd=0.25,type="n", main.cex = 1) #the coordinates correspond to the lat and long of the areas (adjusted manually)
plot(ICES_rect$geometry, col="green",add=T)

### Conversion of the trawl coords as spatial points
Trawls.info.2$Long <- as.numeric(Trawls.info.2$Long)
Trawls.info.2$Lat <- as.numeric(Trawls.info.2$La)

STpt_trawls <- st_multipoint(as.matrix(cbind(Trawls.info.2$Long,Trawls.info.2$Lat))) 
STpt_trawls <- st_sfc(STpt_trawls,crs=Proj) 

################################
#ALREADY DONE: 
ICES_rect_to_keep <- st_read("data/ICES_rectangles/ShapeFiles/Unvert/ICES_rect_BOB_EVHOE_Unvert.shp")
#ICES_rect_to_keep <- st_intersection(ICES_rect,ICES_rect_to_keep)

#plot(ICES_rect_to_keep$geometry,col="blue")
#plot(STpt_trawls,add=T)
#plot(ICES_rect_to_keep[which(ICES_rect_to_keep$ID.1=="168"),],add=TRUE,col="red")

#Removing the ICES recatngles that seem to be too far to be in the area 
ICES_rect_sf_to_keep <- ICES_rect[which((-11.6< ICES_rect$x) & (ICES_rect$x< 0)),] #x longitude
nrow(ICES_rect_sf_to_keep) #328 remaining

ICES_rect_sf_to_keep=ICES_rect_sf_to_keep[which((43<ICES_rect_sf_to_keep$y) & (ICES_rect_sf_to_keep$y<52)),] #y latitude
nrow(ICES_rect_sf_to_keep) #197 remaining

#WARNING: the ICES is divided into multiple polygons so we have to check the intersection between our area and all of them 
list_ICES_no_station <- list() #in case some ices rect do not have stations sampled

for (x in 1:nrow(ICES_rect_sf_to_keep)){
  
  message(paste0("x=",x))
  divided_poly <- st_cast(ICES_rect_sf_to_keep[x,"geometry"], "POLYGON") #dividing the multipol into several polygons
  test_intersec <- st_intersection(divided_poly,STpt_trawls) #test if there is an intersection 
  
  #trsfo the geometry into a normal df
  test_intersec_df<- sf_to_df(test_intersec,fill = TRUE)
  
  if(length(grep("NA",test_intersec_df))==0){ #case where there is an intersection
    for (i in 1:nrow(test_intersec_df)){
      row_number <- grep(test_intersec_df$x[i],Trawls.info.2$Long) #search for the station based on the value of the longitude (could also be done with latitude)
      Trawls.info.2[row_number,"ICES_ID"]=ICES_rect_sf_to_keep$ID[x] #giving the ID to the station 
    }#end for2 
  } else{
    list_ICES_no_station <- c(list_ICES_no_station,ICES_rect_sf_to_keep$ID[x])
    message(paste0("list_ICES_no_station=",list_ICES_no_station))
 } 
}#end for1

table(Trawls.info.2$ICES_ID)
min(table(Trawls.info.2$ICES_ID)) #1: number of station to resample per rectangle 
#saveRDS(Trawls.info.2, file="Results/Trawls.info.Unvert_73_ICES_EVHOE.rds")
saveRDS(Trawls.info.2, file="Results/Trawls.info.Unvert_73_ICES_EVHOE2022.rds")


#saving the sf object with the 73 rectangles 
ICES_rect_sf_evhoe <- ICES_rect_sf_to_keep[which(!ICES_rect_sf_to_keep$ID %in% list_ICES_no_station),]
ICES_rect_evhoe <- ICES_rect[which(ICES_rect$ID %in% ICES_rect_sf_evhoe$ID),] #save as spatial polygon

ICES_rect_to_keep2 <- st_crop(ICES.rect.shp,ICES_rect_evhoe)
ICES_rect_to_keep2.2 <- st_intersection(ICES.rect.shp,ICES_rect_evhoe)

#writeOGR(ICES_rect_evhoe,"data/ICES_rectangles/ShapeFiles/Unvert/","ICES_rect_BOB_EVHOE_Unvert",
 #        driver="ESRI Shapefile",overwrite_layer=T)

################################################################################
###      PART 4:Resampling the station with 1 trawl per rectangle            ###
################################################################################

#checking if all the stations are present in the species table 
table(rownames(mat_sp_stations) %in% Trawls.info.2$Trait) #25, correspond to 2017 trawls that we removed (V0....)
table(Trawls.info.2$Trait %in% rownames(mat_sp_stations)) #25 as well 
Missing_trawls <- Trawls.info.2[which(!Trawls.info.2$Trait %in% rownames(mat_sp_stations)),"Station"]
Trawls.info.2  <- Trawls.info.2[which(Trawls.info.2$Trait %in% rownames(mat_sp_stations)),]#keep only the ones in the mat_sp

Trawls.info.2 <- cbind(Trawls.info.2,Ecoreg=rep("Bay of Biscay and the Iberian Coast",nrow(Trawls.info.2)))


CS <- ICES_rect_to_keep2.2[which(ICES_rect_to_keep2.2$Ecoregion=="Celtic Seas"),]$ID.1

for (i in 1 :length(CS)){
  Trawls.info.2[which(Trawls.info.2$ICES_ID==CS[i]),"Ecoreg"] <- "Celtic Seas"
}

points(x=Trawls.info.2$Longitude,y=Trawls.info.2$Latitude,col=as.factor(Trawls.info.2$Ecoreg))

#creation col 
Trawls.info.2$Poly.Annee <- paste0(Trawls.info.2$ICES_ID,"_",Trawls.info.2$Annee)
NB.Trawl.PolyID.PerYear <-  table(Trawls.info.2$Poly.Annee)
length(which(NB.Trawl.PolyID.PerYear == 1)) ### 237 ICES rectangle over all the year wwith only one trawl sample.

#So we resample 1 trawl in every ICES rectangle over the years we do that 100 times.
Uni.PoliID.Annee <-  unique(Trawls.info.2$Poly.Annee)

List.100.resamp.Com.Dat.ICES.rect_all <-  lapply(seq(1, 100), function(i){
  message(paste0("######## i=",i))
  Species.list.resamp <- do.call(rbind,lapply(seq(1,length(Uni.PoliID.Annee)), function(x){ ### for each year and  for each ICES.rectangle.
    message(paste0("x=",x))
    posi <- sample(which(as.character(Trawls.info.2$Poly.Annee) == as.character(Uni.PoliID.Annee[x])),1)
    station <- Trawls.info.2[posi,"Trait"]
    mat_sp_stations[station,]
  }))
  row.names(Species.list.resamp) <- Uni.PoliID.Annee
  Species.list.resamp
})


CS_dataset <- Trawls.info.2[which(Trawls.info.2$Ecoreg=="Celtic Seas"),]
Uni.PoliID.Annee <-  unique(CS_dataset$Poly.Annee)

List.100.resamp.Com.Dat.ICES.rect_CS <-  lapply(seq(1, 100), function(i){
  message(paste0("######## i=",i))
  Species.list.resamp <- do.call(rbind,lapply(seq(1,length(Uni.PoliID.Annee)), function(x){ ### for each year and  for each ICES.rectangle.
    message(paste0("x=",x))
    posi <- sample(which(as.character(CS_dataset$Poly.Annee) == as.character(Uni.PoliID.Annee[x])),1)
    station <- CS_dataset[posi,"Trait"]
    mat_sp_stations[station,]
  }))
  row.names(Species.list.resamp) <- Uni.PoliID.Annee
  Species.list.resamp
})


Bob_dataset <- Trawls.info.2[which(Trawls.info.2$Ecoreg=="Bay of Biscay and the Iberian Coast"),]
Uni.PoliID.Annee <-  unique(Bob_dataset$Poly.Annee)

List.100.resamp.Com.Dat.ICES.rect_Bob <-  lapply(seq(1, 100), function(i){
  message(paste0("######## i=",i))
  Species.list.resamp <- do.call(rbind,lapply(seq(1,length(Uni.PoliID.Annee)), function(x){ ### for each year and  for each ICES.rectangle.
    message(paste0("x=",x))
    posi <- sample(which(as.character(Bob_dataset$Poly.Annee) == as.character(Uni.PoliID.Annee[x])),1)
    station <- Bob_dataset[posi,"Trait"]
    mat_sp_stations[station,]
  }))
  row.names(Species.list.resamp) <- Uni.PoliID.Annee
  Species.list.resamp
})


#change the aphiaID in the scientific name of the species 
list_missing_aphiaID <- c() #there are some that don't have attributed species

List.100.resamp.Com.Dat.ICES.rect_all.sp_name <- lapply(seq(1:length(List.100.resamp.Com.Dat.ICES.rect_all)), function(x){
  for (i in 1:ncol(List.100.resamp.Com.Dat.ICES.rect_all[[x]])){
    message(paste0("######## i=",i))
    sp_tax=Worrmtax[which(Worrmtax$valid_AphiaID==colnames(List.100.resamp.Com.Dat.ICES.rect_all[[x]])[i]),"SOURCE_Nom_Scientifique"]
    
    if (length(sp_tax)==1){
      colnames(List.100.resamp.Com.Dat.ICES.rect_all[[x]])[i]=sp_tax
    } else{ #AphiaID does not correspond to anything 
      colnames(List.100.resamp.Com.Dat.ICES.rect_all[[x]])[i] <- paste0("NA_",colnames(List.100.resamp.Com.Dat.ICES.rect_all[[x]])[i])
      list_missing_aphiaID <- c(list_missing_aphiaID,colnames(List.100.resamp.Com.Dat.ICES.rect_all[[x]])[i])
    } # end of else
  } # end of i
  List.100.resamp.Com.Dat.ICES.rect_all[[x]]
})


#### Mat abves species name for the celtic sea 
list_missing_aphiaID <- c() #there are some that don't have attributed species

List.100.resamp.Com.Dat.ICES.rect_CS.sp_name <- lapply(seq(1:length(List.100.resamp.Com.Dat.ICES.rect_CS)), function(x){
  for (i in 1:ncol(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])){
    message(paste0("######## i=",i))
    sp_tax=Worrmtax[which(Worrmtax$valid_AphiaID==colnames(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])[i]),"SOURCE_Nom_Scientifique"]
    
    if (length(sp_tax)==1){
      colnames(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])[i]=sp_tax
    } else{ #AphiaID does not correspond to anything 
      colnames(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])[i] <- paste0("NA_",colnames(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])[i])
      list_missing_aphiaID <- c(list_missing_aphiaID,colnames(List.100.resamp.Com.Dat.ICES.rect_CS[[x]])[i])
    } # end of else
  } # end of i
  List.100.resamp.Com.Dat.ICES.rect_CS[[x]]
})

#### Mat abves species name for the Bob 
list_missing_aphiaID <- c() #there are some that don't have attributed species

List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name <- lapply(seq(1:length(List.100.resamp.Com.Dat.ICES.rect_Bob)), function(x){
  for (i in 1:ncol(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])){
    message(paste0("######## i=",i))
    sp_tax=Worrmtax[which(Worrmtax$valid_AphiaID==colnames(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])[i]),"SOURCE_Nom_Scientifique"]
    
    if (length(sp_tax)==1){
      colnames(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])[i]=sp_tax
    } else{ #AphiaID does not correspond to anything 
      colnames(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])[i] <- paste0("NA_",colnames(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])[i])
      list_missing_aphiaID <- c(list_missing_aphiaID,colnames(List.100.resamp.Com.Dat.ICES.rect_Bob[[x]])[i])
    } # end of else
  } # end of i
  List.100.resamp.Com.Dat.ICES.rect_Bob[[x]]
})

#search for the NA aphiaID
test_NA <- List.100.resamp.Com.Dat.ICES.rect_all.sp_name[[1]]
list_NA <- grep("NA", colnames(test_NA))
colnames(test_NA)[c(list_NA)]  #""NA_181504" "NA_515738"
saveRDS(List.100.resamp.Com.Dat.ICES.rect_all.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_all_Unvert.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_Bob_Unvert.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_CS_Unvert.rds")


saveRDS(List.100.resamp.Com.Dat.ICES.rect_all.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_all_Unvert_2022.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_Bob_Unvert_2022.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name, file = "Results/List.100.resamp.Com.Dat.ICES.rect_CS_Unvert_2022.rds")


################################################################################

List.100.resamp.Com.Dat.ICES.rect_all.sp_name <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_all_Unvert_2022.rds")
List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_Bob_Unvert_2022.rds")
List.100.resamp.Com.Dat.ICES.rect_CS.sp_name <-  readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_CS_Unvert_2022.rds")

#conversion into pres/abs
List.100.resamp.Com.Dat.ICES.rect_all.sp_name.PA <- 
  lapply(List.100.resamp.Com.Dat.ICES.rect_all.sp_name, function(x){x[x>0]=1; x})

List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name.PA <- 
  lapply(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name, function(x){x[x>0]=1; x})

List.100.resamp.Com.Dat.ICES.rect_CS.sp_name.PA <- 
  lapply(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name, function(x){x[x>0]=1; x})

saveRDS(List.100.resamp.Com.Dat.ICES.rect_all.sp_name.PA, file = "Results/List.100.resamp.Com.Dat.ICES.rect_all.sp_name.PA_2022.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name.PA, file = "Results/List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name.PA_2022.rds")
saveRDS(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name.PA, file = "Results/List.100.resamp.Com.Dat.ICES.rect_CS.sp_name.PA_2022.rds")

################################################################################
##### extract 100 list of data frame with the community maTrix lumped at the years level.
### Abundance all
List.100.resamp.Com.Dat.Glob.PerYear <- lapply(List.100.resamp.Com.Dat.ICES.rect_all.sp_name, function(x){
  res <-  do.call(rbind,lapply(seq(1,length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  row.names(res) <- Years; res
})
### Abundance Bob
List.100.resamp.Com.Dat.Bob.PerYear <- lapply(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name, function(x){
  res <-  do.call(rbind,lapply(seq(1,length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  row.names(res) <- Years; res
})
### Abundance CS
List.100.resamp.Com.Dat.CS.PerYear <- lapply(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name, function(x){
  res <-  do.call(rbind,lapply(seq(1,length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  row.names(res) <- Years; res
})

### Presence Absence All
List.100.resamp.Com.Dat.Glob.PerYear.PA <-  lapply(List.100.resamp.Com.Dat.ICES.rect_all.sp_name, function(x){
  resPA <- do.call(rbind, lapply(seq(1, length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  resPA[resPA>0] <- 1
  row.names(resPA) <- Years; resPA
})
### Presence Absence Bob
List.100.resamp.Com.Dat.Bob.PerYear.PA <-  lapply(List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name, function(x){
  resPA <- do.call(rbind, lapply(seq(1, length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  resPA[resPA>0] <- 1
  row.names(resPA) <- Years; resPA
})
### Presence Absence Cs
List.100.resamp.Com.Dat.CS.PerYear.PA <- lapply(List.100.resamp.Com.Dat.ICES.rect_CS.sp_name, function(x){
  resPA <- do.call(rbind,lapply(seq(1,length(Years)), function(i){colSums(x[grep(Years[i],rownames(x)),])}))
  resPA[resPA>0] <- 1
  row.names(resPA) <- Years; resPA
})

saveRDS(List.100.resamp.Com.Dat.Glob.PerYear, file = "Results/Year_level/List.100.resamp.Com.Dat.Glob.PerYear.Biomass_2022.rds")
saveRDS(List.100.resamp.Com.Dat.Glob.PerYear.PA, file = "Results/Year_level/List.100.resamp.Com.Dat.Glob.PerYear.PA_2022.rds")

saveRDS(List.100.resamp.Com.Dat.Bob.PerYear,file="Results/List.100.resamp.Com.Dat.Bob.PerYear.Biomass_2022.rds")
saveRDS(List.100.resamp.Com.Dat.Bob.PerYear.PA, file = "Results/List.100.resamp.Com.Dat.Bob.PerYear.PA_2022.rds")

saveRDS(List.100.resamp.Com.Dat.CS.PerYear,file="Results/List.100.resamp.Com.Dat.CS.PerYear.Biomass_2022.rds")
saveRDS(List.100.resamp.Com.Dat.CS.PerYear.PA, file = "Results/List.100.resamp.Com.Dat.CS.PerYear.PA_2022.rds")

################################################################################
###
### Estimate diversity indicators for the invertebrates of the BoB
###
################################################################################

List.MatPA_All <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_all.sp_name.PA_2022.rds")
List.MatPA_Bob <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_Bob.sp_name.PA_2022.rds")
List.MatPA_CS <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_all.sp_name.PA_2022.rds")

#### Estimate the average SR for each ICES rectangle for each year.
all.SR <- do.call(cbind,lapply(List.MatPA_All, function(x){x[x>0] <- 1 ; rowSums(x)}))
row.names(all.SR) <- row.names(List.MatPA_All[[1]])
Av.SR.ICES.perYear.resamp <- apply(all.SR,1,mean,na.rm=T)

Bob.SR <- do.call(cbind,lapply(List.MatPA_Bob, function(x){x[x>0] <- 1 ; rowSums(x)}))
row.names(Bob.SR) <- row.names(List.MatPA_Bob[[1]])
Bob_Av.SR.ICES.perYear.resamp <- apply(Bob.SR ,1,mean,na.rm=T)

CS.SR <- do.call(cbind,lapply(List.MatPA_CS, function(x){x[x>0] <- 1 ; rowSums(x)}))
row.names(CS.SR ) <- row.names(List.MatPA_CS[[1]])
CS_Av.SR.ICES.perYear.resamp <- apply(CS.SR,1,mean,na.rm=T)

### Estimate the average Abundance for each ICES rectangle for each year
all.Abund <- do.call(cbind, lapply(List.MatPA_All, function(x){rowSums(x)}))
row.names(all.Abund) <- row.names(all.Abund[[1]])
Av.Abund.ICES.perYear.resamp <- apply(all.Abund,1,mean,na.rm=T)

Bob.SR <- do.call(cbind, lapply(List.MatPA_Bob, function(x){rowSums(x)}))
row.names(Bob.SR) <- row.names(List.MatPA_Bob[[1]])
Av.Abund.ICES.perYear.resamp_Bob <- apply(Bob.SR,1,mean,na.rm=T)

CS.SR <- do.call(cbind, lapply(List.MatPA_CS, function(x){rowSums(x)}))
row.names(CS.SR) <- row.names(List.MatPA_CS[[1]])
Av.Abund.ICES.perYear.resamp_CS <- apply(CS.SR,1,mean,na.rm=T)

### Estimate the average Shannon index for each ICES rectangle for each year
### Global
all.Shannon <- do.call(cbind, lapply(List.MatPA_All, function(x){diversity(x)}))
row.names(all.Shannon) <- row.names(List.MatPA_All[[1]])
Av.Shannon.ICES.perYear.resamp <- apply(all.Shannon,1,mean,na.rm=T)

### Bob
Bob.Shannon <- do.call(cbind, lapply(List.MatPA_Bob, function(x){diversity(x)}))
row.names(Bob.Shannon) <- row.names(List.MatPA_Bob[[1]])
Av.Shannon.ICES.perYear.resamp_Bob <- apply(Bob.Shannon,1,mean,na.rm=T)

### CS
CS.Shannon <- do.call(cbind, lapply(List.MatPA_CS, function(x){diversity(x)}))
row.names(CS.Shannon) <- row.names(List.MatPA_CS[[1]])
Av.Shannon.ICES.perYear.resamp_CS <- apply(CS.Shannon,1,mean,na.rm=T)

### Estimate the average Simpson index (1-D) for each ICES rectangle for each year, 
all.Simpson <- do.call(cbind, lapply(List.MatPA_All, function(x){diversity(x,index="simpson")}))
row.names(all.Simpson) <- row.names(List.MatPA_All[[1]])
Av.Simpson.ICES.perYear.resamp <- apply(all.Simpson,1,mean,na.rm=T)

### Bob
Bob.Simpson <- do.call(cbind, lapply(List.MatPA_Bob, function(x){diversity(x,index="simpson")}))
row.names(Bob.Simpson) <- row.names(List.MatPA_Bob[[1]])
Av.Simpson.ICES.perYear.resamp_Bob <- apply(Bob.Simpson,1,mean,na.rm=T)

### CS
CS.Simpson <- do.call(cbind,lapply(List.MatPA_CS, function(x){diversity(x,index="simpson")}))
row.names(CS.Simpson) <- row.names(List.MatPA_CS[[1]])
Av.Simpson.ICES.perYear.resamp_CS <- apply(CS.Simpson,1,mean,na.rm=T)

### Compile the final data set
### To match the trawl information witrh the SR
Trawls.info.2 <- Trawls.info.2[which(duplicated(Trawls.info.2$Poly.Annee)==F),]
rownames(Trawls.info.2) <- Trawls.info.2$Poly.Annee
Trawls.info.2 <- Trawls.info.2[names(Av.SR.ICES.perYear.resamp),]


Av.Resamp.ICES.perYear.DF <- data.frame(Trawls.info.2,Av.SR.resamp = Av.SR.ICES.perYear.resamp, 
                                        log10.Av.Abund.resamp = log10(Av.Abund.ICES.perYear.resamp), 
                                        Av.Shannon.resamp = Av.Shannon.ICES.perYear.resamp, 
                                        Av.Simpson.resamp = Av.Simpson.ICES.perYear.resamp) 

head(Av.Resamp.ICES.perYear.DF)
saveRDS(Av.Resamp.ICES.perYear.DF, file = "Data/AV.resamp.Alpha.Div_2022.rds")

##### Average over the time series.
PolyID.uni <-  unique(Av.Resamp.ICES.perYear.DF$ICES_ID)

Av.Resamp.ICES.AvTimeSeries.DF <- do.call(rbind,lapply(seq(1,length(PolyID.uni)),function(x){
  DFtemp <- Av.Resamp.ICES.perYear.DF[which(Av.Resamp.ICES.perYear.DF$ICES_ID == PolyID.uni[x]),]
  apply(DFtemp[,c(14:17)],2,mean,na.rm=T)}))

Av.Resamp.ICES.AvTimeSeries.DF <- data.frame(PolyID=PolyID.uni,Av.Resamp.ICES.AvTimeSeries.DF)
saveRDS(Av.Resamp.ICES.AvTimeSeries.DF, file="Results/Av.Resamp.ICES.AvTimeSeries.DF_2022.rds")

#### Compute the temporal trend of SR, abudnance and evenness indices.
source("Scripts/Time.trend.R")
Resp.var = c(14:17)
family = rep("gaussian", 4) #, rep("betareg", 1))
ylab="";xlab=""

Temporal.Trend.SR.etc.DF <- lapply(seq(1,length(Resp.var)), function(i){
  res <-  do.call(rbind, lapply(seq(1, length(PolyID.uni)), function(x){
    DFtemp <- Av.Resamp.ICES.perYear.DF[which(Av.Resamp.ICES.perYear.DF$ICES_ID==PolyID.uni[x]),]
    DFtemp$Year <- as.numeric(as.character(DFtemp$Annee))
   
     if(dim(DFtemp)[1] < 4){
      c(as.numeric(as.character(PolyID.uni[x])), rep(NA, 14))
    } else {
      res1 <- Time.trend(input=DFtemp,time.s.nb=3,Resp.Var.nb=Resp.var[i],family=family[i],k.basis =3,
                         plot=T,xlab=xlab[i],ylab=ylab[i],cex.axis=0.8,tck=-0.02,mgp=c(1.5,0.5,0))
      c(as.numeric(as.character(PolyID.uni[x])),res1)
    }
  }))
  colnames(res) = c("ID", "Slope.lm", "p.val.slope", "Category", "edf.gam", "p.val.edf", 
                    "Periodicity.obsStat", "Periodicity.pvalue", "Periodicity.binary", 
                    "Mann-Kendall.S.stat", "Var.S.stat", "Tau", "Mann-Kendal.p.value", 
                    "Number.of.year", "Intercept.lm")
  res
})

names(Temporal.Trend.SR.etc.DF) <- names(Av.Resamp.ICES.perYear.DF)[Resp.var]
saveRDS(Temporal.Trend.SR.etc.DF, file = "Data/Temporal.Trend.SR.etc.DF_2022.rds")

################################################################################
#### 
####       Beta diversity analyses
####
################################################################################

#Compute the beta.jac, beta.jutu, beta.jne and beta.ratio = beta.jtu/beta.jac.
# Converte the 100 community data matrix in PA.
List.MatPA_All_PA <- lapply(List.MatPA_All,function(x){x[x>0] = 1; x})

Ave.Resamp.Beta.Div.Year <- do.call(rbind, lapply(List.MatPA_All_PA, function(i){
  Ave.Beta.Div.Year <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
    DFtemp <- i[which(Trawls.info.2[,3]== Years[x]),]
    a <- beta.pair(DFtemp, index.family = "jaccard") [c(3,1,2)]
    a[[4]] <- a[[2]]/a[[1]] #beta ratio beta.jtu/beta.jac
    names(a) <- c("beta.jac", "beta.jtu", "beta.jne", "beta.ratio")
    r <- unlist(lapply(a, mean, na.rm = T)) ## average 
    c(r, length(which(colSums(DFtemp) > 0)))
  }))
  Ave.Beta.Div.Year = data.frame(Ave.Beta.Div.Year[,5], Ave.Beta.Div.Year[,1:4])
  names(Ave.Beta.Div.Year) = c("SR", "beta.jac", "beta.jtu", "beta.jne", "beta.ratio")
  Ave.Beta.Div.Year = data.frame(Years = Years, Ave.Beta.Div.Year)
}))

Ave.Resamp.Beta.Div.Year <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  DFtemp = Ave.Resamp.Beta.Div.Year[which(Ave.Resamp.Beta.Div.Year[,1] == Years[x]),]
  apply(DFtemp[,-1], 2, mean)
}))

##### Compute the Ruzicka index (beta.ruz) including abundance and its decomposition in balance variation in abundnace (beta.ruz.bal) and abundance gradient (beta.ruz.gra). Following Baselga 2013 MEE. Ruzicka is the equivalent of the Jaccard index when presnece absence are used.
List.MatBIOM_All <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_all_Unvert_2022.rds")
Resamp.Beta.Div.Abund.Year <- do.call(rbind,mclapply(List.MatBIOM_All,mc.cores=5,function(i){
  Ave.Beta.Div.Year <- do.call(rbind,lapply(seq(1,length(Years)),function(x){
    DFtemp <- i[which(Trawls.info.2[,3]==Years[x]),]
    a <- beta.pair.abund(DFtemp, index.family = "ruzicka") [c(3,1,2)]
    a[[4]] <- a[[2]]/a[[1]] #beta ratio beta.ruz.bal/beta.ruz.gra
    names(a) <- c("beta.ruz", "beta.ruz.bal", "beta.ruz.gra", "beta.ruz.ratio")
    
    ### also performed the beta.jac, beta.jtu and beta.jne. and contrast beta.ruz.bal beta.jtu
    DFtempPA <- DFtemp
    DFtempPA[DFtempPA>0] <- 1
    b <- beta.pair(DFtempPA, index.family = "jaccard") 
    b <- b[c(3,1,2)] # change the order start by beta.jac
    b[[4]] <- b[[2]]/b[[1]] #beta ratio beta.jtu/beta.jac
    names(b) <- c("beta.jac", "beta.jtu", "beta.jne", "beta.ratio")
    
    ### contribution of abundance to each compornent of beta diversity ,
    c <- list(beta.ruz.only = (a[[1]]-b[[1]]), beta.ruz.bal.only = (a[[2]]-b[[2]]), 
              beta.ruz.gra.only = (a[[3]]-b[[3]]), beta.ratio.ruz.ratio.jac = (a[[4]]-b[[4]]))
    
    ra <- unlist(lapply(a, mean, na.rm = T)) ## average 
    rb <- unlist(lapply(b, mean, na.rm = T)) ## average
    rc <- unlist(lapply(c, mean, na.rm = T)) ## average
    
    c(length(which(colSums(DFtemp) > 0)), ra, rb, rc)
  }))
  colnames(Ave.Beta.Div.Year) = c("SR", "beta.ruz", "beta.ruz.bal", "beta.ruz.gra", 
                                  "beta.ruz.ratio", "beta.jac", "beta.jtu", "beta.jne", 
                                  "beta.ratio", "beta.ruz.only", "beta.ruz.bal.only", 
                                  "beta.ruz.gra.only", "beta.ratio.ruz.ratio.jac")
  Ave.Beta.Div.Year <- data.frame(Years=Years, Ave.Beta.Div.Year)
}))

Ave.Resamp.Beta.Div.Abund.Year <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  DFtemp <- Resamp.Beta.Div.Abund.Year[which(Resamp.Beta.Div.Abund.Year[,1] == Years[x]),]
  apply(DFtemp[,-1], 2, mean)
}))

saveRDS(Ave.Resamp.Beta.Div.Abund.Year, file="Results/Ave.Resamp.Beta.Div.Abund.Year_2022.RDS")

### for the CS
List.MatBIOM_CS <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_CS_Unvert_2022.rds")
Trawls.info.Cs <- Trawls.info.2[which(Trawls.info.2$Ecoreg=="Celtic Seas"),]

Resamp.Beta.Div.Abund.Year_CS <- do.call(rbind,mclapply(List.MatBIOM_CS,mc.cores=5,function(i){
  Ave.Beta.Div.Year <- do.call(rbind,lapply(seq(1,length(Years)),function(x){
    cat(x)
    DFtemp <- i[which(Trawls.info.Cs[,3]==Years[x]),]
    a <- beta.pair.abund(DFtemp, index.family = "ruzicka") [c(3,1,2)]
    a[[4]] <- a[[2]]/a[[1]] #beta ratio beta.ruz.bal/beta.ruz.gra
    names(a) <- c("beta.ruz", "beta.ruz.bal", "beta.ruz.gra", "beta.ruz.ratio")
    
    ### also performed the beta.jac, beta.jtu and beta.jne. and contrast beta.ruz.bal beta.jtu
    DFtempPA <- DFtemp
    DFtempPA[DFtempPA>0] <- 1
    b <- beta.pair(DFtempPA, index.family = "jaccard") 
    b <- b[c(3,1,2)] # change the order start by beta.jac
    b[[4]] <- b[[2]]/b[[1]] #beta ratio beta.jtu/beta.jac
    names(b) <- c("beta.jac", "beta.jtu", "beta.jne", "beta.ratio")
    
    ### contribution of abundance to each compornent of beta diversity ,
    c <- list(beta.ruz.only = (a[[1]]-b[[1]]), beta.ruz.bal.only = (a[[2]]-b[[2]]), 
              beta.ruz.gra.only = (a[[3]]-b[[3]]), beta.ratio.ruz.ratio.jac = (a[[4]]-b[[4]]))
    
    ra <- unlist(lapply(a, mean, na.rm = T)) ## average 
    rb <- unlist(lapply(b, mean, na.rm = T)) ## average
    rc <- unlist(lapply(c, mean, na.rm = T)) ## average
    
    c(length(which(colSums(DFtemp) > 0)), ra, rb, rc)
  }))
  colnames(Ave.Beta.Div.Year) = c("SR", "beta.ruz", "beta.ruz.bal", "beta.ruz.gra", 
                                  "beta.ruz.ratio", "beta.jac", "beta.jtu", "beta.jne", 
                                  "beta.ratio", "beta.ruz.only", "beta.ruz.bal.only", 
                                  "beta.ruz.gra.only", "beta.ratio.ruz.ratio.jac")
  Ave.Beta.Div.Year_CS <- data.frame(Years=Years, Ave.Beta.Div.Year)
}))

Ave.Resamp.Beta.Div.Abund.Year_CS <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  DFtemp <- Resamp.Beta.Div.Abund.Year_CS[which(Resamp.Beta.Div.Abund.Year_CS[,1] == Years[x]),]
  apply(DFtemp[,-1], 2, mean)
}))

saveRDS(Ave.Resamp.Beta.Div.Abund.Year_CS, file="Results/Ave.Resamp.Beta.Div.Abund.Year_CS_2002.RDS")

###############################################################################
### for the CS
List.MatBIOM_Bob <- readRDS("Results/List.100.resamp.Com.Dat.ICES.rect_Bob_Unvert_2022.rds")
Trawls.info.Bob <- Trawls.info.2[which(Trawls.info.2$Ecoreg=="Bay of Biscay and the Iberian Coast"),]

Resamp.Beta.Div.Abund.Year_Bob <- do.call(rbind,mclapply(List.MatBIOM_Bob,mc.cores=5,function(i){
  Ave.Beta.Div.Year <- do.call(rbind,lapply(seq(1,length(Years)),function(x){
    cat(x)
    DFtemp <- i[which(Trawls.info.Bob[,3]==Years[x]),]
    a <- beta.pair.abund(DFtemp, index.family = "ruzicka") [c(3,1,2)]
    a[[4]] <- a[[2]]/a[[1]] #beta ratio beta.ruz.bal/beta.ruz.gra
    names(a) <- c("beta.ruz", "beta.ruz.bal", "beta.ruz.gra", "beta.ruz.ratio")
    
    ### also performed the beta.jac, beta.jtu and beta.jne. and contrast beta.ruz.bal beta.jtu
    DFtempPA <- DFtemp
    DFtempPA[DFtempPA>0] <- 1
    b <- beta.pair(DFtempPA, index.family = "jaccard") 
    b <- b[c(3,1,2)] # change the order start by beta.jac
    b[[4]] <- b[[2]]/b[[1]] #beta ratio beta.jtu/beta.jac
    names(b) <- c("beta.jac", "beta.jtu", "beta.jne", "beta.ratio")
    
    ### contribution of abundance to each compornent of beta diversity ,
    c <- list(beta.ruz.only = (a[[1]]-b[[1]]), beta.ruz.bal.only = (a[[2]]-b[[2]]), 
              beta.ruz.gra.only = (a[[3]]-b[[3]]), beta.ratio.ruz.ratio.jac = (a[[4]]-b[[4]]))
    
    ra <- unlist(lapply(a, mean, na.rm = T)) ## average 
    rb <- unlist(lapply(b, mean, na.rm = T)) ## average
    rc <- unlist(lapply(c, mean, na.rm = T)) ## average
    
    c(length(which(colSums(DFtemp) > 0)), ra, rb, rc)
  }))
  colnames(Ave.Beta.Div.Year) = c("SR", "beta.ruz", "beta.ruz.bal", "beta.ruz.gra", 
                                  "beta.ruz.ratio", "beta.jac", "beta.jtu", "beta.jne", 
                                  "beta.ratio", "beta.ruz.only", "beta.ruz.bal.only", 
                                  "beta.ruz.gra.only", "beta.ratio.ruz.ratio.jac")
  Ave.Beta.Div.Year_Bob <- data.frame(Years=Years, Ave.Beta.Div.Year)
}))

Ave.Resamp.Beta.Div.Abund.Year_Bob <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  DFtemp <- Resamp.Beta.Div.Abund.Year_Bob[which(Resamp.Beta.Div.Abund.Year_Bob[,1] == Years[x]),]
  apply(DFtemp[,-1], 2, mean)
}))

saveRDS(Ave.Resamp.Beta.Div.Abund.Year_Bob, file="Results/Ave.Resamp.Beta.Div.Abund.Year_Bob_2002.RDS")

###############################################################################
####### Global scale

### Average SR of the ICES rectangle over the whole area (resampled 100 times) 
Av.Resamp.Global <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  a <- Av.Resamp.ICES.perYear.DF[which(Av.Resamp.ICES.perYear.DF$Annee == Years[x]),]
  apply(a[,c(14:17)], 2, mean, na.rm =T)
}))

Av.Resamp.CS <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  a<- Av.Resamp.ICES.perYear.DF[which(Av.Resamp.ICES.perYear.DF$Ecoreg=="Celtic Seas"),]
  a <- a[which(a$Annee == Years[x]),]
  apply(a[,c(14:17)], 2, mean, na.rm =T)
}))

Av.Resamp.Bob <- do.call(rbind, lapply(seq(1, length(Years)), function(x){
  a <- Av.Resamp.ICES.perYear.DF[which(Av.Resamp.ICES.perYear.DF$Ecoreg=="Bay of Biscay and the Iberian Coast"),]
  a <- a[which(a$Annee == Years[x]),]
  apply(a[,c(14:17)], 2, mean, na.rm =T)
}))


### Build a single table with all the indicators at the global scale for each year.
Av.Resamp.Global.DF <- data.frame(Years, Av.Resamp.Global, Ave.Resamp.Beta.Div.Year,
                                  Ave.Resamp.Beta.Div.Abund.Year[,-c(1, 6:9)])
Av.Resamp.CS.DF <- data.frame(Years, Av.Resamp.CS, Ave.Resamp.Beta.Div.Abund.Year_CS[,c(6:9)],
                              Ave.Resamp.Beta.Div.Abund.Year_CS[,c(2:5,10:13)])
Av.Resamp.Bob.DF <- data.frame(Years, Av.Resamp.Bob, Ave.Resamp.Beta.Div.Abund.Year_Bob[,c(6:9)],
                              Ave.Resamp.Beta.Div.Abund.Year_Bob[,c(2:5,10:13)])

#Av.Resamp.Global.DF <- data.frame(Years, Av.Resamp.Global, Ave.Resamp.Beta.Div.Year)

saveRDS(Av.Resamp.Global.DF, file = "Data/Av.Resamp.Global.DF_2022.rds")
saveRDS(Av.Resamp.CS.DF, file="Data/Av.Resamp.CS.DF_2022.rds")
saveRDS(Av.Resamp.Bob.DF, file="Data/Av.Resamp.Bob.DF_2022.rds")

################################################################################@
Resp.Var.nNB = c(2:5, 7:10)
familly = rep("gaussian", 9)#, rep("betareg", 5))
Ylab = names(Av.Resamp.Global.DF)[Resp.Var.nNB]

#source("../Data_evhoe/Scripts/Time.trend.R")
par(mfrow=c(3,3))
Result.Global.Eff.Fishing.DF = do.call(rbind, lapply(seq(1, length(Resp.Var.nNB)), function(x){
  Time.trend(input = Av.Resamp.Global.DF, time.s.nb = 1,
             Resp.Var.nb = Resp.Var.nNB[x], family = familly[x], plot = TRUE, xlab = "Year", 
             ylab = Ylab[x], cex.axis = 0.8, tck = -0.02, mgp = c(1.5, 0.5, 0), k.basis = 3)
}))

dev.copy2pdf(device = cairo,  file = "Figures/Resampled.Metrics.GlobalScale_2022.pdf")

################################################################################@

Resp.Var.nNB = c(2:9)
familly = rep("gaussian", 8)#, rep("betareg", 5))
Ylab = names(Av.Resamp.CS.DF)[Resp.Var.nNB]

par(mfrow=c(3,3))
Result.CS.Eff.Fishing.DF = do.call(rbind, lapply(seq(1, length(Resp.Var.nNB)), function(x){
  Time.trend(input = Av.Resamp.CS.DF, time.s.nb = 1,
             Resp.Var.nb = Resp.Var.nNB[x], family = familly[x], plot = TRUE, xlab = "Year", 
             ylab = Ylab[x], cex.axis = 0.8, tck = -0.02, mgp = c(1.5, 0.5, 0), k.basis = 3)
}))

dev.copy2pdf(device = cairo,  file = "Figures/Resampled.Metrics.CS.pdf")
###############################################################################@
Resp.Var.nNB = c(2:9)
familly = rep("gaussian", 8)#, rep("betareg", 5))
Ylab = names(Av.Resamp.CS.DF)[Resp.Var.nNB]

par(mfrow=c(3,3))
Result.Bob.Eff.Fishing.DF = do.call(rbind, lapply(seq(1, length(Resp.Var.nNB)), function(x){
  Time.trend(input = Av.Resamp.Bob.DF, time.s.nb = 1,
             Resp.Var.nb = Resp.Var.nNB[x], family = familly[x], plot = TRUE, xlab = "Year", 
             ylab = Ylab[x], cex.axis = 0.8, tck = -0.02, mgp = c(1.5, 0.5, 0), k.basis = 3)
}))

dev.copy2pdf(device = cairo,  file = "Figures/Resampled.Metrics.Bob.pdf")

##############################################################################@


#### including beta diversity indices taking into account the abundance.
Resp.Var.nNB <- c(11:18)
familly <- rep("betareg", 11)
Ylab <- names(Av.Resamp.Global.DF)[Resp.Var.nNB]

#source("../Data_evhoe/Scripts/Time.trend.R")
par(mfrow=c(2,4))
Result.Global.Eff.Fishing.DF = do.call(rbind, lapply(seq(1, length(Resp.Var.nNB)), function(x){
  Time.trend(input = Av.Resamp.Global.DF, time.s.nb = 1,
             Resp.Var.nb = Resp.Var.nNB[x], family = familly[x], plot = TRUE, xlab = "Year", 
             ylab = Ylab[x], cex.axis = 0.8, tck = -0.02, mgp = c(1.5, 0.5, 0), k.basis = 3)
}))


dev.copy2pdf(device = cairo,  file = "Figures/Resampled.Metrics.GlobalScale.DEta.Abund.NewFrameWork.pdf")





