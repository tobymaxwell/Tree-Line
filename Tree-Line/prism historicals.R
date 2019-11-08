library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(plyr)
##Turn points into polygons
library("sp")
library("rgdal")
library("raster")
gps<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.csv")
oregon<-readOGR("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/oregon shape/tl_2017_41_cousub.shp")

#loop it all up
str(gps)
Sites<-levels(gps$Site)
Aspect<-levels(gps$Aspect)
Zone<-levels(gps$Zone)
Poly<-SpatialPolygons(list())
names.df<-NA
for(i in Sites){
  for(j in Aspect){
    for(k in Zone){
      loc<-gps[gps$Site==i&gps$Aspect==j&gps$Zone==k,]
      dat <- loc[,8:9]
      ch <- chull(dat)
      coords <- dat[c(ch, ch[1]), ]  # closed polygon
      loc.poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=paste0(i,k,j))), proj4string = oregon@proj4string)
      Poly<-bind(loc.poly, Poly)
    }
  }
}
str(gps)
gps$SAZ<-as.factor(paste0(gps$Site, gps$Aspect, gps$Zone))
names<-data.frame(sort(levels(gps$SAZ), decreasing=T))
colnames(names)<-c("Name")
Poly<-SpatialPolygonsDataFrame(Poly, data=names)
library(maptools)
data.frame(Poly)

Poly
plot(oregon)
plot(Poly, add=T,col='red', lwd=5)
library(prism)
#options(prism.path="/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/")
#get_prism_monthlys(type="tmean", years=1895:2014, mon=c(1:12))
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/Tmean/")
folders<-list.files(pattern="*bil")
folders

rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  list<-list.files(i, pattern="*bil.bil")[1]
  temp<-readGDAL(paste0(i,"/",list))
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAT.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(dd)<-Poly$Name
dd$Year<-rep(1895:2014, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "MAT", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
ggplot(dd.long[dd.long$ID=="YNB",], aes(y=MAT, x=date, color=ID))+geom_line()




####Tmin#####

#options(prism.path="/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin")
#get_prism_monthlys(type="tmin", years=1895:2014, mon=c(1:12))
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin/")
folders<-list.files(pattern="*bil")
folders

rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  list<-list.files(i, pattern="*bil.bil")[1]
  temp<-readGDAL(paste0(i,"/",list))
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAT.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
tmin  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(tmin)<-Poly$Name
tmin$Year<-rep(1895:2014, each=12)
tmin$Month<-rep(1:12)
str(tmin)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
tmin$monthname<-rep(months)
tmin$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)

tmin.long<-gather(tmin, key="ID", value = "tmin", YSF:DNB)
str(tmin.long)
tmin.long$date<-paste0(tmin.long$Year, "-", tmin.long$Month, "-01")
tmin.long$date<-as.Date(tmin.long$date)
tmin.long$Site<-substr(tmin.long$ID,1,1)
tmin.long$Aspect<-substr(tmin.long$ID,2,2)
tmin.long$Zone<-substr(tmin.long$ID,3,3)
tmin.long

######TMax#####
#options(prism.path="/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax")
#get_prism_monthlys(type="tmax", years=1895:2014, mon=c(1:12), keepZip = FALSE)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax/")
folders<-list.files(pattern="*bil")
folders

rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  list<-list.files(i, pattern="*bil.bil")[1]
  temp<-readGDAL(paste0(i,"/",list))
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("Tmax.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
tmax  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(tmax)<-Poly$Name
tmax$Year<-rep(1895:2014, each=12)
tmax$Month<-rep(1:12)
str(tmax)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
tmax$monthname<-rep(months)
tmax$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))


tmax.long<-gather(tmax, key="ID", value = "tmax", YSF:DNB)
str(tmax.long)
tmax.long$date<-paste0(tmax.long$Year, "-", tmax.long$Month, "-01")
tmax.long$date<-as.Date(tmax.long$date)
tmax.long$Site<-substr(tmax.long$ID,1,1)
tmax.long$Aspect<-substr(tmax.long$ID,2,2)
tmax.long$Zone<-substr(tmax.long$ID,3,3)
tmax.long
str(tmax.long)

######PPT#####
#options(prism.path="/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt")
#get_prism_monthlys(type="ppt", years=1895:2014, mon=c(1:12), keepZip = FALSE)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt/")
folders<-list.files(pattern="*bil")
folders

rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  list<-list.files(i, pattern="*bil.bil")[1]
  temp<-readGDAL(paste0(i,"/",list))
  temp.raster<-raster(temp)
  assign(paste0("PPT.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("PPT.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
ppt  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(ppt)<-Poly$Name
ppt$Year<-rep(1895:2014, each=12)
ppt$Month<-rep(1:12)
str(ppt)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
ppt$monthname<-rep(months)
ppt$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)

ppt.long<-gather(ppt, key="ID", value = "MAP", YSF:DNB)
str(ppt.long)
ppt.long$date<-paste0(ppt.long$Year, "-", ppt.long$Month, "-01")
ppt.long$date<-as.Date(ppt.long$date)
ppt.long$Site<-substr(ppt.long$ID,1,1)
ppt.long$Aspect<-substr(ppt.long$ID,2,2)
ppt.long$Zone<-substr(ppt.long$ID,3,3)
ppt.long

clim.hist<-merge(ppt.long, dd.long)
clim.hist<-merge(clim.hist, tmax.long)
clim.hist<-merge(clim.hist, tmin.long)
str(clim.hist)
#write.csv(clim.hist, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")
