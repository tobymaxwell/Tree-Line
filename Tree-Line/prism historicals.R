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
write.csv(dd.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.long.csv")



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
write.csv(tmin.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.long.csv")
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

library(tidyr)
tmax.long<-gather(tmax, key="ID", value = "tmax", YSF:DNB)
str(tmax.long)
tmax.long$date<-paste0(tmax.long$Year, "-", tmax.long$Month, "-01")
tmax.long$date<-as.Date(tmax.long$date)
tmax.long$Site<-substr(tmax.long$ID,1,1)
tmax.long$Aspect<-substr(tmax.long$ID,2,2)
tmax.long$Zone<-substr(tmax.long$ID,3,3)
tmax.long
str(tmax.long)
write.csv(tmax.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.long.csv")

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
#write.csv(ppt.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.long.csv")


######VPD######

library(prism)
options(prism.path="/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/")
#get_prism_monthlys(type="vpdmax", years=1895:2014, mon=c(1:12))
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/vpdmax/")
folders<-list.files(pattern="*bil.bil")
folders<-folders[c(TRUE, FALSE)]
tail(folders,200)
rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  temp<-readGDAL(paste0(i,"/"))
  temp.raster<-raster(temp)
  assign(paste0("vpdmax.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("vpdmax.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
dd<-dd[-1489:-1494,]
colnames(dd)<-Poly$Name
dd$Year<-rep(1895:2018, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "vpdmax", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
library(ggplot2)
ggplot(dd.long, aes(y=vpdmax, x=date, color=ID))+geom_line()
#write.csv(dd.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/vpdmax.long.csv")


######ppt recent years######
library(prism)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.modern/")
folders<-list.files(pattern="*bil.bil")
folders<-folders[c(TRUE, FALSE)]
tail(folders,200)
rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  temp<-readGDAL(paste0(i,"/"))
  temp.raster<-raster(temp)
  assign(paste0("ppt.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("ppt.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(dd)<-Poly$Name
dd$Year<-rep(2015:2018, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "ppt", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
#write.csv(dd.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.modern.csv" )
#ppt.modern<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.modern.csv")[-1]
#ppt.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.long.csv")[-1]
#ppt.long<-rbind(ppt.long, ppt.modern)
#write.csv(ppt.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.long.csv" )

######tmax recent years######
library(prism)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.modern/")
folders<-list.files(pattern="*bil.bil")
folders<-folders[c(TRUE, FALSE)]
tail(folders,200)
rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  temp<-readGDAL(paste0(i,"/"))
  temp.raster<-raster(temp)
  assign(paste0("tmax.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("tmax.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(dd)<-Poly$Name
dd$Year<-rep(2015:2018, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "tmax", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
#write.csv(dd.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.modern.csv" )
#tmax.modern<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.modern.csv")[-1]
#tmax.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.long.csv")[-1]
#tmax.long<-rbind(tmax.long, tmax.modern)
str(tmax.long)
tail(ppt.long)
#write.csv(tmax.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.long.csv" )


######tmin recent years######
library(prism)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.modern/")
folders<-list.files(pattern="*bil.bil")
folders<-c(folders[c(TRUE, FALSE)])
tail(folders,200)
rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  temp<-readGDAL(paste0(i,"/"))
  temp.raster<-raster(temp)
  assign(paste0("tmin.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("tmin.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(dd)<-Poly$Name
dd$Year<-rep(2015:2018, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "tmin", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
write.csv(dd.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.modern.csv" )
#tmin.modern<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.modern.csv")[-1]
#tmin.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.long.csv")[-1]
#tmin.long<-rbind(tmin.long, tmin.modern)
str(tmin.long)
tail(tmin.long)
#write.csv(tmin.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.long.csv" )




######tmean recent years######
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.modern/")
folders<-list.files(pattern="*bil.bil")
folders<-folders[c(TRUE, FALSE)]
tail(folders,200)
rasterlist<-NULL
results<-NULL
extracts<-NULL

for(i in folders){
  temp<-readGDAL(paste0(i,"/"))
  temp.raster<-raster(temp)
  assign(paste0("tmean.", data.frame(strsplit(list,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("tmean.", data.frame(strsplit(list,"_"))[5,]))
  extracts<-raster::extract(temp.raster, Poly)
  extracts<-lapply(extracts,FUN=mean)
  results<-rbind(results, extracts)
}

results<-data.frame(results)
dd  <-  as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1]))))
colnames(dd)<-Poly$Name
dd$Year<-rep(2015:2018, each=12)
dd$Month<-rep(1:12)
str(dd)
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
dd$monthname<-rep(months)
dd$monthname<-factor(dd$monthname, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
str(dd)
library(tidyr)
dd.long<-gather(dd, key="ID", value = "tmean", YSF:DNB)
str(dd.long)
dd.long$date<-paste0(dd.long$Year, "-", dd.long$Month, "-01")
dd.long$date<-as.Date(dd.long$date)
dd.long$Site<-substr(dd.long$ID,1,1)
dd.long$Aspect<-substr(dd.long$ID,2,2)
dd.long$Zone<-substr(dd.long$ID,3,3)
dd.long
#write.csv(dd.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.modern.csv" )
#tmean.modern<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.modern.csv")[-1]
#tmean.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.long.csv")[-1]
#tmean.long<-rbind(tmean.long, tmean.modern)
str(tmean.long)
str(tmean.modern)
tail(tmax.long)
#write.csv(tmean.long,"/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.long.csv" )


####aggregate all sites#####

ppt.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/ppt.long.csv")[-1]
tmean.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmean.long.csv")[-1]
tmax.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmax.long.csv")[-1]
tmin.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/tmin.long.csv")[-1]
vpdmax.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/vpdmax.long.csv")[-1]


clim.hist<-merge(ppt.long, tmean.long, by=c("ID", "Site", "Aspect", "Zone", "Year", "Month", "monthname", 'date'))
str(clim.hist)
clim.hist<-merge(clim.hist, tmax.long, by=c("ID", "Site", "Aspect", "Zone", "Year", "Month", "monthname"), all=T)
clim.hist<-merge(clim.hist, tmin.long, by=c("ID", "Site", "Aspect", "Zone", "Year", "Month", "monthname"))
library(tidyr)
clim.hist<-clim.hist[c(-8, -12)]
clim.hist<-merge(clim.hist, vpdmax.long, by=c("ID", "Site", "Aspect", "Zone", "Year", "Month", "monthname"))
str(clim.hist)
tail(clim.hist)
#write.csv(clim.hist, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")
clim.hist<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")
str(clim.hist)
