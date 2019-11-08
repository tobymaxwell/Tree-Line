library(raster)
library(sp)
library(rgdal)
library(rgeos)
oregon<-readOGR("/Users/tobymaxwell/Downloads/tl_2017_41_place/tl_2017_41_place.shp")
##Turn points into polygons
library("sp")
library("rgdal")
library("raster")
gps<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.csv")

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
coords(Poly)

setwd("/Users/tobymaxwell/Downloads/PRISM_ppt_30yr_normal_800mM2_annual_bil/")
list<-list.files(pattern="bil.bil*")
list<-list[c(TRUE,FALSE)]
rasterlist<-NULL
results<-NULL

ppt<-readGDAL(list)
ppt.raster<-raster(ppt)
ppt.data<-extract(ppt.raster, Poly)
ppt.df<-as.data.frame(sapply(ppt.data, FUN=mean))
ppt.df
ppt.df$Site<-rownames(ppt.df)<-Poly$Name
colnames(ppt.df)<-c("ppt", "Site")

########################MAT#########################
setwd("/Users/tobymaxwell/Downloads/PRISM_tmean_30yr_normal_800mM2_annual_bil/")
list<-list.files(pattern="bil.bil*")
list<-list[c(TRUE,FALSE)]

mat<-readGDAL(list)
mat.raster<-raster(mat)
mat.data<-extract(mat.raster, Poly)
mat.df<-as.data.frame(sapply(mat.data, FUN=mean))
mat.df
mat.df$Site<-rownames(mat.df)<-Poly$Name
colnames(mat.df)<-c("mat","Site")
mat.df
clim<-merge(mat.df, ppt.df)


###########################MinT#######################
setwd("/Users/tobymaxwell/Downloads/PRISM_tmin_30yr_normal_800mM2_annual_bil/")
list<-list.files(pattern="bil.bil*")
list<-list[c(TRUE,FALSE)]

tmin<-readGDAL(list)
tmin.raster<-raster(tmin)
tmin.data<-extract(tmin.raster, Poly)
tmin.df<-as.data.frame(sapply(tmin.data, FUN=mean))
tmin.df
tmin.df$Site<-rownames(tmin.df)<-Poly$Name
colnames(tmin.df)<-c("tmin","Site")
tmin.df

clim<-merge(clim, tmin.df)



####################tmax##############################
setwd("/Users/tobymaxwell/Downloads/PRISM_tmax_30yr_normal_800mM2_annual_bil/")
list<-list.files(pattern="bil.bil*")
list<-list[c(TRUE,FALSE)]

tmax<-readGDAL(list)
tmax.raster<-raster(tmax)
tmax.data<-extract(tmax.raster, Poly)
tmax.df<-as.data.frame(sapply(tmax.data, FUN=mean))
tmax.df
tmax.df$Site<-rownames(tmax.df)<-Poly$Name
colnames(tmax.df)<-c("tmax","Site")
tmax.df

clim<-merge(clim, tmax.df)


####################vpd##############################
setwd("/Users/tobymaxwell/Downloads/PRISM_vpdmax_30yr_normal_800mM2_annual_bil/")
list<-list.files(pattern="bil.bil*")
list<-list[c(TRUE,FALSE)]

vpdmax<-readGDAL(list)
vpdmax.raster<-raster(vpdmax)
vpdmax.data<-extract(vpdmax.raster, Poly)
vpdmax.df<-as.data.frame(sapply(vpdmax.data, FUN=mean))
vpdmax.df
vpdmax.df$Site<-rownames(vpdmax.df)<-Poly$Name
colnames(vpdmax.df)<-c("vpdmax","Site")
vpdmax.df

clim<-merge(clim, vpdmax.df)

write.csv(clim, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/clim.norms.csv")
