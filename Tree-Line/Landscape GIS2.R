#Methods
library(raster)
library(rgeos)
library(rgdal)
library(maptools)
library(sp)
library(latticeExtra)
library(GISTools)
oregon<-readOGR("/Users/Maxwell/Downloads/tl_2017_41_cousub/tl_2017_41_cousub.shp")
str(oregon)
dem<-raster("/Users/Maxwell/OneDrive - University Of Oregon/Oregon/DEM/DEM_10m_OR_TIFF/DEM_10m_OR", )
plot(dem)

plot(Hood)

OR<-spTransform(oregon, crs(dem))
Hood<-OR[OR$NAME=="Parkdale",]
plot(Hood)
dem.hood<-crop(dem, Hood)
plot(dem.hood)
plot(Hood, add=T)
slope.Hood<-terrain(dem.hood)
plot(slope.Hood)
aspect.Hood<-terrain(dem.hood, opt='aspect')
plot(aspect.Hood)
twi<-build_chans(atb=dem.hood, drn=NULL, atb.thresh=.95)
plot(twi[[1]])
ua<-upslope.area(dem.hood)

library(dynatopmodel)
plot(slope)
slope<-terrain(dem)
aspect<-terrain(dem, opt='aspect')
tpi<-terrain(dem, opt='TPI')
tri<-terrain(dem, opt='TRI')
plot(aspect)
plot(slope)
plot(tri)
plot(tpi)

##Turn points into polygons
library("sp")
library("rgdal")
library("raster")
gps<-read.csv("/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.csv")

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
      dat <- loc[,9:10]
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
setwd("/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Site Shapes/")
shapenames<-levels(names$Name)
for(i in shapenames){
writeSpatialShape(Poly[Poly$Name==i,], i)
}
#dem.latlong<-projectRaster(dem, crs=crs(Poly))
library(tmap)
plot(dem.latlong)
writeSpatialShape(Poly)
plot(oregon)
data.frame(Poly)
plot(Poly, add=T, lwd=5, col="red")

slope
gps

plot(dem.latlong)
drake.utm<-drawExtent()
drake.utm.dem<-crop(dem, drake.utm)
plot(drake.utm.dem)
drake.latlong<-projectRaster(drake.utm.dem, crs=crs(Poly))
drake<-drawExtent()
dem.drake<-crop(drake.latlong, drake)
plot(dem.drake)
plot(Poly, add=T, lwd=5)
terrain.drake<-terrain(dem.drake, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
plot(terrain.drake)
drake.df<-data.frame(extract(terrain.drake, Poly[25:28,], fun=mean))
drake.df
aspect.drake<-terrain(dem.drake, opt='aspect')
plot(aspect.drake)
plot(Poly, add=T)
drake.df<-data.frame(extract(terrain.drake, Poly[25:28,], fun=mean))
rownames(drake.df)<-Poly[25:28,]$Name

drake.elevations<-extract(dem.drake, gps[gps$Site=='D',9:10])*.3048 #extract elevations and turn to meters
gps[gps$Site=='D',]$elevation<-drake.elevations

plot(dem.latlong)
hood<-drawExtent()
elkhorn<-drawExtent()
tumalo<-drawExtent()
wallowa<-drawExtent()
dem.hood<-crop(dem.latlong, hood)
dem.elkhorn<-crop(dem.latlong, elkhorn)
dem.tumalo<-crop(dem.latlong, tumalo)
dem.wallowa<-crop(dem.latlong, wallowa)

terrain.drake<-terrain(dem.drake, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
terrain.hood<-terrain(dem.hood, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
terrain.elkhorn<-terrain(dem.elkhorn, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
terrain.tumalo<-terrain(dem.tumalo, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
terrain.wallowa<-terrain(dem.wallowa, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')
plot(terrain.elkhorn)

gps$slope<-NA
gps$aspect.deg<-NA
gps$tpi<-NA
gps$tri<-NA

df.drake<-data.frame((extract(terrain.drake, gps[gps$Site=="D",9:10])),gps[gps$Site=="D"] )
gps[gps$Site=="D",12:15]
gps[gps$Site=="H",12:15]<-(extract(terrain.hood, gps[gps$Site=="H",9:10]))
gps[gps$Site=="D",12:15]<-(extract(terrain.drake, gps[gps$Site=="D",9:10]))
gps[gps$Site=="D",12:15]<-(extract(terrain.drake, gps[gps$Site=="D",9:10]))
gps[gps$Site=="D",12:15]<-(extract(terrain.drake, gps[gps$Site=="D",9:10]))
gps[gps$Site=="D",12:15]<-(extract(terrain.drake, gps[gps$Site=="D",9:10]))
gps[gps$Site=="D",12:15]<-(extract(terrain.drake, gps[gps$Site=="D",9:10]))

sites<-list(hood, elkhorn, tumalo, wallowa)
sites.names<-c("hood", "elkhorn", "tumalo", "wallowa")
sapply(sites.names, paste, collapse=":")
for(i in sites){
  assign(x=paste0("dem.", paste0(i)), crop(dem.latlong, i))
  }
}

dem.drake<-crop(drake.latlong, drake)
plot(dem.drake)
plot(Poly, add=T, lwd=5)
terrain.drake<-terrain(dem.drake, opt=c("slope","aspect", "tri", "tpi"), unit='degrees')


plot(dem.hood)
plot(Hood, add=T)
slope.Hood<-terrain(dem.hood)
plot(slope.Hood)
aspect.Hood<-terrain(dem.hood, opt='aspect')
plot(aspect.Hood)
plot(Poly[9:12,], lwd=5, add=T)
data.frame(Poly)
rad2deg(extract(aspect.Hood, Poly[17:20,]))
twi<-build_chans(atb=dem.hood, drn=NULL, atb.thresh=.95)
plot(twi[[1]])
ua<-upslope.area(dem)

str(Poly)
dem.drake
drake.dat<-extract(terrain.drake, Poly, fun=mean)
drake.sd<-extract(terrain.drake, Poly, fun=sd)
drake.elev<-extract(dem.drake, Poly, fun=mean)
drake.dat<-cbind(drake.dat, drake.sd, drake.elev)
rownames(drake.dat)<-Poly$Name
colnames(drake.dat)<-c("tri", "tpi", "slope", "aspect", "trisd", "tpisd", "slopesd", "aspectsd", "elevation")
drake.dat<-data.frame(na.omit(drake.dat))

dem.hood
hood.dat<-extract(terrain.hood, Poly, fun=mean)
hood.sd<-extract(terrain.hood, Poly, fun=sd)
hood.elev<-extract(dem.hood, Poly, fun=mean)
hood.dat<-cbind(hood.dat, hood.sd, hood.elev)
rownames(hood.dat)<-Poly$Name
colnames(hood.dat)<-c("tri", "tpi", "slope", "aspect", "trisd", "tpisd", "slopesd", "aspectsd", "elevation")
hood.dat<-data.frame(na.omit(hood.dat))

dem.elkhorn
elkhorn.dat<-extract(terrain.elkhorn, Poly, fun=mean)
elkhorn.sd<-extract(terrain.elkhorn, Poly, fun=sd)
elkhorn.elev<-extract(dem.elkhorn, Poly, fun=mean)
elkhorn.dat<-cbind(elkhorn.dat, elkhorn.sd, elkhorn.elev)
rownames(elkhorn.dat)<-Poly$Name
colnames(elkhorn.dat)<-c("tri", "tpi", "slope", "aspect", "trisd", "tpisd", "slopesd", "aspectsd", "elevation")
elkhorn.dat<-data.frame(na.omit(elkhorn.dat))

dem.wallowa
wallowa.dat<-extract(terrain.wallowa, Poly, fun=mean)
wallowa.sd<-extract(terrain.wallowa, Poly, fun=sd)
wallowa.elev<-extract(dem.wallowa, Poly, fun=mean)
wallowa.dat<-cbind(wallowa.dat, wallowa.sd, wallowa.elev)
rownames(wallowa.dat)<-Poly$Name
colnames(wallowa.dat)<-c("tri", "tpi", "slope", "aspect", "trisd", "tpisd", "slopesd", "aspectsd", "elevation")
wallowa.dat<-data.frame(na.omit(wallowa.dat))

dem.tumalo
tumalo.dat<-extract(terrain.tumalo, Poly, fun=mean)
tumalo.sd<-extract(terrain.tumalo, Poly, fun=sd)
tumalo.elev<-extract(dem.tumalo, Poly, fun=mean)
tumalo.dat<-cbind(tumalo.dat, tumalo.sd, tumalo.elev)
rownames(tumalo.dat)<-Poly$Name
colnames(tumalo.dat)<-c("tri", "tpi", "slope", "aspect", "trisd", "tpisd", "slopesd", "aspectsd", "elevation")
tumalo.dat<-data.frame(na.omit(tumalo.dat))

Poly$Name
rbind(tumalo.dat, elkhorn.dat, wallowa.dat, hood.dat, drake.dat)

landscape.dat<-rbind(tumalo.dat, elkhorn.dat, wallowa.dat, hood.dat, drake.dat)
write.csv(landscape.dat, "/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/landscapeGIS.csv")
