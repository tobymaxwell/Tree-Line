#Methods
library(raster)
library(rgeos)
library(rgdal)
library(maptools)
library(sp)
library(latticeExtra)
library(GISTools)
NTC<-raster("/Users/Maxwell/Documents/MRT/MOD44B_2016-03-05.Percent_Tree_Cover.tif")
oregon<-readOGR("/Users/Maxwell/Documents/geospatial/orcnty2015/orcntyline.shp")

NTC.crs<- crs(NTC)
OR<-spTransform(oregon, crs(NTC))
plot(OR)

NTC.mask<-crop(NTC, OR)
plot(NTC.mask)

extent(oregon)
bio.rasters<-getData(name = 'worldclim', var = 'bio', res=2.5)
plot(bio.rasters[[12]])
class(NTC.crs)

or.ppt<-projectRaster(bio.rasters[[12]], crs=NTC.crs)
plot(or.ppt)

ppt.crop<-crop(or.ppt, OR)
plot(ppt.crop)
plot(OR, add=T)

dem<-getData(name = 'alt', country="USA")
or.dem<-dem$`/Users/Maxwell/Documents/GitHub/Tree-Line/Tree-Line/USA1_msk_alt.grd`
plot(or.dem)

or.dem<-projectRaster(or.dem, crs=NTC.crs)
plot(or.dem)

dem.crop<-crop(or.dem, OR)
plot(dem.crop)
plot(OR, add=T)

slope<-terrain(dem.crop)
aspect<-terrain(dem.crop, opt='aspect')
tpi<-terrain(dem.crop, opt='TPI')
tri<-terrain(dem.crop, opt='TRI')
plot(aspect)
plot(slope)
plot(tri)
plot(tpi)
library(dynatopmodel)
res(dem.crop)<-c(900,900)
plot(dem.crop)
upslope.area(dem.crop)

par(mfrow=c(1,4))
plot(ppt.crop)
plot(NTC.mask)
plot(slope)
plot("/Users/Maxwell/Documents/geospatial/")

writeRaster(ppt.crop,"/Users/Maxwell/Desktop/ppt.crop.tif")
writeRaster(slope,"/Users/Maxwell/Desktop/slope.tif")
writeRaster(aspect,"/Users/Maxwell/Desktop/aspect.tif")
