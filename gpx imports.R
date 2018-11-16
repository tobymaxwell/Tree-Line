#kml to csv
library(rgdal)
library(raster)
library(maptools)
library(tools)
setwd("/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/GPX files/")
files<-list.files(pattern=".gpx")
files
k<-NA
m<-NA
l<-NA
for(i in files){
  k<-readOGR(i)
  l<-data.frame(subset(k, select=c('ele', "name")))
  m<-merge(l,m, all=T)
}
m
tail(m,100)
write.csv