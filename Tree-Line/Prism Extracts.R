library(raster)
library(sp)
library(rgdal)
library(rgeos)

##Turn points into polygons
library("sp")
library("rgdal")
library("raster")

setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Prism/")
list<-list.files(pattern="*bil.bil")
list
list<-list[c(TRUE,FALSE)]
list
rasterlist<-NULL
results<-NULL

for(i in list){
  temp<-readGDAL(i)
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(i,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAT.", data.frame(strsplit(i,"_"))[5,]))
  results<-rbind(results, extract(temp.raster, Poly))
}

plot(temp.raster)
data.frame(results)
results<-results[1:444,]
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
monthrep<-c(rep(months,37))
monthrep
years<-rep(1981:2017, each = 12)
years
rows<-paste0(monthrep,years)
results.df<-data.frame(results)
rownames(results.df)<-rows
colnames(results.df)<-Poly$Name
results.df$Year<-as.vector(years)
results.df$Month<-as.vector(monthrep)
tmean.df<-results.df
df <- data.frame(matrix(unlist(tmean.df), nrow=444, byrow=F),stringsAsFactors=FALSE)
write.csv(tmean.df, "/Users/Maxwell/Desktop/tmean.csv")
str(YSF)
df[15]<-df[14]
df[18]<-df[20]
df[25]<-df[26]
df<-df[-24]
df<-df[-19]
df<-df[-16]
colnames(df)<-c(as.character(Poly$Name), "Year", "Month")
df[1]<-as.vector(df[1])
str(df)
write.csv(df, "/Users/Maxwell/Desktop/df.csv")
Temp<-read.csv("/Users/Maxwell/Desktop/df.csv")
str(Temp)
Temp<-Temp[-1]
MAT<-Temp[-30] %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
MAT<-data.frame(MAT)
MAT
library(reshape2)
MAT.long<-reshape(MAT, 
                         idvar="Year", ids = "Year",
                         times=names(MAT[-1]), timevar = "Site",
                         varying=list(names(MAT[-1])), v.names="MAT", 
                         direction = "long")

ggplot(MAT.long, aes(y=MAT, x=Year, col=Location))+geom_line()
write.csv(MAT, "/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climate.csv")

###############  MAP  ##################
setwd("/Users/Maxwell/Documents/geospatial/PRISM_ppt_stable_4kmM3_198101_201807_bil/")
list<-list.files(pattern="bil.bil*")
list
list<-list[c(TRUE,FALSE)]
list
rasterlist<-NULL
results<-NULL

for(i in list){
  ppt<-readGDAL(i)
  ppt.raster<-raster(temp)
  assign(paste0("MAP.", data.frame(strsplit(i,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAP.", data.frame(strsplit(i,"_"))[5,]))
  results<-rbind(results, extract(ppt.raster, Poly))
}

plot(temp.raster)
data.frame(results)
results<-results[1:444,]
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
monthrep<-c(rep(months,37))
monthrep
years<-rep(1981:2017, each = 12)
years
rows<-paste0(monthrep,years)
results.df<-data.frame(results)
rownames(results.df)<-rows
colnames(results.df)<-Poly$Name
results.df$Year<-as.vector(years)
results.df$Month<-as.vector(monthrep)
ppt.df<-results.df
df <- data.frame(matrix(unlist(ppt.df), nrow=444, byrow=F),stringsAsFactors=FALSE)
head(results.df)
head(df)
str(YSF)
df[18]<-df[20]
df<-df[-24]
df<-df[-19]
df<-df[-15]
colnames(df)<-c(as.character(Poly$Name), "Year", "Month")
df$ESB<-df$ESF
df$MNF<-df$MNB
head(df)
str(df)
write.csv(df, "/Users/Maxwell/Desktop/df.csv")
ppt<-read.csv("/Users/Maxwell/Desktop/df.csv")
str(ppt)
ppt<-ppt[-1]
MAP<-ppt[-30] %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
MAP<-data.frame(MAP)
MAP
library(reshape2)
MAP.long<-reshape(MAP, 
                  idvar="Year", ids = "Year",
                  times=names(MAP[-1]), timevar = "Site",
                  varying=list(names(MAP[-1])), v.names="MAP", 
                  direction = "long")

ggplot(MAP.long, aes(y=MAP, x=Year, col=Location))+geom_line()
write.csv(MAP, "/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climate.csv")

MAP.long
climate<-merge(MAP.long, MAT.long)
write.csv(climate, "/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climate.csv")
