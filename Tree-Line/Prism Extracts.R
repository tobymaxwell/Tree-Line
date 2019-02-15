library(raster)
library(sp)
library(rgdal)
library(rgeos)
Poly
sites<-Poly

setwd("/Users/Maxwell/Documents/geospatial/PRISM_tmean_stable_4kmM2_198101_201807_bil/")
list<-list.files(pattern="bil.bil*")
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
library(plyr)
results.average<-ddply(df, .(Year), summarise,
                       DNB=mean(DNB),
                       DNF=mean(DNF))
str(df)

results.reshape<-reshape(results.average, 
                         idvar="Year", ids = "Year",
                         times=names(results.avg[-1]), timevar = "Location",
                         varying=list(names(results.avg[-1])), v.names="MAT", 
                         direction = "long")
str(results.avg)
str(results.reshape)
world.mean<-results.reshape


####Tmin####
list<-list.files(path = "/Users/maxwell1/Documents/PRISM/tmin/bils only/")
setwd("/Users/maxwell1/Documents//PRISM/tmin/")
rasterlist<-NULL
results<-NULL
sitecoords<-data.frame(sites)[3:4]
sitecoords

for(i in list){
  temp<-readGDAL(i)
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(i,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAT.", data.frame(strsplit(i,"_"))[5,]))
  results<-rbind(results,extract(temp.raster, sitecoords))
}

data.frame(results)
results<-results[1:432,]
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
monthrep<-c(rep(months,36))
monthrep
years<-(c(rep(1981,12), rep(1982,12),rep(1983,12),rep(1984,12),rep(1985,12),rep(1986,12),rep(1987,12),rep(1988,12),rep(1989,12),rep(1990,12), rep(1991,12), rep(1992,12), rep(1993,12), rep(1994,12), rep(1995,12), rep(1996,12), rep(1997,12), rep(1998,12), rep(1999,12), rep(2000,12),rep(2001,12), rep(2002,12), rep(2003,12), rep(2004,12), rep(2005,12), rep(2006,12), rep(2007,12), rep(2008,12), rep(2009,12), rep(2010,12),rep(2011,12), rep(2012,12), rep(2013,12), rep(2014,12), rep(2015,12), rep(2016,12)))
years
rows<-paste0(monthrep,years)
results.df<-data.frame(results)
rownames(results.df)<-rows
colnames(results.df)<-sites$Name
str(results.df)
results.df$Year<-as.vector(years)
results.df$Month<-as.vector(monthrep)
str(results.df)
results.df$Year<-years
results.df$Month<-monthrep
str(results.df)
tmin.df<-results.df
library(plyr)
results.avg<-ddply(results.df[results.df$Month=="Jan"|results.df$Month=="Feb"|results.df$Month=="Mar",], .(Year),summarise,
      Bieber = mean(Bieber),
      Butte = mean(Butte),
      Clarksburg = mean(Clarksburg),
      Colusa = mean(Colusa),
      Davis = mean (Davis),
      Delta = mean(Delta),
      Delta.Mello = mean(Delta.Mello),
      Delta.Tyler = mean(Delta.Tyler),
      Dunnigan = mean (Dunnigan),
      Etna = mean(Etna),
      Five.Points = mean(Five.Points),
      Fort.Jones = mean(Fort.Jones),
      Fresno = mean(Fresno),
      Glenn = mean(Glenn),
      Glenn.Artois =mean(Glenn.Artois),
      Imperial = mean(Imperial),
      Kern = mean(Kern),
      Kern.Button = mean(Kern.Button),
      Kern.Lake = mean(Kern.Lake),
      Kern.Shaft = mean(Kern.Shaft),
      Kings = mean(Kings),
      Kings.Hansen=mean(Kings.Hansen),
      Kings.McCarthy = mean(Kings.McCarthy),
      Little.Valley = mean(Little.Valley),
      Madera = mean(Madera),
      Merced.Atwater = mean(Merced.Atwater),
      Merced.Dos = mean(Merced.Dos),
      Merced = mean (Merced),
      Montague = mean(Montague),
      Santa.Clara =mean(Santa.Clara),
      SLO = mean(SLO),
      Susanville = mean(Susanville),
      Sutter = mean(Sutter),
      Tehama = mean(Tehama),
      Tulare = mean(Tulare),
      Tulare.Porter = mean(Tulare.Porter),
      Tulelake = mean(Tulelake),
      Yolo = mean(Yolo),
      Delta.Staten = mean(Delta.Staten),
      Kings.Stratford = mean(Kings.Stratford),
      Solano = mean(Solano),
      Stanislaus.Delgran = mean(Stanislaus.Delgran),
      Madera.Chowchilla = mean(Madera.Chowchilla),
      Butte.Dutro = mean(Butte.Dutro),
      SantaB=mean(SantaB))

results.reshape<-reshape(results.avg, 
        idvar="Year", ids = "Year",
        times=names(results.avg[-1]), timevar = "Location",
        varying=list(names(results.avg[-1])), v.names="TMIN", 
        direction = "long")
str(results.avg)
str(results.reshape)
world.tmin<-results.reshape

####TMAX####
list<-list.files(path = "/Users/maxwell1/Documents/PRISM/tmax/bils/")
setwd("/Users/maxwell1/Documents/PRISM/tmax/")
rasterlist<-NULL
results<-NULL
sitecoords<-data.frame(sites)[3:4]
sitecoords

for(i in list){
  temp<-readGDAL(i)
  temp.raster<-raster(temp)
  assign(paste0("MAT.", data.frame(strsplit(i,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("MAT.", data.frame(strsplit(i,"_"))[5,]))
  results<-rbind(results,extract(temp.raster, sitecoords))
}

data.frame(results)
results<-results[1:432,]
months<-(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
monthrep<-c(rep(months,36))
monthrep
years<-(c(rep(1981,12), rep(1982,12),rep(1983,12),rep(1984,12),rep(1985,12),rep(1986,12),rep(1987,12),rep(1988,12),rep(1989,12),rep(1990,12), rep(1991,12), rep(1992,12), rep(1993,12), rep(1994,12), rep(1995,12), rep(1996,12), rep(1997,12), rep(1998,12), rep(1999,12), rep(2000,12),rep(2001,12), rep(2002,12), rep(2003,12), rep(2004,12), rep(2005,12), rep(2006,12), rep(2007,12), rep(2008,12), rep(2009,12), rep(2010,12),rep(2011,12), rep(2012,12), rep(2013,12), rep(2014,12), rep(2015,12), rep(2016,12)))
years
rows<-paste0(monthrep,years)
results.df<-data.frame(results)
rownames(results.df)<-rows
colnames(results.df)<-sites$Name
str(results.df)
results.df$Year<-as.vector(years)
results.df$Month<-as.vector(monthrep)
str(results.df)
results.df$Year<-years
results.df$Month<-monthrep
str(results.df)
tmax.df<-results.df
library(plyr)
results.avg<-ddply(results.df[results.df$Month=="Jan"|results.df$Month=="Feb"|results.df$Month=="Mar",], .(Year),summarise,
                   Bieber = mean(Bieber),
                   Butte = mean(Butte),
                   Clarksburg = mean(Clarksburg),
                   Colusa = mean(Colusa),
                   Davis = mean (Davis),
                   Delta = mean(Delta),
                   Delta.Mello = mean(Delta.Mello),
                   Delta.Tyler = mean(Delta.Tyler),
                   Dunnigan = mean (Dunnigan),
                   Etna = mean(Etna),
                   Five.Points = mean(Five.Points),
                   Fort.Jones = mean(Fort.Jones),
                   Fresno = mean(Fresno),
                   Glenn = mean(Glenn),
                   Glenn.Artois =mean(Glenn.Artois),
                   Imperial = mean(Imperial),
                   Kern = mean(Kern),
                   Kern.Button = mean(Kern.Button),
                   Kern.Lake = mean(Kern.Lake),
                   Kern.Shaft = mean(Kern.Shaft),
                   Kings = mean(Kings),
                   Kings.Hansen=mean(Kings.Hansen),
                   Kings.McCarthy = mean(Kings.McCarthy),
                   Little.Valley = mean(Little.Valley),
                   Madera = mean(Madera),
                   Merced.Atwater = mean(Merced.Atwater),
                   Merced.Dos = mean(Merced.Dos),
                   Merced = mean (Merced),
                   Montague = mean(Montague),
                   Santa.Clara =mean(Santa.Clara),
                   SLO = mean(SLO),
                   Susanville = mean(Susanville),
                   Sutter = mean(Sutter),
                   Tehama = mean(Tehama),
                   Tulare = mean(Tulare),
                   Tulare.Porter = mean(Tulare.Porter),
                   Tulelake = mean(Tulelake),
                   Yolo = mean(Yolo),
                   Delta.Staten = mean(Delta.Staten),
                   Kings.Stratford = mean(Kings.Stratford),
                   Solano = mean(Solano),
                   Stanislaus.Delgran = mean(Stanislaus.Delgran),
                   Madera.Chowchilla = mean(Madera.Chowchilla),
                   Butte.Dutro = mean(Butte.Dutro),
                   SantaB=mean(SantaB))

results.reshape<-reshape(results.avg, 
                         idvar="Year", ids = "Year",
                         times=names(results.avg[-1]), timevar = "Location",
                         varying=list(names(results.avg[-1])), v.names="TMAX", 
                         direction = "long")
str(results.avg)
str(results.reshape)
world.tmax<-results.reshape


#####PPT#####

list<-list.files(path = "/Users/maxwell1/Documents/PRISM/ppt/bils/")
setwd("/Users/maxwell1/Documents/PRISM/ppt/")

rasterlist<-NULL
results<-NULL
sitecoords<-data.frame(sites)[3:4]
sitecoords
for(i in list){
  temp<-readGDAL(i)
  temp.raster<-raster(temp)
  assign(paste0("PPT.", data.frame(strsplit(i,"_"))[5,]), temp.raster)
  rasterlist<-c(rasterlist,paste0("PPT.", data.frame(strsplit(i,"_"))[5,]))
  results<-rbind(results,extract(temp.raster, sitecoords))
}

PPT.df<-data.frame(results[1:432,])
rownames(PPT.df)<-rows
colnames(PPT.df)<-sites$Name
str(PPT.df)
PPT.df$Year<-years
PPT.df$Month<-monthrep
str(PPT.df)
library(plyr)
PPT.avg<-ddply(PPT.df, .(Year),summarise,
                   Bieber = mean(Bieber),
                   Butte = mean(Butte),
                   Clarksburg = mean(Clarksburg),
                   Colusa = mean(Colusa),
                   Davis = mean (Davis),
                   Delta = mean(Delta),
                   Delta.Mello = mean(Delta.Mello),
                   Delta.Tyler = mean(Delta.Tyler),
                   Dunnigan = mean (Dunnigan),
                   Etna = mean(Etna),
                   Five.Points = mean(Five.Points),
                   Fort.Jones = mean(Fort.Jones),
                   Fresno = mean(Fresno),
                   Glenn = mean(Glenn),
                   Glenn.Artois =mean(Glenn.Artois),
                   Imperial = mean(Imperial),
                   Kern = mean(Kern),
                   Kern.Button = mean(Kern.Button),
                   Kern.Lake = mean(Kern.Lake),
                   Kern.Shaft = mean(Kern.Shaft),
                   Kings = mean(Kings),
                   Kings.Hansen=mean(Kings.Hansen),
                   Kings.McCarthy = mean(Kings.McCarthy),
                   Little.Valley = mean(Little.Valley),
                   Madera = mean(Madera),
                   Merced.Atwater = mean(Merced.Atwater),
                   Merced.Dos = mean(Merced.Dos),
                   Merced = mean (Merced),
                   Montague = mean(Montague),
                   Santa.Clara =mean(Santa.Clara),
                   SLO = mean(SLO),
                   Susanville = mean(Susanville),
                   Sutter = mean(Sutter),
                   Tehama = mean(Tehama),
                   Tulare = mean(Tulare),
                   Tulare.Porter = mean(Tulare.Porter),
                   Tulelake = mean(Tulelake),
                   Yolo = mean(Yolo),
                   Delta.Staten = mean(Delta.Staten),
                   Kings.Stratford = mean(Kings.Stratford),
                   Solano = mean(Solano),
                   Stanislaus.Delgran = mean(Stanislaus.Delgran),
                   Madera.Chowchilla = mean(Madera.Chowchilla),
                   Butte.Dutro = mean(Butte.Dutro),
                   SantaB=mean(SantaB))

PPT.reshape<-reshape(PPT.avg, 
                         idvar="Year", ids = "Year",
                         times=names(PPT.avg[-1]), timevar = "Location",
                         varying=list(names(PPT.avg[-1])), v.names="MAP", 
                         direction = "long")
str(PPT.avg)
str(PPT.reshape)
pdsi.bieber
pdsi.bieber<-data.frame(PPT.df$Year[-433], rep(1:12, 36), results.df$Bieber[-433],PPT.df$Bieber[-433])
colnames(pdsi.bieber)<-c("year", "month", "temp", "prec")
result.bieber<-pdsi(awc[1,4], awc[1,3], pdsi.bieber,1982, 2016, mode = "both")


