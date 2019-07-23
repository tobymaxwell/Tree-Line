library(TTR)
library(plyr)
library(ggplot2)

setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Landsat/")

DNB.NDVI<-read.csv("DNB.csv")
DNB.NDVI<-ddply(na.omit(DNB.NDVI[DNB.NDVI$NDVI.DNB>0.1,]), .(Year), summarise,
      NDVI.DNB = max(NDVI.DNB))
plot(NDVI.DNB~Year, DNB.NDVI)
summary(lm(NDVI.DNB~Year, DNB.NDVI))

DNF.NDVI<-read.csv("DNF.csv")
DNF.NDVI<-ddply(na.omit(DNF.NDVI[DNF.NDVI$NDVI.DNF>0.1,]), .(Year), summarise,
                NDVI.DNF = max(NDVI.DNF))
plot(NDVI.DNF~Year, DNF.NDVI)
summary(lm(NDVI.DNF~Year, DNF.NDVI))

NDVI.DN<-merge(DNF.NDVI, DNB.NDVI)
NDVI.DN$Diff.DN<-NDVI.DN$NDVI.DNB/NDVI.DN$NDVI.DNF
plot(Diff.DN~Year, NDVI.DN)
summary(lm(Diff.DN~Year, NDVI.DN))

###DS

DSB.NDVI<-read.csv("DSB.csv")
DSB.NDVI<-ddply(na.omit(DSB.NDVI[DSB.NDVI$NDVI.DSB>0.1,]), .(Year), summarise,
                NDVI.DSB = max(NDVI.DSB))
plot(NDVI.DSB~Year, DSB.NDVI)
summary(lm(NDVI.DSB~Year, DSB.NDVI))

DSF.NDVI<-read.csv("DSF.csv")
DSF.NDVI<-ddply(na.omit(DSF.NDVI[DSF.NDVI$NDVI.DSF>0.1,]), .(Year), summarise,
                NDVI.DSF = max(NDVI.DSF))
plot(NDVI.DSF~Year, DSF.NDVI)
summary(lm(NDVI.DSF~Year, DSF.NDVI))

NDVI.DS<-merge(DSF.NDVI, DSB.NDVI)
NDVI.DS$Diff.DS<-NDVI.DS$NDVI.DSB/NDVI.DS$NDVI.DSF
plot(Diff.DS~Year, NDVI.DS[NDVI.DS$Year>1984,])
summary(lm(Diff.DS~Year, NDVI.DS[NDVI.DS$Year>1984,]))

Drake<-merge(NDVI.DS,NDVI.DN)
ggplot(Drake[Drake$Year>1984,], aes(y=Diff.DS, x=Year))+
  geom_line()+
  geom_line(aes(y=Diff.DN, x=Year), color="forest green")+
  theme_bw()+ylab("NDVI Treeline/NDVI Forest")+
  theme(legend.position = 'none')

ggplot(Drake, aes(y=NDVI.DSF, x=Year))+
  geom_line()+
  geom_line(aes(y=NDVI.DSB, x=Year), color="black", linetype="dashed")+
  geom_line(aes(y=NDVI.DNF, x=Year), color="forest green")+
  geom_line(aes(y=NDVI.DNB, x=Year), color="forest green", linetype="dashed")+
  theme_bw()+ylab("NDVI")+
  theme(legend.position = 'none')

Drake
####Powder NDVI TS

PNB.NDVI<-read.csv("PNB.csv")
PNB.NDVI<-ddply(na.omit(PNB.NDVI[PNB.NDVI$NDVI.PNB>0.1,]), .(Year), summarise,
                NDVI.PNB = max(NDVI.PNB))
plot(NDVI.PNB~Year, PNB.NDVI)
summary(lm(NDVI.PNB~Year, PNB.NDVI))

PNF.NDVI<-read.csv("PNF.csv")
PNF.NDVI<-ddply(na.omit(PNF.NDVI[PNF.NDVI$NDVI.PNF>0.1,]), .(Year), summarise,
                NDVI.PNF = max(NDVI.PNF))
plot(NDVI.PNF~Year, PNF.NDVI)
summary(lm(NDVI.PNF~Year, PNF.NDVI))

NDVI.PN<-merge(PNF.NDVI, PNB.NDVI)
NDVI.PN$Diff.PN<-NDVI.PN$NDVI.PNB/NDVI.PN$NDVI.PNF
plot(Diff.PN~Year, NDVI.PN)
summary(lm(Diff.PN~Year, NDVI.PN))

###PS

PSB.NDVI<-read.csv("PSB.csv")
PSB.NDVI<-ddply(na.omit(PSB.NDVI[PSB.NDVI$NDVI.PSB>0.1,]), .(Year), summarise,
                NDVI.PSB = max(NDVI.PSB))
plot(NDVI.PSB~Year, PSB.NDVI)
summary(lm(NDVI.PSB~Year, PSB.NDVI))

PSF.NDVI<-read.csv("PSF.csv")
PSF.NDVI<-ddply(na.omit(PSF.NDVI[PSF.NDVI$NDVI.PSF>0.1,]), .(Year), summarise,
                NDVI.PSF = max(NDVI.PSF))
plot(NDVI.PSF~Year, PSF.NDVI)
summary(lm(NDVI.PSF~Year, PSF.NDVI))

NDVI.PS<-merge(PSF.NDVI, PSB.NDVI)
NDVI.PS$Diff.PS<-NDVI.PS$NDVI.PSB/NDVI.PS$NDVI.PSF
plot(Diff.PS~Year, NDVI.PS)
summary(lm(Diff.PS~Year, NDVI.PS))

Picea<-merge(NDVI.PS,NDVI.PN)
ggplot(Picea, aes(y=Diff.PS, x=Year))+
  geom_line()+
  geom_line(aes(y=Diff.PN, x=Year), color="forest green")+
  theme_bw()+ylab("NDVI Treeline/NDVI Forest")+
  theme(legend.position = 'none')

ggplot(Picea, aes(y=NDVI.PSF, x=Year))+
  geom_line()+
  geom_line(aes(y=NDVI.PSB, x=Year), color="black", linetype="dashed")+
  geom_line(aes(y=NDVI.PNF, x=Year), color="forest green")+
  geom_line(aes(y=NDVI.PNB, x=Year), color="forest green", linetype="dashed")+
  theme_bw()+ylab("NDVI")+
  theme(legend.position = 'none')

####Twin NDVI TS

TNB.NDVI<-read.csv("TNB.csv")
TNB.NDVI<-ddply(na.omit(TNB.NDVI[TNB.NDVI$NDVI.TNB>0.1&TNB.NDVI$Year!=2008,]), .(Year), summarise,
                NDVI.TNB = max(NDVI.TNB))
plot(NDVI.TNB~Year, TNB.NDVI)
summary(lm(NDVI.TNB~Year, TNB.NDVI))

TNF.NDVI<-read.csv("TNF.csv")
TNF.NDVI<-ddply(na.omit(TNF.NDVI[TNF.NDVI$NDVI.TNF>0.1&TNF.NDVI$Year!=2008,]), .(Year), summarise,
                NDVI.TNF = max(NDVI.TNF))
plot(NDVI.TNF~Year, TNF.NDVI)
summary(lm(NDVI.TNF~Year, TNF.NDVI))

NDVI.TN<-merge(TNF.NDVI, TNB.NDVI)
NDVI.TN$Diff.TN<-NDVI.TN$NDVI.TNB/NDVI.TN$NDVI.TNF
plot(Diff.TN~Year, NDVI.TN)
summary(lm(Diff.TN~Year, NDVI.TN))

###TS

TSB.NDVI<-read.csv("TSB.csv")
TSB.NDVI<-ddply(na.omit(TSB.NDVI[TSB.NDVI$NDVI.TSB>0.1,]), .(Year), summarise,
                NDVI.TSB = max(NDVI.TSB))
plot(NDVI.TSB~Year, TSB.NDVI)
summary(lm(NDVI.TSB~Year, TSB.NDVI))

TSF.NDVI<-read.csv("TSF.csv")
TSF.NDVI<-ddply(na.omit(TSF.NDVI[TSF.NDVI$NDVI.TSF>0.1,]), .(Year), summarise,
                NDVI.TSF = max(NDVI.TSF))
plot(NDVI.TSF~Year, TSF.NDVI)
summary(lm(NDVI.TSF~Year, TSF.NDVI))

NDVI.TS<-merge(TSF.NDVI, TSB.NDVI)
NDVI.TS$Diff.TS<-NDVI.TS$NDVI.TSB/NDVI.TS$NDVI.TSF
plot(Diff.TS~Year, NDVI.TS[NDVI.TS$Year>1984,])
summary(lm(Diff.TS~Year, NDVI.TS[NDVI.TS$Year>1984,]))

Twin<-merge(NDVI.TS,NDVI.TN)
ggplot(Twin, aes(y=Diff.TS, x=Year))+
  geom_line()+
  geom_line(aes(y=Diff.TN, x=Year), color="forest green")+
  theme_bw()+ylab("NDVI Treeline/NDVI Forest")+
  theme(legend.position = 'none')

ggplot(Twin, aes(y=NDVI.TSF, x=Year))+
  geom_line()+
  geom_line(aes(y=NDVI.TSB, x=Year), color="black", linetype="dashed")+
  geom_line(aes(y=NDVI.TNF, x=Year), color="forest green")+
  geom_line(aes(y=NDVI.TNB, x=Year), color="forest green", linetype="dashed")+
  theme_bw()+ylab("NDVI")+
  theme(legend.position = 'none')

####Hood NDVI TS

HNB.NDVI<-read.csv("HNB.csv")
HNB.NDVI<-ddply(na.omit(HNB.NDVI[HNB.NDVI$NDVI.HNB>0.1,]), .(Year), summarise,
                NDVI.HNB = max(NDVI.HNB))
plot(NDVI.HNB~Year, HNB.NDVI)
summary(lm(NDVI.HNB~Year, HNB.NDVI))

HNF.NDVI<-read.csv("HNF.csv")
HNF.NDVI<-ddply(na.omit(HNF.NDVI[HNF.NDVI$NDVI.HNF>0.1,]), .(Year), summarise,
                NDVI.HNF = max(NDVI.HNF))
plot(NDVI.HNF~Year, HNF.NDVI)
summary(lm(NDVI.HNF~Year, HNF.NDVI))

NDVI.HN<-merge(HNF.NDVI, HNB.NDVI)
NDVI.HN$Diff.HN<-NDVI.HN$NDVI.HNB/NDVI.HN$NDVI.HNF
plot(Diff.HN~Year, NDVI.HN)
summary(lm(Diff.HN~Year, NDVI.HN))

###HS

HSB.NDVI<-read.csv("HSB.csv")
HSB.NDVI<-ddply(na.omit(HSB.NDVI[HSB.NDVI$NDVI.HSB>0.1,]), .(Year), summarise,
                NDVI.HSB = max(NDVI.HSB))
plot(NDVI.HSB~Year, HSB.NDVI)
summary(lm(NDVI.HSB~Year, HSB.NDVI))

HSF.NDVI<-read.csv("HSF.csv")
HSF.NDVI<-ddply(na.omit(HSF.NDVI[HSF.NDVI$NDVI.HSF>0.1,]), .(Year), summarise,
                NDVI.HSF = max(NDVI.HSF))
plot(NDVI.HSF~Year, HSF.NDVI)
summary(lm(NDVI.HSF~Year, HSF.NDVI))

NDVI.HS<-merge(HSF.NDVI, HSB.NDVI)
NDVI.HS$Diff.HS<-NDVI.HS$NDVI.HSB/NDVI.HS$NDVI.HSF
plot(Diff.HS~Year, NDVI.HS[NDVI.HS$Year>1984,])
summary(lm(Diff.HS~Year, NDVI.HS[NDVI.HS$Year>1984,]))

Hood<-merge(NDVI.HS,NDVI.HN)
ggplot(Hood, aes(y=Diff.HS, x=Year))+
  geom_line()+
  geom_line(aes(y=Diff.HN, x=Year), color="forest green")+
  theme_bw()+ylab("NDVI Treeline/NDVI Forest")+
  theme(legend.position = 'none')


ggplot(Hood, aes(y=NDVI.HSF, x=Year))+
  geom_line()+
  geom_line(aes(y=NDVI.HSB, x=Year), color="black", linetype="dashed")+
  geom_line(aes(y=NDVI.HNF, x=Year), color="forest green")+
  geom_line(aes(y=NDVI.HNB, x=Year), color="forest green", linetype="dashed")+
  theme_bw()+ylab("NDVI")+
  theme(legend.position = 'none')

ggplot(Hood, aes(y=NDVI.HSB, x=Year))+
  geom_line(color="Black", linetype="dashed")+
  geom_line(aes(y=NDVI.HNB, x=Year), color="Forest Green", linetype="dashed")+
  theme_bw()+ylab("NDVI Treeline/NDVI Forest")+
  theme(legend.position = 'none')

Drake
colnames(Drake)<-c("Year", "DSF", "DSB", "DS", "DNF", "DNB", "DN")
Drake
Drake.long<-reshape(Drake, 
                   idvar="Year", ids = "Year",
                   times=names(Drake[-1]), timevar = "Site",
                   varying=list(names(Drake[-1])), v.names="NDVI", 
                   direction = "long")


Hood
colnames(Hood)<-c("Year", "HSF", "HSB", "HS", "HNF", "HNB", "HN")
Hood.long<-reshape(Hood, 
                    idvar="Year", ids = "Year",
                    times=names(Hood[-1]), timevar = "Site",
                    varying=list(names(Hood[-1])), v.names="NDVI", 
                    direction = "long")
Hood.long

Twin
colnames(Twin)<-c("Year", "TSF", "TSB", "TS", "TNF", "TNB", "TN")
Twin.long<-reshape(Twin, 
                   idvar="Year", ids = "Year",
                   times=names(Twin[-1]), timevar = "Site",
                   varying=list(names(Twin[-1])), v.names="NDVI", 
                   direction = "long")
Twin.long

Picea
colnames(Picea)<-c("Year", "PSF", "PSB", "PS", "PNF", "PNB", "PN")
Picea.long<-reshape(Picea, 
                   idvar="Year", ids = "Year",
                   times=names(Picea[-1]), timevar = "Site",
                   varying=list(names(Picea[-1])), v.names="NDVI", 
                   direction = "long")
Picea.long

NDVI.long<-rbind(Picea.long, Twin.long, Hood.long, Drake.long)

NDVI.clim<-merge(climate, NDVI.long)
NDVI.clim<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/NDVI.Clim.csv")
ggplot(NDVI.clim, aes(y=scale(NDVI), x=Year))+geom_line()+facet_wrap(~Site)+
  geom_line(data=NDVI.clim, aes(y=scale(MAT), x=Year), color="Red")+
  geom_line(data=NDVI.clim, aes(y=scale(MAP), x=Year), color="Blue")+
  theme_classic()

ddply(NDVI.clim, .(Site), summarise,
      MAP=mean(MAP),
      MAT=mean(MAT)/10)
