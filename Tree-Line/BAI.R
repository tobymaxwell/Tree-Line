#install.packages("dplR")
library("dplR")
library(plyr)
# Bring in master file of RWs
Master=read.csv("/Users/tobymaxwell//OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Hall.csv", header=TRUE)
str(Master)
#Basic Ring Width Index calculations
n<-min(na.omit(2017-length(Master$HSFAbLa1)*.5))
Master$Year<-as.numeric(rep(2017:(n+1), each=2))
library(dplyr)
Master.T<-Master %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
Master.T<-data.frame(Master.T)
str(Master.T)
revdf<-function(df)df[seq(dim(df)[1],1),]
#Master.rev<-revdf(Master.T)

RWI<-detrend(Master.T,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rownames(RWI)<-Master.T$Year
spag.plot(RWI[15:22], zfac=.3)
spag.plot(RWI, zfac=.3)
rwl.report(RWI[15:22])

##### Calculate Master RWI using the C-Method (recomended method, see Biondi and Qeadan, 2008) 
library(graphics)
library(utils)
str(Master.T)
HSFAbLa<-Master.T[1:6]
#HSFAbLa<-revdf(HSFAbLa)
rownames(HSFAbLa)<-HSFAbLa$Year
PO<-NULL
PO$series<-c("T1", "T2", "T3", "T4", "T5")
PO<-data.frame(PO)
cols<-colnames(HSFAbLa[-1])
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HSFAbLa, i)[1])))
lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HSFAbLa<-HSFAbLa[-1]
colnames(HSFAbLa)<-c("T1", "T2", "T3", "T4", "T5")
#write.csv(PO, "/Users/tobymaxwell/Desktop/PO.csv")
write.rwl(HSFAbLa,"/Users/tobymaxwell/Desktop/HSFAbLa2.csv")
HSFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/HSFAbLa2.csv", header=TRUE)
rwl.report(HSFAbLa)
RWI_C<-cms(HSFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index", line = 3)
title(sub = "RWI - AR Model", line = -8.5, font.sub = 2)

####HNF AbLa #####################
str(Master.T)
HNFAbLa<-Master.T[31:35]
RWI<-detrend(HNFAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
#HSFAbLa<-revdf(HSFAbLa)
rownames(HNFAbLa)<-1535:2017
PO<-NULL
PO$series<-c("T1", "T2", "T3", "T4", "T5")
PO<-data.frame(PO)
cols<-colnames(HNFAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HNFAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HNFAbLa
colnames(HNFAbLa)<-c("T1", "T2", "T3", "T4", "T5")
write.rwl(HNFAbLa,"/Users/tobymaxwell/Desktop/HNFAbLa.csv")
HNFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/HNFAbLa.csv", header=TRUE)

rwl.report(HNFAbLa)
RWI_C<-cms(HNFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index", line = 3)
title(sub = "RWI - AR Model", line = -28.5, font.sub = 2)


####HSF TSME #####################
str(Master.T)
HSFTsMe<-Master.T[7:15]
RWI<-detrend(HSFTsMe,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
#HSFAbLa<-revdf(HSFAbLa)
rownames(HSFTsMe)<-1535:2017
PO<-NULL
PO$series<-c("T1", "T2", "T3", "T4", "T5", "T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(HSFTsMe)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HSFTsMe, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HSFTsMe
colnames(HSFTsMe)<-c("T1", "T2", "T3", "T4", "T5", "T7", "T8", "T9", "T10")
write.rwl(HSFTsMe,"/Users/tobymaxwell/Desktop/HSFTsMe.csv")
HSFTsMe<-read.rwl("/Users/tobymaxwell/Desktop/HSFTsMe.csv", header=TRUE)

rwl.report(HSFTsMe)
RWI_C<-cms(HSFTsMe, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron.HSFTsMe<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron.HSFTsMe
crn.plot(MeanChron.HSFTsMe, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index", line = 3)
title(sub = "RWI - AR Model", line = -28.5, font.sub = 2)

####HSB TSME #####################
str(Master.T)
HSBTsMe<-Master.T[16:23]
RWI<-detrend(HSBTsMe,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
#HSFAbLa<-revdf(HSFAbLa)
rownames(HSBTsMe)<-1535:2017
PO<-NULL
PO$series<-c("T1", "T2", "T3", "T4", "T5", "T7", "T8", "T9")
PO<-data.frame(PO)
cols<-colnames(HSBTsMe)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HSBTsMe, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HSBTsMe
colnames(HSBTsMe)<-c("T1", "T2", "T3", "T4", "T5", "T7", "T8", "T9")
write.rwl(HSBTsMe,"/Users/tobymaxwell/Desktop/HSBTsMe.csv")
HSBTsMe<-read.rwl("/Users/tobymaxwell/Desktop/HSBTsMe.csv", header=TRUE)

rwl.report(HSBTsMe)
RWI_C<-cms(HSBTsMe, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron.HSBTsMe<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron.HSBTsMe
crn.plot(MeanChron.HSBTsMe, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index - Hood South Border", line = 3)
title(sub = "RWI - AR Model", line = -28.5, font.sub = 2)

########HSB TsMe Normalized Plot#######

MeanChron.HSBTsMe
MeanChron.HSFTsMe
HSmean<-NULL
HSmean$TsMe<-MeanChron.HSBTsMe[1:50,1]/MeanChron.HSFTsMe[1:50,1]
HSmean$Year<-2017:(2017-49)
HSmean<-data.frame(HSmean)
ggplot(HSmean, aes(TsMe, x=Year))+geom_line()
summary(lm(TsMe~Year, HSmean))
####HNF TSME #####################
str(Master.T)
HNFTsMe<-Master.T[24:30]
RWI<-detrend(HNFTsMe,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
#HSFAbLa<-revdf(HSFAbLa)
rownames(HNFTsMe)<-1535:2017
PO<-NULL
PO$series<-c("T2", "T4", "T5", "T6", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(HNFTsMe)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HNFTsMe, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HNFTsMe
colnames(HNFTsMe)<-c("T2", "T4", "T5", "T6", "T8", "T9", "T10")
write.rwl(HNFTsMe,"/Users/tobymaxwell/Desktop/HNFTsMe.csv")
HNFTsMe<-read.rwl("/Users/tobymaxwell/Desktop/HNFTsMe.csv", header=TRUE)

rwl.report(HNFTsMe)
RWI_C<-cms(HNFTsMe, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index", line = 3)
title(sub = "RWI - AR Model", line=-8.5, font.sub = 2)

####HNB TSME #####################
str(Master.T)
HNBTsMe<-Master.T[36:45]
RWI<-detrend(HNBTsMe,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
rownames(HNBTsMe)<-1535:2017
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6","T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(HNBTsMe)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HNBTsMe, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
HNBTsMe
colnames(HNBTsMe)<-c("T1", "T2","T3", "T4", "T5", "T6","T7", "T8", "T9", "T10")
write.rwl(HNBTsMe,"/Users/tobymaxwell/Desktop/HNBTsMe.csv")
HNBTsMe<-read.rwl("/Users/tobymaxwell/Desktop/HNBTsMe.csv", header=TRUE)

rwl.report(HNBTsMe)
RWI_C<-cms(HNBTsMe, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -8.5, font.sub = 2)

####TNB AbLa #####################
TNBAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/TNBAbLa.csv")
str(TNBAbLa)
n<-(length(TNBAbLa[,1])/2)-1
TNBAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
TNBAbLa<-TNBAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
TNBAbLa<-data.frame(TNBAbLa)
str(TNBAbLa)
RWI<-detrend(TNBAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
spag.plot(RWI, zfac = .3)
rownames(TNBAbLa)<-TNBAbLa$Year
TNBAbLa<-TNBAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6","T7", "T8", "T9", "T10", "T11")
PO<-data.frame(PO)
cols<-colnames(TNBAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(TNBAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
TNBAbLa
colnames(TNBAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6","T7", "T8", "T9", "T10", "T11")
write.rwl(TNBAbLa,"/Users/tobymaxwell/Desktop/TNBAbLa.csv")
TNBAbLa<-read.rwl("/Users/tobymaxwell/Desktop/TNBAbLa.csv", header=TRUE)

rwl.report(TNBAbLa)
RWI_C<-cms(TNBAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -8.5, font.sub = 2)

ggplot(MeanChron, aes(y=scale(IZTstd), x=1616:2018, lwd=3))+geom_line()+
  geom_line(data=Twin, aes(y=scale(NDVI.TNB), x=Year), color="Green", lwd=3)+
  #geom_line(data=NDVI.clim[NDVI.clim$Site=="TNB",], aes(y=scale(MAT), x=Year), lwd=2, col="Red")+
  #geom_line(data=NDVI.clim[NDVI.clim$Site=="TNB",], aes(y=scale(MAP), x=Year), lwd=2, col="Blue")+
  xlim(1970,2018)+
  ylab("Standardized NDVI (Green), RWI (Black)")+
  xlab("")


###########   TNFAbLA ##########
TNFAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/TNFAbLa.csv")
str(TNFAbLa)
n<-(length(TNFAbLa$TNFAbLa1)/2)-1
TNFAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
TNFAbLa<-TNFAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
TNFAbLa<-data.frame(TNFAbLa)
str(TNFAbLa)
RWI<-detrend(TNFAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
spag.plot(RWI, zfac = .3)
rownames(TNFAbLa)<-TNFAbLa$Year
TNFAbLa<-TNFAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(TNFAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(TNFAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
TNFAbLa
colnames(TNFAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
write.rwl(TNFAbLa,"/Users/tobymaxwell/Desktop/TNFAbLa.csv")
TNFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/TNFAbLa.csv", header=TRUE)

rwl.report(TNFAbLa)
rwl.stats(TNFAbLa)
RWI_C<-cms(TNFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -28.5, font.sub = 2)

###########   TSBAbLA ##########
TSBAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/TSBAbLa.csv")
str(TSBAbLa)

n<-(length(TSBAbLa$TSBAbLa1)/2)-1
TSBAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))

library(dplyr)
TSBAbLa<-TSBAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
TSBAbLa<-data.frame(TSBAbLa)
str(TSBAbLa)
RWI<-detrend(TSBAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rwi.stats(RWI)
spag.plot(RWI, zfac = .3)
rownames(TSBAbLa)<-TSBAbLa$Year
TSBAbLa<-TSBAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(TSBAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(TSBAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
TSBAbLa
colnames(TSBAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
write.rwl(TSBAbLa,"/Users/tobymaxwell/Desktop/TSBAbLa.csv")
TSBAbLa<-read.rwl("/Users/tobymaxwell/Desktop/TSBAbLa.csv", header=TRUE)

rwl.report(TSBAbLa)
RWI_C<-cms(TSBAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=10, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -8.5, font.sub = 2)

ggplot(MeanChron, aes(y=scale(IZTstd), x=1842:2018))+geom_line(lwd=3)+
 geom_line(data=Twin, aes(y=scale(NDVI.TSB), x=Year), color="Green", lwd=3)+
 geom_line(data=NDVI.clim[NDVI.clim$Site=="TSB",], aes(y=scale(MAT), x=Year), lwd=2, col="Red")+
 geom_line(data=NDVI.clim[NDVI.clim$Site=="TSB",], aes(y=scale(MAP), x=Year), lwd=2, col="Blue")+
  xlim(1970,2018)+
  ylab("Standardized NDVI (Green), RWI (Black)")+
  xlab("")


###########   TSFAbLA ##########
TSFAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/TSFAbLa.csv")
str(TSFAbLa)

n<-(length(TSFAbLa$TSFAbLa1)/2)-1
TSFAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))

library(dplyr)
TSFAbLa<-TSFAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
TSFAbLa<-data.frame(TSFAbLa)
str(TSFAbLa)
RWI<-detrend(TSFAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
spag.plot(RWI, zfac = .3)
rownames(TSFAbLa)<-TSFAbLa$Year
TSFAbLa<-TSFAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(TSFAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(TSFAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
TSFAbLa
colnames(TSFAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
write.rwl(TSFAbLa,"/Users/tobymaxwell/Desktop/TSFAbLa.csv")
TSFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/TSFAbLa.csv", header=TRUE)

rwl.report(TSFAbLa)
RWI_C<-cms(TSFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -9.5, font.sub = 2)

###########   PSFAbLA ##########
PSFAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/PSFAbLa.csv")
str(PSFAbLa)

n<-(length(PSFAbLa$PSFAbLa1)/2)-1
PSFAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
PSFAbLa<-PSFAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
PSFAbLa<-data.frame(PSFAbLa)
str(PSFAbLa)
RWI<-detrend(PSFAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
spag.plot(RWI, zfac = .3)
rownames(PSFAbLa)<-PSFAbLa$Year
PSFAbLa<-PSFAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5")
PO<-data.frame(PO)
cols<-colnames(PSFAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(PSFAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
PSFAbLa
colnames(PSFAbLa)<-c("T1", "T2","T3", "T4", "T5")
write.rwl(PSFAbLa,"/Users/tobymaxwell/Desktop/PSFAbLa.csv")
PSFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/PSFAbLa.csv", header=TRUE)

rwl.report(PSFAbLa)
RWI_C<-cms(PSFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Picea South Forest", line = 3)
title(sub = "RWI - AR Model", line = -9.5, font.sub = 2)

###########   PSBAbLA ##########
PSBAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/PSBAbLa.csv")
str(PSBAbLa)
n<-(length(PSBAbLa$PSBAbLa1)/2)-1
PSBAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
PSBAbLa<-PSBAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
PSBAbLa<-data.frame(PSBAbLa)
str(PSBAbLa)
RWI<-detrend(PSBAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
spag.plot(RWI, zfac = .3)
rownames(PSBAbLa)<-PSBAbLa$Year
PSBAbLa<-PSBAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5")
PO<-data.frame(PO)
cols<-colnames(PSBAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(PSBAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
PSBAbLa
colnames(PSBAbLa)<-c("T1", "T2","T3", "T4", "T5")
write.rwl(PSBAbLa,"/Users/tobymaxwell/Desktop/PSBAbLa.csv")
PSBAbLa<-read.rwl("/Users/tobymaxwell/Desktop/PSBAbLa.csv", header=TRUE)

rwl.report(PSBAbLa)
RWI_C<-cms(PSBAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Hood North Border", line = 3)
title(sub = "RWI - AR Model", line = -28.5, font.sub = 2)

###########   PNFAbLA ##########
PNFAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/PNFAbLa.csv")
str(PNFAbLa)
n<-(length(PNFAbLa$PNFAbLa1)/2)-1
PNFAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
PNFAbLa<-PNFAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
PNFAbLa<-data.frame(PNFAbLa)
str(PNFAbLa)
RWI<-detrend(PNFAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
spag.plot(RWI, zfac = .3)
rownames(PNFAbLa)<-PNFAbLa$Year
PNFAbLa<-PNFAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
PO<-data.frame(PO)
cols<-colnames(PNFAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(PNFAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
PNFAbLa
colnames(PNFAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
write.rwl(PNFAbLa,"/Users/tobymaxwell/Desktop/PNFAbLa.csv")
PNFAbLa<-read.rwl("/Users/tobymaxwell/Desktop/PNFAbLa.csv", header=TRUE)

rwl.report(PNFAbLa)
RWI_C<-cms(PNFAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Picea North Forest", line = 3)
title(sub = "RWI - AR Model", line = -9.5, font.sub = 2)

###########   PNBAbLA ##########
PNBAbLa<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/PNBAbLa.csv")
str(PNBAbLa)
n<-(length(PNBAbLa$PNBAbLa1)/2)-1
PNBAbLa$Year<-as.numeric(rep(2018:(2018-n), each=2))
library(dplyr)
PNBAbLa<-PNBAbLa %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
PNBAbLa<-data.frame(PNBAbLa)
str(PNBAbLa)
RWI<-detrend(PNBAbLa,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
spag.plot(RWI, zfac = .3)
rownames(PNBAbLa)<-PNBAbLa$Year
PNBAbLa<-PNBAbLa[-1]
PO<-NULL
PO$series<-c("T1", "T2","T3", "T4", "T5", "T6")
PO<-data.frame(PO)
cols<-colnames(PNBAbLa)
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(PNBAbLa, i)[1])))
  lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
PNBAbLa
colnames(PNBAbLa)<-c("T1", "T2","T3", "T4", "T5", "T6")
write.rwl(PNBAbLa,"/Users/tobymaxwell/Desktop/PNBAbLa.csv")
PNBAbLa<-read.rwl("/Users/tobymaxwell/Desktop/PNBAbLa.csv", header=TRUE)

rwl.report(PNBAbLa)
RWI_C<-cms(PNBAbLa, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Picea North Border", line = 3)
title(sub = "RWI - AR Model", line = -9.5, font.sub = 2)


ggplot(MeanChron, aes(y=scale(IZTstd), x=1842:2018))+geom_line(lwd=3)+
  geom_line(data=Twin, aes(y=scale(NDVI.TSF), x=Year), color="Green", lwd=3)+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TSF",], aes(y=scale(MAT), x=Year), lwd=2, col="Red")+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TSF",], aes(y=scale(MAP), x=Year), lwd=2, col="Blue")+
  xlim(1970,2018)+
  ylab("Standardized NDVI (Green), RWI (Black)")+
  xlab("")


### Basal Area Increment calculation. DBS's specified by RW length ####
str(PNBAbLa)
PNB_rwl<-data.frame(PNBAbLa)
#PNBAbLa<-bai.out(rwl=PNBAbLa)
PNBAbLa
PNB_rwl<-na.omit(data.frame(apply(PNBAbLa, 1, FUN="mean"), apply(PNBAbLa, 1, FUN="sd")/5))
PNB_rwl$Year<-as.numeric(rownames(PNB_rwl))
colnames(PNB_rwl)<-c("rwlB", "sdB", "Year")

PNF_rwl<-data.frame(PNFAbLa)
#PNFAbLa<-bai.out(rwl=PNFAbLa)
PNF_rwl<-na.omit(data.frame(apply(PNFAbLa, 1, FUN="mean"), apply(PNFAbLa, 1, FUN="sd")/10))
PNF_rwl$Year<-as.numeric(rownames(PNF_rwl))
colnames(PNF_rwl)<-c("rwlF", "sdF", "Year")
str(PNF_rwl)
PN<-merge(PNF_rwl[PNF_rwl$Year>1983,], PNB_rwl[PNB_rwl$Year>1983,])

str(PSBAbLa)
PSB_rwl<-data.frame(PSBAbLa)
#PSBAbLa<-bai.out(rwl=PSBAbLa)
PSBAbLa
PSB_rwl<-na.omit(data.frame(apply(PSBAbLa[-2], 1, FUN="mean"), apply(PSBAbLa[-2], 1, FUN="sd")/5))
PSB_rwl$Year<-as.numeric(rownames(PSB_rwl))
colnames(PSB_rwl)<-c("rwlB", "sdB", "Year")

PSF_rwl<-data.frame(PSFAbLa)
#PSFAbLa<-bai.out(rwl=PSFAbLa)
PSF_rwl<-na.omit(data.frame(apply(PSFAbLa, 1, FUN="mean"), apply(PSFAbLa, 1, FUN="sd")/4))
PSF_rwl$Year<-as.numeric(rownames(PSF_rwl))
colnames(PSF_rwl)<-c("rwlF", "sdF", "Year")
str(PSF_rwl)
PS<-merge(PSF_rwl[PSF_rwl$Year>1983,], PSB_rwl[PSB_rwl$Year>1983,])

####HOOD combined sets North######
str(HNBTsMe)
HNB_rwl<-data.frame(HNBTsMe)
#HNBTsMe<-bai.out(rwl=HNBTsMe)
HNBTsMe
HNB_rwl<-na.omit(data.frame(apply(HNBTsMe, 1, FUN="mean"), apply(HNBTsMe, 1, FUN="sd")/10))
HNB_rwl$Year<-as.numeric(rownames(HNB_rwl))
colnames(HNB_rwl)<-c("rwlB", "sdB", "Year")

HNF_rwl<-data.frame(HNFTsMe)
#HNFTsMe<-bai.out(rwl=HNFTsMe)
HNF_rwl<-na.omit(data.frame(apply(HNFTsMe, 1, FUN="mean"), apply(HNFTsMe, 1, FUN="sd")/7))
HNF_rwl$Year<-as.numeric(rownames(HNF_rwl))
colnames(HNF_rwl)<-c("rwlF", "sdF", "Year")
str(HNF_rwl)

HN<-merge(HNF_rwl, HNB_rwl)
HN<-merge(HNF_rwl, HNB_rwl)
####HOOD combined sets South######
str(HSBTsMe)
HSB_rwl<-data.frame(HSBTsMe)
#HSBTsMe<-bai.out(rwl=HSBTsMe)
HSBTsMe
tail(HSBTsMe,150)
HSB_rwl<-na.omit(data.frame(apply(HSBTsMe, 1, FUN="mean"), apply(HSBTsMe, 1, FUN="sd")/8))
HSB_rwl$Year<-as.numeric(rownames(HSB_rwl))
colnames(HSB_rwl)<-c("rwlB", "sdB", "Year")
HSB_rwl2<-na.omit(data.frame(apply(HSBTsMe[-1:-2], 1, FUN="mean"), apply(HSBTsMe[-1:-2], 1, FUN="sd")/8))
HSB_rwl2$Year<-as.numeric(rownames(HSB_rwl2))
colnames(HSB_rwl2)<-c("rwlB", "sdB", "Year")
HSB_rwl.merge<-rbind(HSB_rwl2[HSB_rwl2$Year<1940,],HSB_rwl)

HSF_rwl<-data.frame(HSFTsMe)
#HSFTsMe<-bai.out(rwl=HSFTsMe)
tail(HSFTsMe, 150)
HSF_rwl<-na.omit(data.frame(apply(HSFTsMe, 1, FUN="mean"), apply(HSFTsMe, 1, FUN="sd")/9))
HSF_rwl$Year<-as.numeric(rownames(HSF_rwl))
HSF.1900<-HSFTsMe[-6]
HSF.1900<-HSF.1900[-2]
HSF_rwl2<-na.omit(data.frame(apply(HSF.1900, 1, FUN="mean"), apply(HSF.1900, 1, FUN="sd")/8))
HSF_rwl2$Year<-as.numeric(rownames(HSF_rwl2))
colnames(HSF_rwl2)<-c("rwlF", "sdF", "Year")
colnames(HSF_rwl)<-c("rwlF", "sdF", "Year")
HSF_rwl.merge<-rbind(HSF_rwl2[HSF_rwl2$Year<1915,],HSF_rwl)

str(HSF_rwl)

HS<-merge(HSF_rwl.merge, HSB_rwl.merge)
####Twin combined sets North######
str(TNBAbLa)
TNB_rwl<-data.frame(TNBAbLa)
#TNBAbLa<-bai.out(rwl=TNBAbLa)
TNBAbLa
TNB_rwl<-na.omit(data.frame(apply(TNBAbLa, 1, FUN="mean"), apply(TNBAbLa, 1, FUN="sd")/11))
TNB_rwl$Year<-as.numeric(rownames(TNB_rwl))
colnames(TNB_rwl)<-c("rwlB", "sdB", "Year")

TNF_rwl<-data.frame(TNFAbLa)
#TNFAbLa<-bai.out(rwl=TNFAbLa)
TNF_rwl<-na.omit(data.frame(apply(TNFAbLa, 1, FUN="mean"), apply(TNFAbLa, 1, FUN="sd")/10))
TNF_rwl$Year<-as.numeric(rownames(TNF_rwl))
colnames(TNF_rwl)<-c("rwlF", "sdF", "Year")
str(TNF_rwl)

TN<-merge(TNF_rwl[TNF_rwl$Year,], TNB_rwl[TNB_rwl$Year,])

####Twin combined sets SOuth######
str(TSBAbLa)
TSB_rwl<-data.frame(TSBAbLa)
#TSBAbLa<-bai.out(rwl=TSBAbLa)
TSBAbLa
TSB_rwl<-na.omit(data.frame(apply(TSBAbLa, 1, FUN="mean"), apply(TSBAbLa, 1, FUN="sd")/10))
TSB_rwl$Year<-as.numeric(rownames(TSB_rwl))
colnames(TSB_rwl)<-c("rwlB", "sdB", "Year")

TSF_rwl<-data.frame(TSFAbLa)
#TSFAbLa<-bai.out(rwl=TSFAbLa)
TSF_rwl<-na.omit(data.frame(apply(TSFAbLa, 1, FUN="mean"), apply(TSFAbLa, 1, FUN="sd")/9))
TSF_rwl$Year<-as.numeric(rownames(TSF_rwl))
colnames(TSF_rwl)<-c("rwlF", "sdF", "Year")
str(TSF_rwl)
TS<-merge(TSF_rwl[TSF_rwl$Year>1983,], TSB_rwl[TSB_rwl$Year>1983,])

#BAI graphical analysis ##
#install.packages("plotly")
#install.packages("plyr")   
library(ggplot2)
library(reshape)

### The script below is to plot each time series by "treatment" with standard errors (calculated in Excel) ###

library(ggplot2)
ggplot(BAI, aes(x = Year, y = North_High))+
  geom_line(aes(y=North_High, colour="North_High"),lwd=1) +                        
  geom_ribbon(aes(ymin = North_High - SE_NH, 
                  ymax = North_High + SE_NH),alpha = 0.2, fill="cyan")+
  geom_line(aes(y=North_Low, colour="North_Low"),lwd=1) +
  geom_ribbon(aes(ymin = North_Low - SE_NL, 
                  ymax = North_Low + SE_NL), alpha =0.2, fill="green")+
  geom_line(aes(y=South_High, colour="South_High"),lwd=1) +
  geom_ribbon(aes(ymin = South_High - SE_SH, 
                  ymax = South_High + SE_SH),lwd=.2, alpha = 0.2, fill="red")+
  geom_line(aes(y=South_Low, colour="South_Low"),lwd=1) +
  geom_ribbon(aes(ymin = South_Low - SE_SL, 
                  ymax = South_Low + SE_SL),lwd=.2, alpha = 0.2, fill="yellow")+
  scale_color_manual(values = c("cyan", "green", "red", "orange")) +
  labs(color="Locations") +
  xlab('Year') +
  ylab("Basal Area Increment (mm2)")+
  ggtitle("Pinus hartwegii Growth by Elevation and Aspect")+
  theme(plot.title = element_text(lineheight=.8, size=20, face="bold"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.title=element_text(size=14, face="bold"))

#### OK, the following script gets more into the analysis, though I haven't really run my own models yet
#### Here I created first a data frame (in the code is called "PO") for all my independent variables which relates to the response variable (BAI) through a "Tree ID" column
#### This is needed to you can then "melt" both data frames by "Tree ID" and run the analysis on your BAI data.

## You need Tidyverse 

library(Tidyverse)

##### Tranformations ####


logBAI_North.m <- melt(logBAI_North, id = "Cambial_Age")
logBAI_South.m <- melt(logBAI_South, id = "Cambial_Age")
ggplot(data = logBAI_North.m, aes(x = Cambial_Age, y = value, color="blue"))+
  geom_point(data=logBAI_South.m, aes(x = Cambial_Age, y = value, color = "red"))

### Bringing the two data frames together

allDat <- logBAI %>% #this says, take the first dataframe (tsDat) and the pipe %>% indicates that you want to perform a function to that dataframe
  gather(-NumRings, key = TreeID, value = log_BAI) %>%  #this changes the first dataframe (tsDat) from wide to long format
  full_join(PO, by = c("TreeID" = "ID")) #this joins the long version of tsDat with rDat and joins them by the RingNum column


### BAI Wide format to Long format

M_BAI_Long <- M_BAI_Wide %>% 
  gather(-NumRings, key = TreeID, value = BAI) %>%  
  full_join(PO, by = c("TreeID" = "ID")) 

### RWI Wide format to long format

RWI_Long <- RWI_Wide %>% 
  gather(-RingNum, key = TreeID, value = RWI) %>%  
  full_join(PO, by = c("TreeID" = "ID")) 

#Drop NAs

RWI_Long <- RWI_Long[!is.na(RWI_Long$RWI), ]

# Generate Cambial Age Column

library(dplyr)
CA <- Dat %>% group_by(TreeID) %>% mutate(id = row_number())

#Plot it
# By aspecta and altitude

ggplot(data = BAI_Cambial, aes(x = Cambial_Age, y = log_BAI, color = CA)) +
  geom_point(size = 2, shape = 2) +
  scale_colour_gradient(low = "green", high = "red")+
  facet_grid(.~ Aspect)+
  theme(legend.position = "bottom") + 
  xlab("Cambial Age") +
  ylab("Log BAI")

# Aspect in color

ggplot(data = BAI_Cambial, aes(x = Cambial_Age, y = log_BAI, color = Aspect)) +
  geom_point(size = 2, shape = 2) +
  scale_color_manual(values=c("green", "red"))+
  facet_grid(.~ Alt_Cat)+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  xlab("Cambial Age") +
  ylab("Log BAI")

# This creates the group means

gm <- Dat %>% 
  group_by("Aspect", "Cambial_Age") %>% 
  summarise(log_BAI = mean(log_BAI))

ggplot(data = Dat, aes(x = Cambial_Age, y = log_BAI, group=Alt_Cat, color = Aspect)) +
  geom_point(aes(group = Aspect), size = 2, shape = 1) +
  geom_line(data=gm)+
  theme(legend.position = "bottom") + 
  xlab("Cambial Age") +
  ylab("Log BAI")

# Plotting by fixed diameters #
library(tidyr)
# Move from wide to long format
RW_Long <- gather(RW, TreeID, RW, A19P3:B18P5, factor_key=FALSE)
RW_Long <- RW_Long[!is.na(RW_Long$RW), ]
# Calculate Running Totals (Radius)
CumRW <- RW_Long %>% group_by(TreeID) %>% mutate(cumsum = cumsum(RW))
# Create cm column
CumRW$Radius <- CumRW$cumsum / 10
# Merge CumRW and BAI_Cambial
BAI_Cambial <- merge(BAI_Cambial, CumRW, by = TreeID)
#Add Cambial Age as sequence # After some stupid error when I deleted CA
BAI_Cambial %>% group_by(TreeID) %>% mutate(id = row_number())
BAI_Cambial$X <- NULL
BAI_Cambial$CA <- sequence(rle(as.character(BAI_Cambial$TreeID))$lengths)
# Finally -> Master Long Format Data Set
Master_Long = read.csv("C:/Users/Paulo/Documents/R/CH2/Master_Long.csv")

# Check Diam ~ CA plots

ggplot(data = Master_Long, aes(x = CambialAge, y = Diam, color = Altitude)) +
  geom_point(size = 2, shape = 2)+
  scale_color_gradient(low="red", high="green")+
  facet_grid(.~ Aspect)
theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  xlab("Cambial Age") +
  ylab("Diameter")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))


# Join Fix Diameter Values to Long Formate Master

FixDiam45 = read.csv("C:/Users/Paulo/Documents/R/CH2/FixDiam45cm.csv")
require(plyr)
Master_Long<-join(Master_Long, FixDiam45, by = "TreeID", type="left")

# Plot BAi and RWI using Fixed Diameters
library(ggplot2)
ggplot(data = Master_Long, aes(x = Year_FD10, y = FD10_RWI, color = Aspect)) +
  geom_point(size = 3, shape = 2, stroke = 1.5) +
  geom_smooth(method="lm", show.legend = FALSE)+
  scale_color_manual(breaks= c("N", "S"),
                     values=c("green", "red"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle("RWI when DBH = 10 - 11 cm")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Year") +
  ylab("RWI")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))+
  theme(strip.text.x = element_text(size = 8, face="bold"))

# Plot growth at fixed diameter by cambia age

ggplot(data = Master_Long, aes(x = Year_FD10, y = FD10_BAI, color = Aspect)) +
  geom_point(size = 3, shape = 2, stroke = 1.5) +
  geom_smooth(method="lm", show.legend = FALSE)+
  scale_color_manual(breaks= c("N", "S"),
                     values=c("green", "red"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle("BAI when DBH = 45 - 46 cm")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Year") +
  ylab("BAI (cm2)")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))

