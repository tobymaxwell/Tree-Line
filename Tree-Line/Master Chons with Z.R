#RWI All
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Final Chrons/")
files<-list.files()

files
library(plyr)
library(dplyr)
library(sjPlot)
mergedrings<-ldply(files, read.csv)
str(mergedrings)

##########################all rings RWI ###########################
library(dplyr)
library(dplR)
PO<-NULL
temp<-NULL
rwi.all<-NULL
allchrons<-NULL
rwi.all$Year<-1535:2018
rwi.all<-as.data.frame(rwi.all)
rownames(rwi.all)<-rwi.all$Year
for(i in files){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  
  RWI<-detrend(temp,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
  rownames(temp)<-temp$Year
  temp<-temp[-1]
  PO<-NULL
  PO$series<-rep(paste0("T",1:length(colnames(temp))))
  PO<-data.frame(PO)
  cols<-colnames(temp)
  lengths<-NULL
  for(j in cols){
    length<-length(rownames(na.omit(select(temp, j)[1])))
    lengths<-c(lengths, length)
  }
  PO$pith.offset<-lengths
  temp
  colnames(temp)<-PO$series
  write.rwl(temp, paste0("/Users/tobymaxwell/Desktop/",i))
  temp<-read.rwl(paste0("/Users/tobymaxwell/Desktop/",i), header=TRUE)
  RWI_C<-cms(temp, PO, c.hat.t = FALSE, c.hat.i = FALSE)
  colnames(RWI_C)<-cols
  MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
  MeanChron<-data.frame(MeanChron$IZTstd, rownames(MeanChron))
  colnames(MeanChron)<-c(paste0("rwi",i), "Year")
  rwi.all<-merge(rwi.all, MeanChron, all=T)
  
}
colnames(rwi.all)<-c("Year",substr(files, 1,7))
str(rwi.all)
tail(rwi.all$MNBAbLa)
library(tidyr)
rwi.long<-gather(rwi.all, key="ID", value = "rwi", DNBPiAl:YSFTsMe)

rwi.long$Site<-substr(rwi.long$ID, 0,1)
rwi.long$Aspect<-substr(rwi.long$ID, 2,2)
rwi.long$Zone<-substr(rwi.long$ID, 3,3)
rwi.long$Species<-substr(rwi.long$ID, 4,7)
str(rwi.long)

library(wesanderson)
ggplot(rwi.long[rwi.long$Species=="AbLa",], aes(rwi,x=Year, color=Species))+
  geom_line()+
  stat_smooth(span=0.2)+
  facet_grid(Aspect~Site)+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)


#####Resiliance Stats#####
library(pointRes)

library(dplyr)
library(dplR)

setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Final Chrons/")
files<-list.files()


PO<-NULL
temp<-NULL
rwi.all<-NULL
rwi.all$Year<-1535:2018
res.all<-NULL
rwi.all<-as.data.frame(rwi.all)
rownames(rwi.all)<-rwi.all$Year
for(i in files){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  
  RWI<-detrend(temp,make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
  rownames(temp)<-temp$Year
  temp<-temp[-1]
  PO<-NULL
  PO$series<-rep(paste0("T",1:length(colnames(temp))))
  PO<-data.frame(PO)
  cols<-colnames(temp)
  lengths<-NULL
  for(j in cols){
    length<-length(rownames(na.omit(select(temp, j)[1])))
    lengths<-c(lengths, length)
  }
  PO$pith.offset<-lengths
  temp
  colnames(temp)<-PO$series
  write.rwl(temp, paste0("/Users/tobymaxwell/Desktop/rings/",i))
  temp<-read.rwl(paste0("/Users/tobymaxwell/Desktop/rings/",i), header=TRUE)
  res<-res.comp(temp, nb.yrs=4, series.thresh = 50, res.thresh.neg = 25)$out
  ###
  
  res$ID<-substr(cols[1],1,7)
  res.all<-rbind(res.all, res)
  
}
str(res.all)

res.all<-res.all[res.all$nature!=0,]
str(res.all)
res.all$Site<-substr(res.all$ID, 0,1)
res.all$Aspect<-substr(res.all$ID, 2,2)
res.all$Zone<-substr(res.all$ID, 3,3)
res.all$Species<-substr(res.all$ID, 4,7)
res.all$Site<-as.factor(res.all$Site)
res.all<-as.data.frame(res.all)
str(res.all)
res.all<-do.call(data.frame,lapply(res.all, function(x) replace(x, is.infinite(x),NA)))
res.all[res.all$Species=="AbLA",]$Species<-"AbLa"
res.all[res.all$Species=="ABLa",]$Species<-"AbLa"
str(res.all)

res.all$nature<-(0)
ggplot(rwi.long, aes(rwi,x=Year, color=Species))+
  geom_line()+
  stat_smooth(span=0.2)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1984,2018)+ylim(0,4)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))


######NDVI ######

setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Landsat/")

library(plyr)
list<-list.files()
IDs<-NULL
ndvi.list<-NULL
for(i in list){
  ndvi<-na.omit(read.csv(i))
  ndvi <- ndvi[order(ndvi[4]), ]
  top3<-by(ndvi, ndvi["Year"], tail, n=3)
  ndvi.list<-list(ndvi.list,top3)
  assign(substr(i, 1,3), top3)
  IDs<-c(IDs, substr(i, 1,3))
}
ndvi.list<-list(DNB, DNF, DSB,DSF,ENB,ENF,ESB,ESF, HNB,HNF,HSB,HSF,MNB,MNF,MSB,MNF,PNB,PNF,PSB,PSF,TNB,TNF,TSB,TSF,YNB,YNF,YSB,YSF)

years<-as.factor(1984:2018)
meanday<-NULL
days<-NULL
NDVImax<-NULL
for(h in ndvi.list){
  days<-NULL
  for(i in 1:35){
    meanday<-data.frame(h[i])[,4]
    days<-rbind(days, meanday)
  }
  rownames(days)<-1984:2018
  NDVImax<-cbind(NDVImax,apply(days, 1, mean))
}
NDVImax<-data.frame(NDVImax)
colnames(NDVImax)<-IDs
NDVImax$Year<-1984:2018
NDVI.long<-gather(NDVImax, key="ID", value="NDVI", DNB:YSF)
NDVI.long$Site<-substr(NDVI.long$ID, 1,1)
NDVI.long$Aspect<-substr(NDVI.long$ID, 2,2)
NDVI.long$Zone<-substr(NDVI.long$ID, 3,3)
NDVI.ag<-NDVI.long%>%
  group_by(Site, Aspect, Zone, Year)%>%
  summarize(NDVI=mean(NDVI))
NDVI.ag<-as.data.frame(NDVI.ag)
NDVI.ag
####masterplot#####
library(dplyr)
climhist<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")[-1]
str(climhist)
climhist.yr<-climhist%>%
  group_by(Site, Aspect, Year)%>%
  summarise(MAP=sum(MAP), MAT=mean(MAT), tmin=mean(tmin), tmax=mean(tmax))
climhist.yr<-as.data.frame(climhist.yr)
climhist.yr

library(wesanderson)
ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Zone~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1985,2018)+ylim(-2.5,2.5)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=NDVI.ag,span=0.2,se=F,color='Black',lwd=0.5,aes(y=scale(NDVI), x=Year))

ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)+ylim(-2,2)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=NDVI.ag,span=0.2,se=F,aes(y=scale(NDVI), x=Year, color='Black'))

ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)+ylim(-2,2)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=climhist.yr,span=0.2,se=F, color='Black',lwd=0.5, aes(y=scale(MAP), x=Year))

ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)+ylim(-2,2)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=climhist.yr,span=0.2,se=F,aes(y=scale(MAT), x=Year, color='Black'))

ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)+ylim(-2,2)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=climhist.yr,span=0.2,se=F,aes(y=scale(tmin), x=Year, color='Black'))

ggplot(rwi.long, aes(scale(rwi),x=Year, color=Species))+
  stat_smooth(span=0.2, se=F)+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()+xlim(1895,2018)+ylim(-2,2)+
  geom_point(data=res.all, aes(nature,x=year, color=Species))+
  stat_smooth(data=climhist.yr,span=0.2,se=F,aes(y=scale(tmax), x=Year, color='Black'))

#####time series breakpoints#####

library(TTR)
library(bfast)

ts.slopes<-function(ts){
  breaks<-bfast(ts, season = "none", max.iter=2,breaks=3)
  niter <- length(breaks$output)
  slopes<-coef(breaks$output[[niter]]$bp.Vt)[,2]
  intercepts<-coef(breaks$output[[niter]]$bp.Vt)[,1]
  ci<-(breaks$output[[niter]]$ci.Vt)
  nm <-deparse(substitute(ts))
  return(c(slopes, intercepts, plot(breaks, type="trend", main = paste(nm))))
}
######all data breaks
d<-na.omit(rwi.long)%>%
  group_by(Year, Species)%>%
  summarize(rwi=mean(rwi))
d<-as.data.frame(d[d$Year>1839,])
ts.d<-ts(d$rwi, start=min(d$Year))
slopes<-ts.slopes(ts.d)
print(slopes)

######

sites<-levels(as.factor(rwi.long$ID))
slopes<-NULL
result<-NULL
breakyears<-NULL

for (i in sites[c(-1,-3,-12:-13,-23,-25:-26,-43, -11,-20, -24,-27,-29, -39,-40, -41,-42)]) {
  d<-na.omit(rwi.long[rwi.long$ID==i&rwi.long$Year>1800,])
  ts.d<-ts(d$rwi, start=min(d$Year))
  slopes<-ts.slopes(ts.d)
  print(slopes)
  
  bks<-bfast(ts.d,season = "none", max.iter=2,breaks=2)
  bks<-bks$output[[2]]$ci.Vt$confint[,2]+min(d$Year)-1
  
  breakyears<-rbind(breakyears, bks)
  result<-rbind(result,slopes)
}
colnames(result)<-c("S1", "S2", "S3", "I1", "I2","I3")
result<-(as.data.frame(result))
rownames(result)<-sites[c(-1,-3,-12:-13,-23,-25:-26,-43, -11, -20, -24,-27,-29, -39,-40, -41,-42)]
rownames(breakyears)<-sites[c(-1,-3,-12:-13,-23,-25:-26,-43, -11, -20, -24,-27,-29, -39,-40, -41,-42)]
length(result)

slopes.1<-as.data.frame(result)
slopes.1$ID<-rownames(slopes.1)
slopes.1$Site<-substr(slopes.1$ID,1,1)
slopes.1$Aspect<-substr(slopes.1$ID,2,2)
slopes.1$Zone<-substr(slopes.1$ID,3,3)
slopes.1$Species<-substr(slopes.1$ID,4,7)
slopes.1$Change2.3<-slopes.1$S3-slopes.1$S2

breakyears<-as.data.frame(breakyears)
breakyears$ID<-rownames(breakyears)
breakyears$Site<-substr(breakyears$ID,1,1)
breakyears$Aspect<-substr(breakyears$ID,2,2)
breakyears$Zone<-substr(breakyears$ID,3,3)
breakyears$Species<-substr(breakyears$ID,4,7)
colnames(breakyears)<-c("B1", "B2", "ID", "Site", "Aspect", "Zone", "Species")
rwi.minyear<-rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1800,]%>%
  group_by(ID)%>%
  summarize(xmin=min(Year))
rwi.minyear
rwi.minyear<-merge(rwi.minyear, breakyears)
rwi.minyear
rwi.minyear$xmax<-2018

rwi.minyear$Y1<-rwi.minyear$xmin*slopes.1$S1+slopes.1$I1
rwi.minyear$Y2<-rwi.minyear$B1*slopes.1$S1+slopes.1$I1
rwi.minyear$Y3<-rwi.minyear$B1*slopes.1$S2+slopes.1$I2
rwi.minyear$Y4<-rwi.minyear$B2*slopes.1$S2+slopes.1$I2
rwi.minyear$Y5<-rwi.minyear$B2*slopes.1$S3+slopes.1$I3
rwi.minyear$Y6<-rwi.minyear$xmax*slopes.1$S3+slopes.1$I3
rwi.minyear


######breaks for sites with 1 point#####

sites<-levels(as.factor(rwi.long$ID))
slopes<-NULL
result<-NULL
breakyear<-NULL


for (i in sites[c(1, 11, 20, 24,27,29, 39,40, 41,42)]) {
  d<-na.omit(rwi.long[rwi.long$ID==i&rwi.long$Year>1800,])
  ts.d<-ts(d$rwi, start=min(d$Year))
  slopes<-ts.slopes(ts.d)
  print(slopes)
  
  bks<-bfast(ts.d,season = "none", max.iter=2,breaks=2)
  bks<-bks$output[[2]]$ci.Vt$confint[,2]+min(d$Year)-1
  
  breakyear<-rbind(breakyear, bks)
  result<-rbind(result,slopes)
}
colnames(result)<-c("S1", "S2", "I1", "I2")
rownames(result)<-sites[c(1,11,20, 24,27,29, 39,40, 41,42)]
rownames(breakyear)<-sites[c(1,11,20, 24,27,29, 39,40, 41,42)]

slopes.2<-as.data.frame(result)
slopes.2$ID<-rownames(slopes.2)
slopes.2$Site<-substr(slopes.2$ID,1,1)
slopes.2$Aspect<-substr(slopes.2$ID,2,2)
slopes.2$Zone<-substr(slopes.2$ID,3,3)
slopes.2$Species<-substr(slopes.2$ID,4,7)

breakyear<-as.data.frame(breakyear)
breakyear$ID<-rownames(breakyear)
breakyear$Site<-substr(breakyear$ID,1,1)
breakyear$Aspect<-substr(breakyear$ID,2,2)
breakyear$Zone<-substr(breakyear$ID,3,3)
breakyear$Species<-substr(breakyear$ID,4,7)
colnames(breakyear)<-c("B1", "ID", "Site", "Aspect", "Zone", "Species")
rwi.minyear.2<-rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1800,]%>%
  group_by(ID)%>%
  summarize(xmin=min(Year))
rwi.minyear.2
rwi.minyear.2<-merge(rwi.minyear.2, breakyear)
rwi.minyear.2
rwi.minyear.2$xmax<-2018

rwi.minyear.2$Y1<-rwi.minyear.2$xmin*slopes.2$S1+slopes.2$I1
rwi.minyear.2$Y2<-rwi.minyear.2$B1*slopes.2$S1+slopes.2$I1
rwi.minyear.2$Y3<-rwi.minyear.2$B1*slopes.2$S2+slopes.2$I2
rwi.minyear.2$Y4<-rwi.minyear.2$xmax*slopes.2$S2+slopes.2$I2
rwi.minyear.2

#####sites with no breaks######
coefs<-NULL
res<-NULL
for(i in sites[c(3,12:13,23,25:26,43)]){
res<-summary(lm(rwi~Year, rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1800&rwi.long$ID==i,]))
coefs<-rbind(coefs,res$coefficients[-5:-6])

}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs<-as.data.frame(coefs)
coefs$ID<-sites[c(3,12:13,23,25:26,43)]
coefs[coefs$pyear<0.05,]
coefs$Site<-substr(coefs$ID,1,1)
coefs$Aspect<-substr(coefs$ID,2,2)
coefs$Zone<-substr(coefs$ID,3,3)
coefs$Species<-substr(coefs$ID,4,7)

ggplot(data=rwi.minyear, aes(x=xmin, xend=B1, y=Y1, yend=Y2))+geom_segment()+
  geom_segment(aes(x=B1, xend=B2, y=Y3, yend=Y4),color='Black', rwi.minyear)+
  geom_segment(aes(x=B2, xend=xmax, y=Y5, yend=Y6),color='Black',rwi.minyear)+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2),color='Black',rwi.minyear.2)+
  geom_segment(aes(x=B1, xend=xmax, y=Y3, yend=Y4),color='Black',rwi.minyear.2)+
  geom_abline(data=coefs[coefs$Species=="PiAl",], aes(slope=Year, intercept=int))+
  ylim(0,3)


ggplot(rwi.long[rwi.long$Species=="PiAl",], aes(rwi, x=Year, color=Aspect))+
  geom_line()+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyears[breakyears$Species=="PiAl",])+
  geom_vline(aes(xintercept=B2, color=Aspect), breakyears[breakyears$Species=="PiAl",])+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyear[breakyear$Species=="PiAl",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2),color='Black',rwi.minyear[rwi.minyear$Species=="PiAl",])+
  geom_segment(aes(x=B1, xend=B2, y=Y3, yend=Y4),color='Black', rwi.minyear[rwi.minyear$Species=="PiAl",])+
  geom_segment(aes(x=B2, xend=xmax, y=Y5, yend=Y6),color='Black',rwi.minyear[rwi.minyear$Species=="PiAl",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2),color='Black',rwi.minyear.2[rwi.minyear.2$Species=="PiAl",])+
  geom_segment(aes(x=B1, xend=xmax, y=Y3, yend=Y4),color='Black',rwi.minyear.2[rwi.minyear.2$Species=="PiAl",])+xlim(1801,2018)+
  geom_abline(data=coefs[coefs$Species=="PiAl",], aes(slope=Year, intercept=int))+
  facet_grid(Site~Zone, scales='free')



ggplot(rwi.long[rwi.long$Species=="AbLa",], aes(rwi, x=Year, color=Aspect))+
  geom_line()+facet_grid(Site~Zone, scales='free')+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyears[breakyears$Species=="AbLa",])+
  geom_vline(aes(xintercept=B2, color=Aspect), breakyears[breakyears$Species=="AbLa",])+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyear[breakyear$Species=="AbLa",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2),color='Black',rwi.minyear[rwi.minyear$Species=="AbLa",])+
  geom_segment(aes(x=B1, xend=B2, y=Y3, yend=Y4),color='Black', rwi.minyear[rwi.minyear$Species=="AbLa",])+
  geom_segment(aes(x=B2, xend=xmax, y=Y5, yend=Y6),color='Black',rwi.minyear[rwi.minyear$Species=="AbLa",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2),color='Black',rwi.minyear.2[rwi.minyear.2$Species=="AbLa",])+
  geom_segment(aes(x=B1, xend=xmax, y=Y3, yend=Y4),color='Black',rwi.minyear.2[rwi.minyear.2$Species=="AbLa",])+xlim(1801,2018)

ggplot(rwi.long[rwi.long$Species=="TsMe",], aes(rwi, x=Year, color=Aspect))+
  geom_line()+facet_grid(Site~Zone, scales='free')+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyears[breakyears$Species=="TsMe",])+
  geom_vline(aes(xintercept=B2, color=Aspect), breakyears[breakyears$Species=="TsMe",])+
  geom_vline(aes(xintercept=B1, color=Aspect), breakyear[breakyear$Species=="TsMe",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2), color='black', rwi.minyear[rwi.minyear$Species=="TsMe",])+
  geom_segment(aes(x=B1, xend=B2, y=Y3, yend=Y4),   color='black',rwi.minyear[rwi.minyear$Species=="TsMe",])+
  geom_segment(aes(x=B2, xend=xmax, y=Y5, yend=Y6), color='black',rwi.minyear[rwi.minyear$Species=="TsMe",])+
  geom_segment(aes(x=xmin, xend=B1, y=Y1, yend=Y2), color='black',rwi.minyear.2[rwi.minyear.2$Species=="TsMe",])+
  geom_segment(aes(x=B1, xend=xmax, y=Y3, yend=Y4), color='black',rwi.minyear.2[rwi.minyear.2$Species=="TsMe",])+theme_bw()+xlim(1801,2018)

modslopes<-c(slopes.1$S3, slopes.2$S2)
modslopes<-as.data.frame(modslopes)
modslopes$ID<-c(slopes.1$ID, slopes.2$ID)
modslopes$Site<-substr(modslopes$ID, 1,1)
modslopes$Aspect<-substr(modslopes$ID,2,2)
modslopes$Zone<-substr(modslopes$ID,3,3)

clim.norms<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/clim.norms.csv")
clim.norms<-clim.norms%>%
  group_by(Site, Aspect, Zone)%>%
  summarize(mat=mean(mat), ppt=mean(ppt), tmin=mean(tmin), tmax=mean(tmax), vpdmax=mean(vpdmax))
clim.norms<-as.data.frame(clim.norms)
str(gps.trees)
gps.trees<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.trees.csv")
gis.SA<-gps.trees%>%
  group_by(Site, Aspect, Zone, Species)%>%
  summarize(tri=mean(tri), tpi=mean(tpi), slope=mean(slope), curv=mean(curv))
climtop<-merge(clim.norms, gis.SA, by=c("Site", "Aspect", "Zone"))
str(soillit)
data.frame(gis.SA)
allsoil.SA<-na.omit(soillit)%>%
  group_by(Site, Aspect, Zone)%>%
  summarize(Sand=mean(Sand), Silt=mean(Silt), Clay=mean(Clay), C.N=mean(C.N), C.pct=mean(C.pct),C.N.L=mean(C.N.L))
climtopsoil<-merge(climtop,allsoil.SA)
str(leafdat)
leafdat.cut<-leafdat%>%
  group_by(Site, Aspect, Zone, Species)%>%
  summarize(DBH=mean(DBH), SLA=mean(SLA), latitude=mean(latitude))
data.frame(leafdat.cut)

Slopedata<-merge(climtopsoil, leafdat.cut, all=T)
str(Slopedata)



#####Allsites breaks defined by master chron#####

ts.slopes<-function(ts){
  breaks<-bfast(ts, season = "none", max.iter=3,breaks=2)
  niter <- length(breaks$output)
  slopes<-coef(breaks$output[[niter]]$bp.Vt)[,2]
  intercepts<-coef(breaks$output[[niter]]$bp.Vt)[,1]
  ci<-(breaks$output[[niter]]$ci.Vt)
  nm <-deparse(substitute(ts))
  return(c(slopes, intercepts, plot(breaks, type="trend", main = paste(nm))))
}
######all data breaks
d<-na.omit(rwi.long)%>%
  group_by(Year)%>%
  summarize(rwi=mean(rwi))
d<-as.data.frame(d[d$Year>1836,])
ts.d<-ts(d$rwi, start=min(d$Year))
slopes<-ts.slopes(ts.d)
print(slopes)

ggplot(rwi.long, aes(y=log(rwi), x=Year))+geom_point()+ylim(0,4)
#####1948-2018#####
coefs<-NULL
res<-NULL
sites<-levels(factor(rwi.long$ID))
for(i in sites){
  res<-summary(lm(rwi~Year, rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1948&rwi.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1948<-as.data.frame(coefs)
coefs.1948$ID<-sites
coefs.1948$sig.05<-ifelse(coefs.1948$pyear<0.05, yes=1, no=0)
coefs.1948$dir<-ifelse(coefs.1948$Year>0, yes=1, no=-1)
coefs.1948$Site<-substr(coefs.1948$ID,1,1)
coefs.1948$Aspect<-substr(coefs.1948$ID,2,2)
coefs.1948$Zone<-substr(coefs.1948$ID,3,3)
coefs.1948$Species<-substr(coefs.1948$ID,4,7)
coefs.1948[coefs.1948$sig.05==1&coefs.1948$dir==-1,]
length(coefs.1948[coefs.1948$dir==1&coefs.1948$sig.05==1,]$dir)

coefs.1948$ParentMaterial<-NA
coefs.1948[coefs.1948$Site=='H'|coefs.1948$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1948[coefs.1948$Site=='E'|coefs.1948$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1948[coefs.1948$Site=='M'|coefs.1948$Site=='P',]$ParentMaterial<-'Granite'
coefs.1948[coefs.1948$Site=='D',]$ParentMaterial<-'Mixed Volvanics'
coefs.1948.dat<-merge(coefs.1948[coefs.1948$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone", 'Species'))
ggplot(rwi.long, aes(rwi, x=Year))+geom_blank()+
  geom_abline(data=coefs.1948.dat, aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.1948[coefs.1948$sig.05==0,],color='grey', size=.5, aes(slope=Year, intercept=int))+
  xlim(1948,2018)+ylim(0,3)+theme_bw()

summary(lm(Year~C.N.L+C.pct, coefs.1948.dat))
#pctc, C.N.L
str(rwi.long)
library(corrplot)
str(coefs.1948.dat)
M<-cor(na.omit(coefs.1948.dat[c(-1:-4,-11:-14, -30)]))
res1 <- cor.mtest(na.omit(coefs.1948.dat[c(-1:-4,-11:-14, -30)], conf.level = .90))

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


coefs<-NULL
res<-NULL
for(i in sites){
  res<-summary(lm(rwi~Year, rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year<1949&rwi.long$Year>1890&rwi.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1890<-as.data.frame(coefs)
coefs.1890$ID<-sites
coefs.1890$sig.05<-ifelse(coefs.1890$pyear<0.05, yes=1, no=0)
coefs.1890$dir<-ifelse(coefs.1890$Year>0, yes=1, no=-1)
coefs.1890$Site<-substr(coefs.1890$ID,1,1)
coefs.1890$Aspect<-substr(coefs.1890$ID,2,2)
coefs.1890$Zone<-substr(coefs.1890$ID,3,3)
coefs.1890$Species<-substr(coefs.1890$ID,4,7)
length(coefs.1890[coefs.1890$sig.05==1&coefs.1890$dir==-1,]$Year)

coefs.1890$ParentMaterial<-NA
coefs.1890[coefs.1890$Site=='H'|coefs.1890$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1890[coefs.1890$Site=='E'|coefs.1890$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1890[coefs.1890$Site=='M'|coefs.1890$Site=='P',]$ParentMaterial<-'Granite'
coefs.1890[coefs.1890$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1890.dat<-merge(coefs.1890[coefs.1890$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone", 'Species'))
ggplot(rwi.long, aes(rwi, x=Year))+geom_blank()+
  geom_abline(data=coefs.1890.dat, aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.1890[coefs.1890$sig.05==0,], color='grey', aes(slope=Year, intercept=int))+
  xlim(1890,1948)+ylim(0,3)+theme_bw()

anova(lm(Year~tmax+tmax:Species, coefs.1890.dat))


str(coefs.1890.dat)

M<-cor(coefs.1890.dat[c(-1:-4,-7:-14, -30)])
res1 <- cor.mtest(coefs.1890.dat[c(-1:-4,-7:-14, -30)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

#########

coefs<-NULL
res<-NULL
sites.1840<-levels(as.factor(rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1836&rwi.long$Year<1891,]$ID))
for(i in sites.1840){
  res<-summary(lm(rwi~Year, rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1836&rwi.long$Year<1891&rwi.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1840<-as.data.frame(coefs)
coefs.1840$ID<-sites.1840
coefs.1840$sig.05<-ifelse(coefs.1840$pyear<0.05, yes=1, no=0)
coefs.1840$dir<-ifelse(coefs.1840$Year>0, yes=1, no=-1)
coefs.1840$Site<-substr(coefs.1840$ID,1,1)
coefs.1840$Aspect<-substr(coefs.1840$ID,2,2)
coefs.1840$Zone<-substr(coefs.1840$ID,3,3)
coefs.1840$Species<-substr(coefs.1840$ID,4,7)
length(coefs.1840[coefs.1840$sig.05==1&coefs.1840$dir==-1,]$Year)


coefs.1840$ParentMaterial<-NA
coefs.1840[coefs.1840$Site=='H'|coefs.1840$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1840[coefs.1840$Site=='E'|coefs.1840$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1840[coefs.1840$Site=='M'|coefs.1840$Site=='P',]$ParentMaterial<-'Granite'
coefs.1840[coefs.1840$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1840.dat<-merge(coefs.1840[coefs.1840$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone", 'Species'))
ggplot(rwi.long, aes(rwi, x=Year))+geom_blank()+
  geom_abline(data=coefs.1840.dat, aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.1840[coefs.1840$sig.05==0,], color='grey', aes(slope=Year, intercept=int))+
  xlim(1840,1890)+ylim(0,3)+theme_bw()

summary(lm(Year~ParentMaterial, coefs.1840.dat))

str(coefs.1840.dat)

M<-cor(coefs.1840.dat[c(-1:-4,-7:-14, -30)])
res1 <- cor.mtest(coefs.1840.dat[c(-1:-4,-7:-14, -30)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
