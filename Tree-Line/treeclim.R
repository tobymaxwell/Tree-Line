#treeclim

library(treeclim)
library(ggplot2)
library(dplyr)
library(tidyr)
bai.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.long.csv")[-1]
bai.spgl<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.spgl.csv")[-1]
climhist<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")[-1]

#BAI All

bai.long$ID<-substr(bai.long$ID,1,3)
bai.tree<-bai.long%>%
  group_by(ID, Species, Year)%>%
  summarize(bai=mean(bai))
bai.tree<-as.data.frame(bai.tree)
bai.tree



climseas<-climhist
climseas$Season<-NULL
climseas$Season[climseas$Month==10|climseas$Month==11|climseas$Month==12]<-"Fall"
climseas$Season[climseas$Month==1|climseas$Month==2|climseas$Month==3]<-"Winter"
climseas$Season[climseas$Month==4|climseas$Month==5|climseas$Month==6]<-"Spring"
climseas$Season[climseas$Month==7|climseas$Month==8|climseas$Month==9]<-"Summer"
str(climseas)

coefvar<-function(x){(sd(x)/mean(x))}

tail(climseas)
climseas.mean<-climseas%>%
  group_by(Site, Aspect, Zone, Year)%>%
  summarise(cvppt=coefvar(ppt),ppt=sum(ppt), tmean=mean(tmean), tmin=min(tmin), tmax=max(tmax), vpdmax=max(vpdmax))

climseas.mean<-climseas.mean%>%
  rename(ppt.yr=ppt, tmean.yr=tmean, tmin.yr=tmin, tmax.yr=tmax, vpdmax.yr=vpdmax)

climseas.mean<-as.data.frame(climseas.mean)
climseas.mean$ppt.lag1<-lag(climseas.mean$ppt.yr, 1)
climseas.mean$ppt.lag2<-lag(climseas.mean$ppt.yr, 2)
climseas.mean$ppt.lag3<-lag(climseas.mean$ppt.yr, 3)
climseas.mean<-climseas.mean[climseas.mean$Year>1897,]

ggplot(climseas.mean, aes(cvppt, x=Site,fill=Aspect))+geom_boxplot()
         
climseas<-climseas%>%
  group_by(Site, Aspect, Zone, Year, Season)%>%
  summarise(ppt=sum(ppt), tmean=mean(tmean), tmin=min(tmin), tmax=max(tmax), vpdmax=max(vpdmax))
climseas<-as.data.frame(climseas)
tail(climseas)

climseas<-climseas %>% 
  gather("climate", "ID", ppt:vpdmax) %>%
  unite( "climseas", Season:climate)%>%
  spread(climseas, ID)

tail(climseas)
str(climseas.mean)
str(bai.spgl)
str(climseas)
bai.clim<-merge(bai.spgl, climseas, by=c("Year", "Site", "Aspect", "Zone"))
bai.clim<-merge(bai.clim, climseas.mean, by=c("Year", "Site", "Aspect", "Zone"))
str(bai.clim)
max(bai.clim$Year)

summary(lm(lnbai~ppt.lag3, bai.clim))

cor()
library(treeclim)


str(climhist)
cordat.month<-c("SEP", "AUG", "JUL", "JUN", "MAY", "APR", "MAR", "FEB", "JAN", "Dec.prev", "Nov.prev", "Oct.prev", "Sep.prev", "Aug.prev" )


climhist <- climhist[order(climhist[2]),]
str(climhist)
bai.prism<-bai.tree[bai.tree$Year>1894&bai.tree$Year<2015,]
bai.prism.abla<-bai.prism[bai.prism$Species=="AbLa",]
Sites<-levels(as.factor(bai.prism.abla$ID))

clim<-NULL
tmean<-NULL
ppt<-NULL
cordat<-NULL
ppt.all<-NULL
tmean.all<-as.data.frame(NULL)

for(i in Sites){
bai<-bai.prism.abla[bai.prism.abla$ID==i,]$bai
bai.abla<-as.data.frame(bai)
n<-2014-length(bai)+1
rownames(bai.abla)<-n:2014

clim<-climhist[climhist$ID==i,][c(1,2,9,10)]

cordat<-seascorr(bai.abla, clim, season_lengths = c(1:12))
plot(cordat)

ppt<-cordat$coef[[1]]$primary
ppt$month<-cordat.month
ppt$Site<-paste0(i)
ppt$Species<-"AbLa"
ppt<-ppt[ppt$significant==TRUE,]


tmean<-cordat$coef[[1]]$secondary
tmean$month<-cordat.month
tmean$Site<-paste0(i)
tmean$Species<-"AbLa"
tmean<-tmean[tmean$significant==TRUE,]

ppt.all<-rbind(ppt.all, ppt)
tmean.all<-rbind(tmean.all, tmean)

  }
}


clim<-NULL
tmin<-NULL
tmax<-NULL
cordat<-NULL
tmin.all<-NULL
tmax.all<-NULL

for(i in Sites){
  bai<-bai.prism.abla[bai.prism.abla$ID==i&bai.prism.abla$Year>2004,]$bai
  bai.abla<-as.data.frame(bai)
  n<-2014-length(bai)+1
  rownames(bai.abla)<-n:2014
  
  clim<-climhist[climhist$ID==i,][c(1,2,11,12)]
  
  cordat<-seascorr(bai.abla, clim, season_lengths = c(3))
  plot(cordat)
  
  tmax<-cordat$coef[[1]]$primary
  tmax$month<-cordat.month
  tmax$Site<-paste0(i)
  tmax$Species<-"AbLa"
  tmax<-tmax[tmax$significant==TRUE,]
  
  
  tmin<-cordat$coef[[1]]$secondary
  tmin$month<-cordat.month
  tmin$Site<-paste0(i)
  tmin$Species<-"AbLa"
  tmin<-tmin[tmin$significant==TRUE,]
  
  tmax.all<-rbind(ppt.all, ppt)
  tmin.all<-rbind(tmean.all, tmean)
  
}




######PiAl
climhist <- climhist[order(climhist[3]),]
rwi.prism<-rwi.long[rwi.long$Year>1894&rwi.long$Year<2015,]
rwi.prism.pial<-rwi.prism[rwi.prism$Species=="PiAl",]
rwi.prism.pial$SA<-paste0(rwi.prism.pial$Site, rwi.prism.pial$Aspect)
Sites<-levels(as.factor(rwi.prism.pial$SA))
climhist$SA<-paste0(climhist$Site, climhist$Aspect)

clim<-NULL
tmean<-NULL
ppt<-NULL
cordat<-NULL

for(i in Sites){
  rwi<-rwi.prism.pial[rwi.prism.pial$SA==i,]$rwi
  rwi.pial<-as.data.frame(rwi)
  rownames(rwi.pial)<-1895:2014
  
  clim<-climhist[climhist$SA==i,][3:8][-5:-6]
  
  cordat<-seascorr(rwi.pial, clim, season_lengths = 3)
  
  ppt<-cordat$coef[[1]]$primary
  ppt$month<-cordat.month
  ppt$Site<-paste0(i)
  ppt$Species<-"PiAl"
  ppt<-ppt[ppt$significant==TRUE,]
  
  
  tmean<-cordat$coef[[1]]$secondary
  tmean$month<-cordat.month
  tmean$Site<-paste0(i)
  tmean$Species<-"PiAl"
  tmean<-tmean[tmean$significant==TRUE,]
  
  ppt.all<-rbind(ppt.all, ppt)
  tmean.all<-rbind(tmean.all, tmean)
  
}
}

climhist <- climhist[order(climhist[3]),]
rwi.prism<-rwi.long[rwi.long$Year>1894&rwi.long$Year<2015,]
rwi.prism.tsme<-rwi.prism[rwi.prism$Species=="TsMe",]
rwi.prism.tsme$SA<-paste0(rwi.prism.tsme$Site, rwi.prism.tsme$Aspect)
Sites<-levels(as.factor(rwi.prism.tsme$SA))
climhist$SA<-paste0(climhist$Site, climhist$Aspect)

clim<-NULL
tmean<-NULL
ppt<-NULL
cordat<-NULL

for(i in Sites){
  rwi<-rwi.prism.tsme[rwi.prism.tsme$SA==i,]$rwi
  rwi.tsme<-as.data.frame(rwi)
  rownames(rwi.tsme)<-1895:2014
  
  clim<-climhist[climhist$SA==i,][3:8][-5:-6]
  
  cordat<-seascorr(rwi.tsme, clim, season_lengths = 3)
  
  ppt<-cordat$coef[[1]]$primary
  ppt$month<-cordat.month
  ppt$Site<-paste0(i)
  ppt$Species<-"TsMe"
  ppt<-ppt[ppt$significant==TRUE,]
  
  
  tmean<-cordat$coef[[1]]$secondary
  tmean$month<-cordat.month
  tmean$Site<-paste0(i)
  tmean$Species<-"TsMe"
  tmean<-tmean[tmean$significant==TRUE,]
  
  ppt.all<-rbind(ppt.all, ppt)
  tmean.all<-rbind(tmean.all, tmean)
  
}
}
