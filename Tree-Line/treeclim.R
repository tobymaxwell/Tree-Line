#treeclim

library(treeclim)
library(ggplot2)


#BAI All
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Master Chrons/")
files<-list.files()

files
library(plyr)
library(dplyr)
library(sjPlot)
library(tidyr)
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
colnames(rwi.all)<-c("Year",substr(files, 1,6))
str(rwi.all)
rwi.long<-gather(rwi.all, key="ID", value = "rwi", DNPiAl:YSTsMe)

rwi.long$Site<-substr(rwi.long$ID, 0,1)
rwi.long$Aspect<-substr(rwi.long$ID, 2,2)
rwi.long$Species<-substr(rwi.long$ID, 3,6)
str(rwi.long)



climhist<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.csv")[-1]
str(climhist)
str(na.omit(climhist))
climseas<-climhist
climseas$Season<-NULL
climseas$Season[climseas$Month==10|climseas$Month==11|climseas$Month==12]<-"Fall"
climseas$Season[climseas$Month==1|climseas$Month==2|climseas$Month==3]<-"Winter"
climseas$Season[climseas$Month==4|climseas$Month==5|climseas$Month==6]<-"Spring"
climseas$Season[climseas$Month==7|climseas$Month==8|climseas$Month==9]<-"Summer"
str(climseas)
climseas.mean<-climseas%>%
  group_by(Site, Aspect, Year, Season)%>%
  summarise(ppt=sum(MAP), tmean=mean(MAT), tmin=mean(tmin), tmax=mean(tmax))
climseas.mean<-as.data.frame(climseas.mean)
str(climseas.mean)
ggplot(climseas.mean[climseas.mean$Season=="Winter",], aes(ppt, x=Year,color=Aspect))+geom_line()+stat_smooth(span=.4)+facet_wrap(~Site, scales='free')


library(treeclim)


str(climhist)
cordat.month<-c("SEP", "AUG", "JUL", "JUN", "MAY", "APR", "MAR", "FEB", "JAN", "Dec.prev", "Nov.prev", "Oct.prev", "Sep.prev", "Aug.prev" )


climhist <- climhist[order(climhist[3]),]
rwi.prism<-rwi.long[rwi.long$Year>1894&rwi.long$Year<2015,]
rwi.prism.abla<-rwi.prism[rwi.prism$Species=="AbLa",]
rwi.prism.abla$SA<-paste0(rwi.prism.abla$Site, rwi.prism.abla$Aspect)
Sites<-levels(as.factor(rwi.prism.abla$SA))
climhist$SA<-paste0(climhist$Site, climhist$Aspect)

clim<-NULL
tmean<-NULL
ppt<-NULL
cordat<-NULL
ppt.all<-NULL
tmean.all<-as.data.frame(NULL)

for(i in Sites){
rwi<-rwi.prism.abla[rwi.prism.abla$SA==i,]$rwi
rwi.abla<-as.data.frame(rwi)
rownames(rwi.abla)<-1895:2014

clim<-climhist[climhist$SA==i,][3:8][-5:-6]

cordat<-seascorr(rwi.abla, clim, season_lengths = 3)

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
