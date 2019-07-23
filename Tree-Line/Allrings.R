#BAI All
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/")
files<-list.files()

files
library(plyr)
library(dplyr)
mergedrings<-ldply(files, read.csv)
str(mergedrings)

##########################all rings RWI ###########################
library(dplyr)
library(dplR)
PO<-NULL
temp<-NULL
rwi.all<-NULL
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
RWI_C<-data.frame(RWI_C)
RWI_C$Year<-rownames(RWI_C)
rwi.all<-merge(rwi.all, RWI_C, all=T, by="Year")
}
str(rwi.all)
tail(rwi.all)
library(tidyr)
str(rwi.long)
rwi.long<-gather(rwi.all, key="ID", value = "rwi", DNBPiAl1:YSFTsMe8)
rwi.long<-na.omit(rwi.long)
str(rwi.long)
rwi.long$Site<-substr(rwi.long$ID, 0,1)
rwi.long$Aspect<-substr(rwi.long$ID, 2,2)
rwi.long$Zone<-substr(rwi.long$ID, 3,3)
rwi.long$Species<-substr(rwi.long$ID, 4,7)
rwi.long[rwi.long$Species=="AbLA",]$Species<-"AbLa"
rwi.long[rwi.long$Species=="ABLa",]$Species<-"AbLa"
rwi.long[rwi.long$Species=="TsME",]$Species<-"TsMe"
rwi.long$number<-substr(rwi.long$ID, 8,11)
head(rwi.long, 1000)
ggplot(rwi.long, aes(y=rwi, x=Year))+geom_point()+facet_grid(~Species)

#########################All Rings BAI#############################
library(dplyr)
library(dplR)
PO<-NULL
temp<-NULL
bai.all<-NULL
bai.all$Year<-1535:2018
bai.all<-as.data.frame(bai.all)
rownames(bai.all)<-bai.all$Year
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
  rwl<-data.frame(temp)
  temp<-bai.out(rwl=temp)
  colnames(temp)<-cols
  bai<-data.frame(temp)
  bai$Year<-rownames(bai)
  bai.all<-merge(bai.all, bai, all=T, by="Year")
}
tail(bai.all)
str(bai.all)

library(tidyr)
bai.long<-gather(bai.all, key="ID", value = "bai", DNBPiAl1:YSFTsMe8)
bai.long<-na.omit(bai.long)
str(bai.long)
bai.long$Site<-substr(bai.long$ID, 0,1)
bai.long$Aspect<-substr(bai.long$ID, 2,2)
bai.long$Zone<-substr(bai.long$ID, 3,3)
bai.long$Species<-substr(bai.long$ID, 4,7)
bai.long[bai.long$Species=="AbLA",]$Species<-"AbLa"
bai.long[bai.long$Species=="ABLa",]$Species<-"AbLa"
bai.long[bai.long$Species=="TsME",]$Species<-"TsMe"
bai.long$number<-substr(bai.long$ID, 8,11)
head(bai.long, 1000)
ggplot(bai.long, aes(y=bai, x=Year))+geom_point()+facet_grid(Zone~Species)

#######################All rings Raw ################################
mergedata<-NA
mergedata$Year<-1535:2017
for(i in files[1:3]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2017:(2017-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata<-merge(mergedata, temp, by=c("Year"), all=T)
  
}
mergedata<-mergedata[-2]
head(mergedata)

mergedata2<-NULL
mergedata2$Year<-1535:2018
for(i in files[4:9]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata2<-merge(mergedata2, temp, by=c("Year"), all=T)
  
}
head(mergedata2)

mergedata3<-NULL
mergedata3$Year<-1535:2017
for(i in files[10:12]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2017:(2017-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata3<-merge(mergedata3, temp, by=c("Year"), all=T)
  
}

mergedata4<-NULL
mergedata4$Year<-1535:2018
for(i in files[13:25]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata4<-merge(mergedata4, temp, by=c("Year"), all=T)
  
}
head(mergedata4)
str(mergedata4)
mergedata5<-NULL
mergedata5$Year<-1535:2018
for(i in files[25:34]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata5<-merge(mergedata5, temp, by=c("Year"), all=T)
  
}
mergedata6<-NULL
mergedata6$Year<-1535:2018
for(i in files[35:39]){
  temp<-read.csv(i)
  n<-(length(temp[,1])/2)-1
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  temp$Year<-as.numeric(rep(2018:(2018-n), each=2))
  temp<-temp %>%
    group_by(Year) %>% 
    summarise_all(funs(sum))
  temp<-data.frame(temp)
  mergedata6<-merge(mergedata6, temp, by=c("Year"), all=T)
  
}
str(mergedata)
head(mergedata6, 25)
mergedata$DNBPiAl11
allrings<-merge(mergedata, mergedata2, by=c("Year"), all=T)
allrings<-merge(allrings, mergedata3, all=T)
allrings<-merge(allrings, mergedata4, all=T)
allrings<-merge(allrings, mergedata5, all=T)
allrings<-merge(allrings, mergedata6, all=T)
str(allrings)

library(tidyr)
rings.long<-gather(allrings, key="ID", value = "mm", PSFAbLa1:YSFTsMe8)
ringslong<-na.omit(rings.long)
str(ringslong)
