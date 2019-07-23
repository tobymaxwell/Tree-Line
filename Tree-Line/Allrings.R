#BAI All
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/")
files<-list.files()

files
library(plyr)
library(dplyr)
mergedrings<-ldply(files, read.csv)
str(mergedrings)

library(dplyr)
PO<-NULL
temp<-NULL
rwi.all<-NULL
rwi.all$Year<-1534:2018
rwi.all<-as.data.frame(rwi.all)
rownames(rwi.all)<-rwi.all$Year
for(i in files[1:2]){
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
rwi.all<-merge(RWI_C, rwi.all)
}


RWI_C<-cms(temp, PO, c.hat.t = FALSE, c.hat.i = FALSE)
RWI_C
rwi.stats(RWI_C)
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=20, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index Picea North Border", line = 3)
title(sub = "RWI - AR Model", line = -9.5, font.sub = 2)




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
