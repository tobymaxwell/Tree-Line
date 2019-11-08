#BAI All
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
RWI_C$Year<-as.numeric(RWI_C$Year)
rwi.all<-merge(rwi.all, RWI_C, all=T, by="Year")
}


str(rwi.all)
head(rwi.all)
library(tidyr)

rwi.long<-gather(rwi.all, key="ID", value = "rwi", DNBPiAl1:YSFTsMe8)
rwi.long<-na.omit(rwi.long)
str(rwi.long)
write.csv(rwi.long, "/Users/tobymaxwell/Desktop/rwi.long.csv")
rwi.long$Site<-substr(rwi.long$ID, 0,1)
rwi.long$Aspect<-substr(rwi.long$ID, 2,2)
rwi.long$Zone<-substr(rwi.long$ID, 3,3)
rwi.long$Species<-substr(rwi.long$ID, 4,7)
rwi.long[rwi.long$Species=="AbLA",]$Species<-"AbLa"
rwi.long[rwi.long$Species=="ABLa",]$Species<-"AbLa"
rwi.long[rwi.long$Species=="TsME",]$Species<-"TsMe"
rwi.long$number<-substr(rwi.long$ID, 8,11)
tail(rwi.long, 1000)
library(lubridate)
library(eeptools)

IDs<-levels(as.factor(rwi.long$ID))
minyear<-NULL
refyears<-NULL
for(i in IDs){
  minyear<-rbind(minyear, min(rwi.long[rwi.long$ID==i,]$Year))
  refyears<-rbind(refyears, paste0(i))
}
minyear<-data.frame(minyear, refyears)
colnames(minyear)<-c("minyear", "ID")
rwi.long<-merge(rwi.long, minyear, by="ID")
rwi.long$age<-(rwi.long$Year+1)-rwi.long$minyear

library(ggplot2)
ggplot(rwi.long[rwi.long$age==100,], aes(y=log(rwi), x=Year))+geom_point()
write.csv(rwi.long, "/Users/tobymaxwell/Desktop/rwi.long.csv")
str(rwi.long)
co2<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/spline_merged_ice_core_yearly.csv", skip=27)
rwi.long<-merge(rwi.long, co2, by="Year")
ggplot(rwi.long[rwi.long$age==10,], aes(y=rwi, x=Year, color = Species))+geom_point()
str(rwi.long)
summary(lm(rwi~co2, rwi.long))
library(lme4)
mer1<-lmer(rwi~Year+(1|ID), rwi.long)
sjt.lmer(mer1, p.kr=F)
summary(mer1)
#########################All Rings BAI#############################
library(dplyr)
library(dplR)
PO<-NULL
temp<-NULL
bai.all<-NULL
bai.all$Year<-1535:2018
bai.all<-as.data.frame(bai.all)
rownames(bai.all)<-bai.all$Year
cumbai<-NULL
cumbai$Year<-1534:2018
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
  
  temp2<-NULL
cumbai<-as.data.frame(cumbai)
  for(k in cols){
    temp2<-cumsum(na.omit(bai[k]))
    temp2$Year<-rownames(temp2)
    cumbai<-merge(cumbai,temp2,all=T)
  }
    
  
  bai$Year<-rownames(bai)
  bai.all<-merge(bai.all, bai, all=T, by="Year")
}
ggplot(bai.all[bai.all$Year>1900,], aes(lnbai, x=Year))+geom_line()+facet_grid(Aspect~Site)
head(bai.all)
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

########wide to long for cumbai########
cumbai.long<-gather(cumbai, key="ID", value = "bai", DNBPiAl1:YSFTsMe8)
cumbai.long<-na.omit(cumbai.long)
str(cumbai.long)
colnames(cumbai.long)<-c("Year", "ID", "cumbai")
cumbai.long$Site<-substr(cumbai.long$ID, 0,1)
cumbai.long$Aspect<-substr(cumbai.long$ID, 2,2)
cumbai.long$Zone<-substr(cumbai.long$ID, 3,3)
cumbai.long$Species<-substr(cumbai.long$ID, 4,7)
cumbai.long[cumbai.long$Species=="AbLA",]$Species<-"AbLa"
cumbai.long[cumbai.long$Species=="ABLa",]$Species<-"AbLa"
cumbai.long[cumbai.long$Species=="TsME",]$Species<-"TsMe"
cumbai.long$number<-substr(cumbai.long$ID, 8,11)
str(cumbai.long)
bai.all<-merge(bai.long, cumbai.long, by=c("Year", "ID", "Site", "Aspect", "Zone", "Species", "number"))
str(bai.all)
bai.all<-merge(bai.all, minyear, by="ID")
bai.all$age<-(bai.all$Year+1)-bai.all$minyear
bai.all$lnbai<-log(bai.all$bai)
bai.all$ln.cumbai<-log(bai.all$cumbai)


summary(lm(lnbai~ln.cumbai, bai.all[bai.all$ln.cumbai>-10&bai.all$lnbai>-10,]))
plot(lnbai~log(age), bai.all)
library(lme4)
bai.all$ln.age<-log(bai.all$age)
head(bai.all)


bai.all.cut<-bai.all[bai.all$ln.cumbai>-10&bai.all$lnbai>-10&bai.all$ln.age>-1,]
str(bai.all)
colnames(bai.all.cut)<-c("ID", "Year", "Site","Aspect","Zone", "Species", "Number", "bai", 'cumbai', 'minyear', 'age', 'lnbai', 'ln.cumbai', 'ln.age')
str(bai.all.cut)

CO2<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/spline_merged_ice_core_yearly.csv", skip=27)
str(bai.all.cut)
bai.all.cut<-merge(bai.all.cut, CO2, by="Year")
str(bai.all.cut)
summary(lm(C.pct~Depth, allsoil))
str(allsoil)
allsoil.ag<-na.omit(allsoil) %>%
  group_by(Site, Aspect, Zone, PM)%>%
  summarise(C.pct=mean(C.pct),
            N.pct=mean(N.pct),
            C.N=mean(C.N),
            Sand = mean(Sand),
            Silt = mean(Silt),
            Clay = mean(Clay),
            Cstock = mean(Cstock),
            Nstock = mean(Nstock))

allsoil.ag<-data.frame(allsoil.ag)
str(allsoil.ag)
str(bai.all.cut)
bai.soil<-merge(na.omit(bai.all.cut), allsoil.ag, by = c("Site", "Aspect", "Zone"))
str(bai.soil)
plot(lnbai~log(age),bai.soil)
summary(lm(lnbai~ln.age, bai.all.cut))
str(bai.soil)
str(leafdat)
leafdat.cut<-na.omit(leafdat) %>%
  select(-ID, -green.pixels, -red.pixels,-area.cm2, -PSA.HSAcoef, -HSA, -wt, -elevation, -longitude, -name)
str(leafdat.cut)

bai.sp<-merge(bai.soil, leafdat.cut, by = c("Site", "Aspect", "Zone", "Species", "Number", "PM"), all=T)
gis<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/landscapeGIS.csv")
str(gis)
gps.trees<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.trees.csv")
str(gps.trees)
gps.trees<-gps.trees%>%
  select(-X, -latitude, -longitude) %>%
  rename(tri.tree=tri, tpi.tree=tpi, slope.tree=slope, aspect.tree=aspect)
str(gps.trees)

summary(lm(tri.tree~Species, gps.trees))
str(gis)
gis.merge<-merge(gis, gps.trees, by=c("Site", "Aspect", "Zone"))
str(gis.merge)

str(CNlit)
CNlit.cut<-na.omit(CNlit)%>%
  group_by(Site, Aspect, Zone)%>%
    summarize(C.lit=mean(C.pct.L),
              N.lit=mean(N.pct.L),
              sdCN.lit=sd(C.N.L),
              CN.lit=mean(C.N.L))
CNlit.cut<-data.frame(CNlit.cut)  


bai.spg<-merge(bai.sp, gis.merge, by=c("Site", "Aspect", "Zone", "Species", "Number"))

bai.spgl<-merge(bai.spg, CNlit.cut, by=c("Site", "Aspect", "Zone"))
str(bai.spgl)
bai.spgl<-bai.spgl[!is.na(bai.spgl$cumbai),]
lai.int<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/LAI.int.csv")

bai.spgl<-merge(bai.spgl, lai.int[-1:-2], by=c("Site", "Aspect", "Zone"))
write.csv(bai.spgl, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.spgl.csv")

str(bai.spgl)
plot(mat~LAIint, bai.spgl)
plot(lnbai~log(age), bai.spgl)
summary(lm(lnbai~log(age), bai.spgl))
plot(lnbai~ln.cumbai, bai.spgl)
summary(lm(lnbai~ln.cumbai, bai.spgl))
plot(lnbai~Year, bai.spgl[bai.spgl$age==200,])
summary(lm(lnbai~Year, bai.spgl[bai.spgl$age==100,]))

summary(lm(C.pct~slope, bai.spgl))
mer1<-lmer(lnbai~ln.cumbai+(1|age)+(1|Year), bai.spgl)
sjt.lmer(mer1, p.kr=F)
plot(predict(mer1)~bai.spgl$lnbai)
bai.spgl$mer1resid<-residuals(mer1)
bai.spgl$mer1pred<-predict(mer1)
summary(lm(mer1pred~lnbai, bai.spgl))
str(bai.spgl)
bai.plot<-bai.spgl%>%
  group_by(Site, PM, Aspect, Zone)%>%
  summarise(sdSLA=sd(SLA),
            sdsilt=sd(Silt),
            sdsand=sd(Sand),
            sdclay=sd(Clay),
            sdcurv=sd(curv),
            sdtpi=sd(tpi.tree),
            sdtri=sd(tri.tree),
            sdC.N=sd(C.N),
            sdC.pct=sd(C.pct),
            sdCstock=sd(Cstock),
            sdNstock=sd(Nstock),
            sdCN.lit=sd(CN.lit),
            latitude=mean(latitude),
            mat=mean(mat),
            ppt = mean(ppt),
            tmin=mean(tmin),
            tmax=mean(tmax),
            vpdmax=mean(vpdmax),
            co2=mean(co2),
            SLA=mean(SLA),
            silt=mean(Silt),
            sand=mean(Sand),
            clay=mean(Clay),
            curv=mean(curv),
            tpi=mean(tpi.tree),
            tri=mean(tri.tree),
            slope=mean(slope),
            C.N=mean(C.N),
            C.pct=mean(C.pct),
            Cstock=mean(Cstock),
            Nstock=mean(Nstock),
            CN.lit=mean(CN.lit))

bai.plot<-as.data.frame(bai.plot)
str(bai.plot)
ggplot(bai.spgl, aes(mer1resid, x=Year))+geom_smooth(color='black')+theme_bw()+theme(legend.position='none')+
  geom_line(data=bai.spgl, aes(y=scale(co2), x=Year, color='red'))

ggplot(bai.spgl, aes(mer1resid, x=co2))+geom_smooth(color='red')+theme_classic()+theme(legend.position='none')
bai.spgl
str(bai.spgl)
bai.spgl<-do.call(data.frame,lapply(bai.spgl, function(x) replace(x, is.infinite(x),NA)))
pca<-prcomp(na.omit(bai.spgl[c(-1:-15, -33)]), center=T, scale.=T)
tail(pca$rotation)
autoplot(pca,data = na.omit(bai.spgl),
         x=3, y=2, colour ='Site',
         loadings = TRUE, loadings.colour = 'Red',
         loadings.label = TRUE, loadings.label.size = 5)+theme_bw()
bai.scaled<-bai.spgl
bai.scaled[c(-1:-15, -33)]<-scale(bai.scaled[c(-1:-15, -33)])
resmer<-lmer(mer1resid~co2+N.lit+C.pct+SLA+curv+vpdmax+slope+(1|Zone)+(1|Species)+(1|Site)+(1|Aspect), bai.scaled)
sjt.lmer(resmer, p.kr = F)
######Multiinference modeling #####
library(MuMIn)
str(bai.scaled)
dredge.dat<-na.omit(bai.spgl)%>%
  select(Aspect, Species, Zone, C.pct, N.lit, C.N, Sand, Silt, Clay, DBH, latitude, ppt, mat, vpdmax, SLA, tri.tree, tpi.tree, slope.tree, curv, elevation,CN.lit, mer1resid)
str(dredge.dat)
model<-lm(mer1resid~., dredge.dat) #the . after the ~ means include the whole dataset
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi<-dredge(model, trace=2, beta="sd", m.lim=c(1,6)) #find all models for mpg, trace setting makes it show progress bar
length(multi$logLik)
multi.cut<-get.models(multi, subset= delta<2.01) #delta is the difference in AIC from the best model, so subsetting the models into those within 2 of best model.
imp.multi<-importance(multi.cut) # importance is sum of akaike weights for each variable. Weights are 
imp.multi #shoes importance of variables from models in multi.cut
barplot(t(imp.multi)) #display importance values
multi.mod<-model.avg(multi.cut)
multi.mod
plot(mer1resid~predict(multi.mod), dredge.dat)
summary(lm(mer1resid~predict(multi.mod), dredge.dat))
summary(multi.mod) #check out model coefficients 

multimer<-lmer(mer1resid~C.N+curv+DBH+elevation+N.lit+SLA+slope.tree+tri.tree+PM+Zone*Aspect+(1|Site)+(1|Species), bai.scaled)
sjt.lmer(mer1, multimer, p.kr = F)


#######################All rings Raw ###############################
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

