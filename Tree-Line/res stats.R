###########pointRes resiliance info############
library(pointRes)

library(dplyr)
library(dplR)
library(tidyr)

bai.spgl<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.spgl.csv")

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
  write.rwl(temp, paste0("/Users/tobymaxwell/Desktop/",i))
  temp<-read.rwl(paste0("/Users/tobymaxwell/Desktop/",i), header=TRUE)
  res<-res.comp(temp, nb.yrs=5, series.thresh = 50)$out
  ###, res.thresh.neg = 25
  
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
str(allsoil)
allsoil.ag<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/allsoil.ag.csv")[-1]
res.s<-merge(res.all, allsoil.ag, by=c("Site", "Aspect", "Zone"))
str(res.s)

anova(lm(rel.resil_mean~PM, res.s))
leafdat.cut<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/leafdat.cut.csv")[-1]

leafdat.ag<-leafdat.cut %>%
  group_by(Site, Aspect, Zone) %>%
  summarise(latitude=mean(latitude),
            mat=mean(mat),
            ppt = mean(ppt),
            tmin=mean(tmin),
            tmax=mean(tmax),
            vpdmax=mean(vpdmax),
            sdSLA=sd(SLA),
            SLA=mean(SLA))
leafdat.ag<-data.frame(leafdat.ag)

res.sp<-merge(res.s, leafdat.ag)
gps.trees<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.trees.csv")
str(gps.trees)
gps.trees<-gps.trees%>%
  select(-X, -latitude, -longitude) %>%
  rename(tri.tree=tri, tpi.tree=tpi, slope.tree=slope, aspect.tree=aspect)
str(gps.trees)
gps.trees.ag<-gps.trees%>%
  group_by(Site, Species)%>%
  summarise(curv=mean(curv),
            tri=mean(tri.tree),
            tpi=mean(tpi.tree),
            slope=mean(slope.tree))
gps.tree.ag<-as.data.frame(gps.trees.ag)
str(res.sp)
str(gps.trees.ag)
res.spg<-merge(res.sp, gps.trees.ag)
str(res.spg)

res.spglai<-merge(res.spg, lai.int[-1:-2], by=c("Site","Aspect", "Zone"))
str(res.spglai)
climhist.yr<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/climhist.yr.csv")[-1]
colnames(climhist.yr)<-c("Site", "Aspect", "year", "MAP.yr", "MAT.yr", "tmin.yr", "tmax.yr")
str(climhist.yr)
res.spglaic<-merge(res.spglai, climhist.yr, by=c("Site", "Aspect", "year"))
str(res.spglaic)

###### Resilience Stats ######
library(corrplot)
str(res.spglaic)
res.num<-res.spglaic%>%
  select(-nature, -nb.series, -perc.neg, -resist_sd,-recov_sd, -resil_sd, -rel.resil_sd)
str(res.num)
res.num<-data.frame(na.omit(res.num))
str(res.num)

dredge.dat<-na.omit(res.num)
str(dredge.dat)
dredge.dat<-na.omit(res.num)%>%
  select(-resist_mean, -recov_mean, -rel.resil_mean,-ID, -latitude, -year, -mat, -ppt,-tmin,-tmax, -sdSLA, -N.pct, -Cstock, -Nstock)
str(dredge.dat)
dredge.dat[c(5,7:22)]<-scale(dredge.dat[c(5, 7:22)])
str(dredge.dat)
model<-lm(resil_mean~., dredge.dat)
options(na.action = "na.fail") 
#options(na.action = "na.omit")
library(MuMIn)
multi<-dredge(model, trace=2, beta="sd", m.lim=c(1,6))
resil.multi.all<-multi
resil.multi.cut.all<-get.models(resil.multi.all, subset= delta<2.01)
resil.imp.multi.all<-importance(resil.multi.cut.all) # importance is sum of akaike weights for each variable.
resil.imp.multi.all #shows importance of variables from models in multi.cut
str(resil.imp.multi.all)
data.frame(resil.imp.multi.all)
barplot(t(resil.imp.multi.all))
resil.multi.mod.all<-model.avg(resil.multi.cut.all)
summary(resil.multi.mod.all)
plot(resil_mean~predict(resil.multi.mod.all), dredge.dat)
summary(lm(resil_mean~predict(resil.multi.mod.all), dredge.dat))
summary(resil.multi.mod.all) #check out model coefficients 

######Ordinations######
library(ggplot2)
str(res.num)
str(res.num[c(-1:-5, -10:-11,-37:-39)])
res.pca<-prcomp(res.num[c(-1:-5, -10:-11,-37:-39)], center=T, scale.=T)
res.num$PlotID<-paste0(res.num$Site, res.num$Aspect, res.num$Zone)
res.num$AspectID<-paste0(res.num$Site, res.num$Aspect)
res.num$ZoneID<-paste0(res.num$Site, res.num$Zone)
library(ggfortify)
autoplot(res.pca, data = res.num, x=1, y=2, colour = 'Site', shape='Species',
         loadings = TRUE, loadings.colour = 'black', frame=T,
         loadings.label = TRUE, loadings.label.size = 3)+theme_bw()
summary(lm(resist_mean~Clay*PM, res.num))
resil.mer<-lmer(resil_mean~scale(ppt)+scale(SLA)+scale(mat)+Zone+(1|Site)+(1|Species), res.num)
sjt.lmer(resil.mer)

str(res.spglaic)
res.sp.ag<-na.omit(res.spglaic)%>%
  group_by(Site,Aspect,Zone, PM, Species, year)%>%
  summarise(latitude=mean(latitude),
            resist=mean(resist_mean),
            recov=mean(recov_mean),
            resil=mean(resil_mean),
            rel.resil=mean(rel.resil_mean),
            MAT=mean(MAT.yr),
            MAP = mean(MAP.yr),
            tmin=mean(tmin.yr),
            tmax=mean(tmax.yr),
            vpdmax=mean(vpdmax),
            SLA=mean(SLA),
            silt=mean(Silt),
            sand=mean(Sand),
            clay=mean(Clay),
            curv=mean(curv),
            tpi=mean(tpi),
            tri=mean(tri),
            C.N=mean(C.N),
            C.pct=mean(C.pct),
            LAIint=mean(LAIint))
res.sp.ag<-data.frame(res.sp.ag)
str(res.sp.ag)

ggplot(res.sp.ag, aes(MAT, x=LAIint))+geom_point()
summary(lm(LAIint~MAT,res.sp.ag))

###
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
###
str(res.sp.ag)
res.cor<-cor(res.sp.ag[-1:-5])
str(res.sp.ag)
res1<-cor.mtest(res.cor, conf.level = .9)
corrplot(res.cor, method="color", sig.level = 0.1,  
         type="upper",
         p.mat = res1,
         diag=FALSE,
         
)

str(res.sp.ag)

autoplot(prcomp(res.sp.ag[c(-1:-5)], center=T, scale.=T), data = res.sp.ag,
         x=1, y=2, colour ='PM', shape='Species',
         loadings = TRUE, loadings.colour = 'black', frame=T,
         loadings.label = TRUE, loadings.label.size = 3)

library(MuMIn)
dredge.dat<-na.omit(res.sp.ag)
str(dredge.dat)
dredge.dat<-na.omit(res.sp.ag)%>%
  select(-resist, -recov, -rel.resil)
str(dredge.dat)
dredge.dat[7:23]<-scale(dredge.dat[7:23])
str(dredge.dat)
model<-lm(resil~., dredge.dat)
options(na.action = "na.fail") 
#options(na.action = "na.omit")
multi<-dredge(model, trace=2, beta="sd", m.lim=c(1,6))
resil.multi<-multi
resil.multi.cut<-get.models(resil.multi, subset= delta<2.01)
resil.imp.multi<-importance(resil.multi.cut) # importance is sum of akaike weights for each variable.
resil.imp.multi #shows importance of variables from models in multi.cut
str(resil.imp.multi)
data.frame(resil.imp.multi)
barplot(t(resil.imp.multi))
resil.multi.mod<-model.avg(resil.multi.cut)
summary(resil.multi.mod)
plot(resil~predict(resil.multi.mod), dredge.dat)
summary(lm(resil~predict(resil.multi.mod), dredge.dat))
summary(resil.multi.mod) #check out model coefficients 

resmod1<-lm(resil~scale(tpi)+scale(clay)+scale(vpdmax)+Site*scale(SLA), res.sp.ag) #no zone
anova(resmod1)
sjt.lm(resmod1)
summary(resmod1)

resmod2<-lm(resil~scale(mat)+scale(tmax)+scale(tmin)+scale(ppt)+scale(tpi)+scale(tri)+scale(tpi)+scale(sand)+scale(silt)+scale(curv)+Zone+Aspect, res.sp.ag)#yes zone
anova(resmod2)
sjt.lm(resmod2)
summary(resmod2)

dredge.dat<-na.omit(res.sp.ag)
str(dredge.dat)
dredge.dat<-na.omit(res.sp.ag)%>%
  select(-resist, -resil, -rel.resil, -latitude)
dredge.dat[7:22]<-scale(dredge.dat[7:22])
str(dredge.dat)
model<-lm(recov~., dredge.dat)
options(na.action = "na.fail") 
#options(na.action = "na.omit")
multi<-dredge(model, trace=2, beta="sd", m.lim=c(1,6))
rec.multi<-multi
rec.multi.cut<-get.models(multi, subset= delta<2.01)
rec.imp.multi<-importance(rec.multi.cut) # importance is sum of akaike weights for each variable.
rec.imp.multi #shows importance of variables from models in multi.cut
barplot(t(rec.imp.multi))
rec.multi.mod<-model.avg(rec.multi.cut)
rec.multi.mod
plot(recov~predict(rec.multi.mod), dredge.dat)
summary(lm(recov~predict(rec.multi.mod), dredge.dat))
summary(rec.multi.mod) #check out model coefficients 

recmod1<-lm(recov~scale(SLA)+scale(C.pct)+scale(vpdmax)+scale(tpi)+Aspect, res.sp.ag)
anova(recmod1)
summary(recmod1)

res.sp.ag$predres1<-predict(resmod1)
ggplot(res.sp.ag, aes(resil, x=predres1))+geom_point()
resmer<-lmer(resil~scale(tpi)+scale(vpdmax)+scale(clay)+scale(SLA)+(scale(SLA)||Site), res.sp.ag)
sjt.lmer(resmer)
res.sp.ag$pred<-predict(resmer)
plot(resil_mean~pred, res.all)
#######dredge resist
dredge.dat<-na.omit(res.sp.ag)
str(dredge.dat)
dredge.dat<-na.omit(res.sp.ag)%>%
  select(-resil, -recov, -rel.resil, -latitude)
str(dredge.dat)
dredge.dat[7:22]<-scale(dredge.dat[7:22])
str(dredge.dat)
model<-lm(resist~., dredge.dat)
options(na.action = "na.fail") 
#options(na.action = "na.omit")
multi<-dredge(model, trace=2, beta="sd", m.lim=c(1,6))
resist.multi<-multi
resist.multi.cut<-get.models(resist.multi, subset= delta<2.01)
resist.imp.multi<-importance(resist.multi.cut) # importance is sum of akaike weights for each variable.
resist.imp.multi #shows importance of variables from models in multi.cut
str(resist.imp.multi)
data.frame(resist.imp.multi)
barplot(t(resist.imp.multi))
resist.multi.mod<-model.avg(resist.multi.cut)
summary(resist.multi.mod)
plot(resist~predict(resist.multi.mod), dredge.dat)
summary(lm(resist~predict(resist.multi.mod), dredge.dat))


summary(resil.multi.mod) #check out model coefficients 
summary(resist.multi.mod)
summary(rec.multi.mod)

levels(point.all$Site)

ptrg<-pointer.rgc(temp, series.thresh=50) #pointer years and relative growth changers
rgc.plot(ptrg)

ptwin<-pointer.norm(temp)#pointer years and cropper values
tail(ptwin$out,200)
norm.plot(ptwin)
event.plot(ptrg)
pointer.plot(ptrg)

res<-res.comp(temp, nb.yrs=5, series.thresh = 50)
res
res.plot(res)