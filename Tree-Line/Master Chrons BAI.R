#BAI All
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Final Chrons/")
files<-list.files()

files
library(plyr)
library(dplyr)
library(sjPlot)
mergedrings<-ldply(files, read.csv)
str(mergedrings)

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

#write.csv(bai.long, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.long.csv")
library(ggplot2)
library(wesanderson)
ggplot(bai.long, aes(bai,x=Year, color=Species, shape=Zone))+geom_point()+
  facet_grid(Aspect~Site, scales="free")+
  scale_color_manual(values=wes_palette(name="GrandBudapest1"))+
  theme_bw()


bai.long<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.long.csv")[-1]
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

######all data breaks#####
d<-na.omit(bai.long)%>%
  group_by(Year)%>%
  summarize(sdbai=sd(bai),n=length(bai),bai=mean(bai))
d<-as.data.frame(d)
d.20th<-as.data.frame(d[d$Year>1916,])
ts.d<-ts(d.20th$bai, start=min(d.20th$Year))
slopes<-ts.slopes(ts.d)
print(slopes)


ggplot(data=bai.long,aes(y=bai, x=Year))+geom_point()+
  geom_ribbon(data=d, aes(ymin=bai-sdbai,
                          ymax=bai+sdbai))+
  geom_line(data=d, color='Red',aes(y=bai, x=Year))+
  theme_bw()

ggplot(data=d, color='Red',aes(y=bai, x=Year))+
  geom_ribbon(data=d, aes(ymin=bai-sdbai,
                          ymax=bai+sdbai))+
  geom_line(color='Red')+
  geom_segment(aes(x=1536, xend=1757, y=0.8819282*1536-1349.8, yend=0.8819282*1757-1349.8), color='Cyan')+
  geom_segment(aes(x=1758, xend=1916, y=0.4577999*1758-635.4, yend=0.4577999*1916-635.4), color="Cyan")+
  geom_segment(aes(x=1917, xend=1952, y=4.660652*1917-8740.44, yend=4.660652*1952-8740.44), color="Cyan")+
  geom_segment(aes(x=1953, xend=2003, y=3.3834814*1953-6275.22, yend=3.3834814*2003-6275), color='Cyan')+
  geom_segment(aes(x=2004, xend=2018, y=3.029650*2004-5502, yend=3.029650*2018-5502), color='cyan')+
  theme_bw()

######soil/plant/climatedata aggregation#

clim.norms<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/clim.norms.csv")
clim.norms<-clim.norms%>%
  group_by(Site, Aspect, Zone)%>%
  summarize(mat=mean(mat), ppt=mean(ppt), tmin=mean(tmin), tmax=mean(tmax), vpdmax=mean(vpdmax))
clim<-as.data.frame(clim.norms)

gps.trees<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/gps.trees.csv")[-1]
gis.SA<-gps.trees%>%
  group_by(Site, Aspect, Zone, Species)%>%
  summarize(tri=mean(tri), tpi=mean(tpi), slope=mean(slope), curv=mean(curv))
climtop<-merge(clim, gis.SA, by=c("Site", "Aspect", "Zone"))

soillit<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/soillit.csv")
allsoil.SA<-na.omit(soillit)%>%
  group_by(Site, Aspect, Zone)%>%
  summarize(Sand=mean(Sand), Silt=mean(Silt), Clay=mean(Clay), C.N=mean(C.N), C.pct=mean(C.pct),C.N.L=mean(C.N.L))
climtopsoil<-merge(climtop,allsoil.SA)

leafdat<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/leafdat.csv")
leafdat.cut<-leafdat%>%
  group_by(Site, Aspect, Zone, Species)%>%
  summarize(DBH=mean(DBH),SLA=mean(SLA), latitude=mean(latitude))
data.frame(leafdat.cut)
Slopedata<-merge(climtopsoil, leafdat.cut, all=T)
str(Slopedata)
Slopedata[Slopedata$Site=='H',]
#####Slope cahracterizations#####
bai.long$ID<-substr(bai.long$ID, 1, 7)
bai.long[bai.long$ID=="TNBAbLA",]$ID<-"TNBAbLa"


######2004-2018#####
levels(factor(bai.long$ID))
str(bai.long)
coefs<-NULL
res<-NULL
sites<-levels(factor(bai.long$ID))
for(i in sites){
  res<-summary(lm(log(bai)~Year, bai.long[is.finite(bai.long$bai)&bai.long$Year>2003&bai.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.2004<-as.data.frame(coefs)
coefs.2004$ID<-sites
coefs.2004$sig.1<-ifelse(coefs.2004$pyear<0.1, yes=1, no=0)
coefs.2004$sig.05<-ifelse(coefs.2004$pyear<0.05, yes=1, no=0)
coefs.2004$dir<-ifelse(coefs.2004$Year>0, yes=1, no=0)
coefs.2004$Site<-substr(coefs.2004$ID,1,1)
coefs.2004$Aspect<-substr(coefs.2004$ID,2,2)
coefs.2004$Zone<-substr(coefs.2004$ID,3,3)
coefs.2004$Species<-substr(coefs.2004$ID,4,7)
coefs.2004[coefs.2004$sig.1==1,]
length(coefs.2004[coefs.2004$sig.05==1,]$dir)

coefs.2004$ParentMaterial<-NA
coefs.2004[coefs.2004$Site=='H'|coefs.2004$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.2004[coefs.2004$Site=='E'|coefs.2004$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.2004[coefs.2004$Site=='M'|coefs.2004$Site=='P',]$ParentMaterial<-'Granite'
coefs.2004[coefs.2004$Site=='D',]$ParentMaterial<-'Mixed Volvanics'
coefs.2004.dat<-merge(coefs.2004, Slopedata, by=c("Site", "Aspect", "Zone", "Species"))
coefs.2004.dat<-coefs.2004.dat[coefs.2004.dat$sig.05==1&is.finite(coefs.2004.dat$Year),]
ggplot(bai.long, aes(bai, x=Year))+geom_blank()+
  geom_abline(data=coefs.2004.dat, size=1,aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.2004[coefs.2004$sig.05==0,], color='grey',size=.5,aes(slope=Year, intercept=int))+
  xlim(2004,2018)+ylim(0,8)+theme_bw()
library(lme4)
mer<-lmer(Year~(1|Species), coefs.2004.dat)
library(MuMIn)
r.squaredGLMM(mer)
#ppt+vpdmax+tri+slope+Clay+C.N.L
str(coefs.2004.dat)
M<-cor(coefs.2004.dat[c(-1:-4,-7:-15, -31)])
res1 <- cor.mtest(coefs.2004.dat[c(-1:-4,-7:-15, -31)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

#####1953-2003#####

coefs<-NULL
res<-NULL
sites<-levels(factor(bai.long[is.finite(bai.long$bai)&bai.long$bai>0&bai.long$Year<2004&bai.long$Year>1952,]$ID))
for(i in sites){
  res<-summary(lm(log(bai)~Year, bai.long[is.finite(bai.long$bai)&bai.long$bai>0&bai.long$Year<2004&bai.long$Year>1952&bai.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1952<-as.data.frame(coefs)
coefs.1952$ID<-sites
coefs.1952$sig.05<-ifelse(coefs.1952$pyear<0.05, yes=1, no=0)
coefs.1952$dir<-ifelse(coefs.1952$Year>0, yes=1, no=0)
coefs.1952$Site<-substr(coefs.1952$ID,1,1)
coefs.1952$Aspect<-substr(coefs.1952$ID,2,2)
coefs.1952$Zone<-substr(coefs.1952$ID,3,3)
coefs.1952$Species<-substr(coefs.1952$ID,4,7)
length(coefs.1952[coefs.1952$sig.05==1,]$Year)

coefs.1952$ParentMaterial<-NA
coefs.1952[coefs.1952$Site=='H'|coefs.1952$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1952[coefs.1952$Site=='E'|coefs.1952$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1952[coefs.1952$Site=='M'|coefs.1952$Site=='P',]$ParentMaterial<-'Granite'
coefs.1952[coefs.1952$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1952.dat<-merge(coefs.1952[coefs.1952$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone"))
ggplot(bai.long, aes(log(bai), x=Year))+geom_blank()+
  geom_abline(data=coefs.1952.dat, size=1,aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.1952[coefs.1952$sig.05==-0,], color='grey',size=.5,aes(slope=Year, intercept=int))+
  xlim(1952,2003)+ylim(0,8)+theme_bw()

anova(lm(Year~curv, coefs.1952.dat))

str(coefs.1952.dat)

M<-cor(coefs.1952.dat[c(-1:-4,-7:-14, -30)])
res1 <- cor.mtest(coefs.1952.dat[c(-1:-4,-7:-14, -30)], conf.level = .90)

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
#sites.1878<-levels(as.factor(rwi.long[is.finite(rwi.long$rwi)&rwi.long$Year>1840&rwi.long$Year<1891,]$ID))
for(i in sites){
  res<-summary(lm(log(bai)~Year, bai.long[is.finite(bai.long$bai)&bai.long$bai>0&bai.long$Year>1912&bai.long$Year<1952&bai.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1918<-as.data.frame(coefs)
coefs.1918$ID<-sites
coefs.1918$sig.05<-ifelse(coefs.1918$pyear<0.05, yes=1, no=0)
coefs.1918$dir<-ifelse(coefs.1918$Year>0, yes=1, no=-1)
coefs.1918$Site<-substr(coefs.1918$ID,1,1)
coefs.1918$Aspect<-substr(coefs.1918$ID,2,2)
coefs.1918$Zone<-substr(coefs.1918$ID,3,3)
coefs.1918$Species<-substr(coefs.1918$ID,4,7)
length(coefs.1918[coefs.1918$sig.05==1,]$Year)

coefs.1918$ParentMaterial<-NA
coefs.1918[coefs.1918$Site=='H'|coefs.1918$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1918[coefs.1918$Site=='E'|coefs.1918$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1918[coefs.1918$Site=='M'|coefs.1918$Site=='P',]$ParentMaterial<-'Granite'
coefs.1918[coefs.1918$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1918.dat<-merge(coefs.1918[coefs.1918$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone"))
ggplot(bai.long, aes(log(bai), x=Year))+geom_blank()+
  geom_abline(data=coefs.1918.dat, size=1, aes(slope=Year, intercept=int, color=factor(dir)))+
  geom_abline(data=coefs.1918[coefs.1918$sig.05==-0,], color='grey',size=.5,aes(slope=Year, intercept=int))+
  xlim(1918,1952)+ylim(0,8)+theme_bw()

summary(lm(Year~tmin, coefs.1918.dat))


M<-cor(coefs.1918.dat[c(-1:-4,-10:-14, -30)])
res1 <- cor.mtest(coefs.1918.dat[c(-1:-4,-10:-14, -30)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


#########1758-1916######

coefs<-NULL
res<-NULL
sites.1758<-levels(as.factor(bai.long[is.finite(bai.long$bai)&bai.long$Year>1757&bai.long$Year<1917,]$ID))
for(i in sites.1758){
  res<-summary(lm(log(bai)~Year, bai.long[is.finite(bai.long$bai)&bai.long$bai>0&bai.long$Year>1757&bai.long$Year<1917&bai.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1758<-as.data.frame(coefs)
coefs.1758$ID<-sites.1758
coefs.1758$sig.05<-ifelse(coefs.1758$pyear<0.05, yes=1, no=0)
coefs.1758$dir<-ifelse(coefs.1758$Year>0, yes=1, no=-1)
coefs.1758$Site<-substr(coefs.1758$ID,1,1)
coefs.1758$Aspect<-substr(coefs.1758$ID,2,2)
coefs.1758$Zone<-substr(coefs.1758$ID,3,3)
coefs.1758$Species<-substr(coefs.1758$ID,4,7)
length(coefs.1758[coefs.1758$sig.05==1&coefs.1758$dir==-1,]$Year)

coefs.1758$ParentMaterial<-NA
coefs.1758[coefs.1758$Site=='H'|coefs.1758$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1758[coefs.1758$Site=='E'|coefs.1758$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1758[coefs.1758$Site=='M'|coefs.1758$Site=='P',]$ParentMaterial<-'Granite'
coefs.1758[coefs.1758$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1758.dat<-merge(coefs.1758[coefs.1758$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone"))
ggplot(bai.long, aes(log(bai), x=Year))+geom_blank()+
  geom_abline(data=coefs.1758.dat, aes(slope=Year, intercept=int, color=tmax))+
  xlim(1758,1918)+theme_bw()

summary(lm(Year~Species, coefs.1758.dat))


M<-cor(coefs.1758.dat[c(-1:-4,-10:-14, -30)])
res1 <- cor.mtest(coefs.1758.dat[c(-1:-4,-10:-14, -30)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

#########1917-2018######

coefs<-NULL
res<-NULL
for(i in sites){
  res<-summary(lm(log(bai)~Year, bai.long[is.finite(bai.long$bai)&bai.long$bai>0&bai.long$Year>1916&bai.long$ID==i,]))
  coefs<-rbind(coefs,res$coefficients[-5:-6])
  
}
colnames(coefs)<-c("int", "Year", "seint", "seyear", "pint", "pyear")
coefs.1917<-as.data.frame(coefs)
coefs.1917$ID<-sites
coefs.1917$sig.05<-ifelse(coefs.1917$pyear<0.05, yes=1, no=0)
coefs.1917$dir<-ifelse(coefs.1917$Year>0, yes=1, no=-1)
coefs.1917$Site<-substr(coefs.1917$ID,1,1)
coefs.1917$Aspect<-substr(coefs.1917$ID,2,2)
coefs.1917$Zone<-substr(coefs.1917$ID,3,3)
coefs.1917$Species<-substr(coefs.1917$ID,4,7)
length(coefs.1917[coefs.1917$sig.05==1&coefs.1917$dir==1,]$Year)

coefs.1917$ParentMaterial<-NA
coefs.1917[coefs.1917$Site=='H'|coefs.1917$Site=='Y',]$ParentMaterial<-'Basalt'
coefs.1917[coefs.1917$Site=='E'|coefs.1917$Site=='T',]$ParentMaterial<-'MetaSedimentary'
coefs.1917[coefs.1917$Site=='M'|coefs.1917$Site=='P',]$ParentMaterial<-'Granite'
coefs.1917[coefs.1917$Site=='D',]$ParentMaterial<-'Mixed Volcanics'
coefs.1917.dat<-merge(coefs.1917[coefs.1917$sig.05==1,], Slopedata, by=c("Site", "Aspect", "Zone"))
ggplot(bai.long, aes(log(bai), x=Year))+geom_blank()+
  geom_abline(data=coefs.1917.dat, aes(slope=Year, intercept=int, color=tmax))+
  xlim(1917,1918)+ylim(-1,6)+theme_bw()

summary(lm(Year~curv, coefs.1917.dat))


M<-cor(coefs.1917.dat[c(-1:-4,-10:-14, -30)])
res1 <- cor.mtest(coefs.1917.dat[c(-1:-4,-10:-14, -30)], conf.level = .90)

corrplot(M, method="color", sig.level = 0.1,  
         type="upper",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = res1, 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)




#######Cor.mtest function######

library(corrplot)
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



######multimodel inference ######
library(MuMIn)
library(lme4)
library(dplyr)
bai.clim<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/bai.clim.csv")
str(bai.clim)
#bai.clim$aspect.tree
#bai.clim$northness<-cos(bai.clim$aspect.tree*pi/180)
dredge.dat<-na.omit(bai.clim)%>%
  select(ID, Site, Aspect, Species, PM, Zone, Year, lnbai, age, cumbai, co2, C.pct, N.lit, C.N, Clay, latitude, ppt, mat, tmean.yr,Summer_vpdmax, Winter_ppt, Winter_tmean, ppt.yr, ppt.lag1, ppt.lag2, ppt.lag3, cvppt, SLA, tri.tree, tpi.tree, slope.tree, curv, elevation, CN.lit, northness)
str(dredge.dat)
dredge.dat[c(-1:-7)]<-scale(dredge.dat[c(-1:-7)])
str(dredge.dat)
dredge.dat<-dredge.dat[is.finite(dredge.dat$lnbai),]
dredge.2004<-dredge.dat[dredge.dat$Year>2003,]
globmod<-lmer(lnbai~C.pct+C.N+CN.lit+Clay+curv+SLA+ppt+cvppt+ppt.yr+ppt.lag1+mat+tmean.yr+Aspect+Zone+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.2004) #the . after the ~ means include the whole dataset
summary(globmod)
r.squaredGLMM(globmod)
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi<-dredge(globmod, trace=2, beta="sd") #find all models for mpg, trace setting makes it show progress bar
length(multi$logLik)
multi.cut<-get.models(multi, subset= delta<6) #delta is the difference in AIC from the best model, so subsetting the models into those within 2 of best model.
imp.multi<-importance(multi.cut) # importance is sum of akaike weights for each variable. Weights are 
imp.multi #shows importance of variables from models in multi.cut
barplot(t(imp.multi)) #display importance values
multi.mod<-model.avg(multi.cut)
multi.mod
summary(lm(lnbai~predict(multi.mod), dredge.2004))
summary(multi.mod) #check out model coefficients 
r.squaredGLMM(multi.mod)

#####m2004 2.0#####

dredge.2004$Zone<-factor(dredge.2004$Zone, levels=c("F", "B"))
globmod<-lmer(lnbai~0+C.N+C.pct+mat+northness+ppt+Zone+slope.tree+PM+tpi.tree+Zone+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.2004) #the . after the ~ means include the whole dataset
summary(globmod)
r.squaredGLMM(globmod)
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi.2<-dredge(globmod, trace=2, beta="sd") #find all models for mpg, trace setting makes it show progress bar

multi.cut.2<-get.models(multi.2, subset= delta<6)
imp.multi.2<-importance(multi.cut.2) # importance is sum of akaike weights for each variable. Weights are 
multi.mod.2<-model.avg(multi.cut.2)
summary(lm(lnbai~predict(multi.mod.2), dredge.2004))
summary(multi.mod.2)

####1953-2003#####

dredge.1953<-dredge.dat[dredge.dat$Year<2003&dredge.dat$Year>1952,]
globmod<-lmer(lnbai~C.pct+C.N+CN.lit+Clay+curv+SLA+ppt+cvppt+ppt.yr+ppt.lag1+mat+tmean.yr+Aspect+Zone+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1953) #the . after the ~ means include the whole dataset
summary(globmod)
r.squaredGLMM(globmod)
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi.1953<-dredge(globmod, trace=2, beta="sd") #find all models, trace setting makes it show progress bar
length(multi$logLik)
multi.cut.1953<-get.models(multi.1953, subset= delta<6) #delta is the difference in AIC from the best model, so subsetting the models into those within 2 of best model.
imp.multi.1953<-importance(multi.cut.1953) # importance is sum of akaike weights for each variable. Weights are 
imp.multi.1953 #shows importance of variables from models in multi.cut
barplot(t(imp.multi.1953)) #display importance values
multi.mod.1953<-model.avg(multi.cut.1953)
multi.mod.1953
summary(lm(lnbai~predict(multi.mod.1953), dredge.1953))
summary(multi.mod.1953) #check out model coefficients 
r.squaredGLMM(multi.mod.1953)

#####1953mod 2.0#####
globmod<-lmer(lnbai~0+C.N+CN.lit+mat+ppt+SLA+northness+Winter_ppt+Winter_tmean+slope.tree+tpi.tree+PM+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1953)
summary(globmod)
r.squaredGLMM(globmod)
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi.1953.2<-dredge(globmod, trace=2, beta="sd") #find all models, trace setting makes it show progress bar
length(multi.1953.2$logLik)
multi.cut.1953.2<-get.models(multi.1953.2, subset= delta<6) #delta is the difference in AIC from the best model, so subsetting the models into those within 2 of best model.
imp.multi.1953.2<-importance(multi.cut.1953.2) # importance is sum of akaike weights for each variable. Weights are 
imp.multi.1953.2 #shows importance of variables from models in multi.cut
barplot(t(imp.multi.1953.2)) #display importance values
multi.mod.1953.2<-model.avg(multi.cut.1953.2)
multi.mod.1953.2
summary(lm(lnbai~predict(multi.mod.1953.2), dredge.1953))
summary(multi.mod.1953.2)



testmod<-lm(lnbai~0+C.N+CN.lit+mat+ppt+SLA+northness+Winter_ppt+Winter_tmean+slope.tree+tpi.tree+PM, dredge.1953)
testdredge<-dredge(testmod, trace=2, beta="sd")
testcut<-get.models(testdredge, subset= delta<6)
testavg<-model.avg(testcut)
summary(testavg)

dredge.1953%>%
  summarize(C.N=mean(C.N))

summary(lm(lnbai~C.N, dredge.1953))


####1917-1952#####

dredge.1917<-dredge.dat[dredge.dat$Year<1953&dredge.dat$Year>1916,]
globmod<-lmer(lnbai~C.pct+C.N+CN.lit+Clay+curv+SLA+ppt+cvppt+ppt.yr+ppt.lag1+mat+tmean.yr+Aspect+Zone+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1917) #the . after the ~ means include the whole dataset
summary(globmod)
r.squaredGLMM(globmod)
options(na.action = "na.fail") #needs this option set to run model tool
#options(na.action = "na.omit") #needs this option set to run model tool
multi.1917<-dredge(globmod, trace=2, beta="sd") #find all models, trace setting makes it show progress bar
length(multi.1917$logLik)
multi.cut.1917<-get.models(multi.1917, subset= delta<6) #delta is the difference in AIC from the best model, so subsetting the models into those within 2 of best model.
imp.multi.1917<-importance(multi.cut.1917) # importance is sum of akaike weights for each variable. Weights are 
imp.multi.1917 #shows importance of variables from models in multi.cut
barplot(t(imp.multi.1917)) #display importance values
multi.mod.1917<-model.avg(multi.cut.1917)
multi.mod.1917
summary(lm(lnbai~predict(multi.mod.1917), dredge.1917))
summary(multi.mod.1917) #check out model coefficients 
r.squaredGLMM(multi.mod.1917)


###dredge 1917 2.0####
dredge.1917<-dredge.dat[dredge.dat$Year<1953&dredge.dat$Year>1916,]
globmod<-lmer(lnbai~C.N+Clay+ppt+CN.lit+ppt.yr+tmean.yr+ppt.lag1+slope.tree+Summer_vpdmax+Winter_tmean+Zone+PM+tpi.tree+northness+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1917)
summary(globmod)
r.squaredGLMM(globmod)

multi.1917.2<-dredge(globmod, trace=2, beta="sd") #find all models, trace setting makes it show progress bar
length(multi.1917.2$logLik)
multi.cut.1917.2<-get.models(multi.1917.2, subset= delta<6)
imp.multi.1917.2<-importance(multi.cut.1917.2) 
barplot(t(imp.multi.1917.2)) #display importance values
multi.mod.1917.2<-model.avg(multi.cut.1917.2)
summary(lm(lnbai~predict(multi.mod.1917.2), dredge.1917))
summary(multi.mod.1917.2) #check out model coefficients

m2004<-lmer(lnbai~Aspect+C.N+C.pct+mat+ppt+Zone+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.2004)
m2004.2<-lmer(lnbai~C.N+C.pct+mat+northness+ppt+Zone+slope.tree+PM+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.2004)
m1953<-lmer(lnbai~Aspect+C.N+CN.lit+curv+mat+ppt+SLA+C.pct+Clay+tmean.yr+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1953)
m1953.2<-lmer(lnbai~Aspect+C.N+CN.lit+mat+ppt+SLA+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1953)
m1917<-lmer(lnbai~C.N+Clay+ppt+Zone+Aspect+CN.lit+ppt.yr+tmean.yr+ppt.lag1+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1917)
m1917.2<-lmer(lnbai~C.N+Clay+ppt+(1|Site/Species)+(1|cumbai)+(1|Year), dredge.1917)
r.squaredGLMM(m2004)
r.squaredGLMM(m2004.2)
library(sjPlot)
sjt.lmer(m2004, m1953.2, m1917, p.kr=F)
