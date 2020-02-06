setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/")
area<-read.csv("LeafArea.All.csv")
weight<-read.csv("needleweights.csv")
clim.norms<-read.csv("clim.norms.csv")
gps<-read.csv("gps.csv")
gps<-gps[gps$Species!="Litter",]
str(gps)
str(area)
str(weight)

leafdat<-merge(area, weight, by=c("Site", "Aspect", "Zone", "Species","Number"))
str(leafdat)
leafdat<-leafdat[-15:-16]
str(gps)
leafdat<-merge(leafdat, gps, by=c("Site", "Aspect", "Zone", "Species","Number"))
str(leafdat)
leafdat<-merge(leafdat, clim.norms[-1], by=c("Site", "Aspect", "Zone"))
leafdat$SLA<-leafdat$wt/leafdat$HSA
str(leafdat)
anova(lm(SLA~Aspect*Site*Zone, leafdat))
summary(lm(SLA~Site*Zone+Species*PM, leafdat))
write.csv(leafdat, "/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/leafdat.csv")

library(ggplot2)
library(wesanderson)
leafdat$Site<-factor(leafdat$Site, levels=c("H", "Y", "P", "M", "T", "E", "D"))
ggplot(leafdat[leafdat$Species=="PiAl",], aes(SLA, x=Site))+geom_boxplot(data=leafdat, aes(fill=PM), position = position_dodge2(preserve = "single")) +
  geom_point(data=leafdat, aes(fill=PM), shape=21, position = position_jitterdodge())+
  ylab("SLA (g Leaf/ cm2 HSA)")+facet_grid(Zone~Aspect)+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling1"))

  

summary(lm(SLA~Zone+PM+Species+Zone:PM+Zone:PM:Species, leafdat))
summary(lm(SLA~Aspect*Species, leafdat))
plot(SLA~DBH, leafdat)
library(agricolae)
str(leafdat)
leafdat$AZ<-paste0(leafdat$Aspect, leafdat$Zone)
mod<-lm(SLA~AZ, leafdat[leafdat$Species=="TsMe",])
anova(mod)
tuk<-HSD.test(mod,trt="AZ", DFerror=69, MSerror=.0000052158,group=T )
tuk

t.test(SLA~Zone, leafdat[leafdat$Species=="TsMe"&leafdat$Aspect=="S",])

leaftex<-merge(leafdat, allsoil, by=c("Site", "Aspect", "Zone", "Number"))

summary(lm(SLA~C.N, leaftex))
