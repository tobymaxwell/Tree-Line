library(ggplot2)
setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/")
soil<-read.csv("Soil Data.csv")
CN<-read.csv("CNsoil.csv")
CNlit<-read.csv("CNlit.csv")
texture<-read.csv("texture/Final Data.csv")
clim.norms<-read.csv("clim.norms.csv")

str(soil)
tail(soil)
allsoil<-merge(CN, texture, by=c("Site", "Aspect", "Zone", "Number", "Depth", "PM"), all=T)
str(allsoil)
soil<-soil%>%
  select(Site, Aspect, Zone, Number, Depth, PM, soildens_gcm3, litdens_sqin)
str(soil)
tail(allsoil,30)
allsoil<-merge(allsoil, soil, by =c("Site", "Aspect", "Zone", "Number", "Depth", "PM"), all=T)
str(texture)
allsoil$AZ<-paste0(allsoil$Aspect, allsoil$Zone)
write.csv(allsoil, "/Users/tobymaxwell/Desktop/allsoils.csv")
allsoil$SiltClay<-allsoil$Silt+allsoil$Clay
allsoil$SandSilt<-allsoil$Sand+allsoil$Silt
ggplot(allsoil, aes(y=C.pct, x=Sand, col=Site))+geom_point()
summary(lm(C.pct~Sand, allsoil))
#to /100 turns to g/g, then multiply by 1000g/kg for units of g/kg
allsoil$C_g_kg<-((allsoil$C.pct/100)*1000) 
allsoil$N_g_kg<-((allsoil$N.pct/100)*1000)

#calculate carbon per area by multiplying C (g/kg)*BD (g/cm3)*depth(6cm)
allsoil$C_bd<-(allsoil$C_g_kg/1000)*allsoil$soildens_gcm3*6
allsoil$N_bd<-(allsoil$N_g_kg/1000)*allsoil$soildens_gcm3*6

#carbon stocks
allsoil$Cstock<-allsoil$C_bd*100000000*(1/1000000)
allsoil$Nstock<-allsoil$N_bd*100000000*(1/1000000)
m1<-lmer(Cstock~Clay+(1|PM)+(1|Aspect)+(1|Zone), allsoil)
sjt.lmer(m1, p.kr=F)
ggplot(allsoil, aes(y=Cstock, x=Silt, col=Site))+geom_point()
CNlit$AZ<-paste0(CNlit$Aspect, CNlit$Zone)
library(ggplot2)
library(wesanderson)
allsoil$Site<-factor(allsoil$Site, levels=c("D","H", "Y", "P", "M", "T", "E"))
ggplot(allsoil, aes(y=Nstock, x=Site, fill = PM))+geom_boxplot()+
  theme_bw()+
  ylab("Carbon Stocks (Mg/Ha)")+ xlab("")+
  scale_fill_manual(values=wes_palette(name="Darjeeling1"))+
  facet_grid(Zone~Aspect)

anova((lm(C.N~CN.lit, allsoil)))

sjt.lmer(lmer(Cstock~Aspect*PM+(1|Site), allsoil), p.kr=F)

ggplot(allsoil, aes(y=N.pct, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Nitrogen (g/kg)")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  facet_wrap(~Depth)+
  theme(legend.position="none")+
  xlab("")

ggplot(allsoil[allsoil$Site=="H"|allsoil$Site=="T"|allsoil$Site=="P",], aes(y=C_bd, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Carbon (g/cm^3)")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  theme(legend.position="none")+
  xlab("")

ggplot(allsoil[allsoil$Site=="H"|allsoil$Site=="T"|allsoil$Site=="P",], aes(y=C.N, x=Site, fill=AZ))+geom_boxplot()+theme_bw()+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+ylab("C:N")+
  facet_wrap(~Depth)+
  theme(legend.position="none")+
  xlab("")

ggplot(CNlit, aes(y=C.N, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Litter C:N")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  theme(legend.position="none")+
  xlab("")

ggplot(CNlit, aes(y=N.pct, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Litter N")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  theme(legend.position="none")+
  xlab("")

mod<-lm(C.N~Site, CN)
anova(mod)
m2<-lm(C.N~Site*Zone, allsoil[allsoil$Depth=="six",])
anova(m2)

ggplot(CN, aes(y=C.N, x=Site, color=Aspect))+geom_boxplot()+facet_grid(Zone~Depth)+theme_bw()

library(agricolae)
mod<-lm(litdens_sqin~PM, soil[soil$Depth=="tw",])
anova(mod)
tuk<-HSD.test(mod,trt="PM", DFerror=117, MSerror=3.7605,group=T )
tuk

mod<-lm(C.N~Site, CN[CN$Type=="Soil"&CN$Depth=="tw",])
anova(mod)
tuk<-HSD.test(mod,trt="Site", DFerror=57, MSerror=18.79,group=T )
tuk

t.test(C.N~Zone, CN[CN$Type=="Soil"&CN$Depth=="tw"&CN$Site=="T",])

soil$Type
m1<-lm(rockwt~Zone*PM, soil[soil$Depth=="tw",])
anova(m1)

summary(lm(C.N~Site,CN[CN$Type=="Soil"&CN$Depth=="six",] ))

#################Texture#######################
str(texture)
texture$Site<-factor(texture$Site, levels=c("H", "Y", "P", "M", "T", "E", "D"))
allsoil$Site<-factor(allsoil$Site, levels=c("H", "Y", "P", "M", "T", "E", "D"))

anova(lm(Silt~PM*Aspect, texture))
ggplot(texture, aes(Sand, x=Site))+geom_boxplot(data=texture, aes(fill=PM), position = position_dodge2(preserve = "single")) +
  geom_point(data=texture, aes(fill=PM), shape=21, position = position_jitterdodge())+
  ylab("Sand %")+facet_grid(Zone~Aspect)+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling1"))

ggplot(allsoil, aes(N.pct, x=Site))+geom_boxplot(data=allsoil, aes(fill=PM), position = position_dodge2(preserve = "single")) +
  geom_point(data=allsoil, aes(fill=PM), shape=21, position = position_jitterdodge())+
  ylab("C %")+facet_grid(Zone~Aspect)+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling1"))

str(allsoil)
soillit<-merge(allsoil, CNlit)
sjt.lmer(lmer(C.N.L~PM*Aspect+(1|Site), soillit), p.kr=F)
str(allsoil)
allsoil<-(merge(allsoil, clim.norms))
summary(lm(Clay~vpdmax, allsoil))
clay<-lmer(Clay~PM+Zone+mat+(1|Site), allsoil)
silt<-lmer(Silt~Zone+Aspect+mat+(1|Site), allsoil)
sand<-lmer(Sand~Zone+Aspect+mat+(1|Site), allsoil)
sjt.lmer(clay, silt, sand, p.kr=F)
