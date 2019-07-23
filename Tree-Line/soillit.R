setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/")
soil<-read.csv("Soil Data.csv")
CN<-read.csv("CNsoil.csv")
CNlit<-read.csv("CNlit.csv")

allsoil<-merge(CN,soil, by=c("Site", "Aspect", "Zone", "Number", "Depth"))

str(allsoil)
allsoil$AZ<-paste0(allsoil$Aspect, allsoil$Zone)
allsoil$C_g_kg<-((allsoil$C.pct/100)*1000)
allsoil$N_g_kg<-((allsoil$C.pct/100)*1000)
allsoil$C_bd<-(allsoil$N.pct/100)/allsoil$soilvol
m1<-lm(C_g_kg~Site*Aspect, allsoil)
anova(m1)
#allsoil$Cdens<-C_g_kg
CNlit$AZ<-paste0(CNlit$Aspect, CNlit$Zone)
library(ggplot2)
allsoil.tw<-allsoil[allsoil$Depth=="tw",]
ggplot(allsoil.tw[allsoil.tw$Site=="H"|allsoil.tw$Site=="T"|allsoil.tw$Site=="P",], aes(y=C.N, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("C:N")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  theme(legend.position="none")+
  xlab("")

ggplot(allsoil[allsoil$Site=="H"|allsoil$Site=="T"|allsoil$Site=="P",], aes(y=N.pct, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Nitrogen (g/kg)")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
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

ggplot(CNlit[CNlit$Site=="H"|CNlit$Site=="T"|CNlit$Site=="P",], aes(y=C.N, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Litter C:N")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
  theme(legend.position="none")+
  xlab("")

ggplot(CNlit[CNlit$Site=="H"|CNlit$Site=="T"|CNlit$Site=="P",], aes(y=N.pct, x=Site, fill = AZ))+geom_boxplot()+theme_bw()+ylab("Litter N")+scale_fill_manual(values=c("Forest Green", "Green", "Black", "Grey"))+
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
