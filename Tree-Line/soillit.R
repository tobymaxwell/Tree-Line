setwd("/Users/Maxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/")
soil<-read.csv("Soil Data.csv")
str(soil)

library(agricolae)
mod<-lm(litdens_sqin~PM, soil[soil$Depth=="tw",])
anova(mod)
tuk<-HSD.test(mod,trt="PM", DFerror=117, MSerror=3.7605,group=T )
tuk

t.test(soilwt~Site, soil[soil$Depth=="tw",])

soil$Type
m1<-lm(rockwt~Zone*PM, soil[soil$Depth=="tw",])
anova(m1)
