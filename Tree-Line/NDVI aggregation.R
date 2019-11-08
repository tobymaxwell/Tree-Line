setwd("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/HOPS")
sisk<-read.csv("HOPS South NDVI.csv")
library(plyr)


sisk <- sisk[order(sisk$NDVI), ]
top3<-by(sisk, sisk["Year"], tail, n=3)
top3

years<-as.factor(1984:2018)
meanday<-NULL
days<-NULL
for(i in years){
  meanday<-data.frame(top3[i])[1,2]
  days<-rbind(days, meanday)
}

NDVImax<-data.frame(apply(days, 1, mean))
NDVImax<-data.frame(days)
NDVImax$Year<-1984:2018
colnames(NDVImax)<-c("DOY.avg", "Year")
NDVImax$DOY.avg<-round(NDVImax$DOY.avg, 0)
plot(DOY.avg~Year, NDVImax[NDVImax$DOY.avg<182,])


NDVImax$Date<-substr(as.Date(NDVImax$DOY.avg, format = "%d", origin = "1.1.2019"), 6, 10)

NDVImax
mean(NDVImax[NDVImax$DOY.avg<182,]$DOY.avg)

####NDVI avg max3
years<-as.factor(1984:2018)
meanNDVI<-NULL
days<-NULL
for(i in years){
  meanNDVI<-data.frame(top3[i])[,4]
  days<-rbind(days, meanNDVI)
}

NDVImax<-data.frame(apply(days, 1, mean))

