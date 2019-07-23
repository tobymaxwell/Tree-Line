Picea
PS
PN
ggplot(PN, aes(y=scale(rwlB), x=Year))+ geom_line(color = 'thistle3')+
  geom_line(data=Picea, aes(y=scale(PNB), x=Year), color="Forest Green")
ggplot(PN, aes(y=scale(rwlF), x=Year))+ geom_line(color="thistle3")+
  geom_line(data=Picea, aes(y=scale(PNF), x=Year), color="Forest Green")
ggplot(PN, aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color="Thistle3")+
  geom_line(data=Picea, aes(y=scale(PN), x=Year), color="Forest Green")+
  geom_smooth(data=PN, aes(y=scale(rwlB/rwlF), x=Year), color = 'thistle3',method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")

summary(lm(rwlB/rwlF~Year, data =merge(Picea,PN)))

ggplot(PS, aes(y=scale(rwlB), x=Year))+ geom_line(color = "thistle3")+
  geom_line(data=Picea, aes(y=scale(PSB), x=Year), color="Black")
ggplot(PS, aes(y=scale(rwlF), x=Year))+ geom_line(color="thistle3")+
  geom_line(data=Picea, aes(y=scale(PSF), x=Year), color="Black")
ggplot(PS, aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color="Thistle3")+
  geom_line(data=Picea, aes(y=scale(PS), x=Year), color="Black")+
  geom_smooth(data=PS, aes(y=scale(rwlB/rwlF), x=Year), color = 'thistle3',method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")+
  geom_smooth(data=Picea, aes(y=scale(PS), x=Year), color = 'Black',method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")

summary(lm(PS~Year, data =merge(Picea,PS)))


####HOOO####

ggplot(HN[HN$Year>1900,], aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color="Forest Green")+
  #geom_line(data=PS, aes(y=scale(rwlF), x=Year))+         
  #geom_line(data=Picea, aes(y=scale(PSF), x=Year), color="red")+
  geom_line(data=HS, aes(y=scale(rwlB/rwlF), x=Year), color="Black")+
  ylim(-2,4)+
  theme_bw()
ggplot(HN, aes(y=scale(rwlF), x=Year))+ geom_line()+
  #geom_line(data=PN, aes(y=scale(rwlF), x=Year))+         
  #geom_line(data=Picea, aes(y=scale(PNF), x=Year), color="red")+
  geom_line(data=Hood, aes(y=scale(HNF), x=Year), color="Forest Green")

ggplot(HN[HN$Year>1983,], aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color='Forest Green')+geom_line(data=Hood, aes(y=scale(HN), x=Year), color="Thistle3")+
#+geom_smooth(data=Hood, aes(y=scale(HN), x=Year), color="Forest Green", method='lm', se=FALSE)+
#geom_smooth(data=HN[HN$Year>1925,], aes(y=scale(rwlB/rwlF), x=Year), color="thistle3", method='lm', se=FALSE)+
 # geom_smooth(data=HN[HN$Year<1925&HN$Year>1875,], aes(y=scale(rwlB/rwlF), x=Year), color="thistle3", method='lm', se=FALSE)+
  theme_bw()+ylab("Standardized Response")
summary(lm(HN[HN$Year>1983,]$rwlB/HN[HN$Year>1983,]$rwlF~Hood[Hood$Year<2018,]$HN))
summary(lm(HS[HS$Year>1983,]$rwlB/HS[HS$Year>1983,]$rwlF~Hood[Hood$Year<2018,]$HS))
plot(HS[HS$Year>1983,]$rwlB/HS[HS$Year>1983,]$rwlF~Hood[Hood$Year<2018,]$HS)

ggplot(HS, aes(y=scale(rwlB/rwlF), x=Year))+ geom_line()+
  geom_line(data=Hood, aes(y=scale(HSB), x=Year), color="Black", lty='dashed')
ggplot(HS, aes(y=scale(rwlF), x=Year))+ geom_line()+
  geom_line(data=Hood, aes(y=scale(HSF), x=Year), color="thistle3", lty='dashed')
ggplot(HS[HS$Year>1983,], aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color="Black")+
  #geom_smooth(data=HS, aes(y=scale(rwlB/rwlF), x=Year), method='lm', color='thistle3', se=F)+
  geom_line(data=Hood, aes(y=scale(HS), x=Year), color="Thistle3")+
  #geom_smooth(data=Hood, aes(y=scale(HS), x=Year), color="Black", method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")

summary(lm(scale(rwlB/rwlF)~Year, merge(HN, Hood)))
summary(lm(scale(rwlB/rwlF)~Year, merge(HS, Hood)))

####Twin Plots####
ggplot(TN, aes(y=scale(rwlB), x=Year))+ geom_line(color='thistle3')+
  geom_line(data=Twin, aes(y=scale(TNB), x=Year), color="Forest Green")+
theme_bw()+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TNF",], aes(y=scale(MAT), x=Year), color="Blue")

ggplot(TN, aes(y=scale(rwlF), x=Year))+ geom_line(color="thistle3")+
  geom_line(data=Twin, aes(y=scale(TNF), x=Year), color="Forest Green")+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TNF",], aes(y=scale(MAP), x=Year), color="Blue")+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TNF",], aes(y=scale(MAT), x=Year), color="Red")

ggplot(TN, aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color='thistle3')+
 # geom_smooth(data=TN, aes(y=scale(rwlB/rwlF), x=Year), color = 'thistle3',method="lm", se=F)+
  geom_line(data=Twin, aes(y=scale(TN), x=Year), color="Forest Green")+
  #geom_smooth(data=Twin, aes(y=scale(TN), x=Year), color="Forest Green", method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")+
  geom_line(data=NDVI.clim[NDVI.clim$Site=="TNF",], aes(y=scale(MAT), x=Year), color="Blue")


ggplot(TS, aes(y=scale(rwlB), x=Year))+ geom_line()+
  geom_line(data=Twin, aes(y=scale(TSB), x=Year), color="red", lty='dashed')
ggplot(TS, aes(y=scale(rwlF), x=Year))+ geom_line()+
  geom_line(data=Twin, aes(y=scale(TSF), x=Year), color="red", lty='dashed')
ggplot(TS, aes(y=scale(rwlB/rwlF), x=Year))+ geom_line(color='thistle3')+
  geom_line(data=Twin, aes(y=scale(TS), x=Year), color="Black")+
  geom_smooth(data=Twin, aes(y=scale(TS), x=Year), color="Black", method="lm", se=F)+
  theme_bw()+ylab("Standardized Response")+
  geom_smooth(data=TS, aes(y=scale(rwlB/rwlF), x=Year), color = 'thistle3',method="lm", se=F)

summary(lm(scale(rwlB/rwlF)~Year, merge(TN, Twin)))
summary(lm(scale(rwlB/rwlF)~Year, merge(TS, Twin)))
