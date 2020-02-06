str(climhist.lags)
library(lme4)
library(lavaan)
library(semPlot)
mer1<-lmer(lnbai~(1|cumbai), climhist.lags)


plot(residuals(mer1)~climhist.lags$lnbai)

climhist.lags$resid<-residuals(mer1)

str(bai.clim)

summary(lm(resid~co2, climhist.lags))

str(bai.clim)
bai.clim.sem<-climhist.lags%>%
  group_by(SiteID, Aspect, PM,Year, Species)%>%
  select(lnbai, resid, mat, ppt, tmin, SLA, Clay, curv, lnbai, elevation, latitude, tmin.yr, tmean.yr, ppt.yr, LAIint, C.N, tmean.yr, co2, ppt.lag1, ppt.prev3, Winter_ppt, Spring_ppt, Spring_vpdmax, Summer_vpdmax)
bai.clim.sem<-as.data.frame(bai.clim.sem)
str(bai.clim.sem)

cor(na.omit(bai.clim.sem[-1:-4]))

model= '
lnbai~a9*Clay+b9*SLA+tmin.yr+tmean.yr+ppt.yr+d1*tmin+d2*mat+d3*ppt+d4*elevation+d5*latitude+d6*LAIint+d7*curv
LAIint~l1*mat+l2*ppt+l3*tmin+curv+elevation+latitude+Clay+SLA+l4*C.N
Clay~a1*mat+a2*ppt+a3*tmin+a4*curv+a5*elevation+a6*latitude+a7*LAIint
indirectpptclay:=a2*a9
indirecttminclay:=a3*a9
indirectmatlai:=l1*d6
indirectpptlai:=l2*d6
indirectCNLAI:=l4*d6
'

y = c(-0.5, 0, 0, 1, 1, 1, 1, 1,1,1,1,1,1 )
x = c(0,-0.7, .7, .8, -1, -.8, -.6, -.4, -.2, 0, .2, .4, .6)
ly = matrix(c(x, y), ncol=2)

model.fit=sem(model, data=scale(bai.clim.sem[-1:-3][bai.clim.sem$Year>2003&bai.clim.sem$Aspect=='S',]))
summary(model.fit, standardized=TRUE, rsquare=T, fit.measures=T)
semPaths(model.fit, what ='est', whatLabels='est', layout=ly, intercepts = F)

model.fit=sem(model, data=scale(bai.clim.sem[-1:-3][bai.clim.sem$Year<2004&bai.clim.sem$Year>1952&bai.clim.sem$Aspect=='S',]))
summary(model.fit, standardized=TRUE, rsquare=T, fit.measures=T)
semPaths(model.fit, what ='est', whatLabels='est', layout=ly)

model.fit=sem(model, data=scale(bai.clim.sem[-1:-3][bai.clim.sem$Year<1953&bai.clim.sem$Year>1912&bai.clim.sem$Aspect=='N',]))
summary(model.fit, standardized=TRUE, rsquare=T, fit.measures=T)
semPaths(model.fit, what ='est', whatLabels='est', layout=ly)

results.sem<-read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/semNSindresults.csv")
str(results.sem)
ggplot(results.sem, aes(abs(Value), x=Year, color=ID, linetype=Aspect))+geom_line()+theme_bw()
ggplot(na.omit(results.sem), aes(abs(ValuevsPPT), x=Year, color=ID, linetype=Aspect))+geom_line()+theme_bw()
