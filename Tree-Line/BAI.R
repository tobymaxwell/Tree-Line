install.packages("dplR")
library("dplR")
library(plyr)
# Bring in master file of RWs
Master=read.csv("/Users/tobymaxwell/OneDrive - University Of Oregon/Oregon/Nat Geo/Data/Rings/Hall.csv", header=TRUE)
#Basic Ring Width Index calculations
n<-min(na.omit(Master$Year))
Master$Year<-as.numeric(rep(2017:n, each=2))
library(dplyr)
Master.T<-Master[-1] %>%
  group_by(Year) %>% 
  summarise_all(funs(sum))
Master.T<-data.frame(Master.T)
str(Master.T)
write.csv(Master.T, "/Users/tobymaxwell/Desktop/Master.csv")
H.rwl<-read.csv("/Users/tobymaxwell/Desktop/Master.csv")
#revdf<-function(df)df[seq(dim(df)[1],1),]
#Master.T<- Master.T[seq(dim(Master.T)[1],1),]
RWI<-detrend(H.rwl[-1],make.plot=TRUE,method=("Spline"),nyrs=NULL, verbose=TRUE)
rownames(RWI)<-H.rwl$Year
spag.plot(RWI[15:22], zfac=.3)
rwl.report(RWI[15:22])

##### Calculate Master RWI using the C-Method (recomended method, see Biondi and Qeadan, 2008) 
library(graphics)
library(utils)
data(gp.rwl)
str(gp.rwl)
gp.rwl$`50B`
data(gp.po)
gp.rwi <- cms(rwl = gp.rwl, po = gp.po)
gp.crn <- chron(gp.rwi)
crn.plot(gp.crn, add.spline = TRUE)
str(Master.T)
HSFAbLa<-Master.T[1:6]
HSFAbLa<-revdf(HSFAbLa)
rownames(HSFAbLa)<-HSFAbLa$Year
PO<-NULL
PO$series<-colnames(HSFAbLa[-1])
PO<-data.frame(PO)
cols<-colnames(HSFAbLa[-1])
lengths<-NULL
for(i in cols){
  length<-length(rownames(na.omit(select(HSFAbLa, i)[1])))
lengths<-c(lengths, length)
}
PO$pith.offset<-lengths
colnames(HSFAbLa[-1])
RWI_C<-cms(HSFAbLa[-1], PO, c.hat.t = FALSE, c.hat.i = FALSE)
#write.csv(RWI_C, "C:/Users/Paulo/Documents/R/CH2/RWI_C.csv")

### Calculate and plot a mean RWI
MeanChron<-chron(RWI_C, prefix = "IZT", biweight = TRUE, prewhiten = TRUE)
MeanChron
crn.plot(MeanChron, add.spline=TRUE, nyrs=10, f=0.5, xlab="Years", ylab="Ring Width Index", 
         crn.lwd= 2, spline.lwd = 2, abline.pos = 1)
title("Ring Width Index", line = 3)
title(sub = "RWI - AR Model", line = -11, font.sub = 2)

### Basal Area Increment calculation. DBS's specified by RW length ####
str(Master.T)
Master_BAI<-Master.T
Master_BAI[2:45]<-bai.out(rwl=Master.T[2:45])
Master_BAI[1:6]

write.csv(Master_BAI, "/Users/tobymaxwell/Desktop/Master_BAI.csv")

#BAI graphical analysis ##
install.packages("plotly")
install.packages("plyr")
library(ggplot2)
library(reshape)

### The script below is to plot each time series by "treatment" with standard errors (calculated in Excel) ###

library(ggplot2)
ggplot(BAI, aes(x = Year, y = North_High))+
  geom_line(aes(y=North_High, colour="North_High"),lwd=1) +                        
  geom_ribbon(aes(ymin = North_High - SE_NH, 
                  ymax = North_High + SE_NH),alpha = 0.2, fill="cyan")+
  geom_line(aes(y=North_Low, colour="North_Low"),lwd=1) +
  geom_ribbon(aes(ymin = North_Low - SE_NL, 
                  ymax = North_Low + SE_NL), alpha =0.2, fill="green")+
  geom_line(aes(y=South_High, colour="South_High"),lwd=1) +
  geom_ribbon(aes(ymin = South_High - SE_SH, 
                  ymax = South_High + SE_SH),lwd=.2, alpha = 0.2, fill="red")+
  geom_line(aes(y=South_Low, colour="South_Low"),lwd=1) +
  geom_ribbon(aes(ymin = South_Low - SE_SL, 
                  ymax = South_Low + SE_SL),lwd=.2, alpha = 0.2, fill="yellow")+
  scale_color_manual(values = c("cyan", "green", "red", "orange")) +
  labs(color="Locations") +
  xlab('Year') +
  ylab("Basal Area Increment (mm2)")+
  ggtitle("Pinus hartwegii Growth by Elevation and Aspect")+
  theme(plot.title = element_text(lineheight=.8, size=20, face="bold"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.title=element_text(size=14, face="bold"))

#### OK, the following script gets more into the analysis, though I haven't really run my own models yet
#### Here I created first a data frame (in the code is called "PO") for all my independent variables which relates to the response variable (BAI) through a "Tree ID" column
#### This is needed to you can then "melt" both data frames by "Tree ID" and run the analysis on your BAI data.

## You need Tidyverse 

library(Tidyverse)

##### Tranformations ####


logBAI_North.m <- melt(logBAI_North, id = "Cambial_Age")
logBAI_South.m <- melt(logBAI_South, id = "Cambial_Age")
ggplot(data = logBAI_North.m, aes(x = Cambial_Age, y = value, color="blue"))+
  geom_point(data=logBAI_South.m, aes(x = Cambial_Age, y = value, color = "red"))

### Bringing the two data frames together

allDat <- logBAI %>% #this says, take the first dataframe (tsDat) and the pipe %>% indicates that you want to perform a function to that dataframe
  gather(-NumRings, key = TreeID, value = log_BAI) %>%  #this changes the first dataframe (tsDat) from wide to long format
  full_join(PO, by = c("TreeID" = "ID")) #this joins the long version of tsDat with rDat and joins them by the RingNum column


### BAI Wide format to Long format

M_BAI_Long <- M_BAI_Wide %>% 
  gather(-NumRings, key = TreeID, value = BAI) %>%  
  full_join(PO, by = c("TreeID" = "ID")) 

### RWI Wide format to long format

RWI_Long <- RWI_Wide %>% 
  gather(-RingNum, key = TreeID, value = RWI) %>%  
  full_join(PO, by = c("TreeID" = "ID")) 

#Drop NAs

RWI_Long <- RWI_Long[!is.na(RWI_Long$RWI), ]

# Generate Cambial Age Column

library(dplyr)
CA <- Dat %>% group_by(TreeID) %>% mutate(id = row_number())

#Plot it
# By aspecta and altitude

ggplot(data = BAI_Cambial, aes(x = Cambial_Age, y = log_BAI, color = CA)) +
  geom_point(size = 2, shape = 2) +
  scale_colour_gradient(low = "green", high = "red")+
  facet_grid(.~ Aspect)+
  theme(legend.position = "bottom") + 
  xlab("Cambial Age") +
  ylab("Log BAI")

# Aspect in color

ggplot(data = BAI_Cambial, aes(x = Cambial_Age, y = log_BAI, color = Aspect)) +
  geom_point(size = 2, shape = 2) +
  scale_color_manual(values=c("green", "red"))+
  facet_grid(.~ Alt_Cat)+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  xlab("Cambial Age") +
  ylab("Log BAI")

# This creates the group means

gm <- Dat %>% 
  group_by("Aspect", "Cambial_Age") %>% 
  summarise(log_BAI = mean(log_BAI))

ggplot(data = Dat, aes(x = Cambial_Age, y = log_BAI, group=Alt_Cat, color = Aspect)) +
  geom_point(aes(group = Aspect), size = 2, shape = 1) +
  geom_line(data=gm)+
  theme(legend.position = "bottom") + 
  xlab("Cambial Age") +
  ylab("Log BAI")

# Plotting by fixed diameters #
library(tidyr)
# Move from wide to long format
RW_Long <- gather(RW, TreeID, RW, A19P3:B18P5, factor_key=FALSE)
RW_Long <- RW_Long[!is.na(RW_Long$RW), ]
# Calculate Running Totals (Radius)
CumRW <- RW_Long %>% group_by(TreeID) %>% mutate(cumsum = cumsum(RW))
# Create cm column
CumRW$Radius <- CumRW$cumsum / 10
# Merge CumRW and BAI_Cambial
BAI_Cambial <- merge(BAI_Cambial, CumRW, by = TreeID)
#Add Cambial Age as sequence # After some stupid error when I deleted CA
BAI_Cambial %>% group_by(TreeID) %>% mutate(id = row_number())
BAI_Cambial$X <- NULL
BAI_Cambial$CA <- sequence(rle(as.character(BAI_Cambial$TreeID))$lengths)
# Finally -> Master Long Format Data Set
Master_Long = read.csv("C:/Users/Paulo/Documents/R/CH2/Master_Long.csv")

# Check Diam ~ CA plots

ggplot(data = Master_Long, aes(x = CambialAge, y = Diam, color = Altitude)) +
  geom_point(size = 2, shape = 2)+
  scale_color_gradient(low="red", high="green")+
  facet_grid(.~ Aspect)
theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  xlab("Cambial Age") +
  ylab("Diameter")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))


# Join Fix Diameter Values to Long Formate Master

FixDiam45 = read.csv("C:/Users/Paulo/Documents/R/CH2/FixDiam45cm.csv")
require(plyr)
Master_Long<-join(Master_Long, FixDiam45, by = "TreeID", type="left")

# Plot BAi and RWI using Fixed Diameters
library(ggplot2)
ggplot(data = Master_Long, aes(x = Year_FD10, y = FD10_RWI, color = Aspect)) +
  geom_point(size = 3, shape = 2, stroke = 1.5) +
  geom_smooth(method="lm", show.legend = FALSE)+
  scale_color_manual(breaks= c("N", "S"),
                     values=c("green", "red"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle("RWI when DBH = 10 - 11 cm")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Year") +
  ylab("RWI")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))+
  theme(strip.text.x = element_text(size = 8, face="bold"))

# Plot growth at fixed diameter by cambia age

ggplot(data = Master_Long, aes(x = Year_FD10, y = FD10_BAI, color = Aspect)) +
  geom_point(size = 3, shape = 2, stroke = 1.5) +
  geom_smooth(method="lm", show.legend = FALSE)+
  scale_color_manual(breaks= c("N", "S"),
                     values=c("green", "red"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 17)) + 
  theme(legend.key.size = unit(2,"line"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle("BAI when DBH = 45 - 46 cm")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Year") +
  ylab("BAI (cm2)")+
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))

