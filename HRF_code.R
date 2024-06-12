library(gridExtra)
library(ggplot2)
library(reshape2)
library(MCMCglmm)
library(devtools)
# devtools::install_github("ctkremer/priceTools")
library(priceTools)
library(readr)
library(grid)
library(ggpubr)
library(multcomp)
library(agricolae)
library(lmerTest)
library(MuMIn)
library(HH)

#### -------- 1. Price equation analysis -------- ####
# Load data (2017 as example):
data<-read.csv("N2017.csv",stringsAsFactors = F)

# Grouping data by key treatment/experimental columns:
grouped.data<-data %>% group_by(trt,plot)

# Run pairwise comparisons generating Price components:
res1<- pairwise.price(grouped.data,species="spp",func="anpp")

# Organize output:
group.vars<-c('plot')
treat.vars<-c('trt')

# Save each year result
# write.csv(res1,file = "2017res.csv")

## Vector plot of each N level 
pp1<-group.columns(res1,gps=c(group.vars,treat.vars),drop=F)
datN0<-pp1[pp1$trt %in% c('N0 N2'),]
datN2<-pp1[pp1$trt %in% c('N2 N4'),]
datN4<-pp1[pp1$trt %in% c('N4 N8'),]
datN8<-pp1[pp1$trt %in% c('N8 N16'),]
datN16<-pp1[pp1$trt %in% c('N16 N32'),]

s1<-leap.zig(datN0,type='price',standardize=FALSE,
             xlim=c(15,30),ylim=c(100,550),error.bars=F,
             vectors=T,raw.points = F,legend=TRUE)+
  annotate("text", x = mean(datN0$x.rich)+1, y = mean(datN0$x.func), 
           label = "N0",size=3)+
  annotate("text", x = mean(datN0$x.rich), y = mean(datN0$x.func), 
           label = "*",size=5)

s2<-leap.zig(datN2,type='price',standardize=FALSE,
             xlim=c(15,30),ylim=c(100,550),error.bars=F,
             vectors=T,raw.points = F,legend=TRUE,add=T,old.plot=s1)+
  annotate("text", x = mean(datN2$x.rich)+1, y = mean(datN2$x.func), 
           label = "N2",size=3)+
  annotate("text", x = mean(datN2$x.rich), y = mean(datN2$x.func), 
           label = "*",size=5)

s3<-leap.zig(datN4,type='price',standardize=FALSE,
             xlim=c(15,30),ylim=c(100,550),error.bars=F,
             vectors=T,raw.points = F,legend=TRUE,add=T,old.plot=s2)+
  annotate("text", x = mean(datN4$x.rich)+1, y = mean(datN4$x.func), 
           label = "N4",size=3)+
  annotate("text", x = mean(datN4$x.rich), y = mean(datN4$x.func), 
           label = "*",size=5)

s4<-leap.zig(datN8,type='price',standardize=FALSE,
             xlim=c(15,30),ylim=c(100,550),error.bars=F,
             vectors=T,raw.points = F,legend=TRUE,add=T,old.plot=s3)+
  annotate("text", x = mean(datN8$x.rich)+1, y = mean(datN8$x.func), 
           label = "N8",size=3)+
  annotate("text", x = mean(datN8$x.rich), y = mean(datN8$x.func), 
           label = "*",size=5)

s5<-leap.zig(datN16,type='price',standardize=FALSE,
             xlim=c(5,25),ylim=c(100,500),error.bars=F,
             vectors=T,raw.points = F,legend=T,add=T,old.plot=s4)+
  annotate("text", x = mean(datN16$x.rich)+1, y = mean(datN16$x.func), 
           label = "N16",size=3)+
  annotate("text", x = mean(datN16$x.rich), y = mean(datN16$x.func), 
           label = "*",size=5)+
  annotate("segment", x = mean(datN16$y.rich)-1, xend = mean(datN16$y.rich)+1, 
           y = mean(datN16$y.func), yend = mean(datN16$y.func),colour = "black")+
  annotate("text", x = mean(datN16$y.rich)+1, y = mean(datN16$y.func), 
           label = "N32",size=3)+
  theme_bw()+
  theme(text= element_text(family = "TNM", face = "bold", size = 16))+
  ggtitle("2017")
#s5

# save figures
# png(filename = "2017-each Nlevel.png",width=2400,height=1800,res = 400) 
# grid.arrange(s5 + theme(aspect.ratio=1),nrow=1)
# dev.off()

## Vector plot of overall N effect (N0 → N+,N+ include N2,N4,N8,N16,N32)
data<-read.csv("N2017.csv",stringsAsFactors = F)

# Grouping data by key treatment/experimental columns:
grouped.data<-data %>% group_by(trt2,plot)

# Run pairwise comparisons generating Price components:
res1<- pairwise.price(grouped.data,species="spp",func="anpp")

# Organize output:
group.vars<-c('plot')
treat.vars<-c('trt2')

pp1<-group.columns(res1,gps=c(group.vars,treat.vars),drop=F)
pp1<-pp1[pp1$trt2 %in% c('N0 N0','N0 N+'),]

# Separate out the comparisions isolating treatment effects
dat1<-pp1[pp1$trt2 %in% c('N0 N+'),]

s2017<-leap.zig(dat1,type='price',standardize=F,
                xlim=c(10,25),ylim=c(100,500),error.bars=F,
                vectors=T,raw.points = F,legend=F)+
  scale_y_continuous(name = expression('ANPP  g*m'^'-2'))+
  annotate("text", x = mean(dat1$x.rich), y = mean(dat1$x.func), 
           label = "*",size=8)+
  annotate("segment", x = mean(dat1$y.rich)-1, xend = mean(dat1$y.rich)+1, 
           y = mean(dat1$y.func), yend = mean(dat1$y.func),colour = "black")+
  theme_bw()+ 
  theme(text= element_text(family = "A"))+
  ggtitle("2017")
s2017

# save figure
# png(filename='2017-Total Neffect.png',width=3000,height=2000,res = 600) # save figures
# grid.arrange(s2017+ theme(aspect.ratio=1))
# dev.off()

#### -------- 2. ANOVA analysis -------- ####
data2 <- read.csv("HRF_data2014-2017.csv")
#ANPP
fit_ANPP <- aov(ANPP~N*Year+(1|Block),data=data2)
summary(fit_ANPP)
interaction2wt(ANPP~N*Year,data = data2)
dcbj_ANPP<- LSD.test(fit_ANPP,'N',p.adj='none')
plot(dcbj_ANPP)
dcbj_ANPP$groups
lm_ANPP <- lmer(ANPP ~ N*Year + (1|Block),data = data2)
summary(lm_ANPP)
anova(lm_ANPP)

#height
fit_Hsp <- aov(Height~N*Year+(1|Block),data=data2)
summary(fit_Hsp)
interaction2wt(Height~N*Year,data = data2)
dcbj_Hsp<- LSD.test(fit_Hsp,'N',p.adj='none')
plot(dcbj_Hsp)
dcbj_Hsp$groups
lm_Hsp <- lmer(Height ~ N*Year + (1|Block),data = data2)
summary(lm_Hsp)
anova(lm_Hsp)

#coverage
fit_Cover <- aov(Cover~N*Year+(1|Block),data=data2)
summary(fit_Cover)
interaction2wt(Cover~N*Year,data = data2)
dcbj_Cover<- LSD.test(fit_Cover,'N',p.adj='none')
plot(dcbj_Cover)
dcbj_Cover$groups
lm_Cover <- lmer(Cover ~ N*Year + (1|Block),data = data2)
summary(lm_Cover)
anova(lm_Cover)

#C:N
fit_CNratio <- aov(CNratio~N*Year+(1|Block),data=data2) 
summary(fit_CNratio)
interaction2wt(CNratio~N*Year,data = data2)
dcbj_CNratio <- LSD.test(fit_CNratio,'N',p.adj='none')
plot(dcbj_CNratio)
dcbj_CNratio$groups
lm_CNratio <- lmer(CNratio ~ N*Year + (1|Block),data = data2)
summary(lm_CNratio)
anova(lm_CNratio)

#### -------- 3. SEM analysis -------- ####
# fit model use psem
library(piecewiseSEM) 
SEM <- data2[,3:11]
SEM <- scale(SEM,center=T,scale=T)
SEM <- data.frame(SEM)
SEM <- na.omit(SEM)
# initial sem model
model2 <- psem(CNratio %~~% Cover,Height %~~% Cover,Height %~~% CNratio,
  lm(Height ~ N,data = SEM),
  lm(Cover ~  N,data = SEM),
  lm(CNratio ~  N,data = SEM),
  lm(PCA ~   Height + Cover + CNratio + N, data = SEM),
  lm(Rorder ~   Height + Cover + CNratio + N + PCA, data = SEM),
  lm(ST ~   Height + Cover + CNratio + N + Rorder + PCA, data = SEM),
  lm(ANPP ~  PCA + Rorder + ST + N + Height + Cover + CNratio,data = SEM),
  data = SEM)
summary(model2)

# final optimal model
model1 <- psem(CNratio %~~% Cover,Height %~~% Cover,Height %~~% CNratio,
  lm(Height ~ N,data = SEM),
  lm(Cover ~  N,data = SEM),
  lm(CNratio ~  N,data = SEM),
  lm(PCA ~   Height + Cover + CNratio + N, data = SEM),
  lm(Rorder ~   Height + Cover + CNratio + N, data = SEM),
  lm(ST ~   Height + Cover + CNratio + N , data = SEM),
  lm(ANPP ~  PCA + Rorder + ST + N,data = SEM),
  data = SEM)

dSep(model1,conditioning = T)
fisherC(model1)
summary(model1)
