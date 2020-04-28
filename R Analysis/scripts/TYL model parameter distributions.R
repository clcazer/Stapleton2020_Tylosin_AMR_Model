#title: "Tylosin Model Parameter Distributions"
#author: "Casey Cazer"
#Last update: "April 28, 2020"


#R version: 3.6.1

#distribution fitting for degradation rate
#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf

#required packages
install.packages("checkpoint")
library(checkpoint)
checkpoint("2020-03-01")

library(xlsx)
library(stringr)
library(fitdistrplus)
library(dplyr)
library(tidyr)
library(ggplot2)

#set up sink to record output
sink(file="results/parameter distributions.txt")
"These are the parameter distributions for the TYL model"
sink()

########################################################
#Degradation rate

#load data
deg <- read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Degradation Rate", colIndex=1:5)

#drop un-necessary columns
deg <- select(deg, -Half.Life..days., -Units)

#average estimates thate are not independent, bind with those that are independent, make one vector of parameter estimates
deg<- unlist(rbind(filter(deg, Independent.Experiment=="No") %>% group_by(Study) %>% summarise(Parameter.Value=mean(Parameter.Value)), filter(deg, Independent.Experiment=="Yes") %>% select(-Independent.Experiment)) %>% select(Parameter.Value))

#visualize and describe distribution
plotdist(deg)
descdist(deg) #falls near gamma

#fit distributions
deg.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('deg.',distributions[i], sep="")
  deg.list[[i]]<-fitdist(deg, distributions[i])
  names(deg.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(deg.list, plot) #all except normal ok
sapply(deg.list, gofstat) #AIC/BIC/Chi-sq/AD prefers lnorm. CVM prefers gamma. KS prefers beta.
#smaller KS, CVM, AD statistics indicate a better fit

#use log normal dist.

#if deg rate is log-normal then ln(deg rate) is normally distributed
#check
logdeg=log(deg)
descdist(logdeg, boot=100)
plotdist(logdeg) #looks good

logdeg.norm=fitdist(logdeg, "norm")
summary(logdeg.norm)
plot(logdeg.norm)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"deg follows lognorm distribution with parameters:"
deg.meanlog <- round(deg.list[["deg.lnorm"]]$estimate['meanlog'],1)
deg.sdlog <- round(deg.list[["deg.lnorm"]]$estimate['sdlog'],1)
deg.meanlog
deg.sdlog
sink()
#lnorm cannot be <0, so appropriate

####################################################
#dist fitting for sorption rates

#see "All Parameters + Citations_cc_2-13-19.xlsx" "All Parameters" sheet for these values
sorp=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Sorption", colIndex=1:5)

#all independent estimates
sorp=sorp$Measurement
#convert to decimal representation of percentage to allow beta dist
sorp=sorp/100

plotdist(sorp)
descdist(sorp) 

#fit distributions
sorp.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('sorp.',distributions[i], sep="")
  sorp.list[[i]]<-fitdist(sorp, distributions[i])
  names(sorp.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(sorp.list, plot) #none stand out
sapply(sorp.list, gofstat) #nothing clearly fits best.

sink(file="results/parameter distributions.txt", append=TRUE)
""
"preference for normal if nothing fits well. Sorption normal distribution:"
sorp.mean <- round(sorp.list[["sorp.norm"]]$estimate['mean'],1)
sorp.sd <- round(sorp.list[["sorp.norm"]]$estimate['sd'],1)
sorp.mean
sorp.sd
sink()

##################################################################

#environment intake
env_intake=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Environment Intake", colIndex=1:5)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Environment/Pen Intake, only two estimates, uniform"
env_in.min<-round(min(env_intake$Measurement),1)
env_in.max<-round(max(env_intake$Measurement),1) #note R rounds 5 to even number
env_in.min
env_in.max
sink()

#####################################################################

#feed intake
feed_intake=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Feed Intake", colIndex=1:5)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Feed intake: only two estimates, based on personal communications, uniform dist"
feed_in.min<-round(min(feed_intake$Measurement),0)
feed_in.max<-round(max(feed_intake$Measurement),0)
feed_in.min
feed_in.max
sink()

###################################################################

#water consumption

water_intake=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Water Consumption", colIndex=1:7)

#use mL/hr estimates
water_intake<-water_intake$Measurement2

plotdist(water_intake)
descdist(water_intake) 

#fit distributions (exclude beta because values not 0-1)
water_intake.list <- list()
distributions <- c("norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('water_intake.',distributions[i], sep="")
  water_intake.list[[i]]<-fitdist(water_intake, distributions[i])
  names(water_intake.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(water_intake.list, plot) #none stand out
sapply(water_intake.list, gofstat) #chi-sq/KS favors weibull. CVM/AD/AIC/BIC favors lnorm. 

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Water intake, LogNorm"
water_intake.meanlog <- round(water_intake.list[["water_intake.lnorm"]]$estimate['meanlog'],1)
water_intake.sdlog <- round(water_intake.list[["water_intake.lnorm"]]$estimate['sdlog'],1)
water_intake.meanlog
water_intake.sdlog
sink()

#################################################################

#fitness cost

fit=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Fitness Cost", colIndex=1:7)

#all independent experiments
#turn measurment percent into decimal
fit<-fit$Measurement/100

plotdist(fit)
descdist(fit) 

#fit distributions
fit.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('fit.',distributions[i], sep="")
  fit.list[[i]]<-fitdist(fit, distributions[i])
  names(fit.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(fit.list, plot) #uniform bad fit
sapply(fit.list,gofstat) #chisq/CVM/AD favors weibull. KS favors beta. AIC/BIC favors lnorm.

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Fitness cost, Weibull. convert back to percentage if needed"
fit.weibull.shape<-round(fit.list[["fit.weibull"]]$estimate['shape'],1)
fit.weibull.scale<-round(fit.list[["fit.weibull"]]$estimate['scale'],2)
fit.weibull.scale
fit.weibull.shape
sink()

#####################################################################

#anaerobic mic penalty
anaerobe_mic=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Aerobicity", colIndex=1:6)

#need to put data in long form (use Number of Independent Strains)
anaerobe<-anaerobe_mic[rep(row.names(anaerobe_mic), anaerobe_mic$Number.of.Independent.Strains),2]

#all independent estimates
plotdist(anaerobe)
descdist(anaerobe)

#since there are negative values and 0, can't fit weibull, gamma, or lnorm
#histogram supports normal over uniform
mean(anaerobe)
sd(anaerobe) #sd is large relative to mean--will lead to substantial amount of negative values (fitness gains)

#however we are considering transmittable genes--the one negative value arose from ribosome mutation induced in lab
#use only positive values and 0
anaerobe<-anaerobe[anaerobe>=0]

#most of these dist won't fit a 0 value (weibull, gamma, lnorm)
#folded normal would make logical sense based on histogram
anaerobe_unfold <- c(anaerobe, -anaerobe)

anaerobe.norm<-fitdist(anaerobe_unfold, 'norm')
plot(anaerobe.norm)
gofstat(anaerobe.norm) #KS rejects. but this distribution still makes logical sense if we treat MIC as continuous rather than discrete

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Anaerobic MIC penalty, folded normal"
anaerobe.foldnorm.mean <- round(anaerobe.norm$estimate['mean'],1)
anaerobe.foldnorm.sd <- round(anaerobe.norm$estimate['sd'],1)
anaerobe.foldnorm.mean
anaerobe.foldnorm.sd
sink()

########################################################################

#Enterococcus growth

Egrowth=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Enterococcus Growth", colIndex=1:9)

Egrowth<- unlist(rbind(filter(Egrowth, Experiments.Independent.=="No") %>% group_by(Citation) %>% summarise(Measurement=mean(Measurement)), filter(Egrowth, Experiments.Independent.=="Yes") %>% select(Citation, Measurement)) %>% select(Measurement))

min(Egrowth) #0
max(Egrowth) # <1

#drop the two 0 observations--doesn't make sense
Egrowth<-Egrowth[Egrowth>0]

Egrowth.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('Egrowth.',distributions[i], sep="")
  Egrowth.list[[i]]<-fitdist(Egrowth, distributions[i])
  names(Egrowth.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(Egrowth.list, plot) #norm, beta, unif poor fits
sapply(Egrowth.list, gofstat) #chi-sq/CVM/AD/KS/AIC/BIC all favor lnorm. 

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Enterococcus Growth, LogNormal; truncate 0,1"
Egrowth.meanlog <- round(Egrowth.list[["Egrowth.lnorm"]]$estimate['meanlog'],1)
Egrowth.sdlog <- round(Egrowth.list[["Egrowth.lnorm"]]$estimate['sdlog'],1)
Egrowth.meanlog
Egrowth.sdlog
sink()

#############################################################

#death rate in pen
Edeathpen=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Enterococcus Death Pen", colIndex=1:7)

#all independent. make positive values to fit distribution
Edeathpen <- -Edeathpen$Measurement


Edeathpen.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('Edeathpen.',distributions[i], sep="")
  Edeathpen.list[[i]]<-fitdist(Edeathpen, distributions[i])
  names(Edeathpen.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(Edeathpen.list, plot) #norm, unif poor fits
sapply(Edeathpen.list, gofstat) #chi-sq/CVM/AD/KS/AIC/BIC all favor lnorm. 

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Enterococcus death rate in pen, Log Normal; negate the value after random draw"
Edeathpen.meanlog <-round(Edeathpen.list[["Edeathpen.lnorm"]]$estimate['meanlog'],1)
Edeathpen.sdlog <- round(Edeathpen.list[["Edeathpen.lnorm"]]$estimate['sdlog'],1)
Edeathpen.meanlog
Edeathpen.sdlog
sink()

#enterococcus death in water and feed have just one point each
#death rate in water, one estimate
Edeathwater=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Enterococcus Death Water", colIndex=1:7)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Enterococcus death rate in water, Uniform"
round(Edeathwater[1,]$Measurement*1.25,2)
round(Edeathwater[1,]$Measurement*0.75,2)
sink()

#death rate in feed, one estimate
Edeathfeed=read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Enterococcus Death Feed", colIndex=1:7)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Enterococcus death rate in feed, Uniform"
round(Edeathfeed[1,]$Measurement*1.25,2)
round(Edeathfeed[1,]$Measurement*0.75,2)
sink()

#####################################################################

#starting fraction resistant
frac.res <- read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Starting Fraction Resistant", colIndex=1:7)

#none independent. average within study
frac.res <- unlist(filter(frac.res, Independent.Experiments.=="No") %>% group_by(Citation) %>% summarise(Measurement=mean(Measurement)) %>% select(Measurement))

#express percentage as decimal
frac.res<- frac.res/100

plotdist(frac.res)
descdist(frac.res)

frac.res.list <- list()
distributions <- c("beta", "norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('frac.res.',distributions[i], sep="")
  frac.res.list[[i]]<-fitdist(frac.res, distributions[i])
  names(frac.res.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(frac.res.list, plot) #norm, unif poor fits
sapply(frac.res.list, gofstat) #chi-sq favors norm. CVM/AD favors gamma/beta. KS favors gamma. AIC/BIC favors beta.

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Starting fraction resistant, beta is logical because bounded on [0,1]"
frac.res.beta.shape1 <- round(frac.res.list[["frac.res.beta"]]$estimate['shape1'],1)
frac.res.beta.shape2 <- round(frac.res.list[["frac.res.beta"]]$estimate['shape2'],1)
frac.res.beta.shape1
frac.res.beta.shape2
sink()

###################################################################

#plasmid transfer freq
plasmid.data <- read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Plasmid Transfer Frequency", colIndex=1:10)

plasmid <- plasmid.data
#some independent experiments, some not. For not-independent: average within study-donor-recipient pair
plasmid$pair <- paste(plasmid$Donor, plasmid$Recipient, sep=":")

plasmid.not.indep <- filter(plasmid, Independent.Experiments.=="No") %>% group_by(Citation, pair) %>% summarise(Rate=mean(Calculate.Rate..T.D...1.h.)) %>% select(Citation, Rate, pair)

plasmid.indep <- filter(plasmid, Independent.Experiments.=="Yes") %>% select(Citation, Rate=Calculate.Rate..T.D...1.h., pair)

plasmid <- c(plasmid.not.indep$Rate, plasmid.indep$Rate)

rm(plasmid.not.indep, plasmid.indep)

#use log10 to reduce skew, negate to fit dist
plasmidlog10 <- -log10(plasmid)

#fit dist
plotdist(plasmidlog10)
descdist(plasmidlog10)

plasmidlog10.list <- list()
distributions <- c("norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('plasmidlog10.',distributions[i], sep="")
  plasmidlog10.list[[i]]<-fitdist(plasmidlog10, distributions[i])
  names(plasmidlog10.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(plasmidlog10.list, plot) #lnorm, unif poor fits
sapply(plasmidlog10.list, gofstat) #chi-sq p-value rejects all distributions. CVM/AD rejects gamma/weibull. KS rejects gamma/lNorm. 


#try non-transmformed data
plotdist(plasmid)
descdist(plasmid)

plasmid.list <- list()
distributions <- c("beta", "gamma", "lnorm", "unif") #norm and weibull throws errors. data not at all normal
for (i in 1:length(distributions)){
  name <- paste('plasmid.',distributions[i], sep="")
  plasmid.list[[i]]<-fitdist(plasmid, distributions[i])
  names(plasmid.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(plasmid.list, plot) #beta/gamma fit best and look equivalent. lnorm looks poor. unif wost.
sapply(plasmid.list, gofstat) #all chi-sq p-values reject. KS rejects beta/gamma

#still no good option for this distribution
#eliminate outliers (>3x mean)
plasmid.no.outlier <- plasmid[plasmid<mean(plasmid)*3]

plasmid.no.outlier.list <- list()
distributions <- c("beta", "gamma", "lnorm", "unif") #norm and weibull throws errors. data not at all normal
for (i in 1:length(distributions)){
  name <- paste('plasmid.no.outlier.',distributions[i], sep="")
  plasmid.no.outlier.list[[i]]<-fitdist(plasmid.no.outlier, distributions[i])
  names(plasmid.no.outlier.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(plasmid.no.outlier.list, plot) #beta/gamma fit best and look equivalent. lnorm looks poor. unif wost.
sapply(plasmid.no.outlier.list, gofstat) #chisq rejects all , KS rejects gamma/beta


#average all estimates within each study to attempt to find better fit. studies with large number of experiments
ggplot(plasmid.data, aes(x = Calculate.Rate..T.D...1.h.)) +
  geom_histogram(aes(color = Citation, fill = Citation), 
                 position = "identity", bins = 30, alpha = 0.4) #data clusters by citation--experiments in the same study tend to be similar to one another

#average all estimates within each citation
plasmid.avg <- unlist(group_by(plasmid.data, Citation) %>% summarise(Rate=mean(Calculate.Rate..T.D...1.h.)) %>% select(Rate))

plasmid.avg.list <- list()
distributions <- c("beta", "gamma",  "lnorm", "unif") #norm and weibull throws errors. data not at all normal
for (i in 1:length(distributions)){
  name <- paste('plasmid.avg.',distributions[i], sep="")
  plasmid.avg.list[[i]]<-fitdist(plasmid.avg, distributions[i])
  names(plasmid.avg.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(plasmid.avg.list, plot) #beta/gamma fit best and look equivalent. lnorm looks poor. unif wost.
sapply(plasmid.avg.list, gofstat) #chi-sq/CVM/AD/KS/AIC/BIC select lNorm

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Plasmid transfer rate, LogNormal"
plasmid.meanlog <- round(plasmid.avg.list[["plasmid.avg.lnorm"]]$estimate['meanlog'],1)
plasmid.sdlog <- round(plasmid.avg.list[["plasmid.avg.lnorm"]]$estimate['sdlog'],1)
plasmid.meanlog
plasmid.sdlog
sink()

##############################################################

#cattle carrying capacity
cattle.cc <- read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Cattle Carrying Capacity", colIndex=1:7)

#some independent experiments, some not. For not-independent: average within study across time and groups
cattle.cc.not.indep <- filter(cattle.cc, Independent.Experiments.=="No") %>% group_by(Citation) %>% summarise(Measurement=mean(Measurement)) %>% select(Citation, Measurement)

cattle.cc.indep <- filter(cattle.cc, Independent.Experiments.=="Yes") %>% select(Citation, Measurement)

cattle.cc <- c(cattle.cc.not.indep$Measurement, cattle.cc.indep$Measurement)
rm(cattle.cc.not.indep, cattle.cc.indep)

#fit dist
plotdist(cattle.cc)
descdist(cattle.cc)

cattle.cc.list <- list()
distributions <- c("norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('cattle.cc.',distributions[i], sep="")
  cattle.cc.list[[i]]<-fitdist(cattle.cc, distributions[i])
  names(cattle.cc.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(cattle.cc.list, plot) 
sapply(cattle.cc.list, gofstat) #chi-sq favors norm. CVM/AD/KS/AIC/BIC favors weibull


sink(file="results/parameter distributions.txt", append=TRUE)
""
"Cattle carrying capacity, Weibull"
cattle.cc.weibull.shape <-round(cattle.cc.list[["cattle.cc.weibull"]]$estimate['shape'],1)
cattle.cc.weibull.scale <- round(cattle.cc.list[["cattle.cc.weibull"]]$estimate['scale'],1)
cattle.cc.weibull.shape
cattle.cc.weibull.scale
sink()

##################################################################

#pen carrying capacity
pen.cc <- read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Pen Carrying Capacity", colIndex=1:7)

#some independent experiments, some not. For not-independent: average within study across time and groups
pen.cc.not.indep <- filter(pen.cc, Independent.Experiments=="No") %>% group_by(Citation) %>% summarise(Measurement=mean(Measurement)) %>% select(Citation, Measurement)

pen.cc.indep <- filter(pen.cc, Independent.Experiments=="Yes") %>% select(Citation, Measurement)

pen.cc <- c(pen.cc.not.indep$Measurement, pen.cc.indep$Measurement)
rm(pen.cc.not.indep, pen.cc.indep)

#fit dist
plotdist(pen.cc)
descdist(pen.cc)

pen.cc.list <- list()
distributions <- c("norm", "gamma", "weibull", "lnorm", "unif")
for (i in 1:length(distributions)){
  name <- paste('pen.cc.',distributions[i], sep="")
  pen.cc.list[[i]]<-fitdist(pen.cc, distributions[i])
  names(pen.cc.list)[[i]] <- name
  rm(name)
}

#compare distributions
sapply(pen.cc.list, plot) 
sapply(pen.cc.list, gofstat) #chi-sq/CVM/AD/KS/AIC/BIC all favors weibull. 

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Pen carrying capacity, Weibull"
pen.cc.weibull.shape <-round(pen.cc.list[["pen.cc.weibull"]]$estimate['shape'],1)
pen.cc.weibull.scale <- round(pen.cc.list[["pen.cc.weibull"]]$estimate['scale'],1)
pen.cc.weibull.shape
pen.cc.weibull.scale
sink()

############################################################################

#water carrying capacity

water_cc<-read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Water Carrying Capacity", colIndex=1:6)
#one study, all not independent. average (geometric mean of log10 values), then uniform +/- .25
water_cc<-mean(water_cc$Measurement)

sink(file="results/parameter distributions.txt", append=TRUE)
""
"Water carrying capacity, Uniform"
water_cc.min<-round(water_cc*.75,1)
water_cc.max<-round(water_cc*1.25,1)
water_cc.max
water_cc.min
sink()

###############################################################################

#feed carrying capacity

feed.cc<-read.xlsx("data/All Parameters + Citations_cc_2-13-20.xlsx", sheetName="Feed Carrying Capacity", colIndex=1:7)
#one study has non-independent values
feed.cc.not.indep <- filter(feed.cc, Independent.Experiment=="No") %>% group_by(Citation) %>% summarise(Measurement=mean(Measurement)) %>% select(Citation, Measurement)

feed.cc.indep <- filter(feed.cc, Independent.Experiment=="Yes") %>% select(Citation, Measurement)

feed.cc <- c(feed.cc.not.indep$Measurement, feed.cc.indep$Measurement)
rm(feed.cc.not.indep, feed.cc.indep)

sink(file="results/parameter distributions.txt", append=TRUE)
""
#Feed carrying capacity, 3 estimates: triangle distribution with min, max, mean (geometric)
feed.cc.min<-round(min(feed.cc),1)
feed.cc.mean<-round(mean(feed.cc),1)
feed.cc.max<-round(max(feed.cc),1)
feed.cc.min
feed.cc.mean
feed.cc.max
sink()

#######################################################################

#MIC distributions from NARMS beef cecal

#NARMS data from cecal samples
mic <- read.csv("data/NARMS-Data-Cecal-2013-2017.csv", header=TRUE)

#beef and Enterococcus only
mic <- filter(mic, SOURCE=="Cecal (Beef)", GENUS=="E")

#tylosin
mic <- select(mic, TYL, TYL.Sign)
summary(mic)

#drop NA
mic <- filter(mic, !is.na(TYL))

#pair TYL.sign and TYL
mic$TYL.pair <- as.factor(paste(mic$TYL.Sign, mic$TYL, sep=""))
summary(mic$TYL.pair) #consistent range tested from <=0.25 to >32. Also have =32, so need to distinguish them

#turn >32 into 64
mic$TYL[mic$TYL.pair==">32"] <- 64

#susceptible distribution (using NARMS breakpoints, lumping intermediate and resistant together)
mic_S <- as.numeric(unlist(filter(mic, TYL<16) %>% select(TYL)))

#log2 transform to fit dist
mic_S <- log2(mic_S)

#fit dist. negative values limit to normal and unif
plotdist(mic_S) #appears normal
descdist(mic_S)

mic_S.norm <- fitdist(mic_S, 'norm')

plot(mic_S.norm)
gofstat(mic_S.norm)
gofstat(mic_S.norm)$cvmtest
gofstat(mic_S.norm)$chisqpvalue
gofstat(mic_S.norm)$adtest
gofstat(mic_S.norm)$kstest

#chi-sq p-value is small and KS rejects. however this is discrete data representing a continuous distribution. Biologically it makes sense for this data to be normal and histogram looks ok; and we are favoring normal if nothing fits well.

sink(file="results/parameter distributions.txt", append=TRUE)
""
"MIC Susceptible distribution, Normal, truncate distribution at -3 (half of smallest tested value) and  (MIC = 16; open interval), then reverse log2 transform"
mic_S.mean <- round(mic_S.norm$estimate['mean'],1)
mic_S.sd <- round(mic_S.norm$estimate['sd'],1)
mic_S.mean
mic_S.sd
sink()