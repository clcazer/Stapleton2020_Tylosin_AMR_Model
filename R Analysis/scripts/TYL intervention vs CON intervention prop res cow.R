#title: "TYL Model figures"
#author: "Casey Cazer"
#Last update: "April 28, 2020"
#required packages

install.packages("checkpoint")
library(checkpoint)
checkpoint("2020-03-01")

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(xlsx)
library(gridExtra)

#function for making functional boxplots
source('scripts/TYL func boxplots.R')

#variables used through-out
#days of simulation
days<-read.table("data/NoTYL_NI_Days.txt", sep=",")
days<-as.numeric(days[,1]) - 50 #subtract burn-1n

quants<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) #percentiles

#line colors and styles
ltype=c("dotted", "dotdash", "longdash", "solid", "longdash", "dotdash", "dotted")
lcolor=c("black", "#F8766D", "#7CAE00", "#619CFF", "#7CAE00", "#F8766D", "black")


#TYL concentration in LI
#import data: each timepoint (1 hour) in a row; one column per simulation (1000 simulations per scenario)
tyl.conc.ni <- read.table("data/TYL_NI_TYL_li_conc.txt", sep=",")

#summarize the simulations, get percentiles (quants) at each timepoint
tyl.conc.summary<-data.frame(t(apply(tyl.conc.ni, 1, quantile, probs=quants)))

#functional boxplot of TYL LI concentrations
ggsave('figures/Fig 3 tyl_conc.png',
func.boxplots(quants, tyl.conc.ni, days, ltype, lcolor, "Tylosin Concentration (ug/mL)" )+
  scale_y_continuous(limits=c(0,2), breaks=c(seq(0,2,0.25)), expand=c(0,0)),
device="png",
width=6,
height=5,
units="in",
dpi=320)

sink("results/tyl intestine concentrations.txt")
"descriptive results (median, 5th, 95th percentiles) for peak (max) and end of treatment (index 3433)"
"median, peak and end"
round(max(tyl.conc.summary$X50.),2)
round(tyl.conc.summary$X50.[3433],2)

""
"5th percentile, peak and end"
round(max(tyl.conc.summary$X5.),2)
round(tyl.conc.summary$X5.[3433],2)

""
"95th percentile, peak and end"
round(max(tyl.conc.summary$X95.),2)
round(tyl.conc.summary$X95.[3433],2)

"maximum tyl conc still substantially lowr than intermediate MIC (16 ug/mL)"
max(tyl.conc.ni)
sink()

#############################
#validation plot. 
#load data from 2020 tylosin meta-analysis (Cazer et al)
enc=read.xlsx("data/Extracted Data_long form_meta.xlsx", header=TRUE, sheetName="Enterococcus Phenotypic")

#select only measurements of Perc Isolates Resistant against ERY (not all studies measured TYL) and not-missing day (exclude Molitoris abstract with day=="NR")
enc <- filter(enc, AM.R=="ERY" & Outcome.Type=="PercIsolatesR" & Day!="NR") %>% select("Study", "Day", "AM.R", "Group", "Value")
enc$Day <- as.numeric(as.character(enc$Day))
enc$Value <- enc$Value/100 #Value was in percent format, convert to proportion (decimal) to match schmidt
enc$Study <- as.character(enc$Study)
enc$AM.R <- as.character(enc$AM.R)
enc$Group <- as.character(enc$Group)

#exclude Beukers day 225--occured after end of tylosin administration. No other studies sampled on day 225
enc<-filter(enc, Day!=225)

#add schmidt 2020 data
schmidt<-read.xlsx("data/Schmidt validation data.xlsx", header=TRUE, sheetName="Sheet1")
schmidt <- select(schmidt, "Study", "Day", "AM.R", "Group", "Value")
schmidt$Study <- as.character(schmidt$Study)
schmidt$AM.R <- as.character(schmidt$AM.R)
schmidt$Group <- as.character(schmidt$Group)

#combine
enc <- bind_rows(enc, schmidt)

#split tyl and con groups
tyl=enc[which(enc$Group=="Tylosin"),]
con=enc[which(enc$Group=="Control"),]


#merge and take difference: TYL - CON
enc.diff=merge(tyl, con, by=c("Study", "AM.R", "Day"))
enc.diff$Difference=enc.diff$Value.x-enc.diff$Value.y
enc.diff <- dplyr::rename(enc.diff, days=Day, value=Difference, variable=Study)

#rename for later
names(tyl)<-c("variable", "days", "AM.R", "Group", "value")
tyl$variable<-as.factor(tyl$variable)
names(con)<-c("variable", "days", "AM.R", "Group", "value")
con$variable<-as.factor(con$variable)

#validation model simulation data
valid.notyl.cow <- read.table(paste("data/NoTYL_validation_Prop_cow_res.txt", sep=""), sep=",")
valid.tyl.cow <- read.table(paste("data/TYL_validation_Prop_cow_res.txt", sep=""), sep=",")
valid.diff.cow <- valid.tyl.cow - valid.notyl.cow

#days is longer in validation runs
days.v<-read.table("data/NoTYL_validation_Days.txt", sep=",")
days.v<-as.numeric(days.v[,1]) - 50 #subtract burn-1n

#functional box plots with validation points
A<-func.boxplots(quants, valid.tyl.cow, days.v, ltype, lcolor, "Proportion Resistant")+
  scale_y_continuous(limits=c(-0.2,1), breaks=c(seq(-0.2,1,0.1)), expand=c(0,0))+
  geom_point(data=tyl, aes(x=days, y=value, stroke=2, shape=factor(variable)))+
  scale_shape_manual(values=c(21,22,4,24,23))

B<-func.boxplots(quants, valid.notyl.cow, days.v, ltype, lcolor, "Proportion Resistant")+
  scale_y_continuous(limits=c(-0.2,1), breaks=c(seq(-0.2,1,0.1)), expand=c(0,0))+
  geom_point(data=con, aes(x=days, y=value, stroke=2, shape=factor(variable)))+
  scale_shape_manual(values=c(21,22,4,24,23))

C <- func.boxplots(quants, valid.diff.cow, days.v, ltype, lcolor, "Difference in Proportion Resistant")+
  scale_y_continuous(limits=c(-0.2,1), breaks=c(seq(-0.2,1,0.1)), expand=c(0,0))+
  geom_point(data=enc.diff, aes(x=days, y=value, stroke=2, shape=factor(variable)))+
  scale_shape_manual(values=c(21,22,4,24,23))


ggsave('figures/Fig 2 Prop_cow_res_validation.png',
grid.arrange(plot(A)+ggtitle('A'), plot(B)+ggtitle('B'), plot(C)+ggtitle('C'), ncol=1),
width=6,
height=15,
units="in",
dpi=320)


#number of validation obs that fall within 25-75th percentiles (range1) and 5th-95th percentiles (range2)
#tyl and con values
for (i in 1:nrow(enc)){
  D=enc$Day[i]
  G=enc$Group[i]
  V=enc$Value[i]
  index=match(D,days.v)
  if (G=="Tylosin"){
    p25=quantile(valid.tyl.cow[index,], 0.25)
    p75=quantile(valid.tyl.cow[index,], 0.75)
    p05=quantile(valid.tyl.cow[index,], 0.05)
    p95=quantile(valid.tyl.cow[index,], 0.95)
  }else{
    p25=quantile(valid.notyl.cow[index,], 0.25)
    p75=quantile(valid.notyl.cow[index,], 0.75)
    p05=quantile(valid.notyl.cow[index,], 0.05)
    p95=quantile(valid.notyl.cow[index,], 0.95)
  }
  enc$in.range1[i]<-p25<V & V<p75
  enc$in.range2[i]<-p05<V & V<p95
  rm(D,G,V,index,p25,p75, p05, p95)
}

sink("results/validation.txt")
"Number of validation points for TYL and CON groups that fall in 25-75th percentile range (Range1) and 5-95th percentile range (Range2)"
group_by(enc, Group) %>%  dplyr::summarise(n=n(),n.inrange1=sum(in.range1), inrange2=sum(in.range2), inrange1perc=sum(in.range1)/n())
sink()

#difference (tyl-con) values
for (i in 1:nrow(enc.diff)){
  D=enc.diff$days[i]
  V=enc.diff$value[i]
  index=match(D,days.v)
 
  p25=quantile(valid.diff.cow[index,], 0.25)
  p75=quantile(valid.diff.cow[index,], 0.75)
  p05=quantile(valid.diff.cow[index,], 0.05)
  p95=quantile(valid.diff.cow[index,], 0.95)
  
  enc.diff$in.range1[i]<-p25<V & V<p75
  enc.diff$in.range2[i]<-p05<V & V<p95
  rm(D,V,index,p25,p75, p05, p95)
}

sink("results/validation.txt", append=TRUE)
""
"Number of validation points for TYL - CON (TYL effect) that fall in 25-75th percentile range (Range1) and 5-95th percentile range (Range2)"
dplyr::summarise(enc.diff, n=n(),n.inrange1=sum(in.range1), inrange2=sum(in.range2), inrange2perc=sum(in.range2)/n())
sink()

#compare the one trial not used in model parameterization (Jacob 2008)
#percentile of Jacob results in the simulation data
day.index <- match(enc.diff %>% filter(variable=="Jacob08") %>% select(days), days.v)
sink("results/validation.txt", append=TRUE)
""
"Jacob 2008 (not used in model parameterization), where does it fall in the model results? percentile of the TYL-CON value"
ecdf(valid.diff.cow[day.index,])(enc.diff %>% filter(variable=="Jacob08") %>% select(value)) 
sink()


#plot the difference in proportion resistant between TYL and no TYL control for each scenario (NI, DFM, AFTP, RWT, ALL).
#Take the difference between individual simulations because parameter values are the same for each simulation (due to random number seed)
#check assumption that parameter values are equal for each simulation
NoTYL_NI_param <- read.table("data/NoTYL_NI_MC_parameters.txt", sep=",")
TYL_NI_param <- read.table("data/TYL_NI_MC_parameters.txt", sep=",")
all.equal(NoTYL_NI_param, TYL_NI_param) #correct
rm(NoTYL_NI_param, TYL_NI_param)

#automate loading data and calculating differences (note that No TYL NI scenario == No TYL RWT scenario)
datafiles<-c("NoTYL_NI_", "TYL_NI_", "NoTYL_NI_", "TYL_RWT_", "NoTYL_AFTP_", "TYL_AFTP_", "NoTYL_DFM_", "TYL_DFM_", "NoTYL_ALL_", "TYL_ALL_") #control followed by TYL intervention group. RWT control is the same as NoTYL_NI becuase with no tylosin there can be no withdrawal of tylosin. 

#load in proportion resistant
for (i in seq(1,length(datafiles),2)){
  notyl.cow <- read.table(paste("data/",datafiles[i],"Prop_cow_res.txt", sep=""), sep=",")
  tyl.cow <- read.table(paste("data/",datafiles[i+1],"Prop_cow_res.txt", sep=""), sep=",")
 
  diff.cow <- tyl.cow - notyl.cow
  
  label <- paste(strsplit(datafiles[i+1], "_")[[1]][2], "_Prop_res", sep="")
  
  list <- list(notyl.cow, tyl.cow, diff.cow)
  names(list) <- c("notyl.cow", "tyl.cow", "diff.cow")
  
  assign(label, list)
}
rm(notyl.cow, tyl.cow, diff.cow, list, label)

#descriptive stats: 
tyl.cow.summary<-as.data.frame(t(apply(NI_Prop_res[['tyl.cow']],1, quantile, probs=quants)))

sink("results/TYL effect.txt")
"TYL_NI, median. Min, max over time and end of simulation"
round(min(tyl.cow.summary$`50%`),4)
round(max(tyl.cow.summary$`50%`),4)
round(tyl.cow.summary$`50%`[3433],4)

"TYL_NI, 5th percentile. Min, max over time and end of simulation"
round(min(tyl.cow.summary$`5%`),4)
round(max(tyl.cow.summary$`5%`),4)
round(tyl.cow.summary$`5%`[3433],4)

"TYL_NI, 95th percentile. Min, max over time and end of simulation"
round(min(tyl.cow.summary$`95%`),4)
round(max(tyl.cow.summary$`95%`),4)
round(tyl.cow.summary$`95%`[3433],4)
sink()

#change over time: minimal (<10 pp), moderate (10-50 pp) or substantial (>50pp) changes
#for TYL NI sims
time.PR.change=NI_Prop_res[['tyl.cow']][3433,]-NI_Prop_res[['tyl.cow']][1,]

sink("results/TYL effect.txt", append=TRUE)
""
"change over time in TYL_NI simulations"
"% with decrease"
ecdf(time.PR.change)(0) 
"maximum decrease"
min(time.PR.change)

"% with increase"
ecdf(time.PR.change)(1) - ecdf(time.PR.change)(0)

"% with minimal increase"
ecdf(time.PR.change)(0.1)-ecdf(time.PR.change)(0) #percent of sims with between 10pp increase and no change

"% with moderate increases"
ecdf(time.PR.change)(0.5)-ecdf(time.PR.change)(0.1) #% sims with moderate increase

"% with substantial increase"
ecdf(time.PR.change)(1)-ecdf(time.PR.change)(0.5) #% sims with substantial change

""
"CON NI change over time"
time.PR.change.con=NI_Prop_res[['notyl.cow']][3433,]-NI_Prop_res[['notyl.cow']][1,]
"% with decrease"
ecdf(time.PR.change.con)(0) 
"maximum decrease"
min(time.PR.change.con)

"% with minimal increase"
ecdf(time.PR.change.con)(0.1)-ecdf(time.PR.change.con)(0) #percent of sims with between 10pp increase and no change

"% with moderate increases"
ecdf(time.PR.change.con)(0.5)-ecdf(time.PR.change.con)(0.1) #% sims with moderate increase

"% with substantial increase"
ecdf(time.PR.change.con)(1)-ecdf(time.PR.change.con)(0.5) #% sims with substantial change

""
"difference between TYL and CON at the end of the feeding period"
TYL.effect <- NI_Prop_res[['tyl.cow']][3433,]-NI_Prop_res[['notyl.cow']][3433,]

"% with decrease"
ecdf(TYL.effect)(0)
"largest decrease"
min(TYL.effect) 
"the largest decrease is 0 percentage points. so 5% of simulations had no difference between TYL and CON"

"% with minimal increase"
ecdf(TYL.effect)(0.1) - ecdf(TYL.effect)(0)

"% with minimal increase or decrease"
ecdf(abs(TYL.effect))(0.1)

"% with moderate increase"
ecdf(TYL.effect)(0.5) - ecdf(TYL.effect)(0.1)

"% with substantial increase"
ecdf(TYL.effect)(1) - ecdf(TYL.effect)(0.5)

"median difference"
quantile(TYL.effect, 0.5)
sink()


#plotting
scenarios <- c("NI", "RWT", "DFM", "AFTP", "ALL")
#create figure of Prop res in cow for all scenarios
for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_Prop_res", sep=""))
  
  plot <- func.boxplots(quants, data[['diff.cow']], days, ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))
  
  label <- paste("plot.cow.Pr.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}


#plotting Prop Res in each intervention
ggsave('figures/Fig 4 Prop_cow_res_interventions.png',
grid.arrange(plot.cow.Pr.RWT+
               ylab("Difference in Proportion Resistant")+
               xlab("")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("A"), 
             plot.cow.Pr.DFM+
               ylab("")+
               xlab("")+
               ggtitle("B"), 
             plot.cow.Pr.AFTP+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("C"), 
             plot.cow.Pr.ALL+
               ylab("")+
               xlab("Days")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("D"), ncol=2),
width=6*2,
height=5*2,
units="in",
dpi=320)


#descriptive stats
sink("results/TYL effect.txt", append=TRUE)
""
"difference between TYL and CON at the end of the feeding period for different interventions"
"RWT"
TYL.effect.RWT <- RWT_Prop_res[['tyl.cow']][3433,]-RWT_Prop_res[['notyl.cow']][3433,]

"% with minimal increase or decrease"
ecdf(abs(TYL.effect.RWT))(0.1)

"% with moderate increase"
ecdf(TYL.effect.RWT)(0.5) - ecdf(TYL.effect.RWT)(0.1)

"% with substantial increase"
ecdf(TYL.effect.RWT)(1) - ecdf(TYL.effect.RWT)(0.5)


"DFM"
TYL.effect.DFM <- DFM_Prop_res[['tyl.cow']][3433,]-DFM_Prop_res[['notyl.cow']][3433,]

"% with minimal increase or decrease"
ecdf(abs(TYL.effect.DFM))(0.1) #same as NI

"% with moderate increase"
ecdf(TYL.effect.DFM)(0.5) - ecdf(TYL.effect.DFM)(0.1) #same as NI

"% with substantial increase"
ecdf(TYL.effect.DFM)(1) - ecdf(TYL.effect.DFM)(0.5) #same as NI


"AFTP"
TYL.effect.AFTP <- AFTP_Prop_res[['tyl.cow']][3433,]-AFTP_Prop_res[['notyl.cow']][3433,]

"% with minimal increase or decrease"
ecdf(abs(TYL.effect.AFTP))(0.1)

"% with moderate increase"
ecdf(TYL.effect.AFTP)(0.5) - ecdf(TYL.effect.AFTP)(0.1)

"% with substantial increase"
ecdf(TYL.effect.AFTP)(1) - ecdf(TYL.effect.AFTP)(0.5)


"#ALL"
TYL.effect.ALL <- ALL_Prop_res[['tyl.cow']][3433,]-ALL_Prop_res[['notyl.cow']][3433,]

"% with minimal increase or decrease"
ecdf(abs(TYL.effect.ALL))(0.1)

"% with moderate increase"
ecdf(TYL.effect.ALL)(0.5) - ecdf(TYL.effect.ALL)(0.1)

"% with substantial increase"
ecdf(TYL.effect.ALL)(1) - ecdf(TYL.effect.ALL)(0.5)

""
"time and magnitude of intervention impact in TYL effect"
"AFTP"
filter(plot.cow.Pr.AFTP[['data']],  days>=113 & days<143) %>% group_by(variable) %>% 
  summarise(d113=nth(value,1), max=max(value), min=min(value), d143=last(value),
                day.max=days[which.max(value)], day.min=days[which.min(value)],
            d143.reduction=nth(value,1)-last(value), lrgst.reduction=nth(value,1)-min(value))

"RWT"
filter(plot.cow.Pr.RWT[['data']],  days>=113 & days<143) %>% group_by(variable) %>% 
  summarise(d113=nth(value,1), max=max(value), min=min(value), d143=last(value),
            day.max=days[which.max(value)], day.min=days[which.min(value)],
            d143.reduction=nth(value,1)-last(value), lrgst.reduction=nth(value,1)-min(value))

"ALL"
filter(plot.cow.Pr.ALL[['data']],  days>=113 & days<143) %>% group_by(variable) %>% 
  summarise(d113=nth(value,1), max=max(value), min=min(value), d143=last(value), d120=nth(value,24*7),
            day.max=days[which.max(value)], day.min=days[which.min(value)],
            d143.reduction=nth(value,1)-last(value), week.reduction=nth(value,1)-nth(value,24*7), lrgst.reduction=nth(value,1)-min(value))

""
"DFM vs NI. They are essentially the same for proportion resistant"
"DFM"
filter(plot.cow.Pr.DFM[['data']],  days==143) %>% group_by(variable) %>% summarise(v=value)
"NI"
filter(plot.cow.Pr.NI[['data']],  days==143) %>% group_by(variable) %>% summarise(v=value)
sink()