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
days<-read.table("data/NoTYL_NI_LGpen_Days.txt", sep=",")
days<-as.numeric(days[,1]) - 50 #subtract burn-1n

quants<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99) #percentiles

#line colors and styles
ltype=c("dotted", "dotdash", "longdash", "solid", "longdash", "dotdash", "dotted")
lcolor=c("black", "#F8766D", "#7CAE00", "#619CFF", "#7CAE00", "#F8766D", "black")


#TYL concentration in LI####
#import data: each timepoint (1 hour) in a row; one column per simulation (1000 simulations per scenario)
tyl.conc.ni <- read.table("data/TYL_NI_LGpen_TYL_li_conc.txt", sep=",")

#summarize the simulations, get percentiles (quants) at each timepoint
tyl.conc.summary<-data.frame(t(apply(tyl.conc.ni, 1, quantile, probs=quants)))

#functional boxplot of TYL LI concentrations
ggsave('figures/Fig LGPen3 tyl_conc.png',
func.boxplots(quants, tyl.conc.ni, days, ltype, lcolor, "Tylosin Concentration (ug/mL)" )+
  scale_y_continuous(limits=c(0,2), breaks=c(seq(0,2,0.25)), expand=c(0,0)),
device="png",
width=6,
height=5,
units="in",
dpi=320)

sink("results/tyl intestine concentrations LGpen.txt")
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


#plot the difference in proportion resistant between TYL and no TYL control for each scenario (NI, DFM, AFTP, RWT, ALL) ####
#Take the difference between individual simulations because parameter values are the same for each simulation (due to random number seed)
#check assumption that parameter values are equal for each simulation
NoTYL_NI_LGpen_param <- read.table("data/NoTYL_NI_LGpen_MC_parameters.txt", sep=",")
TYL_NI_LGpen_param <- read.table("data/TYL_NI_LGpen_MC_parameters.txt", sep=",")
all.equal(NoTYL_NI_LGpen_param, TYL_NI_LGpen_param) #correct
rm(NoTYL_NI_LGpen_param, TYL_NI_LGpen_param)

#automate loading data and calculating differences (note that No TYL NI scenario == No TYL RWT scenario)
datafiles<-c("NoTYL_NI_LGpen_", "TYL_NI_LGpen_", "NoTYL_NI_LGpen_", "TYL_RWT_LGpen_", "NoTYL_AFTP_LGpen_", "TYL_AFTP_LGpen_", "NoTYL_DFM_LGpen_", "TYL_DFM_LGpen_") #control followed by TYL intervention group. RWT control is the same as NoTYL_NI becuase with no tylosin there can be no withdrawal of tylosin. 

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

sink("results/TYL effect LGpen.txt")
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

sink("results/TYL effect LGpen.txt", append=TRUE)
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
scenarios <- c("NI", "RWT", "DFM", "AFTP")
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
ggsave('figures/Fig LGPen4 Prop_cow_res_interventions.png',
grid.arrange(plot.cow.Pr.RWT+
               ylab("Difference in Proportion Resistant")+
               xlab("")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("A"), 
             plot.cow.Pr.DFM+
               ylab("Difference in Proportion Resistant")+
               xlab("")+
               ggtitle("B"), 
             plot.cow.Pr.AFTP+
               ylab("Difference in Proportion Resistant")+
               xlab("Days")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("C"), 
             ncol=1),
width=6*1,
height=5*3,
units="in",
dpi=320)


#descriptive stats
sink("results/TYL effect LGpen.txt", append=TRUE)
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

""
"DFM vs NI. They are essentially the same for proportion resistant"
"DFM"
filter(plot.cow.Pr.DFM[['data']],  days==143) %>% group_by(variable) %>% summarise(v=value)
"NI"
filter(plot.cow.Pr.NI[['data']],  days==143) %>% group_by(variable) %>% summarise(v=value)
sink()