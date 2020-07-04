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


#plot the difference in proportion resistant between TYL_intervention and TYL_NI
#Take the difference between individual simulations because parameter values are the same for each simulation (due to random number seed)

#automate loading data and calculating differences
TYLdatafiles<-c("TYL_NI_", "TYL_RWT_", "TYL_NI_", "TYL_AFTP_", "TYL_NI_", "TYL_DFM_") #control (TYL_NI) followed by TYL intervention group.

#load in proportion resistant
for (i in seq(1,length(TYLdatafiles),2)){
  tyl.NI.cow <- read.table(paste("data/",TYLdatafiles[i],"Prop_cow_res.txt", sep=""), sep=",")
  tyl.I.cow <- read.table(paste("data/",TYLdatafiles[i+1],"Prop_cow_res.txt", sep=""), sep=",")
 
  diff.cow <- tyl.I.cow - tyl.NI.cow
  
  label <- paste("TYL_", strsplit(TYLdatafiles[i+1], "_")[[1]][2], "_Prop_res", sep="")
  
  list <- list(tyl.NI.cow, tyl.I.cow, diff.cow)
  names(list) <- c("tyl.NI.cow", "tyl.I.cow", "diff.cow")
  
  assign(label, list)
}
rm(tyl.NI.cow, tyl.I.cow, diff.cow, list, label)

#difference between CON_INT and CON_NI. Note that there is no CON_RWT because TYL can't be withdrawn in CON scenario
#automate loading data and calculating differences
CONdatafiles<-c("NoTYL_NI_", "NoTYL_AFTP_", "NoTYL_NI_", "NoTYL_DFM_")

#load in proportion resistant
for (i in seq(1,length(CONdatafiles),2)){
  con.NI.cow <- read.table(paste("data/",CONdatafiles[i],"Prop_cow_res.txt", sep=""), sep=",")
  con.I.cow <- read.table(paste("data/",CONdatafiles[i+1],"Prop_cow_res.txt", sep=""), sep=",")
  
  diff.cow <- con.I.cow - con.NI.cow
  
  label <- paste("CON_", strsplit(CONdatafiles[i+1], "_")[[1]][2], "_Prop_res", sep="")
  
  list <- list(con.NI.cow, con.I.cow, diff.cow)
  names(list) <- c("con.NI.cow", "con.I.cow", "diff.cow")
  
  assign(label, list)
}
rm(con.NI.cow, con.I.cow, diff.cow, list, label)


#difference between TYL_INT and CON_NI.
#automate loading data and calculating differences
datafiles<-c("NoTYL_NI_", "TYL_RWT_", "NoTYL_NI_", "TYL_AFTP_", "NoTYL_NI_", "TYL_DFM_") #control (NoTYL_NI) followed by TYL intervention group.

#load in proportion resistant
for (i in seq(1,length(datafiles),2)){
  con.NI.cow <- read.table(paste("data/",datafiles[i],"Prop_cow_res.txt", sep=""), sep=",")
  tyl.I.cow <- read.table(paste("data/",datafiles[i+1],"Prop_cow_res.txt", sep=""), sep=",")
  
  diff.cow <- tyl.I.cow - con.NI.cow
  
  label <- paste(strsplit(datafiles[i+1], "_")[[1]][2], "_Prop_res", sep="")
  
  list <- list(con.NI.cow, tyl.I.cow, diff.cow)
  names(list) <- c("con.NI.cow", "tyl.I.cow", "diff.cow")
  
  assign(label, list)
}
rm(con.NI.cow, tyl.I.cow, diff.cow, list, label)



#plotting
scenarios <- c("RWT", "DFM", "AFTP")
#create figure of Prop res in cow for all scenarios
for (i in 1:length(scenarios)){
  TYLdata <- get(paste("TYL_", scenarios[i], "_Prop_res", sep=""))
  
  TYLplot <- func.boxplots(quants, TYLdata[['diff.cow']], days, ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(-1,0.1), breaks=c(seq(-1,0.1,0.1)), expand=c(0,0))
  
  TYLlabel <- paste("TYLplot.cow.Pr.",scenarios[i], sep="")
  assign(TYLlabel, TYLplot)
  
  
  data <- get(paste(scenarios[i], "_Prop_res", sep=""))
  
  plot <- func.boxplots(quants, data[['diff.cow']], days, ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))
  
  label <- paste("plot.cow.Pr.",scenarios[i], sep="")
  assign(label, plot)
  
  
  if (i>1){ #no RWT in CON scenario
  CONdata <- get(paste("CON_", scenarios[i], "_Prop_res", sep=""))
  
  CONplot <- func.boxplots(quants, CONdata[['diff.cow']], days, ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(-1,0.1), breaks=c(seq(-1,0.1,0.1)), expand=c(0,0))
  
  CONlabel <- paste("CONplot.cow.Pr.",scenarios[i], sep="")
  assign(CONlabel, CONplot)
  
  rm(CONplot, CONlabel, CONdata)
  }
  
  rm(TYLplot, TYLlabel, TYLdata, plot, label, data)
}


#plotting Prop Res in each intervention for TYL and CON groups
ggsave('figures/Sup Fig 3 INT_vs_NI_Prop_cow_res.png',
grid.arrange(
grid.arrange(TYLplot.cow.Pr.RWT+
               ylab("Difference in Proportion Resistant")+
               xlab("")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("A"), 
             TYLplot.cow.Pr.DFM+
               ylab("Difference in Proportion Resistant")+
               xlab("")+
               ggtitle("B"), 
             TYLplot.cow.Pr.AFTP+
               ylab("Difference in Proportion Resistant")+
               xlab("Days")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("D"), 
             ncol=1),
grid.arrange(ggplot()+
               theme_bw()+
               theme(panel.border=element_blank())+
               ylab("")+
               ggtitle(""), 
             CONplot.cow.Pr.DFM+
               ylab("")+
               xlab("")+
               ggtitle("C"), 
             CONplot.cow.Pr.AFTP+
               ylab("")+
               xlab("Days")+
               geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
               ggtitle("E"),
             ncol=1),
ncol=2),
width=6*2,
height=5*3,
units="in",
dpi=320)

#########################################

#descriptive stats
sink("results/intervention effects.txt")
"difference between TYL_INT and TYL_NI, CON_INT and CON_NI, at the end of the feeding period"
"RWT, TYL (no CON_RWT)"
RWT.effect.TYL <- unlist(TYL_RWT_Prop_res[['tyl.I.cow']][3433,]-TYL_RWT_Prop_res[['tyl.NI.cow']][3433,])

"% with minimal increase or decrease"
ecdf(abs(RWT.effect.TYL))(0.1)

"% with any increase. maximum difference is 0. none have increase due to RWT"
ecdf(RWT.effect.TYL)(0)
max(RWT.effect.TYL)

"% maximum decrease"
min(RWT.effect.TYL) #32 percentage point decrease


""
"DFM, TYL followed by CON"
DFM.effect.TYL <- unlist(TYL_DFM_Prop_res[['tyl.I.cow']][3433,]-TYL_DFM_Prop_res[['tyl.NI.cow']][3433,])
DFM.effect.CON <- unlist(CON_DFM_Prop_res[['con.I.cow']][3433,]-CON_DFM_Prop_res[['con.NI.cow']][3433,])

"% with minimal increase or decrease"
ecdf(abs(DFM.effect.TYL))(0.1) #99%
ecdf(abs(DFM.effect.CON))(0.1) #99.8%

"% with minimal increase"
ecdf(DFM.effect.TYL)(0.1) - ecdf(DFM.effect.TYL)(0) #46%
ecdf(DFM.effect.CON)(0.1) - ecdf(DFM.effect.CON)(0) #45%

"% with minimal decrease"
ecdf(DFM.effect.TYL)(0) - ecdf(DFM.effect.TYL)(-0.1) #53%
ecdf(DFM.effect.CON)(0) - ecdf(DFM.effect.CON)(-0.1) #55%

""
"AFTP, TYL"
AFTP.effect.TYL <- TYL_AFTP_Prop_res[['tyl.I.cow']][3433,]-TYL_AFTP_Prop_res[['tyl.NI.cow']][3433,]

"% with minimal increase or decrease"
ecdf(abs(AFTP.effect.TYL))(0.1) #94%
sink()