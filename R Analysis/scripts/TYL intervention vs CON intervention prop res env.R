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


#plot the difference in proportion resistant (in pen, feed, water) between TYL and no TYL control for each scenario (NI, DFM, AFTP, RWT, ALL).
#automate loading data and calculating differences
datafiles<-c("NoTYL_NI_", "TYL_NI_", "NoTYL_NI_", "TYL_RWT_", "NoTYL_AFTP_", "TYL_AFTP_", "NoTYL_DFM_", "TYL_DFM_", "NoTYL_ALL_", "TYL_ALL_") #control followed by TYL intervention group. RWT control is the same as NoTYL_NI becuase with no tylosin there can be no withdrawal of tylosin. 

#load in proportion resistant
for (i in seq(1,length(datafiles),2)){
  
  notyl.feed <- read.table(paste("data/",datafiles[i],"Prop_feed_res.txt", sep=""), sep=",")
  tyl.feed <- read.table(paste("data/",datafiles[i+1],"Prop_feed_res.txt", sep=""), sep=",")
  
  notyl.water <- read.table(paste("data/",datafiles[i],"Prop_water_res.txt", sep=""), sep=",")
  tyl.water <- read.table(paste("data/",datafiles[i+1],"Prop_water_res.txt", sep=""), sep=",")
  
  notyl.pen <- read.table(paste("data/",datafiles[i],"Prop_pen_res.txt", sep=""), sep=",")
  tyl.pen <- read.table(paste("data/",datafiles[i+1],"Prop_pen_res.txt", sep=""), sep=",")
  
  diff.feed <- tyl.feed - notyl.feed
  diff.water <- tyl.water - notyl.water
  diff.pen <- tyl.pen - notyl.pen
  
  label <- paste(strsplit(datafiles[i+1], "_")[[1]][2], "_Prop_res", sep="")
  
  list <- list(tyl.feed, notyl.feed, diff.feed, notyl.water, tyl.water, diff.water, notyl.pen, tyl.pen, diff.pen)
  names(list) <- c("tyl.feed", "notyl.feed", "diff.feed", "notyl.water", "tyl.water", "diff.water", "notyl.pen", "tyl.pen", "diff.pen")
  
  assign(label, list)
}
rm(tyl.feed, notyl.feed, diff.feed, notyl.water, tyl.water, diff.water, notyl.pen, tyl.pen, diff.pen, list)


#need days for plotting
days<-read.table("data/NoTYL_NI_Days.txt", sep=",")
days<-as.numeric(days[,1]) - 50 #subtract burn-1n

#parameters needed for functional boxplots
quants<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

ltype=c("dotted", "dotdash", "longdash", "solid", "longdash", "dotdash", "dotted")
lcolor=c("black", "#F8766D", "#7CAE00", "#619CFF", "#7CAE00", "#F8766D", "black")


#plotting Prop Res in each intervention. pen, feed, water
scenarios <- c("NI", "RWT", "DFM", "AFTP", "ALL")

#create figure of Prop res in pen for all scenarios
for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_Prop_res", sep=""))
  
  plot <- func.boxplots(quants, data[['diff.pen']], days, ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))
  
  label <- paste("plot.pen.Pr.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}


#for feed, plot only at the end of every day to capture long-term trends and not focus on fluctations from feed being replaced every day
#create figure of Prop res in feed for all scenarios
for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_Prop_res", sep=""))
  data <- data[['diff.feed']][seq(0,3433,24),]
  
  plot <- func.boxplots(quants, data, days[seq(0,3433,24)], ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))
  
  label <- paste("plot.feed.Pr.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}

#for water, also plot end of every day. Water is cleaned every two weeks
#create figure of Prop res in water for all scenarios
for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_Prop_res", sep=""))
  data <- data[['diff.water']][seq(0,3433,24),]
  
  plot <- func.boxplots(quants, data, days[seq(0,3433,24)], ltype, lcolor, "Difference in Proportion Resistant")+
    scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))
  
  label <- paste("plot.water.Pr.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}


#plot together
ggsave('figures/Sup Fig 1 Prop_env_res_interventions.png',
grid.arrange(
  #pen plots, need scenario titles
  plot.pen.Pr.RWT+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("Pen: Diff. Prop. Resistant")+
     xlab("")+
     ggtitle("RWT")+
     theme(plot.title=element_text(hjust=0.5)),
   plot.pen.Pr.DFM+
     ylab("")+
     xlab("")+
     ggtitle("DFM")+
     theme(plot.title=element_text(hjust=0.5)), 
   plot.pen.Pr.AFTP+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("")+
     ggtitle("AFTP")+
     theme(plot.title=element_text(hjust=0.5)), 
   plot.pen.Pr.ALL+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("")+
     ggtitle("ALL")+
     theme(plot.title=element_text(hjust=0.5)),

  #feed plots  
  plot.feed.Pr.RWT+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("Feed: Diff. Prop. Resistant")+
     xlab("")+
     ggtitle(""), 
   plot.feed.Pr.DFM+
     ylab("")+
     xlab("")+
     ggtitle(""), 
   plot.feed.Pr.AFTP+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("")+
     ggtitle(""), 
   plot.feed.Pr.ALL+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("")+
     ggtitle(""),

  #water plots           
  plot.water.Pr.RWT+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("Water: Diff. Prop. Resistant")+
     xlab("Days")+
     ggtitle(""), 
   plot.water.Pr.DFM+
     ylab("")+
     xlab("Days")+
     ggtitle(""), 
   plot.water.Pr.AFTP+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("Days")+
     ggtitle(""), 
   plot.water.Pr.ALL+
     geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
     ylab("")+
     xlab("Days")+
     ggtitle(""),
  ncol=4),
width=6*4,
height=5*3,
units="in",
dpi=320)