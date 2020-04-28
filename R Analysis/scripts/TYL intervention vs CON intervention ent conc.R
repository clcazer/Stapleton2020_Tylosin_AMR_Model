#title: "TYL model Ent conc figures"
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


#enterococcus concentrations
#only save data from each day (not every hour) to save memory
datafiles<-c("NoTYL_NI_", "TYL_NI_", "NoTYL_NI_", "TYL_RWT_", "NoTYL_AFTP_", "TYL_AFTP_", "NoTYL_DFM_", "TYL_DFM_", "NoTYL_ALL_", "TYL_ALL_") #control followed by TYL intervention group. RWT control is the same as NoTYL_NI becuase with no tylosin there can be no withdrawal of tylosin. 

for (i in seq(1,length(datafiles),2)){
  notyl.cow <- read.table(paste("data/",datafiles[i],"Cow_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  tyl.cow <- read.table(paste("data/",datafiles[i+1],"Cow_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  
  notyl.feed <- read.table(paste("data/",datafiles[i],"Feed_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  tyl.feed <- read.table(paste("data/",datafiles[i+1],"Feed_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  
  notyl.water <- read.table(paste("data/",datafiles[i],"Water_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  tyl.water <- read.table(paste("data/",datafiles[i+1],"Water_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  
  notyl.pen <- read.table(paste("data/",datafiles[i],"Pen_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  tyl.pen <- read.table(paste("data/",datafiles[i+1],"Pen_total_conc.txt", sep=""), sep=",")[seq(0,3433,24),]
  
  
  diff.cow <- log10(tyl.cow) - log10(notyl.cow)
  diff.feed <- log10(tyl.feed) - log10(notyl.feed)
  diff.water <- log10(tyl.water+0.0001) - log10(notyl.water+0.0001) #some water values have 0 CFU/g. add in small amount to avoid infinite/NA values
  diff.pen <- log10(tyl.pen) - log10(notyl.pen)
  
  label <- paste(strsplit(datafiles[i+1], "_")[[1]][2], "_total_conc", sep="")
  
  list <- list(tyl.cow, notyl.cow, diff.cow, diff.feed, diff.water, diff.pen)
  names(list) <- c("tyl.cow", "notyl.cow", "diff.cow", "diff.feed", "diff.water", "diff.pen")
  
  assign(label, list)

}

rm(notyl.cow, tyl.cow, diff.cow, tyl.feed, notyl.feed, diff.feed, notyl.water, tyl.water, diff.water, notyl.pen, tyl.pen, diff.pen, list, label)

#is there a difference between ent conc in tyl and con?
min(NI_total_conc[['diff.cow']])
max(NI_total_conc[['diff.cow']])
mean(colMeans(NI_total_conc[['diff.cow']]))
#very little difference when considered in log10

#parameters for functional boxplots
quants<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
ltype=c("dotted", "dotdash", "longdash", "solid", "longdash", "dotdash", "dotted")
lcolor=c("black", "#F8766D", "#7CAE00", "#619CFF", "#7CAE00", "#F8766D", "black")

#load day variable
days<-read.table("data/NoTYL_NI_Days.txt", sep=",")
days<-as.numeric(days[,1]) - 50 #subtract burn-in
days <- days[seq(0,3433,24)]

#plotting
scenarios <- c("NI", "RWT", "DFM", "AFTP", "ALL")
#create figure of Ent conc in cow/pen/feed/water for all scenarios

for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_total_conc", sep=""))
  plot <- func.boxplots(quants, data[['diff.cow']], days, ltype, lcolor, "Difference in Enterococci Concentrations")
  label <- paste("plot.cow.EntConcdiff.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}

cow <- grid.arrange(plot.cow.EntConcdiff.NI+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("Cattle: Diff. Enterococci Conc.")+
                      xlab("")+
                      ggtitle("NI")+
                      theme(plot.title=element_text(hjust=0.5)),
                    plot.cow.EntConcdiff.RWT+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("")+
                      xlab("")+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ggtitle("RWT")+
                      theme(plot.title=element_text(hjust=0.5)), 
                    plot.cow.EntConcdiff.DFM+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("")+
                      xlab("")+
                      ggtitle("DFM")+
                      theme(plot.title=element_text(hjust=0.5)), 
                    plot.cow.EntConcdiff.AFTP+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("")+
                      xlab("")+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ggtitle("AFTP")+
                      theme(plot.title=element_text(hjust=0.5)), 
                    plot.cow.EntConcdiff.ALL+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("")+
                      xlab("")+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ggtitle("ALL")+
                      theme(plot.title=element_text(hjust=0.5)),
                    ncol=5)

for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_total_conc", sep=""))
  plot <- func.boxplots(quants, data[['diff.pen']], days, ltype, lcolor, "Difference in Enterococci Concentrations")
  label <- paste("plot.pen.EntConcdiff.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}

pen <- grid.arrange(plot.pen.EntConcdiff.NI+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("Pen: Diff. Enterococci Conc.")+
                      xlab("")+
                      ggtitle(""),
                    plot.pen.EntConcdiff.RWT+
                      scale_y_continuous(limits=c(-0.125,0))+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ylab("")+
                      xlab("")+
                      ggtitle(""), 
                    plot.pen.EntConcdiff.DFM+
                      scale_y_continuous(limits=c(-0.125,0))+
                      ylab("")+
                      xlab("")+
                      ggtitle(""), 
                    plot.pen.EntConcdiff.AFTP+
                      scale_y_continuous(limits=c(-0.125,0))+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ylab("")+
                      xlab("")+
                      ggtitle(""), 
                    plot.pen.EntConcdiff.ALL+
                      scale_y_continuous(limits=c(-0.125,0))+
                      geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                      ylab("")+
                      xlab("")+
                      ggtitle(""),
                    ncol=5)

for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_total_conc", sep=""))
  plot <- func.boxplots(quants, data[['diff.feed']], days, ltype, lcolor, "Difference in Enterococci Concentrations")
  label <- paste("plot.feed.EntConcdiff.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}

feed <- grid.arrange(plot.feed.EntConcdiff.NI+
                       scale_y_continuous(limits=c(-0.125,0))+
                       ylab("Feed: Diff. Enterococci Conc.")+
                       xlab("")+
                       ggtitle(""),
                     plot.feed.EntConcdiff.RWT+
                       scale_y_continuous(limits=c(-0.125,0))+
                       geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                       ylab("")+
                       xlab("")+
                       ggtitle(""), 
                     plot.feed.EntConcdiff.DFM+
                       scale_y_continuous(limits=c(-0.125,0))+
                       ylab("")+
                       xlab("")+
                       ggtitle(""), 
                     plot.feed.EntConcdiff.AFTP+
                       scale_y_continuous(limits=c(-0.125,0))+
                       geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                       ylab("")+
                       xlab("")+
                       ggtitle(""), 
                     plot.feed.EntConcdiff.ALL+
                       scale_y_continuous(limits=c(-0.125,0))+
                       geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                       ylab("")+
                       xlab("")+
                       ggtitle(""),
                     ncol=5)

#some water values must have 0 CFU/g--resulting in log10 = Inf
for (i in 1:length(scenarios)){
  data <- get(paste(scenarios[i], "_total_conc", sep=""))
  plot <- func.boxplots(quants, data[['diff.water']], days, ltype, lcolor, "Difference in Enterococci Concentrations")
  label <- paste("plot.water.EntConcdiff.",scenarios[i], sep="")
  assign(label, plot)
  rm(plot, label, data)
}

water <- grid.arrange(plot.water.EntConcdiff.NI+
                        scale_y_continuous(limits=c(-0.125,0))+
                        ylab("Water: Diff. Enterococci Conc.")+
                        xlab("Days")+
                        ggtitle(""),
                      plot.water.EntConcdiff.RWT+
                        scale_y_continuous(limits=c(-0.125,0))+
                        geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                        ylab("")+
                        xlab("Days")+
                        ggtitle(""), 
                      plot.water.EntConcdiff.DFM+
                        scale_y_continuous(limits=c(-0.125,0))+
                        ylab("")+
                        xlab("Days")+
                        ggtitle(""), 
                      plot.water.EntConcdiff.AFTP+
                        scale_y_continuous(limits=c(-0.125,0))+
                        geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                        ylab("")+
                        xlab("Days")+
                        ggtitle(""), 
                      plot.water.EntConcdiff.ALL+
                        scale_y_continuous(limits=c(-0.125,0))+
                        geom_vline(aes(xintercept=113), linetype="longdash", size=1.25)+
                        ylab("")+
                        xlab("Days")+
                        ggtitle(""),
                      ncol=5)


ggsave('figures/Sup Fig 2 Ent_conc_env_interventions.png', 
grid.arrange(cow,
             pen,
             feed,
             water,
      ncol=1),
width=6*5,
height=5*4,
dpi=320)


#descriptive analysis
sink("results/TYL effect on enterococci concentrations.txt")
"effect of TYL (NI scenario) on enterococci concentrations (quantiles) by the end of treatment; TYL-CON, TYL, CON"
quantile(NI_total_conc[['diff.cow']][143,], quants)
quantile(NI_total_conc[['tyl.cow']][143,], quants)
quantile(NI_total_conc[['notyl.cow']][143,], quants)
"very little difference in percentiles between tyl and con. difference between counterfactuals (log10 difference) also small"

""
"the largest difference occurs around day 30 (see figure), but it is still a small difference"
"day 30, NI scenario, enterococci concentration quantiles for TYL-CON, TYL, CON"
min(NI_total_conc[['diff.cow']])
quantile(NI_total_conc[['diff.cow']][30,], quants)
quantile(NI_total_conc[['tyl.cow']][30,], quants)
quantile(NI_total_conc[['notyl.cow']][30,], quants) 

""
"DFM, by design reduces the concentration of enterococci, but did not impact the TYL- CON difference"
"DFM, median for TYL, CON at end of treatment"
quantile(log10(DFM_total_conc[['tyl.cow']][143,]),0.5)
quantile(log10(DFM_total_conc[['notyl.cow']][143,]),0.5)

"ALL, median for TYL, CON at end of treatment"
quantile(log10(ALL_total_conc[['tyl.cow']][143,]),0.5)
quantile(log10(ALL_total_conc[['notyl.cow']][143,]),0.5)

"NI, median for TYL, CON at end of treatment"
quantile(log10(NI_total_conc[['tyl.cow']][143,]),0.5)
quantile(log10(NI_total_conc[['notyl.cow']][143,]),0.5)
sink()