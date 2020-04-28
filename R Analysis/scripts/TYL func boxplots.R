#title: "functional boxplot"
#author: "Casey Cazer"
#Last update: "April 28, 2020"

#returns a ggplot of specific percentiles over time for a given input

#Parameters
##percentiles: vector of percentiles (e.g. 0.95) to calculate
##data: row for each timepoint, column for each simulation. Percentiles of the simulations are calculated at each timepoint
##time: vector of timepoints
##ltype: vector of linetypes, length = percentiles
##lcolor: vector of linecolors, length = percentiles
##ylabel: character string for Y-axis label

#Returns ggplot functional boxplot, graphing the percentiles over time
#custom y-axis using scale_y can be added after

func.boxplots<-function(percentiles, data, time, linetype, linecolor, ylabel){
  require(ggplot2)
  require(reshape2)
  
  summary <- data.frame(t(apply(data, 1, quantile, probs=percentiles)))
  summary$days <- time
  
  plot <- melt(summary, id.vars="days")
  
  p<-ggplot(plot, aes(x=days, y=value, group=variable))+
    geom_line(aes(linetype=variable, color=variable, size=variable))+
    ylab(ylabel)+
    scale_x_continuous(limits=c(0, ceiling(max(plot$days))), breaks=c(seq(0,ceiling(max(plot$days)),20)), labels=c(seq(0,ceiling(max(plot$days)),20)), expand=c(0,0))+
    xlab('Days')+
    scale_linetype_manual(values=linetype)+
    scale_color_manual(values=linecolor)+
    scale_size_manual(values=c(1,1,1,2,1,1,1))+
    theme_bw(base_size = 20)+
    theme(legend.position = "none")
  
  p
}
