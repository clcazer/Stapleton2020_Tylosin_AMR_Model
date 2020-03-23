func.boxplots<-function(quantiles, data, time, linetype, linecolor, log10, proportion, ylabel){
  require(ggplot2)
  require(reshape2)
  
  summary <- data.frame(t(apply(data, 1, quantile, probs=quantiles)))
  summary$days <- time
  
  plot <- melt(summary, id.vars="days")
  
  p<-ggplot(plot, aes(x=days, y=value, group=variable))+
    geom_line(aes(linetype=variable, color=variable, size=variable))+
    
    {if (log10==TRUE){
      scale_y_log10(limits=c(10^0,10^7), breaks=c(10^seq(0,7,1)), labels=c(seq(0,7,1)), expand=c(0,0))} else if (proportion==TRUE){
        scale_y_continuous(limits=c(0,1), breaks=c(seq(0,1,0.1)), expand=c(0,0))} else{
          scale_y_continuous(limits=c(NA,NA), expand=c(0,0))}}+
    
    ylab(ylabel)+
    scale_x_continuous(limits=c(0, ceiling(max(days))), breaks=c(seq(0,ceiling(max(days)),20)), labels=c(seq(0,ceiling(max(days)),20)), expand=c(0,0))+
    xlab('Days')+
    scale_linetype_manual(values=linetype)+
    scale_color_manual(values=linecolor)+
    scale_size_manual(values=c(1,1,1,2,1,1,1))+
    theme_bw(base_size = 20)+
    theme(legend.position = "none")
  
  p
}
