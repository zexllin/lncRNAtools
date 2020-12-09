##' @title Plots for lncRNA cis analysis
##' @description Plots for lncRNA cis analysis,lnc-pc cor pvalue distance data visualization
##' @param NearPlotData a data frame of lncRNA cis analysis
##' @param Output file type,default is not save
##' @param data.dir Specify save path,default is Current working path
##' @export
##' @author zexl
##' @import grid
##' @import gtable

plotcis<-function(NearPlotData,data.dir='.',save=''){
  
  NearPlotData$cor<-round(as.numeric(NearPlotData$cor),5)
  NearPlotData$PValue<-round(as.numeric(NearPlotData$PValue),5)
  lim<-c(round(min(NearPlotData$cor),1)-0.05,round(max(NearPlotData$cor),1)+0.05)
  p1 <- ggplot(data=NearPlotData, mapping=aes(x=distance, y=1:nrow(NearPlotData),color=PValue,size=cor))
  p1 <- p1+geom_point()+xlab('Distance(bp)')+ylab('LncRNA Name')+
    scale_size("COR",limits=lim,range=c(2,8))+
    scale_y_continuous(breaks = c(1:nrow(NearPlotData)), label = NearPlotData$name_Lnc)+
    scale_y_discrete(limits=NearPlotData$name_Lnc)+
    theme(axis.text.y=element_text(size=14), #raw 14
          axis.title=element_text(size=16),
          axis.text.x=element_text(size=14)) + # raw 14
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p2 <- ggplot(data=NearPlotData, mapping=aes(x=distance, y=1:nrow(NearPlotData),color=PValue,size=cor))
  p2 <- p2+geom_point()+xlab('Distance(bp)')+ylab('NearGene Name')+
    scale_size("COR",limits=lim,range=c(2,8))+
    scale_y_continuous(breaks = c(1:nrow(NearPlotData)), label = NearPlotData$name_Lnc)+
    scale_y_discrete(limits=NearPlotData$name_NearGene)+
    theme(axis.text.y=element_text(size=14), #raw 14
          axis.title=element_text(size=16),
          axis.text.x=element_text(size=14)) + # raw 14
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))

  tifname = paste(data.dir,'lncRNA_cis.tif',sep='/')
  pdfname = paste(data.dir,'lncRNA_cis.pdf',sep='/')
  if (save=='tif') {
    tiff(tifname,width=300, height=200, units='mm',res=300)
    ggplot2.two_y_axis(p1,p2)
    dev.off()
  }
  if(outpdf==TRUE){
    pdf(pdfname,width=8, height=6)
    ggplot2.two_y_axis(p1,p2)
    dev.off()
  }
  ggplot2.two_y_axis(p1,p2)
}

