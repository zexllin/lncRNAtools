##' @import ggplot2
##' @author zexl
##' @title PCA plot
##' @description Principal component analysis and visualization
##' @param rnaExpr a matrix ,Expression data after standardization 
##' @param group Sample grouping
##' @param datatype prefix for Save file names
##' @param save Output file type,default is not save
##' @export
##' @author zexl
PCA_GGplot <- function(rnaExpr = rnaExpr, group,datatype='', save=''){

  project.pca <- prcomp(t(rnaExpr)) #prcomp 
  project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
  xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")
  ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")
  data_x <- data.frame(samplenames=rownames(project.pca$x), project.pca$x)
  p <- ggplot(data_x, aes(PC1,PC2))+geom_point(aes(color=group),size=3)+
    #coord_equal(ratio=1)+
    geom_text(label=data_x$samplenames,hjust = -0.3)+
    #geom_text(label=data_x$samplenames,colour="black",size=3,vjust = -0.1, hjust = -0.1)+
    labs(x =xlab,y = ylab,title = "Principal components analysis bi-plot")+
    theme(axis.title.x = element_text(size = 27, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 27, vjust = 0.5, hjust = 0.5)
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank()
    )+
    theme(plot.title = element_text(size = 27,hjust = 0.5))+ 
    theme(axis.text=element_text(size=20),axis.title=element_text(size=27))+
    theme(legend.text = element_text(size = 20),legend.title=element_text(size = 20))

  tifname = paste(datatype,'PCA.tif',sep='_')
  pdfname = paste(datatype,'PCA.pdf',sep='_')
  if(save=='tif'){
    tiff(tifname,width=200, height=150, units='mm',res=300)
    print(p)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=8, height=6)
    print(p)
    dev.off()
  }
  print(p)
}

##' @title Combination of box line diagram and density diagram
##' @description Box line diagram and density diagram before and after standardization
##' @param countDataRaw Expression data
##' @param exprSet data after normlizing
##' @param datatype prefix for Save file names
##' @export
##' @author zexl
normPlot <- function(countDataRaw, exprSet,datatype){


  logcountDataRaw <- log2(countDataRaw+1)
  #tiff(paste(datatype,'Normalization.tif',sep='_'),width=200, height=170, units='mm',res=300)
  pdf(file=paste0(datatype,'_Normalization.pdf'), width=12, height=8)
  par(cex = 0.7)
  n.sample=ncol(countDataRaw)
  samplenames <- colnames(countDataRaw)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2) #rainbow 
  par(mfrow=c(2,2))
  boxplot(logcountDataRaw, col = cols,main="expression value",las=2,axes = T)
  #axis(2,seq(from=0,length=5,by=5),lwd = 1,cex.axis=1)
  #axis(1,seq(from=0,length=66,by=1),lwd = 1,cex.axis=1)

  boxplot(exprSet, col = cols,main="expression value",las=2,axes = T)
  #axis(2,seq(from=0,length=12,by=5),lwd = 1,cex.axis=1)
  #axis(1,seq(from=0,length=66,by=1),lwd = 1,cex.axis=1)
  plot(density(logcountDataRaw[,1]), col=cols[1], lwd=2,  las=2,
       main="", xlab="")
  title(main="Raw data", xlab="Log-Raw data")
  abline(v=1, lty=3)
  for (i in 2:n.sample){
    den <- density(logcountDataRaw[,i])
    lines(den$x, den$y, col=cols[i], lwd=2)
  }
  #legend(0.6*max(logcountDataRaw[,1]),max(density(logcountDataRaw[,1])$y)+0.15,
  #samplenames, text.col=cols, bty="n",cex = 0.8,y.intersp=0.25)
  plot(density(exprSet[,1]), col=cols[1], lwd=2,  las=2,
       main="", xlab="")
  title(main="Filtered data", xlab="Log-cpm")
  abline(v=1, lty=3)
  for (i in 2:n.sample){
    den <- density(exprSet[,i])
    lines(den$x, den$y, col=cols[i], lwd=2)
  }
  #legend(0.6*max(exprSet[,1]),max(density(exprSet[,1])$y)+0.03, samplenames, text.col=cols, bty="n",cex = 0.8,y.intersp=0.25)
  dev.off()
}


##' @title volcano plot
##' @description A volcano plot showing differentially expressed genes/miRNAs
##' @param deg.all a dataframe generated from  \code{\link{DEAnalysis}}
##'     containing all genes of analysis no matter they are differentially
##'     expressed or not
##' @param fc a numeric value specifying the threshold of fold change
##' @param pval a nuemric value specifying the threshold of p value
##'             default is 0.05 
##' @param nlab a nuemric value,Specify the number of genes to display the gene name
##' @param data.type prefix for Save file names
##' @param save Output file type,default is not save
##' @export
##' @author zexl
VolcanoPlot <- function(deg.all, fc=2, pval=0.05, data.type = 'mRNA',nlab=0,save='') {

  #geneList <- na.omit(deg.all)
  geneList <- deg.all[which((!is.na(deg.all$logFC))&(!(is.na(deg.all$PValue)))),]
  geneList$threshold <- c()
  geneList$threshold[geneList$logFC>log(fc,2) & geneList$PValue<pval] <- "up"
  geneList$threshold[abs(geneList$logFC)<=log(fc,2) | geneList$PValue>=pval] <- "nosig"
  geneList$threshold[geneList$logFC < -log(fc,2) & geneList$PValue<pval] <- "down"
  geneList$threshold <- as.factor(geneList$threshold)

  if(!'symbol' %in%colnames(deg.all)) {
    geneList$symbol<-rownames(geneList)
  }
  geneList$logPvalue <- -log10(geneList$PValue)
  texts<-geneList[which(abs(geneList$logFC)>=1),]

  if(sum(is.na(texts[order(texts$PValue)[0:nlab],'symbol']))>0){
    message('nlab:0 to nlab symbol contain NA or no symbol')
  }

  lim <- max(max(geneList$logFC), abs(min(geneList$logFC)))+0.5
  volcano <- ggplot(data=geneList, aes(x=geneList$logFC,
                                       y = -log10(geneList$PValue)))
  p<-volcano+geom_point(aes(color=geneList$threshold), alpha=1,size=2) + #alpha 透明度
    xlab("log2(Fold Change)") + ylab("-log10(PValue)") +
    scale_colour_manual(breaks = geneList$threshold,
                        values = c('red','black','green3')) + xlim(c(-lim,lim)) +
    geom_vline(xintercept = c(-log(fc,2),log(fc,2)),
               color='darkgreen', linetype=3,size=0.6) +
    geom_hline(yintercept = -log(pval,10), color='darkgreen',linetype=3,size=0.6)+
    geom_text(data=texts[order(texts$PValue)[0:nlab],],
              aes(x=logFC,y=logPvalue,label=symbol),hjust=0,vjust=0,size=5)+
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='black'),
          panel.background = element_blank())+
    theme(legend.title=element_blank())+
    #theme(legend.key.size=unit(2,'cm'));
    #theme(legend.key.width=unit(2,'cm'));
    theme(legend.text = element_text(size = 30),legend.position = 'right')+ #raw 12
    theme(axis.text=element_text(size=30), # raw 13
          axis.title=element_text(size=30))  # raw 16
  
  tifname = paste(data.type,"_VolCano",".tif",sep="")
  pdfname = paste(data.type,"_VolCano",".pdf",sep="")
  if(save=='tif'){
    tiff(tifname,width=200, height=150, units='mm',res=300)
    print(p)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=8, height=6)
    print(p)
    dev.off()
  }
  print(p)
}

##' @title heatmap
##' @description A heatmap showing unsupervised hierarchical clustering of
##'     DE genes/miRNAs
##' @import pheatmap
##' @param exprdata a dataframe for plot heatmap
##' @param Group sample group
##' @param datatype prefix for Save file names
##' @param save Output file type,default is not save
##' @export
##' @author zexl
HeatmapPlot <- function(exprdata,Group = group,datatype='',save=''){
  
  degDa <-exprdata
  annotation_col <- data.frame('Sample' = Group)
  rownames(annotation_col) <- colnames(degDa)
  #annotation_row<- data.frame(diff=ifelse(deg.sig[order(deg.sig$logFC),]$logFC > 1,"UP","DOWN" ))
  #rownames(annotation_row)<-rownames(degDa)
  pdfname=paste(datatype,"DegHeatMap.pdf",sep='_')
  tifname=paste(datatype,"DegHeatMap.tif",sep='_')
  #pdf(file=pdfname, width=12, height=8)
  p<-pheatmap(degDa,
              #main = 'heatmap', # title
              scale = 'row', #column row none 
              annotation_col = annotation_col, 
              #annotation_row = annotation_row, 
              #legend_labels = NA,
              cluster_cols = FALSE,          
              #cluster_rows = FALSE,         
              #clustering_method = "complete", # complete average median 
              show_rownames = F, 
              #show_colnames = F, 
              #gaps_row = 1169,
              fontsize = 24,
              border_color = NA,
              angle_col=45)
  if(save=='tif'){
    tiff(tifname,width=200, height=150, units='mm',res=300)
    print(p)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=8, height=6)
    print(p)
    dev.off()
  }
  print(p)
}

