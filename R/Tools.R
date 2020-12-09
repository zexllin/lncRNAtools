##' @title write difference analysis results to csv
##' @description Write the result of differential analysis to 
##'     CSV to avoid excel turning on gene name into month,
##'      and combine the result of differential analysis with counts
##' @param wfile a datafram need to write
##' @param c1name The row name is the gene name, the row name is write to 
##'     the first column, and the column name is specified
##' @param outfileprex prefix for Save file names
##' @param Counts counts matrix,normdata result of DEAnalysis function, 
##'     or raw count matrix
##' @export
##' @author zexl
writeDEdata<-function(wfile,c1name='',outfileprex='',Counts=NULL){
  if('symbol' %in% colnames(wfile)){

    wfile$symbol<-unlist(lapply(wfile$symbol,function(x) return(paste0("=\"",x,"\""))))
  }
  if (!is.null(Counts))
    wfile<-cbind(wfile,Counts[match(rownames(wfile), rownames(Counts)),])
  cbindfile <- cbind(rownames(wfile), wfile)
  colnames(cbindfile)<-c(c1name,colnames(wfile))
  write.csv(cbindfile,file=paste0(outfileprex,'.csv'),row.names = F)
}

##' @import scatterplot3d
##' @title 3D PCA
##' @description PCA analysis and plot 3D PCA
##' @param exprdata matrix for PCA
##' @param Group samples group
##' @param title Title of PCA map
##' @export
##' @author zexl
plot_3dPCA <-function(exprdata,Group,title='3dPCA'){
  project.pca <- prcomp(t(exprdata)) #prcomp 
  project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
  data_x <- data.frame(samplenames=rownames(project.pca$x), project.pca$x)
  
  Level = unique(Group)
  cors<-data.frame(row.names = Level,cor=rainbow(length(Level)),stringsAsFactors = F)
  corlor <- cors$cor[match(Group,rownames(cors))]
  
  plotPCA<-with(data_x,scatterplot3d(data_x[,2:4] ,color=corlor,pch=20,font.axis = 2,
                                     font.lab = 2,cex.lab=2,cex.symbols = 2))
  
  legend(plotPCA$xyz.convert(2,1,3),pch=20,yjust=0,legend=Level,title=title,
         cex = 1.2, bty = 'n', xjust = 0.65, horiz = F,
         col =unique(corlor))
}

##' @title Remove outliers
##' @description separation vector outliers replace outliers
##' @param datavector target vector
##' @param qnum Specify quantile for separation outliers,default is 0.75
##' @return Vectors with outliers process
##' @export
##' @author zexl
outline<-function(datavector,qnum=0.75){
  if(qnum<0.5){return(datavector)}
  Qu<-quantile(datavector,qnum)
  Ql<-quantile(datavector,1-qnum)
  IQR <- Qu-Ql
  Up<-Qu+1.5*IQR
  Low<-Ql-1.5*IQR
  ls<-lapply(datavector,function(x){return(min(max(Low,x),Up))})
  return(unlist(ls))
}

##' @title write DeGene
##' @description write DeGene that sorted with P value top 10
##' @export
##' @author zexl
degene.report <- function(deall.sig,rownumber=10,data.type=''){
  if(substring(rownames(deall.sig)[1],1,3)=='hsa'){
    deall.sig$miRNA<-rownames(deall.sig)
    cols<-c('miRNA','logFC','PValue')
  }
  else{
    cols<-c('symbol','logFC','PValue')
  }
  de.up<-deall.sig[which(deall.sig$PValue<0.05 & deall.sig$logFC>=1),cols]
  cat(paste0('total up gene:',nrow(de.up),'\n'))
  de.up<-de.up[order(abs(de.up$logFC),decreasing=T)[1:rownumber],]
  
  de.down<-deall.sig[which(deall.sig$PValue<0.05 & deall.sig$logFC <= -1),cols]
  cat(paste0('total down gene:',nrow(de.down),'\n'))
  de.down<-de.down[order(abs(de.down$logFC),decreasing=T)[1:rownumber],]
  
  outtable<-cbind(de.up,de.down)
  rownames(outtable)<-c(1:nrow(outtable))
  if(data.type!='')
    write.xlsx(outtable,file = paste0(data.type,'_dereport_top',rownumber,'.xlsx'),row.names = T,showNA=F)
  else
    return(outtable)
  
}

##' @title ROC curve
##' @description ROC curve drawing
##' @param rdata pROC results
##' @param data.type prefix for Save file names
##' @param save Output file type,default is not save
##' @export
##' @author zexl
rocplot<-function(rdata,data.type='',save = ''){
  df<-data.frame(x=1-rdata$specificities/100,y=rdata$sensitivities/100)
  p <- ggplot(df) + 
    geom_path(aes(x = x, y = y), size=1) + 
    labs(x = "1 - Specificities", y = "sensitivities", title ="ROC") +
    theme(plot.title = element_text(face = 'bold',size=15))+
    theme_bw()+ #The background is white
    theme(legend.position="none", 
          axis.text.x=element_text(colour="black",size=25), 
          axis.text.y=element_text(size=25,face="plain"), 
          axis.title.y=element_text(size = 25,face="plain"), 
          axis.title.x=element_text(size = 25,face="plain"),
          #axis.title.x=element_blank(),
          plot.title = element_text(size=30,face="bold",hjust = 0.5), 
          panel.grid.major = element_blank(), #Do not show gridlines
          panel.grid.minor = element_blank())+
    theme(axis.line=element_line(colour="black",size = 1))+
    theme(panel.background=element_blank())+
    geom_abline(intercept=0,slope=1 ,size=1, linetype="dashed")+
    geom_text(aes(x=0.75,y=0.45,label=paste0('AUC ',round(rdata$auc/100,4))),size=8)
  tifname = paste0(data.type,'_ROC.tif')
  pdfname = paste0(data.type,'_ROC.pdf')
  if(save=='tif'){
    tiff(tifname,width=200, height=180, units='mm',res=300)
    print(p)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=10, height=8)
    print(p)
    dev.off()
  }
  print(p)
}

##' @title Target prediction
##' @description predict miRNA target genes
##' @param miRNA miRNAs
##' @param database miRNA target genes prediction data base,
##'     default is 'starBase',or 'mirTarBase' 'miRcode'
##' @param datatype prefix for Save file names
##' @export
##' @author zexl
WriteTargets <- function(miRNA,database='starBase',datatype='Gene'){
  tar <- sapply(miRNA,function(mir){mirTargets[[database]][[mir]]})
  targetout<-matrix(,,2)
  for (i in miRNA) {
    for (j in tar[[i]]) {
      targetout<-rbind(targetout,c(i,j))
    }
  }
  targetout<-targetout[-1,]
  id<-biotype[match(targetout[,2],biotype$ensemblID),3]
  targetout<-data.frame(targetout,geneSymbol=id)
  colnames(targetout)<-c('miRNA','GeneID','geneSymbol')
  write.csv(targetout,file=paste('Targets_',datatype,'.csv',sep=''),
            quote=F,row.names=F)
}

##' @title Target prediction
##' @description Combine ggplot2-based plots into a single plot,copy from Seurat
##' @param plots A list of gg objects
##' @param ncol Number of columns
##' @param legend Combine legends into a single legend choose
##'               from 'right' or 'bottom'; pass 'none' to remove 
##'               legends, or NULL to leave legends as they are
##' @param ...  Extra parameters passed to plot_grid
##' @importFrom cowplot plot_grid
##' @importFrom cowplot get_legend
##' @importFrom cowplot theme_cowplot
##' @export
##' @author zexl
CombinePlot<-function (plots, ncol = NULL, legend = NULL, ...) {
  plots.combined <- if (length(x = plots) > 1) {
    if (!is.null(x = legend)) {
      if (legend != "none") {
        plot.legend <- get_legend(plot = plots[[1]] + 
                                    theme(legend.position = legend))
      }
      plots <- lapply(X = plots, FUN = function(x) {
        return(x + NoLegend())
      })
    }
    plots.combined <- plot_grid(plotlist = plots, ncol = ncol, 
                                align = "hv", ...)
    if (!is.null(x = legend)) {
      plots.combined <- switch(EXPR = legend, bottom = plot_grid(plots.combined, 
                                                                 plot.legend, ncol = 1, rel_heights = c(1, 0.2)), 
                               right = plot_grid(plots.combined, plot.legend, 
                                                 rel_widths = c(3, 0.3)), plots.combined)
    }
    plots.combined
  }
  else {
    plots[[1]]
  }
  return(plots.combined)
}
