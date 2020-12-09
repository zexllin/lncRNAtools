##' @title Plots for enrichment analysis
##' @description Bar plot and bubble plot for GO enrichment analysis
##' @param enrichment a dataframe generated from
##'     \code{\link{GOEnrich}}
##' @param type type of the plot, should be one of \code{'bar'}
##'     and \code{'bubble'}
##' @param num.terms number of terms to be plotted. Default is \code{10}
##' @param axis y axis,GO number \code{term},GO Description \code{description},
##' or GO number~Description \code{all},Default is \code{all}
##' @param data.type prefix for Save file names
##' @param w image weight
##' @param h image height
##' @param save Output file type,default is not save
##' @export
##' @import ggplot2
##' @import stringr

GOplots<-function(enrichment, type='bar',num.terms=10,axis='all',
                  data.type='',w=0,h=0,save='') {



  goBP <- enrichment[enrichment$Category=='GO_BP',]
  goCC <- enrichment[enrichment$Category=='GO_CC',]
  goMF <- enrichment[enrichment$Category=='GO_MF',]

  if (nrow(goBP) > num.terms) {
    goBP <- goBP[seq_len(num.terms),]
  }

  if (nrow(goCC) > num.terms) {
    goCC <- goCC[seq_len(num.terms),]
  }

  if (nrow(goMF) > num.terms) {
    goMF <- goMF[seq_len(num.terms),]
  }
  go <- rbind(goBP, goCC, goMF)
  
  #save plot section
  outpath <- paste0(data.type,'_enrichment/')
  if(!file.exists(outpath))
    dir.create(outpath)
  if(save!=''){
    outw<-go[,c('Counts','pValue','FDR','Category')]
    outw$ID<-substring(go$Terms,1,10)
    outw$Description<-substring(go$Terms,12)
    outw<-outw[,c('ID','Description','Counts','pValue','FDR','Category')]
    write.csv(outw,file = paste0(outpath,data.type,'_EnrichGOplot.csv'),
              row.names = F)
  }
  
  #plot y axis
  if(axis=='description'){
    go$Terms<-substring(go$Terms,12)
  }
  if(axis=='term'){
    go$Terms<-substring(go$Terms,1,10)
  }

  sz <- 16
  if (type=='bubble') {

    go$Terms <- paste(go$Terms, '[', go$Category, ']', sep='')
    goBasic = ggplot(data=go, mapping=aes(x=Terms,
                      y=foldEnrichment,color=FDR,size=Counts))
    GOP<-goBasic+geom_point()+ coord_flip() +
        scale_x_discrete(limits=rev(go$Terms)) +
        scale_colour_gradientn(colors= c("red","green")) +
        xlab('')+ylab('Fold enrichment') +
        guides(shape = guide_legend(order=1),
               colour = guide_colourbar(order=2)) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='black'),
                         panel.background = element_blank()) +
        ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
        theme(axis.text.y=element_text(size=sz,colour = 'black'),
              axis.title=element_text(size=sz),
              axis.text.x=element_text(size=sz,colour = 'black')) +
        theme(legend.text = element_text(size = sz),
              legend.title = element_text(size = sz)) +
        theme(strip.text = element_text(size = sz),
              legend.key.size = unit(0.8,'cm'))

  }

  else if (type=='bar') {

    goBasic = ggplot(data=go, mapping=aes(x=Terms, y=-log(FDR,10),
                                              fill=Category))
    GOP<-goBasic + geom_bar(stat='identity') +
        scale_x_discrete(limits=rev(go$Terms)) +
        ylim(0, max(-log(go$FDR,10))) +
        theme(legend.title=element_blank())+
        ylab('-log10(FDR)')+xlab('') + coord_flip() +
        scale_fill_hue(name='',breaks=go$Category,
                       labels=go$Category) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank()) +
        theme(axis.text=element_text(size=sz,colour = 'black'),
              axis.title=element_text(size=sz)) +
        theme(legend.text = element_text(size=sz))
  }
  
  tifname<-paste(outpath,data.type,"_EnrichGO_",type,".tif",sep="")
  pdfname = gsub('tif','pdf',tifname)
  if(save=='tif'){
    tiff(tifname,width=360+w, height=240+h, units='mm',res=300)
    print(GOP)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=12+w, height=8+h)
    print(GOP)
    dev.off()
  }
  print(GOP)
}

##' @title Plots for enrichment analysis
##' @description Bar plot and bubble plot for KEGG/DO enrichment analysis
##' @param kegg a dataframe generated from
##'     \code{\link{KEGGEnrich}} or \code{\link{DOEnrich}}
##' @param type type of the plot, should be one of \code{'bar'}
##'     and \code{'bubble'}
##' @param num.terms number of terms to be plotted. Default is \code{10}
##' @param data.type prefix for Save file names
##' @param w image weight
##' @param h image height
##' @param save Output file type,default is not save
##' @export

KEGGplots<-function(kegg, type='bar',num.terms=10,w=0,h=0,
                    data.type='',save='') {

  outpath <- paste0(data.type,'_enrichment/')
  if(!file.exists(outpath)) dir.create(outpath)

  if (unique(substr(kegg$Terms,1,4))=="DOID"){
    tifname<-paste(outpath,data.type,"_EnrichDO_",type,".tif",sep="")
    outtable <-paste0(outpath,data.type,"_EnrichDOplot.csv")
  }
  else {
    tifname<-paste(outpath,data.type,"_EnrichKEGG_",type,".tif",sep="")
    outtable <-paste0(outpath,data.type,"_EnrichKEGGplot.csv")
  }

  if (nrow(kegg) > num.terms) {
    kegg <- kegg[seq_len(num.terms),]
  }

  sz = 20 #plot size
  if (type=='bubble') {

    keggBasic = ggplot(data=kegg, mapping=aes(x=Terms,
                                              y=foldEnrichment,color=FDR,size=Counts))
    KEGGP<-keggBasic+geom_point()+ coord_flip() +
        scale_x_discrete(limits=rev(kegg$Terms)) +
        scale_colour_gradientn(colors= c("red","green")) +
        xlab('')+ylab('Fold enrichment') +
        guides(shape = guide_legend(order=1),
               colour = guide_colourbar(order=2)) +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='black'),
                         panel.background = element_blank()) +
        ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=20)) +
        theme(axis.text.y=element_text(size=sz),
              axis.title=element_text(size=sz),
              axis.text.x=element_text(size=sz)) +
        theme(legend.text = element_text(size = sz),
              legend.title = element_text(size = sz)) +
        theme(strip.text = element_text(size = sz),
              legend.key.size = unit(0.8,'cm'))

  }
  else if (type=='bar') {

    keggBasic = ggplot(data=kegg, mapping=aes(x=Terms, y=-log(FDR,10)))
    KEGGP<-keggBasic + geom_bar(stat='identity') +
        scale_x_discrete(limits=rev(kegg$Terms)) +
        ylim(0, max(-log(kegg$FDR,10))) +
        theme(legend.title=element_blank())+ylab('-log10(FDR)')+
        xlab('') + coord_flip() +
        theme_bw()+theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour='white'),
                         panel.background = element_blank()) +
        theme(axis.text.y=element_text(size=sz),
              axis.title=element_text(size=sz),
              axis.text.x=element_text(size=sz)) +
        theme(legend.position = 'none')
  }
  if (save=='tif') {
    tiff(file=tifname,width=300+w, height=200+h, units='mm',res=400)
    print(KEGGP)
    dev.off()
  }
  if(save=='pdf'){
    pdf(gsub('tif','pdf',tifname),width=12+w, height=8+h)
    print(KEGGP)
    dev.off()
  }
  if (save!=''){ 
    outw<-kegg[,c('Counts','pValue','FDR')]
    outw$ID<-substring(kegg$Terms,1,8)
    outw$Description<-substring(kegg$Terms,10)
    outw<-outw[,c('ID','Description','Counts','pValue','FDR')]
    write.csv(outw,file = outtable,row.names = F)
  }
  print(KEGGP)
}
