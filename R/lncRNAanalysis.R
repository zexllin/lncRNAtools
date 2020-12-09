##' @title Analysis of lncRNA co expression with mRNA
##' @description compute correlation between lncRNA with mRNA
##' @param lncname a vector of Ensembl lncRNA id
##' @param pcname a vector of Ensembl protein coding gene id
##' @param normData normlized Expression matrix after log Transformation,by \code{\link{DEAnalysis}}
##' @return a datafram of correlation analysis
##' @export
Coexpression <- function(lncname,pcname,normData){
  i=0
  data<-t(normData[c(lncname,pcname),])
  CoOut<-list()
  for (l in lncname) {
    for (p in pcname) {
      i = i + 1
      cor <- cor.test(data[,l],data[,p])
      CoOut[[i]]<-c(l,p,cor$estimate,cor$p.value)
    }
  }
  CoOut<-do.call(rbind,CoOut)
  colnames(CoOut)<-c('lnc','gene','cor','PValue')
  CoOut <- as.data.frame(as.matrix(CoOut), stringsAsFactors=FALSE)
  CoOut$cor<-as.numeric(CoOut$cor)
  CoOut$PValue<-as.numeric(CoOut$PValue)
  return(CoOut)
}

##' @title lncRNA cis
##' @description Calculate the distance between lncrna and its adjacent coding genes,
##'  and combine the correlation results
##' @param correlationresult a datafram of correlation analysis
##' @return a datafram of lncRNA cis analysis results
##' @export
mergeCorForNeargene<-function(correlationresult){
  NearCodeGene_Lnc$cor=rep(NA,nrow(NearCodeGene_Lnc))
  NearCodeGene_Lnc$PValue=rep(NA,nrow(NearCodeGene_Lnc))
  NearCodeGene_Lnc$distance=rep(NA,nrow(NearCodeGene_Lnc))
  for (i in 1:nrow(NearCodeGene_Lnc)) {
    if ((NearCodeGene_Lnc[i,'start_NearGene'] <= NearCodeGene_Lnc[i,'start_Lnc'])&
        (NearCodeGene_Lnc[i,'end_NearGene'] >= NearCodeGene_Lnc[i,'end_Lnc'])) {
      NearCodeGene_Lnc[i,'distance']=0
    }
    else if ((NearCodeGene_Lnc[i,'start_NearGene'] >= NearCodeGene_Lnc[i,'start_Lnc'])&
             (NearCodeGene_Lnc[i,'end_NearGene'] <= NearCodeGene_Lnc[i,'end_Lnc'])){
      NearCodeGene_Lnc[i,'distance']=0
    }
    else{
      absdistance=min(abs(NearCodeGene_Lnc[i,'start_NearGene']-NearCodeGene_Lnc[i,'start_Lnc']),
                      abs(NearCodeGene_Lnc[i,'end_NearGene']-NearCodeGene_Lnc[i,'end_Lnc']))
      if (NearCodeGene_Lnc[i,'start_NearGene']-NearCodeGene_Lnc[i,'start_Lnc']>0){
        NearCodeGene_Lnc[i,'distance']=absdistance
      }
      else{
        NearCodeGene_Lnc[i,'distance']=-absdistance
      }
    }
    l<-as.character(NearCodeGene_Lnc[i,'geneid_Lnc'])
    p<-as.character(NearCodeGene_Lnc[i,'geneid_NearGene'])
    rownum<-which(correlationresult$lnc==l & correlationresult$gene==p)
    if (length(rownum)==0) next
    NearCodeGene_Lnc[i,'cor']=correlationresult[rownum,'cor']
    NearCodeGene_Lnc[i,'PValue']=correlationresult[rownum,'PValue']
  }
  NearCodeGene_Lnc_pearson<-NearCodeGene_Lnc[which(!is.na(NearCodeGene_Lnc$cor)),]
  #rownames(NearCodeGene_Lnc_pearson)<-NearCodeGene_Lnc_pearson$geneid_Lnc
  return(NearCodeGene_Lnc_pearson)
}


##' @title co expression lncRNA enrichment analysis
##' @description Based on the results of co expression, GO and KEGG enrichment analysis of
##'  single lncrna was carried out, and bar and bubble charts were drawn
##' @param lncname lncRNA Ensembl ID
##' @param Coexpressgene a datafram of correlation analysis by \code{\link{Coexpression}}
##' @return Generate folder named by lncrna ensemble ID, and output result file, including
##'  enrichment analysis table and visualization results
##' @export
##' @author zexl
LncRNAEnrich<-function(lncname,Coexpressgene) {

  EnrichGene_Lnc <- Coexpressgene[which(Coexpressgene$lnc==lncname),'gene']

  ## GO
  goout <-GOEnrich(gene = EnrichGene_Lnc, simplify = F, gene.type = type)
  GOplots(goout,type='bar',num.terms=10,data.type=type,save = T)
  ## KEGG
  KEGGOut<-KEGGEnrich(gene = EnrichGene_Lnc,gene.type = lncname)
  KEGGplots(KEGGOut,type='bubble',num.terms=10,data.type=lncname,save = T)
}

##' @title Correlation analysis of transcription factors
##' @description Based on the results of co expression, Fisher's exact
##'  test was used to analyze the transcription factor association of single lncrna
##' @param test.frame Correlation analysis results of designated lncrna
##' @param prex prefix for Save file names
##' @export
##' @author zexl
##' @import TFEA.ChIP
TFanalysis.Fisher<- function(test.frame, prex){
  #Use built-in functions to transform ID
  test.frame$entrezgene<-biotype[match(test.frame$gene,biotype$ensemblID),'entrezgene']
  DEAlls_table <- test.frame[!is.na(test.frame$entrezgene), ]
  DEAlls_table <- DEAlls_table[order(DEAlls_table$cor,decreasing = TRUE), ]
  #DEAlls_table <- preprocessInputData( test.frame )
  Genes.diff <- as.character(DEAlls_table[which(abs(DEAlls_table$cor)>0.8
                                                & DEAlls_table$PValue <= 0.05),'entrezgene'])
  ## Take the genes with high correlation as the gene set, and take all other genes contained in
  ## the chip SEQ data set as the gene set
  ## Take genes as background, and count 2x2 table of each chip SEQ data set
  CM_list_diff <- contingency_matrix(Genes.diff)
  ## Calculate OR values and P values
  pval_mat_diff <- getCMstats( CM_list_diff )
  #head( pval_mat_diff )
  #colnames(pval_mat_diff)[9]<-'-log10.adj.pVal'
  pval_mat_diff$log10.adj.pVal<- -pval_mat_diff$log10.adj.pVal
  enrFisher <- pval_mat_diff
  ## Increase gene number statistics and gene names
  enrFisher$GenesLength <- numeric(nrow(enrFisher))
  enrFisher$Genes <- character(nrow(enrFisher))
  if (!exists("Mat01")) {
    Mat01 <- NULL
    data("Mat01", package = "TFEA.ChIP", envir = environment())
  }
  for (i in 1:nrow(enrFisher)) {
    met01 <- Mat01[,enrFisher[i,'Accession']]
    targetGenes <- names(met01[which(met01 == 1)])
    enrGenes <- intersect(Genes.diff,targetGenes)
    enrGenes <- biotype[match(enrGenes,biotype$entrezgene),'geneSymbol']
    enrFisher$Genes[i] <- paste(enrGenes,collapse =',')
    enrFisher$GenesLength[i] <- length(enrGenes)
  }
  ## write 
  write.csv(enrFisher,file = paste0(prex,'_TFsFisherTest_enrichment.csv'),row.names = F)
  enrFisher_sig <- enrFisher[which(enrFisher$OR >= 1.5 & enrFisher$adj.p.value <= 0.05),]
  if (length(enrFisher_sig$TF)==0){
    write(c('No significance results'),file=paste0(prex,'readme.txt'))
    return(message('specialTFs 0'))
  }
  write.csv(enrFisher_sig,file = paste0(prex,'_TFsFisherTest_enrichment_Sig.csv'),row.names = F)
  ## The transcription factors with significant differences are shown in the figure
  specialTFs <- as.character(unique(enrFisher[which(enrFisher$OR >= 1.5 & enrFisher$adj.p.value <= 0.05),'TF']))
  names(specialTFs) = as.character(specialTFs)
  p <- plot_CM( pval_mat_diff, specialTF = specialTFs) #plot p-values against ORs highlighting indicated TFs
  htmlwidgets::saveWidget(p, paste0(prex,'_TFsFisherTest_enrichment.html'))
  return(enrFisher)
}
