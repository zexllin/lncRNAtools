##' @importFrom clusterProfiler enrichKEGG
##' @import org.Hs.eg.db
##' @importFrom DOSE enrichDO
##' @importFrom biomaRt useMart
##' @importFrom biomaRt getBM
##' @import edgeR
##' @import DESeq2
##' @importFrom SummarizedExperiment assay
##' @importFrom BiocParallel MulticoreParam
##' @importFrom BiocParallel register
##' @references
##'     Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package
##'     for differential expression analysis of digital gene expression data.
##'     Bioinformatics. 2010 Jan 1;26(1):139-40. \cr
##'     Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK.
##'     limma powers differential expression analyses for RNA-sequencing and
##'     microarray studies. Nucleic acids research. 2015 Jan 20;
##'     43(7):e47-e47. \cr
##'     Love MI, Huber W, Anders S. Moderated estimation of fold change and
##'     dispersion for RNA-seq data with DESeq2. Genome biology. 2014 Dec 5;
##'     15(12):550.
##'     @author zexl
##'     @examples
##'     DeResult<-DEAnalysis(Counts,
##'             group,control='Z')
##'     DEGAll<-DeResult[[1]]
##'     normData<-DeResult[[2]]
##'     normData_log <- assay(
##'        rlogTransformation(DeResult[[3]])
##'        )
##' @title Routine analysis of transcriptome sequencing
##' @description Difference analysis of transcriptome sequencing
##' @param countDataRaw Expression matrix of transcriptome sequencing
##' @param Group a vector of sample group
##' @param control Set up control group
##' @param filter a logical ,filter gene use cpm 
##' @param lognormdata a logical,rlogTransformation for normdata
##' @return a list,data frame of Difference analysis,normlized Expression matrix ,
##'  and normlized Expression matrix after log Transformation
##' @export
##' @author zexl

DEAnalysis <- function(countDataRaw, Group,control,filter=TRUE,lognormdata = TRUE){
  
  if (is.logical(lognormdata) && (length(lognormdata) != 1L || is.na(lognormdata)))
    stop("'lognormdata' must be 'TRUE', 'FALSE'")
  lognormdata <- if (is.logical(lognormdata))
    lognormdata
  else lognormdata
  
  if (is.logical(filter) && (length(filter) != 1L || is.na(filter)))
    stop("'lognormdata' must be 'TRUE', 'FALSE'")
  filter <- if (is.logical(filter))
    filter
  else filter
  
  
  dge = DGEList(counts = countDataRaw)
  keep <- rowSums(cpm(dge) > 1) >= 0.5*length(Group)
  cat (paste('Total Number of genes: ', nrow(countDataRaw), '\n', sep=''))
  cat ('filter condition: rowSums(cpm(dge) > 1) >= 0.5*length(Group)\n')
  cat (paste('Number of genes for downstream analysis: ', length(which(keep)), '\n', sep=''))
  mycountData <- countDataRaw
  condition = factor(Group,levels = c(unique(Group)))
  coldata <- data.frame(row.names=colnames(mycountData), condition)
  coldata$condition <- relevel(coldata$condition, ref = control)
  dds <- DESeqDataSetFromMatrix(countData = mycountData,
                                colData = coldata, design=~condition)
  #
  if (filter==TRUE) {
    dds <- dds[keep, ]
  }
  register(MulticoreParam(10))
  dds2 <- DESeq(dds, parallel=TRUE)
  exprSetdata <- counts(dds2, normalized=TRUE)
  res <- results(dds2)
  DEGAll <- data.frame(res)
  colnames(DEGAll) <- c('baseMean', 'logFC', 'lfcSE', 'stat',
                        'PValue', 'FDR')
  o <- order(DEGAll$FDR)
  DEGAll <- DEGAll[o,]
  if (length(grep('^ENSG',rownames(countDataRaw)))>1) { 
    degList <- biotype[match(rownames(DEGAll), biotype$ensemblID),]
    degOutput <- data.frame(symbol=degList$geneSymbol, biotype=degList$biotype,
                            group=degList$group, DEGAll)

  }
  else{
    degOutput<-DEGAll
  }
  if(lognormdata==TRUE){
    normData_log<-assay(rlogTransformation(dds2))
    return(list(degOutput,exprSetdata,normData_log))
  }
  else{
    return(list(degOutput,exprSetdata,dds2))
  }
}


##' @title GO enrichment analysis
##' @description Performs Gene Ontology (GO) enrichment
##'     analyses by clusterProfiler packages
##' @param gene a vector of Ensembl gene id
##' @param simplify logical, specifying whether to remove redundant GO terms.
##'     Default simplify=TRUE
##' @param level a numeric value, restrict the GO enrichment result at a
##'     specific GO level. Default is \code{0}, which means all terms
##'     should be returned
##' @param data.type prefix for Save file names
##' @references
##'     Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for
##'     comparing biological themes among gene clusters.
##'     Omics: a journal of integrative biology. 2012 May 1;16(5):284-7.
##' @references
##'     Yu G, Wang LG, Yan GR, He QY. DOSE: an R/Bioconductor package for
##'     disease ontology semantic and enrichment analysis. Bioinformatics.
##'     2014 Oct 17;31(4):608-9.
##' @export
##' @author zexl
##' @importFrom clusterProfiler enrichGO
##' @importFrom clusterProfiler gofilter
##' @importFrom clusterProfiler simplify
##' @examples
##'     Upgene<-deall[which(deall$logFC>=1),]
##'     UpgeneGoOut<-GOEnrich(
##'             gene = rownames(Upgene),
##'             data.type = 'RNA_Up'
##'             )
##' @return return a data frame of GO enrichment analysis result, and write it to file
GOEnrich <- function(gene, simplify=TRUE, level=0, data.type = '') {

  message ('### This step may take a few minutes ###\n')
  if (length(gene)<1){
    stop("gene must not empty" )
  }
  goBP <- enrichGO(gene = gene,
                   universe = biotype$ensemblID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   keyType = 'ENSEMBL',
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.01,
                   readable = FALSE)

  if (level != 0) {
    goBP <- gofilter(goBP, level=level)
  }

  if (simplify==TRUE) {
    goBP <- simplify(goBP, cutoff=0.7, by="p.adjust", select_fun=min)
  }

  message ('BP analysis done!')

  goCC <- enrichGO(gene = gene,
                   universe = biotype$ensemblID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   keyType = 'ENSEMBL',
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.01,
                   readable = FALSE)

  if (level != 0) {
    goCC <- gofilter(goCC, level=level)
  }

  if (simplify==TRUE) {
    goCC <- simplify(goCC, cutoff=0.7, by="p.adjust", select_fun=min)
  }

  message ('CC analysis done!')

  goMF <- enrichGO(gene = gene,
                   universe = biotype$ensemblID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   keyType = 'ENSEMBL',
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.01,
                   readable = FALSE)

  if (level != 0) {
    goMF <- gofilter(goMF, level=level)
  }

  if (simplify==TRUE) {
    goMF <- simplify(goMF, cutoff=0.7, by="p.adjust", select_fun=min)
  }


  goBP <- organizeEnrichFun(data.frame(goBP@result))
  goCC <- organizeEnrichFun(data.frame(goCC@result))
  goMF <- organizeEnrichFun(data.frame(goMF@result))
  message ('MF analysis done!')

  Path<-paste0(data.type,'_enrichment/')
  if (!file.exists(Path)){
    dir.create(Path)
  }

  # write.table(goBP, file = paste0(Path,data.type,'_AllEnrichGOBP.txt'),sep="\t",quote=F,row.names=F)
  # write.table(goCC, file = paste0(Path,data.type,'_AllEnrichGOCC.txt'),sep="\t",quote=F,row.names=F)
  # write.table(goMF, file = paste0(Path,data.type,'_AllEnrichGOMF.txt'),sep="\t",quote=F,row.names=F)

  enrichOutput <- data.frame(rbind(goBP, goCC, goMF))
  enrichOutput$Category <- rep(c('GO_BP','GO_CC','GO_MF'),c(nrow(goBP),nrow(goCC),nrow(goMF)))
  write.csv(enrichOutput,file = paste0(Path,data.type,'_AllenrichGO.csv'),row.names = F)

  return (enrichOutput)
}

##' @title KEGG Enrichment
##' @description Performs GeneKyoto Encyclopedia of Genes and Genomes KEGG
##'     pathway enrichment analyses by clusterProfiler packages
##' @param gene a vector of Ensembl gene id
##' @param data.type prefix for Save file names
##' @references
##'     Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for
##'     comparing biological themes among gene clusters.
##'     Omics: a journal of integrative biology. 2012 May 1;16(5):284-7. \cr
##'     Yu G, Wang LG, Yan GR, He QY. DOSE: an R/Bioconductor package for
##'     disease ontology semantic and enrichment analysis. Bioinformatics.
##'     2014 Oct 17;31(4):608-9.
##' @export
##' @author zexl
##' @examples
##'     Upgene<-deall[which(deall$logFC>=1),]
##'     UpgeneKEGGOut<-KEGGEnrich(
##'            gene = rownames(Upgene),
##'            data.type = 'RNA_Up'
##'            )
##' @return return a data frame of KEGG enrichment analysis result, and write it to file
KEGGEnrich <- function(gene,data.type = '') {

  message ('### This step may take a few minutes ###\n')
  genes <- biotype[match(gene, biotype$ensemblID),]
  genes <- genes[! is.na(genes$entrezgene),]
  universe <- biotype[!is.na(biotype$entrezgene),]

  kegg <- enrichKEGG(gene = as.character(genes$entrezgene),
                     organism = 'hsa',
                     universe = as.character(unique(
                       universe$entrezgene[!is.na(universe$entrezgene)])),
                     minGSSize = 10,
                     maxGSSize = 500,
                     pAdjustMethod = 'fdr',
                     pvalueCutoff = 0.01)

  kegg <- data.frame(kegg@result)

  kegg$geneID <- unlist(lapply(kegg$geneID, function(v)
    paste(genes$ensemblID[match(strsplit(v, '/', fixed=TRUE)[[1]],
                                genes$entrezgene)], collapse = '/')))
  kegg <- organizeEnrichFun(kegg)
  message ('KEGG analysis done!')

  Path<-paste0(data.type,'_enrichment/')
  if (!file.exists(Path)){
    dir.create(Path)
  }
  write.csv(kegg, file = paste0(Path,data.type,'_AllEnrichKEGG.csv'),row.names=F)
  return(kegg)
}



##' @title DO Enrichment
##' @description Performs Disease Ontology enrichment analyses by
##'      clusterProfiler packages
##' @param gene a vector of Ensembl gene id
##' @param data.type prefix for Save file names
##' @references
##'     Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for
##'     comparing biological themes among gene clusters.
##'     Omics: a journal of integrative biology. 2012 May 1;16(5):284-7. \cr
##'     Yu G, Wang LG, Yan GR, He QY. DOSE: an R/Bioconductor package for
##'     disease ontology semantic and enrichment analysis. Bioinformatics.
##'     2014 Oct 17;31(4):608-9.
##' @export
##' @author zexl
##' @examples
##'     Upgene<-deall[which(deall$logFC>=1),]
##'     UpgeneDoOut<-DOEnrich(
##'        gene = rownames(Upgene),
##'        data.type = 'RNA_Up'
##'        )
##' @return return a data frame of KEGG enrichment analysis result, and write it to file
##'
DOEnrich <- function(gene,data.type = '') {

  genes <- biotype[match(gene, biotype$ensemblID),]
  genes <- genes[! is.na(genes$entrezgene),]
  universe <- biotype[!is.na(biotype$entrezgene),]
  do <- enrichDO(gene = as.character(genes$entrezgene),
                 universe = as.character(unique(
                   universe$entrezgene[!is.na(universe$entrezgene)])),
                 ont = "DO",
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.01,
                 readable = FALSE)

  do <- data.frame(do@result)

  do$geneID <- unlist(lapply(do$geneID, function(v)
    paste(genes$ensemblID[match(strsplit(v, '/', fixed=TRUE)[[1]],
                                genes$entrezgene)], collapse = '/')))
  do   <- organizeEnrichFun(do)
  message ('DO analysis done!')

  Path<-paste0(data.type,'_enrichment/')
  if (!file.exists(Path)){
    dir.create(Path)
  }

  write.csv(do, file = paste0(Path,data.type,'_AllEnrichDO.csv'),row.names=F)
  return(do)
}
