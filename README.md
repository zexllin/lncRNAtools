# README
## install lncRNAtools

lncRNAtools runs in the R statistical computing environment. You will need R version 3.6.1 or higher, Bioconductor version 3.10. 

install a few Bioconductor dependencies that aren't automatically installed

```
BiocManager::install(c('limma', 'edgeR','DESeq2','BiocParallel','TFEA.ChIP','scatterplot3d','clusterProfiler','DOSE'))
```
install a few Bioconductor dependencies that aren't automatically installed,check is True.

```
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
  
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if (!requireNamespace("lncRNAtools", quietly = TRUE)){
  
  if (!requireNamespace("edgeR", quietly = TRUE))
    BiocManager::install("edgeR")
  if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
  if (!requireNamespace("BiocParallel", quietly = TRUE))
    BiocManager::install("BiocParallel")
  if (!requireNamespace("TFEA.ChIP", quietly = TRUE))
    BiocManager::install("TFEA.ChIP")
  if (!requireNamespace("scatterplot3d", quietly = TRUE))
    BiocManager::install("scatterplot3d")
  if (!requireNamespace("DOSE", quietly = TRUE))
    BiocManager::install("DOSE")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")
  devtools::install_github('zexllin/lncRNAtools')
}
```

Now, install lncRNAtools through the zexllin GitHub, execute:
```
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github('zexllin/lncRNAtools')
```
Testing the installation
To ensure that lncRNAtools was installed correctly, start a new R session and run:
```
library(lncRNAtools)
```
