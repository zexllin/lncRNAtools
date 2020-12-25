# README
## install lncRNAtools

lncRNAtools runs in the R statistical computing environment. You will need R version 3.6.1 or higher, Bioconductor version 3.10. install a few Bioconductor dependencies that aren't automatically installed

```
BiocManager::install(c('limma', 'edgeR','DESeq2','BiocParallel','TFEA.ChIP','scatterplot3d','clusterProfiler','DOSE'))
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
