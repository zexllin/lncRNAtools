% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lncRNAanalysis.R
\name{LncRNAEnrich}
\alias{LncRNAEnrich}
\title{co expression lncRNA enrichment analysis}
\usage{
LncRNAEnrich(lncname, Coexpressgene)
}
\arguments{
\item{lncname}{lncRNA Ensembl ID}

\item{Coexpressgene}{a datafram of correlation analysis by \code{\link{Coexpression}}}
}
\value{
Generate folder named by lncrna ensemble ID, and output result file, including
 enrichment analysis table and visualization results
}
\description{
Based on the results of co expression, GO and KEGG enrichment analysis of
 single lncrna was carried out, and bar and bubble charts were drawn
}
\author{
zexl
}
