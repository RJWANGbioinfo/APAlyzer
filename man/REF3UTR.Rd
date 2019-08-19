\name{REF3UTR}
\alias{REF3UTR}
\title{REF3UTR, build reference regions for 3'UTR PASs}
\usage{
    REF3UTR(refUTR)
}
\arguments{
    \item{refUTR}{a dataframe containing 6 colmuns for 3'UTR PASs:
    'gene_symbol',
    'Chrom',
    'Strand',
    'Proximal',
    'Distal',
    'cdsend'}
}
\description{
    Build 3'UTR PAS Reference for distal and proximal PAS.
}
\value{
    The function \code{REF3UTR()} returns a genomic ranges of aUTR(pPAS to dPAS)
    and cUTR(cdsend to pPAS) regions for each gene
}
\author{Ruijia Wang}
\examples{
## build Reference ranges for 3'UTR PASs in human
    extpath = system.file("extdata", "mm9_REF.RData", package="APAlyzer")
    load(extpath)
    refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
    UTRdbraw=REF3UTR(refUTRraw)
}
