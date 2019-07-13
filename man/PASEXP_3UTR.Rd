\name{PASEXP_3UTR}
\alias{PASEXP_3UTR}
\title{PASEXP_3UTR, calculate relative expression of aUTR and cUTR regions}
\usage{
    PASEXP_3UTR(UTRdb, flS, Strandtype="NONE")
}
\arguments{
    \item{UTRdb}{a genomic ranges of aUTR(pPAS to dPAS)
    and cUTR(cdsend to pPAS) regions for each gene}
    \item{flS}{bamfile lists containing the file and path of bam files}
    \item{Strandtype}{strand type of the bam file;
    "forward" is forwand sequencing,
    "invert" is reverse sequencing and
    "NONE" is non-strand specific, Default is "NONE".}
}
\description{
    Map reads to 3'UTR APA regions and
    calculate relative expression of aUTR and cUTR regions.
}
\value{
    The function \code{PASEXP_3UTR()} return a dataframe
    containning reads count, RPKM and
    relative expression of aUTR and cUTR for each gene
}
\author{Ruijia Wang}
\examples{
## count reads mapped to 3'UTR APA regions and
## calculate relative expression of aUTR and cUTR regions
## using forward sequencing
    library("TBX20BamSubset")
    library("Rsamtools")
    library("GenomicAlignments")
    flsall <- getBamFileList()
    extpath <- system.file("extdata", "mm9_REF.RData", package="APAlyzer")
    load(extpath)
    refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
    UTRdbraw=REF3UTR(refUTRraw)
    DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")

}
