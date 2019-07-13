\name{REFCDS}
\alias{REFCDS}
\title{REFCDS, build reference regions for CDS of protein coding genes}
\usage{
    REFCDS(txdb,IDDB)
}
\arguments{
    \item{txdb}{a TranscriptDb generate using GenomicFeatures}
    \item{IDDB}{Genome annotation of the corresponding species,
    e.g., "org.Hs.eg.db"}
}
\description{
    Build CDS reference for protein coding genes.
}
\value{
    The function \code{REFCDS()}
    returns a genomic ranges of CDS regions for each coding gene
}
\author{Ruijia Wang}
\examples{
## build Reference ranges for CDS in human coding genes
    library("GenomicFeatures")
    library("org.Mm.eg.db")
    txdb <- makeTxDbFromUCSC(genome="mm9", tablename="refGene")
    seqlevels(txdb) <-"chr19"
    IDDB <- org.Mm.eg.db
    CDSdbraw=REFCDS(txdb,IDDB)
}
