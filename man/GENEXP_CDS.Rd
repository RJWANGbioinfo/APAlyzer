\name{GENEXP_CDS}
\alias{GENEXP_CDS}
\title{GENEXP_CDS, count reads mapped to CDS regions and
calculate TPM for coding gene}
\usage{
    GENEXP_CDS(CDSbygene, flS, Strandtype="NONE")
}
\arguments{
    \item{CDSbygene}{a genomic ranges of CDS regions for each coding gene}
    \item{flS}{bamfile lists containing the file and path of bam files}
    \item{Strandtype}{strand type of the bam file;
    "forward" is forwand sequencing,
    "invert" is reverse sequencing and
    "NONE" is non-strand specific, Default is "NONE".}
}
\description{
    Map reads to CDS regions and calculate TPM for each gene.
}
\value{
    The function \code{GENEXP_CDS()}
    return a dataframe containing reads count, TPM for each gene
}
\author{Ruijia Wang}
\examples{
## count reads mapped to CDS regions and calculate TPM for each gene
## using forward sequencing
    library("TBX20BamSubset")
    library("Rsamtools")
    library("GenomicAlignments")
    library("GenomicFeatures")
    library("org.Mm.eg.db")
    flsall = getBamFileList()
    extpath = system.file("extdata", "mm9.chr19.refGene.R.DB", package="APAlyzer")
    txdb = loadDb(extpath, packageName='GenomicFeatures')
    IDDB = org.Mm.eg.db
    CDSdbraw = REFCDS(txdb,IDDB)
    DFGENEraw = GENEXP_CDS(CDSdbraw, flsall, Strandtype="forward")

}
