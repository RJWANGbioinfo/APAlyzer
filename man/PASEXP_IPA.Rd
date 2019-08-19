\name{PASEXP_IPA}
\alias{PASEXP_IPA}
\title{PASEXP_IPA, calculate relative expression of IPA regions}
\usage{
    PASEXP_IPA(dfIPAraw, dfLEraw, flS, Strandtype="NONE", nts=1)
}
\arguments{
    \item{dfIPAraw}{a dataframe containing 8 colmuns for Intronic PASs:
    'PASid',
    'gene_symbol',
    'Chrom',
    'Strand',
    'Pos',
    'upstreamSS',
    'downstreamSS'.
    'upstreamSS' means closest 5'/3' splice site to IPA,
    'downstreamSS' means closest 3' splice site}
    \item{dfLEraw}{a dataframe containing 5 colmuns for 3' least exon:
    'gene_symbol',
    'Chrom',
    'Strand',
    'LEstart',
    'TES'.
    'LEstart' means the start position of last 3' exon.}
    \item{flS}{bamfile lists containing the file and path of bam files}
    \item{Strandtype}{strand type of the bam file;
    "forward" is forwand sequencing,
    "invert" is reverse sequencing and
    "NONE" is non-strand specific, Default is "NONE".}
    \item{nts}{number of threads used for computing, 
    parameter used by \link[Rsubread:featureCounts]{featureCounts}, 
    nthread option, Default is 1}
    }
\description{
    Map reads to IPA regions and
    calculte relative expression of aUTR and cUTR regions.
}
\value{
    The function \code{PASEXP_IPA()} return a dataframe
    containning reads count, RPKM and
    relative expression of aUTR and cUTR for each gene
}
\author{Ruijia Wang}
\examples{
## count reads mapped to IPA regions and
## calculte relative expression of aUTR and cUTR regions
## using forward sequencing
    library("TBX20BamSubset")
    library("Rsamtools")
    library("GenomicAlignments")
    flsall = getBamFileList()
    extpath = system.file("extdata", "mm9_REF.RData", package="APAlyzer")
    load(extpath)
    IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=1)
}
