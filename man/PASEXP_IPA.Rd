\name{PASEXP_IPA}
\alias{PASEXP_IPA}
\title{PASEXP_IPA, calculate relative expression of IPA regions}
\usage{
    PASEXP_IPA(dfIPAraw, dfLEraw, flS, Strandtype="NONE", nts=1, minMQS=0, SeqType = "SingleEnd")
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
    \item{minMQS}{minimum mapping quality score of counted reads, 
    parameter used by \link[Rsubread:featureCounts]{featureCounts}, 
    minMQS option, Default is 0}
    \item{SeqType}{set the sequencing type of reads in bam files can be either 
        'SingleEnd' (default) or 'ThreeMostPairEnd'.}
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
	library("repmis")
    flsall = getBamFileList()
	URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
	file="mm9_REF.RData"
	source_data(paste0(URL,file,"?raw=True"))
    IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=1)
}
