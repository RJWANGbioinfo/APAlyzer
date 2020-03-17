\name{REF4PAS}
\alias{REF4PAS}
\title{REF4PAS, build reference regions for 3'UTR and Intronic PAS using 
		dataframe formated input}
\usage{
    REF4PAS(refUTRraw, dfIPAraw, dfLEraw)
}
\arguments{
    \item{refUTRraw}{a dataframe containing 6 colmuns for 3'UTR PASs:
    'gene_symbol',
    'Chrom',
    'Strand',
    'Proximal',
    'Distal',
    'cdsend'}
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
	}
\description{
    build reference regions for 3'UTR and Intronic PAS using dataframe formated input
}
\value{
    The function \code{REF4PAS()} returns list a genomic ranges of 3'UTR, 
	Intronic PAS and last 3'exon regions for each gene
}
\author{Ruijia Wang}
\examples{
## build Reference ranges for 3'UTR and Intronic PAS in mouse (mm9)
	library(repmis)
	URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
	file="mm9_REF.RData"
	source_data(paste0(URL,file,"?raw=True"))
    refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
	dfIPAraw=dfIPA[which(dfIPA$Chrom=='chr19'),]
	dfLEraw=dfLE[which(dfLE$Chrom=='chr19'),]	
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
	UTRdbraw=PASREF$refUTRraw
    dfIPA=PASREF$dfIPA
	dfLE=PASREF$dfLE	
}