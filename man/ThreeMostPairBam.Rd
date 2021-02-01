\name{ThreeMostPairBam}
\alias{ThreeMostPairBam}
\title{ThreeMostPairBam, extract 3 prime most alignment of a paired-end bam file}
\usage{
    ThreeMostPairBam(BamfilePath, OutDirPath, StrandType="NONE")
}
\arguments{
    \item{BamfilePath}{file path of a bam file}
    \item{OutDirPath}{output folder path}
    \item{StrandType}{strand type of the bam file;
    "forward-reverse": read 1 forward but read 2 is reverse sequencing,
    "reverse-forward": read 2 forward but read 1 is reverse sequencing, and
    "NONE" is non-strand specific, Default is "NONE".}
}
\description{
    extract 3 prime most alignment of a paired-end bam file and saved into a new bam file.
}
\value{
    The function \code{ThreeMostPairBam()} return a single-end bam file
    containning 3 prime most alignment of the input paired-end file
}
\author{Ruijia Wang}
\examples{
## Extract 3 prime most alignment of a paired-end 
## bam file and saved into a new bam file
    library("pasillaBamSubset")
	
    ThreeMostPairBam (BamfilePath=untreated3_chr4(), 
						OutDirPath=getwd(), 
						StrandType='forward-reverse')

}
