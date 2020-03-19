\name{PAS2GEF}
\alias{PAS2GEF}
\title{PAS2GEF, build reference regions for 3'UTR PASs}
\usage{
    PAS2GEF(GTFfile)
}
\arguments{
    \item{GTFfile}{GTF file of gene annotation}
}
\description{
    Build 3'UTR PAS and IPA (IPA and LE) Reference using GTF file.
}
\value{
    The function \code{PAS2GEF()} returns 3 input tables of PAS references:
    PASREF$refUTRraw is for 3'UTR PAS, 
	PASREF$dfIPA and PASREF$dfLE are for IPA references.

}
\author{Ruijia Wang}
\examples{
## build Reference ranges for 3'UTR PASs in mouse
	download.file(url='ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz',
              destfile='Mus_musculus.GRCm38.99.gtf.gz')			  
	GTFfile="Mus_musculus.GRCm38.99.gtf.gz"	
    PASREF=PAS2GEF(GTFfile)	
	refUTRraw=PASREF$refUTRraw
    dfIPA=PASREF$dfIPA
	dfLE=PASREF$dfLE
}