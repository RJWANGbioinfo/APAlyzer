\name{APABox}
\alias{APABox}
\title{APABox, APA RED Box plotting}
\usage{
	APABox (df, xlab = "APAreg", ylab = "RED",
						plot_title = NULL)	
}
\arguments{
    \item{df}{a dataframe of APAdiff output}
	\item{xlab}{lable of x-axis, default is 'APAreg'}
	\item{ylab}{lable of y-axis, default is 'RED'}
	\item{plot_title}{Main title of plot}
    }
\description{
	APA RED Box plotting
}
\value{
    The function \code{APABox} return a Box plot.
}
\author{Ruijia Wang}
\examples{
    library("TBX20BamSubset")
    library("Rsamtools")
    flsall = getBamFileList()
    extpath = system.file("extdata",
    "mm9_TBX20.APAout.RData", package="APAlyzer")
    load(extpath)
    sampleTable1 = data.frame(samplename = c(names(flsall)),
        condition = c(rep("NT",3),rep("KD",3)))
    sampleTable2 = data.frame(samplename = c("SRR316184","SRR316187"),
        condition = c("NT","KD"))
    ## 3'UTR APA plot
    test_3UTRmuti=APAdiff(sampleTable1,DFUTRraw,
    conKET='NT',trtKEY='KD',PAS='3UTR',CUTreads=0)
	UTR_APA_PLOTBOX=APABox(test_3UTRmuti, plot_title='3UTR APA')

    ## IPA plot
    test_IPAmuti=APAdiff(sampleTable1,IPA_OUTraw,
    conKET='NT',trtKEY='KD',PAS='IPA',CUTreads=0)
	IPA_PLOTBOX=APABox(test_IPAmuti, plot_title='IPA')


}