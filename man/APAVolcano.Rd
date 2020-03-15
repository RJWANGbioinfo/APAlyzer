\name{APAVolcano}
\alias{APAVolcano}
\title{APAVolcano, APA Volcano plotting}
\usage{
	APAVolcano (df, Pcol = "pvalue",PAS='3UTR',
						top = 5, markergenes = NULL,
						y_cutoff = 0.05,xlab = "RED", ylab = "-Log10(P-value)",
						PAScolor = c("gray80", "red", "blue"),
						alpha = 0.75, plot_title = NULL,
						width = 4, height = 2.5)	
}
\arguments{
    \item{df}{a dataframe of APAdiff output}
    \item{Pcol}{p-value column used to for y-axis of volcano plot, default is 'pvalue'}
    \item{top}{number of genes/IPA to label in the plot, default is 5}
    \item{markergenes}{a set of genes to label in the plot}
    \item{PAS}{type of PAS analyzed, either '3UTR' or 'IPA', default is '3UTR'}
    \item{y_cutoff}{y cutoff line, default is 0.05}
	\item{xlab}{lable of x-axis, default is 'RED'}
	\item{ylab}{lable of y-axis, default is '-Log10(P-value)'}
	\item{PAScolor}{dot color for 'NC','UP' and 'DN' gene/IPAs, 
	default is "gray80", "red", and "blue"}
	\item{alpha}{alpha of the dot, default is 0.75}
	\item{plot_title}{Main title of plot}
	\item{width}{width of the dot, default is 4}
	\item{height}{height of the dot, default is 2.5}
    }
\description{
	APA Volcano plotting
}
\value{
    The function \code{APAVolcano} return a Volcano plot.
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
	UTR_APA_PLOT=APAVolcano(test_3UTRmuti, PAS='3UTR', Pcol = "pvalue", top=5, plot_title='3UTR APA')

    ## IPA plot
    test_IPAmuti=APAdiff(sampleTable1,IPA_OUTraw,
    conKET='NT',trtKEY='KD',PAS='IPA',CUTreads=0)
	IPA_PLOT=APAVolcano(test_IPAmuti, PAS='IPA', Pcol = "pvalue", top=5, plot_title='IPA')


}