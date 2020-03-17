\name{APAdiff}
\alias{APAdiff}
\title{APAdiff, calculate delta relative expression (RED) and
statistics significance between two sample groups}
\usage{
    APAdiff(sampleTable,mutiraw, conKET='NT',
    trtKEY='KD',PAS='3UTR',CUTreads=0,p_methods="fdr")
}
\arguments{
    \item{sampleTable}{a dataframe of sample table
    containing 8 colmuns for Intronic PASs: 'samplename','condition'}
    \item{mutiraw}{a dataframe output
    obtained using either PASEXP_3UTR or PASEXP_IPA}
    \item{conKET}{the name of control in the sampletable, default is 'NT'}
    \item{trtKEY}{the name of control in the sampletable, default is 'KD'}
    \item{PAS}{type of PAS analyzed, either '3UTR' or 'IPA', default is '3UTR'}
    \item{CUTreads}{reads cutoff used for the analysis, default is 0}
	\item{p_methods}{p value correction method, the method can be 
	"holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none", 
	default is "fdr"}
    }
\description{
    Calculate delta relative expression (RED) and
    statistics significance between two sample groups.
}
\value{
    The function \code{APAdiff} return a dataframe containning RED, pvalue and
    regulation pattern (UP, DN or NC) for
    either each gene (3'UTR APA) or each PAS (IPA).
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
    ## Analysis 3'UTR APA between KD and NT group using muti-replicates
    test_3UTRmuti=APAdiff(sampleTable1,DFUTRraw,
    conKET='NT',trtKEY='KD',PAS='3UTR',CUTreads=0,p_methods="fdr")

    ## Analysis 3'UTR APA between KD and NT group without replicates
    test_3UTRsing=APAdiff(sampleTable2,DFUTRraw,
    conKET='NT',trtKEY='KD',PAS='3UTR',CUTreads=0,p_methods="fdr")

    ## Analysis IPA between KD and NT group
    test_IPAmuti=APAdiff(sampleTable1,IPA_OUTraw,
    conKET='NT',trtKEY='KD',PAS='IPA',CUTreads=0,p_methods="fdr")

    ## Analysis IPA between KD and NT group without replicates
    test_IPAsing=APAdiff(sampleTable2,IPA_OUTraw,
    conKET='NT',trtKEY='KD',PAS='IPA',CUTreads=0,p_methods="fdr")

}
