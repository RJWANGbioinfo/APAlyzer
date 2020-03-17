test_that("APAPlots works properly", {
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
	UTR_APA_PLOT=APAVolcano(test_3UTRmuti, PAS='3UTR', 
							Pcol = "pvalue", top=5, plot_title='3UTR APA')
	expect_type(UTR_APA_PLOT, "list")
})