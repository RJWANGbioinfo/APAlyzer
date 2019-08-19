test_that("IPA_analysis works properly", {
    extpath = system.file("extdata", "mm9_REF.RData", package="APAlyzer")
	load(extpath)
	dfIPA=dfIPA[which(dfIPA$Chrom=='chr19'),]
	dfLE=dfLE[which(dfLE$Chrom=='chr19'),]	
	library("TBX20BamSubset")
	library("Rsamtools")
	flsall = getBamFileList()
	IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall[1:2], Strandtype="forward", nts=1)
	expect_true(median(IPA_OUTraw$SRR316185_IPA_RE,na.rm=TRUE)>=0)	
	expect_true(median(IPA_OUTraw$SRR316184_IPA_RE,na.rm=TRUE)>=0)
})
