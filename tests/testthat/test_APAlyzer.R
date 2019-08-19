test_that("APAlyzer works properly", {
    extpath = system.file("extdata", "mm9_REF.RData", package="APAlyzer")
	load(extpath)
	refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
	UTRdbraw=REF3UTR(refUTRraw)
	expect_true(unique(as.data.frame(UTRdbraw)$seqname)=='chr19')	
	
	library("TBX20BamSubset")
	library("Rsamtools")
	flsall = getBamFileList()
	DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
	expect_true(median(DFUTRraw$SRR316189_3UTR_RE,na.rm=TRUE)<0)	
})
