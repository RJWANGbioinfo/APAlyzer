test_that("APAlyzer works properly", {
	URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
	file="mm9_REF.RData"
	source_data(paste0(URL,file,"?raw=True"))
	refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
	UTRdbraw=REF3UTR(refUTRraw)
	expect_true(unique(as.data.frame(UTRdbraw)$seqname)=='chr19')	
	
	library("TBX20BamSubset")
	library("Rsamtools")
	flsall = getBamFileList()
	DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
	expect_true(median(DFUTRraw$SRR316189_3UTR_RE,na.rm=TRUE)<0)	
})
