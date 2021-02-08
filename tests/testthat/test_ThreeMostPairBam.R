test_that("ThreeMostPairBam works properly", {
	library("pasillaBamSubset")
	BamfilePath=untreated3_chr4()
    ThreeMostPairBam (BamfilePath=untreated3_chr4(), 
						OutDirPath=getwd(), 
						StrandType='forward-reverse')	
	expect_true(endsWith(basename(BamfilePath), '.bam'))
})
