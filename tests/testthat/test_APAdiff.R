test_that("APAdiff works properly", {
	library("TBX20BamSubset")
	library("Rsamtools")
	flsall = getBamFileList()
	extpath = system.file("extdata", "mm9_TBX20.APAout.RData", package="APAlyzer")
	load(extpath)
    sampleTable1 = data.frame(samplename = c(names(flsall)),
                    condition = c(rep("NT",3),rep("KD",3)))
	test_IPAmuti=APAdiff(sampleTable1,
							IPA_OUTraw, 
							conKET='NT',
							trtKEY='KD',
							PAS='IPA',
							CUTreads=0)
	expect_true(median(test_IPAmuti$RED)>0)
	expect_true(table(test_IPAmuti$APAreg)[1]<table(test_IPAmuti$APAreg)[3])
})
