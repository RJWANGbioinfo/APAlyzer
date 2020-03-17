test_that("GTF2PAS works properly", {
	download.file(url='ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz',
              destfile='Mus_musculus.GRCm38.99.gtf.gz')			  
	GTFfile="Mus_musculus.GRCm38.99.gtf.gz"	
    PASREF=GTF2PAS(GTFfile)	
	refUTRraw=PASREF$refUTRraw
    dfIPA=PASREF$dfIPA
	dfLE=PASREF$dfLE
	expect_true(nrow(refUTRraw)>0)
})