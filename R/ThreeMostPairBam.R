ThreeMostPairBam <- function(BamfilePath, OutDirPath, 
							StrandType='NONE'){
strandMode=0
if (StrandType=='forward-reverse'){
strandMode=1
} else if(StrandType=='reverse-forward'){
strandMode=2
}
print(paste0("Reading and extracting (Strand = ",StrandType,"):  ",basename(BamfilePath)))
alnReads <- GenomicAlignments::readGAlignmentPairs(BamfilePath, 
													param = Rsamtools::ScanBamParam(
													flag = Rsamtools::scanBamFlag(isPaired = TRUE, 
													isProperPair = TRUE), what = c("qname","flag", 
																					"rname","strand", 
																					"pos", "qwidth", 
																					"mapq", "cigar", 
																					"mrnm", "mpos", 
																					"isize", "seq", 
																					"qual", "groupid", 
																					"mate_status")),
													strandMode=strandMode)													
if (StrandType=='reverse-forward'){
Reads <- GenomicAlignments::first(alnReads)
} else if(StrandType=='forward-reverse' | StrandType=='NONE'){
Reads <- GenomicAlignments::last(alnReads)
} else {
Print("ERROR: StrandType is NOT Valid")						
}

OutfileName = paste0(OutDirPath,"/",
					gsub('.bam','.3most.bam',
					basename(BamfilePath)))
					
print(paste0("Exporting: ",basename(OutfileName)))				
rtracklayer::export(Reads, 
					Rsamtools::BamFile(OutfileName))
}

