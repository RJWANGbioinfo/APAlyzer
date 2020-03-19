#### 3'UTR APA and gene EXP analysis ####
.getname<-function(fls){
    samplenameraw=gsub(".bam","",names(fls))
    samplename=gsub(".Aligned.sortedByCoord.out","",samplenameraw)
    return(samplename)
}

.getsize<-function(DB){
    REGIONsize=as.data.frame(width(DB))
    names(REGIONsize)[3]="size"
    names(REGIONsize)[2]="gene_symbol"
    REGIONsize=REGIONsize[,c("gene_symbol","size")]
    return(REGIONsize)
}

.getSUMsize<-function(DB){
    REGIONsize=as.data.frame(sum(width(DB)))
    names(REGIONsize)[1]="length"
    REGIONsize$gene_symbol=rownames(REGIONsize)
    REGIONsize=REGIONsize[,c("gene_symbol","length")]
    return(REGIONsize)
}

.getcount<-function(DB,BMAfile,STRINFOR){
    preprocess.reads.function = NULL
    if (STRINFOR == "INVERT")
        {
            preprocess.reads.function = invertStrand
        }
    counttbl= summarizeOverlaps(DB, BMAfile,
    mode="IntersectionNotEmpty",
    singleEnd=TRUE, ignore.strand=TRUE,
    preprocess.reads=preprocess.reads.function)
    return(counttbl)
}

.divxdb<-function(xdb,keystr){
    x = unlist(xdb)
    xp=x[strand(x)==keystr,]
    xpdb = split(xp, xp$gene_name)
    return(xpdb)
}

.combine_count<-function(xdb,fls,STRINFOR){
    xdbP=.divxdb(xdb,'+')
    xtblP=.getcount(xdbP,fls,STRINFOR)
    xdfP=as.data.frame(assays(xtblP)$counts)
    xdfP$gene_symbol=rownames(xdfP)
    xdbN=.divxdb(xdb,'-')
    xtblN=.getcount(xdbN,fls,STRINFOR)
    xdfN=as.data.frame(assays(xtblN)$counts)
    xdfN$gene_symbol=rownames(xdfN)
    xdf=rbind(xdfP,xdfN)
    return(xdf)
}

REF3UTR<-function(refUTR){
    stopifnot(ncol(refUTR) == 6)
    colnames(refUTR)=c('gene_symbol','Chrom','Strand','First','Last','cdsend')
    ##aUTR##
    refaUTR = refUTR[,c('gene_symbol','Chrom','Strand','First','Last')]
    refaUTR$note = "__aUTR"
    refaUTR$gene_name=paste0(refaUTR$gene_symbol, refaUTR$note)
    refaUTR$Start=refaUTR$First
    refaUTR$End=refaUTR$Last
    refaUTR[which(refaUTR$Strand=='-'),]$Start=
        refaUTR[which(refaUTR$Strand=='-'),]$Last
    refaUTR[which(refaUTR$Strand=='-'),]$End=
        refaUTR[which(refaUTR$Strand=='-'),]$First
    refaUTR=refaUTR[,-c(4:6)]
    ##cUTR##
    refcUTR = refUTR[,c('gene_symbol','Chrom','Strand','cdsend','First')]
    refcUTR$note = "__cUTR"
    refcUTR$gene_name=paste0(refcUTR$gene_symbol, refcUTR$note)
    refcUTR$Start=refcUTR$cdsend
    refcUTR$End=refcUTR$First
    refcUTR[which(refcUTR$Strand=='-'),]$Start=
        refcUTR[which(refcUTR$Strand=='-'),]$First
    refcUTR[which(refcUTR$Strand=='-'),]$End=
        refcUTR[which(refcUTR$Strand=='-'),]$cdsend
    refcUTR=refcUTR[,-c(4:6)]
    ref3UTR=rbind(refaUTR, refcUTR)
    UTRdb = makeGRangesFromDataFrame(ref3UTR, keep.extra.columns = TRUE)
    UTRdb = split(UTRdb, UTRdb$gene_name)
    return(UTRdb)
}

PASEXP_3UTR<-function(UTRdb, flS, Strandtype="NONE"){
    stopifnot(is(UTRdb, "GRangesList"))
    for (k in seq_len(length(flS))){
        fls=flS[k]
        STRINFOR=Strandtype
        SPNAME=.getname(fls)
        acSIZE=.getsize(UTRdb)
        UTRdf=.combine_count(UTRdb,fls,STRINFOR)
        names(UTRdf)[1]=paste0(SPNAME,"_reads")
        readscols=names(UTRdf)[1]
        RPMcols = paste0(SPNAME,"_RPKM")
        UTRdf[RPMcols]=UTRdf[readscols]/colSums(UTRdf[readscols])*1000000
        UTRdf=merge(UTRdf,acSIZE, by="gene_symbol",all.x=TRUE)
        UTRdf[RPMcols]=UTRdf[RPMcols]/UTRdf$size*1000
        aUTRdf=UTRdf[grepl("__aUTR", UTRdf$gene_symbol),]
        cUTRdf=UTRdf[grepl("__cUTR", UTRdf$gene_symbol),]
        aUTRdf$gene_symbol = gsub('__aUTR', '', aUTRdf$gene_symbol)
        names(aUTRdf) = gsub('_RPKM', '_aRPKM', names(aUTRdf))
        names(aUTRdf) = gsub('_reads', '_areads', names(aUTRdf))
        cUTRdf$gene_symbol = gsub('__cUTR', '', cUTRdf$gene_symbol)
        names(cUTRdf) = gsub('_RPKM', '_cRPKM', names(cUTRdf))
        names(cUTRdf) = gsub('_reads', '_creads', names(cUTRdf))
        dfSRSutr=merge(x = aUTRdf, y = cUTRdf, by = "gene_symbol")
        dfSRSutr[is.na(dfSRSutr)] = 0
        dfSRSutr=dfSRSutr[,c(1,2,3,5,6)]
        dfSRSutr[,paste0(SPNAME,'_3UTR_RE')]=
            log2(dfSRSutr[,paste0(SPNAME,'_aRPKM')]/
                    dfSRSutr[,paste0(SPNAME,'_cRPKM')])
        if(k==1)
        {
            df3UTR_OUT=dfSRSutr
        } else {
            df3UTR_OUT=merge(df3UTR_OUT,dfSRSutr,
                            by = "gene_symbol", all.x=TRUE)
        }

        print(paste0(SPNAME, ", Strand: ", STRINFOR, ", finished"))
        }
    return(df3UTR_OUT)
}

REFCDS<-function(txdb,IDDB){
    CDSbygene = cdsBy(txdb, by="gene")
    x = unlist(CDSbygene)
    acc2sym = AnnotationDbi::select(IDDB, keys = names(x),
                                    keytype = "ENTREZID", columns = "SYMBOL")
    acc2sym[is.na(acc2sym$SYMBOL),]$SYMBOL =
        acc2sym[is.na(acc2sym$SYMBOL), ]$ENTREZID
    names(x) = acc2sym$SYMBOL
    CDSbygene = split(x, names(x))
    CDSbygene = reduce(CDSbygene)
    return(CDSbygene)
}

GENEXP_CDS<-function(CDSbygene, flS, Strandtype="NONE"){
    for (k in seq_len(length(flS))){
        fls=flS[k]
        STRINFOR=Strandtype
        SPNAME=.getname(fls)
        CDStbl=.getcount(CDSbygene,fls,STRINFOR)
        cdslength=.getSUMsize(CDSbygene)
        tepdf=as.data.frame(assays(CDStbl)$counts)
        tepdf$gene_symbol=rownames(tepdf)
        names(tepdf)[1]=paste0(SPNAME,"_CDSreads")
        readscols=names(tepdf)[1]
        tepdf=merge(tepdf, cdslength, by="gene_symbol")
        RPKcols=paste0(SPNAME,"_RPK")
        tepdf[,RPKcols]=tepdf[,readscols]/(tepdf[,"length"]/1000)
        TPMcols=paste0(SPNAME,"_TPM")
        tepdf[TPMcols]=tepdf[RPKcols]/colSums(tepdf[RPKcols])*1000000
        tepdf=tepdf[,c("gene_symbol", readscols, TPMcols)]
        if(k==1)
        {
            dfEXP_OUT=tepdf
        } else {
            dfEXP_OUT=merge(dfEXP_OUT,tepdf, by = "gene_symbol", all.x=TRUE)
        }

        print(paste0(SPNAME, ", Strand: ", STRINFOR, ", finished"))
        }

    return(dfEXP_OUT)
}

REF4PAS<-function(refUTRraw, dfIPAraw, dfLEraw){
	UTRdbraw=REF3UTR(refUTRraw)
	dfIPA=dfIPAraw
	dfLE=dfLEraw
	PASREF=list(UTRdbraw,dfIPA,dfLE)
	names(PASREF)=c('UTRdbraw','dfIPA','dfLE')
	return(PASREF)
}

download_testbam<-function(){
URL="https://media.githubusercontent.com/media/RJWANGbioinfo/PAS_reference_RData_and_testing_data/master/bam/"	
SAMPLES=c("Heart_rep1",
"Heart_rep2",
"Heart_rep3",
"Heart_rep4",
"Testis_rep1",
"Testis_rep2",
"Testis_rep3",
"Testis_rep4")
for(SAMPLE in SAMPLES)
{
print(paste0("Download  ",SAMPLE))
download.file(url = paste0(URL,SAMPLE,".bam")
                                   , destfile = paste0(SAMPLE,".bam"))
}
}
