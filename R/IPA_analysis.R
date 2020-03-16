#### IPA analysis ####

.getname<-function(fls){
    samplenameraw=gsub(".bam","",names(fls))
    samplename=gsub(".Aligned.sortedByCoord.out","",samplenameraw)
    return(samplename)
}

.get_strand_IPA<-function(STRINFOR){
    if(STRINFOR=="forward")
        {
            IPAstr=1
        } else if(STRINFOR=="invert"){
            IPAstr=2
        } else {
            IPAstr=0
        }
    return(IPAstr)
}

.getIPAtbl<-function(rawtbl, samplename,rawlength,nametail){
    tepdf=as.data.frame(rawtbl$counts)
    tepdf$IPAID=rownames(tepdf)
    names(tepdf)[1]=paste0(samplename,paste0('_',nametail,'reads'))
    readscols=names(tepdf)[1]
    tepdf=merge(tepdf, rawlength, by="IPAID")
    RPKcols=paste0(samplename,"_RPK")
    tepdf[,RPKcols]=tepdf[,readscols]/(tepdf[,"length"]/1000)
    TPMcols=paste0(samplename,paste0('_',nametail,'RPK'))
    tepdf[TPMcols]=tepdf[RPKcols]
    tepdf=tepdf[,c("IPAID","PASid", "gene_symbol",readscols, TPMcols)]
    return(tepdf)
}

.getTEtbl<-function(rawtbl, samplename,rawlength,nametail){
    tepdf=as.data.frame(rawtbl$counts)
    tepdf$gene_symbol=rownames(tepdf)
    names(tepdf)[1]=paste0(samplename,paste0('_',nametail,'reads'))
    readscols=names(tepdf)[1]
    tepdf=merge(tepdf, rawlength, by="gene_symbol")
    RPKcols=paste0(samplename,"_RPK")
    tepdf[,RPKcols]=tepdf[,readscols]/(tepdf[,"length"]/1000)
    TPMcols=paste0(samplename,paste0('_',nametail,'RPK'))
    tepdf[TPMcols]=tepdf[RPKcols]
    tepdf=tepdf[,c("gene_symbol", readscols, TPMcols)]
    return(tepdf)
}

.simple_size<-function(dflength){
    dflength$length=abs(dflength$end-dflength$start)
    return(dflength)
}


.REFIPA<-function(dfIPAraw){
    dfIPAraw$IPA_up_start=dfIPAraw$upstreamSS
    dfIPAraw$IPA_up_end=dfIPAraw$Pos
    dfIPAraw$IPA_dn_start=dfIPAraw$Pos
    dfIPAraw$IPA_dn_end=dfIPAraw$downstreamSS
    dfIPAraw[which(dfIPAraw$Strand=='-'),]$IPA_up_start=
        dfIPAraw[which(dfIPAraw$Strand=='-'),]$Pos
    dfIPAraw[which(dfIPAraw$Strand=='-'),]$IPA_up_end=
        dfIPAraw[which(dfIPAraw$Strand=='-'),]$upstreamSS
    dfIPAraw[which(dfIPAraw$Strand=='-'),]$IPA_dn_start=
        dfIPAraw[which(dfIPAraw$Strand=='-'),]$downstreamSS
    dfIPAraw[which(dfIPAraw$Strand=='-'),]$IPA_dn_end=
        dfIPAraw[which(dfIPAraw$Strand=='-'),]$Pos
    dfIPAout=dfIPAraw
    dfIPAout$IPAID=paste0(dfIPAout$PASid, dfIPAout$gene_symbol)
    dfIPAout$UP_size=abs(dfIPAout$IPA_up_end-dfIPAout$IPA_up_start)
    dfIPAout$DN_size=abs(dfIPAout$IPA_dn_end-dfIPAout$IPA_dn_start)
    return(dfIPAout)
}

.REFLE<-function(dfLEraw){
    dfLEraw$start=dfLEraw$LEstart
    dfLEraw$end=dfLEraw$TES
    dfLEraw[which(dfLEraw$Strand=='-'),]$start=
        dfLEraw[which(dfLEraw$Strand=='-'),]$TES
    dfLEraw[which(dfLEraw$Strand=='-'),]$end=
        dfLEraw[which(dfLEraw$Strand=='-'),]$LEstart
    LEref=dfLEraw[,c('gene_symbol','Chrom','start','end','Strand')]
    colnames(LEref)=c('GeneID','Chr','Start','End','Strand')
    return(LEref)
}
.calIPA_RE<-function(dfREtbl,samplename){
    colnames(dfREtbl)=c('gene_symbol','IPAID','PASid',
                        'IPA_UPreads', 'IPA_UPRPK',
                        'IPA_DNreads', 'IPA_DNRPK','LEreads','LERPK')
    dfsubIPA=dfREtbl[which(dfREtbl$IPA_UPreads>0 & dfREtbl$LEreads>0),]
    dfsubIPA=dfsubIPA[which(dfsubIPA$IPA_UPRPK > dfsubIPA$IPA_DNRPK),]
    dfsubIPA[,'IPA_RE']=
        log2((dfsubIPA$IPA_UPRPK-dfsubIPA$IPA_DNRPK)/dfsubIPA$LERPK)
    dfsubIPAXXX=dfsubIPA[,c('IPAID','IPA_RE')]
    colnames(dfsubIPAXXX)=c('IPAID',paste0(samplename,'_IPA_RE'))
    dfsubIPAXXX=dfsubIPAXXX[!duplicated(dfsubIPAXXX),]
    return(dfsubIPAXXX)

}

PASEXP_IPA<-function(dfIPAraw, dfLEraw, flS, 
						Strandtype="NONE", nts=1, minMQS=0){
    for (k in seq_len(length(flS))){
        fls=flS[k]
        STRINFOR=Strandtype
        SPNAME=.getname(fls)
        dfIPAout=.REFIPA(dfIPAraw)
        LEref=.REFLE(dfLEraw)
        comnames=c('IPAID','PASid','gene_symbol','Chrom',
                    'Strand','start','end')
        dfIPA_UP=dfIPAout[,c('IPAID','PASid',
                            'gene_symbol','Chrom','Strand',
                            'IPA_up_start','IPA_up_end')]
        dfIPA_DN=dfIPAout[,c('IPAID','PASid',
                            'gene_symbol','Chrom','Strand',
                            'IPA_dn_start','IPA_dn_end')]
        colnames(dfIPA_UP)=comnames
        colnames(dfIPA_DN)=comnames

        oldcolnames=c('IPAID','Chrom','start','end','Strand')
        newcolnames=c('GeneID','Chr','Start','End','Strand')
        UPref=dfIPA_UP[,oldcolnames]
        colnames(UPref)=newcolnames
        DNref=dfIPA_DN[,oldcolnames]
        colnames(DNref)=newcolnames

        UPlength=.simple_size(dfIPA_UP)
        DNlength=.simple_size(dfIPA_DN)
        lencols=c('IPAID','PASid','gene_symbol','Chrom','length')
        UPlength=UPlength[,lencols]
        DNlength=DNlength[,lencols]

        dfLEraw$start=dfLEraw$LEstart
        dfLEraw$end=dfLEraw$TES
        dfLEraw[which(dfLEraw$Strand=='-'),]$start=
            dfLEraw[which(dfLEraw$Strand=='-'),]$TES
        dfLEraw[which(dfLEraw$Strand=='-'),]$end=
            dfLEraw[which(dfLEraw$Strand=='-'),]$LEstart
        LElength=.simple_size(dfLEraw)
        LElength=LElength[,c('gene_symbol','Chrom','Strand','length')]

        IPAstr=.get_strand_IPA(STRINFOR)
        UPtblraw = featureCounts(fls,annot.ext=UPref,
                                    strandSpecific=IPAstr,minMQS=minMQS,
                                    allowMultiOverlap=TRUE,nthreads=nts)
        DNtblraw = featureCounts(fls,annot.ext=DNref,
                                    strandSpecific=IPAstr,minMQS=minMQS,
                                    allowMultiOverlap=TRUE,nthreads=nts)
        LEtblraw = featureCounts(fls,annot.ext=LEref,
                                    strandSpecific=IPAstr,minMQS=minMQS,
                                    allowMultiOverlap=TRUE,nthreads=nts)

        UPtbl=.getIPAtbl(UPtblraw, SPNAME,UPlength,'IPA_UP')
        DNtbl=.getIPAtbl(DNtblraw , SPNAME,DNlength,'IPA_DN')
        LEtbl=.getTEtbl(LEtblraw, SPNAME,LElength,'LE')

        dfout=merge(UPtbl,DNtbl,by=c("IPAID", "PASid","gene_symbol"))
        dfout=merge(dfout,LEtbl,by="gene_symbol")
        dfoutRE=.calIPA_RE(dfout,SPNAME)
        dfout=merge(dfout,dfoutRE, by = "IPAID", all.x=TRUE)
        dfout[is.na(dfout)] = 0

        if(k==1)
        {
            dfIPA_OUT=dfout
        } else {
            dfIPA_OUT=merge(dfIPA_OUT,dfout,
                            by = c("gene_symbol","IPAID","PASid"), all.x=TRUE)
        }
        print(paste0(SPNAME, ", Strand: ", STRINFOR, ", finished"))
        }
    return(dfIPA_OUT)

}

