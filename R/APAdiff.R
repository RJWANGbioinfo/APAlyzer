#### APA comparison between groups ####
.absMIN <- function(x) {
    abs(x)[which.min( abs(x) )]
}
.getSPs<-function(sptbl,grpname){
    GRsamples=sptbl[which(sptbl$condition==grpname),]$samplename
    return(GRsamples)
}

.judge_rep<-function(trtsamples,consamples){
    if (length(trtsamples)==1 & length(consamples)==1){
        reptype='single'
        } else if(length(trtsamples)>1 & length(consamples)>1){
            reptype='multi'
        } else {
            reptype='ERROR'
        }
    return(reptype)
}

.calDAbn<-function(tblin,trt,con){
    col1=paste0(trt,'_aRPKM')
    col3=paste0(trt,'_cRPKM')
    col2=paste0(con,'_aRPKM')
    col4=paste0(con,'_cRPKM')
    tblin$AbnT=tblin[,col3]/(tblin[,col1]+tblin[,col3])
    tblin$AbnC=tblin[,col4]/(tblin[,col2]+tblin[,col4])
    tblin$DAbn=tblin$AbnT-tblin$AbnC
    return(tblin)
}

.calIPADAbn<-function(tblin,trt,con){
    col1=paste0(trt,'_IPA_UPRPK')
    col9=paste0(trt,'_LERPK')
    col2=paste0(con,'_IPA_UPRPK')
    col10=paste0(con,'_LERPK')
    tblin$AbnT=tblin[,col1]/(tblin[,col1]+tblin[,col9])
    tblin$AbnC=tblin[,col2]/(tblin[,col2]+tblin[,col10])
    tblin$DAbn=tblin$AbnT-tblin$AbnC
    return(tblin)
}


.caltype<-function(tblin){
    tblin$type='NC'
    if(length(tblin[which(tblin$pv<0.05 & tblin$DAbn>0.05 & tblin$RED<0),]$type)>0){
        tblin[which(tblin$pv<0.05 & tblin$DAbn>0.05 & tblin$RED<0),]$type='DN'
    }

    if(length(tblin[which(tblin$pv<0.05 & tblin$DAbn< -0.05 & tblin$RED>0),]$type)>0){
        tblin[which(tblin$pv<0.05 & tblin$DAbn< -0.05 & tblin$RED>0),]$type='UP'
    }
    return(tblin)
}

.caltype4<-function(tblin){
    tblin$type='NC'
    if(length(tblin[which(tblin$pv<0.05 & tblin$DAbn>0.05 & tblin$RED>0),]$type)>0){
        tblin[which(tblin$pv<0.05 & tblin$DAbn>0.05 & tblin$RED>0),]$type='UP'
    }

    if(length(tblin[which(tblin$pv<0.05 & tblin$DAbn< -0.05 & tblin$RED<0),]$type)>0){
        tblin[which(tblin$pv<0.05 & tblin$DAbn< -0.05 & tblin$RED<0),]$type='DN'
    }
    return(tblin)
}

.caltype2<-function(tblin){
    tblin$APAreg='NC'
    if(length(tblin[which(tblin$pvalue<0.05 & tblin$RED<0),]$APAreg)>0){
        tblin[which(tblin$pvalue<0.05 & tblin$RED<0),]$APAreg='DN'
    }
    if(length(tblin[which(tblin$pvalue<0.05 & tblin$RED>0),]$APAreg)>0){
        tblin[which(tblin$pvalue<0.05 & tblin$RED>0),]$APAreg='UP'
    }
    return(tblin)
}

.caltype3<-function(tblin){
    tblin$APAreg='NC'
    if(length(tblin[which(tblin$pvalue<0.05 & tblin$RED<0),]$APAreg)>0){
        tblin[which(tblin$pvalue<0.05 & tblin$RED<0),]$APAreg='DN'
    }
    if(length(tblin[which(tblin$pvalue<0.05 & tblin$RED>0),]$APAreg)>0){
        tblin[which(tblin$pvalue<0.05 & tblin$RED>0),]$APAreg='UP'
    }
    return(tblin)
}

.calDRUD<-function(dfinraw,trt,con,CUTreads){
    col1=paste0(trt,'_aRPKM')
    col2=paste0(con,'_aRPKM')
    col3=paste0(trt,'_cRPKM')
    col4=paste0(con,'_cRPKM')
    col5=paste0(trt,'_areads')
    col6=paste0(con,'_areads')
    col7=paste0(trt,'_creads')
    col8=paste0(con,'_creads')
    log2col=paste0('indiRUD_',trt,'_',con)
    pcol=paste0('pv_',trt,'_',con)
    tycol=paste0('APAreg_',trt,'_',con)
    dfin=dfinraw[which(dfinraw[,col5]>CUTreads &
                            dfinraw[,col6]>CUTreads &
                            dfinraw[,col7]>CUTreads &
                            dfinraw[,col8]>CUTreads),]
    dfsub=dfin[,c('gene_symbol',col5,col6,col7,col8,col1,col2,col3,col4)]
    dfsub$pv = apply(dfsub[,c(col5,col6,col7,col8)], 1,
                    function(x) fisher.test(matrix(x,nrow=2))$p.value)
    dfsub$log2pA2=log2(dfsub[,col1]/dfsub[,col2])
    dfsub$log2pA1=log2(dfsub[,col3]/dfsub[,col4])
    dfsub$RED=dfsub$log2pA2-dfsub$log2pA1
    dfsub=.calDAbn(dfsub,trt,con)
    dfsub=.caltype(dfsub)
    dfsub=dfsub[c('gene_symbol','RED','pv','type')]
    colnames(dfsub)=c('gene_symbol',log2col,pcol,tycol)
    dfsub=dfsub[!duplicated(dfsub),]
    dfinraw=merge(dfinraw,dfsub,by='gene_symbol',all.x=TRUE)
    return(dfinraw)
}

.calsizeCUT<-function(dfpair,proCUT){
    pairszie=length(dfpair$consamples)
    sizeCUT=pairszie*proCUT
    return(sizeCUT)
}

.calt_p<-function(dfsubXXX,col9S,col10S,adjust_methods){
    trtlen=length(col9S)
    conlen=length(col10S)    
    if(trtlen>1 & conlen>1){
    dfsubXXX$pvalue = apply(dfsubXXX[,c(col9S,col10S)], 
                    1, function (x) {t.test(x[seq_len(trtlen)],
                    x[(1+trtlen):(trtlen+conlen)])$p.value})	
    }

    if(trtlen>1 & conlen==1){
    dfsubXXX$pvalue = apply(dfsubXXX[,c(col9S,col10S)], 
                    1, function (x) {t.test(x[seq_len(trtlen)],
                    mu=x[(1+trtlen):(trtlen+conlen)])$p.value})
    }    
    
    if(trtlen==1 & conlen>1){
    dfsubXXX$pvalue = apply(dfsubXXX[,c(col9S,col10S)], 
                    1, function (x) {t.test(mu=x[seq_len(trtlen)],
                    x[(1+trtlen):(trtlen+conlen)])$p.value})
    }
	dfsubXXX$p_adj	= p.adjust(dfsubXXX$pvalue, method = adjust_methods)
    return(dfsubXXX)

}

.trimCommon<-function(dfinraw,READSCOLS,CUTreads){
    for (colp in READSCOLS){
        dfinraw=dfinraw[which(dfinraw[,colp]>CUTreads),]
    }    
    return(dfinraw)
}

.APA3_muti<-function(dfinput, dfpair, CUTreads){
    for (l in seq_len(length(dfpair$trtsamples))){
            trt=as.character(dfpair$trtsamples[l])
            con=as.character(dfpair$consamples[l])
            dfinput=.calDRUD(dfinput, trt,con,CUTreads)

        }
    return(dfinput)
}

.finaltype<-function(tblin,sizeCUT){
    tblin$APAreg='NC'
    if(length(tblin[which(tblin$coutUP>=sizeCUT & tblin$RED<0),]$APAreg)>0){
        tblin[which(tblin$coutUP>=sizeCUT & tblin$RED<0),]$APAreg='UP'
    }
    if(length(tblin[which(tblin$coutDN>=sizeCUT & tblin$RED>0),]$APAreg)>0){
        tblin[which(tblin$coutDN>=sizeCUT & tblin$RED>0),]$APAreg='DN'
    }
    return(tblin)
}

.final_tbl_3mutil<-function(mutiraw,dfpair,CUTreads,proCUT){
    xxxx=.APA3_muti(mutiraw, dfpair, CUTreads)
    RUDcols=c(grep("indiRUD_", colnames(xxxx)))
    pvcols=c(grep("pv_", colnames(xxxx)))
    typecols=c(grep("APAreg_", colnames(xxxx)))
    xxxx$RED=rowMeans(xxxx[,RUDcols], na.rm = TRUE)
    xxxx$Min_pv=apply(xxxx[, pvcols], 1, .absMIN)
    xxxx$coutUP=rowSums(xxxx[, typecols]=='UP')
    xxxx$coutDN=rowSums(xxxx[, typecols]=='DN')
    sizeCUT=.calsizeCUT(dfpair,proCUT)
    xxxx=.finaltype(xxxx,sizeCUT)
    xxxx=xxxx[,c('gene_symbol','RED','Min_pv','APAreg')]
    xxxx=xxxx[!duplicated(xxxx),]
    return(xxxx)
}

.final_tbl_3mutil2<-function(mutiraw, trtsamples,consamples, CUTreads,p_methods){
    col5S=paste0(trtsamples,'_areads')
    col6S=paste0(consamples,'_areads')
    col7S=paste0(trtsamples,'_creads')
    col8S=paste0(consamples,'_creads')
    col9S=paste0(trtsamples,'_3UTR_RE')
    col10S=paste0(consamples,'_3UTR_RE')
    dfsubXXX=mutiraw[,c('gene_symbol',col5S,col6S,col7S,col8S,col9S,col10S)]
    READSCOLS=c(col5S,col6S,col7S,col8S)    
    dfsubXXX=.trimCommon(dfsubXXX,READSCOLS,CUTreads)    
    dfsubXXX$trt_RE=rowMeans(dfsubXXX[,col9S], na.rm = TRUE)
    dfsubXXX$con_RE=rowMeans(dfsubXXX[,col10S], na.rm = TRUE)
    dfsubXXX$RED=dfsubXXX$trt_RE-dfsubXXX$con_RE
    dfsubXXX=.calt_p(dfsubXXX,col9S,col10S,p_methods)    
    dfsubXXX=.caltype2(dfsubXXX)
    dfsubXXX=dfsubXXX[c('gene_symbol','RED','pvalue','p_adj','APAreg')]
    dfsubXXX=dfsubXXX[!duplicated(dfsubXXX),]
    return(dfsubXXX)    
}

.final_tbl_3singl<-function(mutiraw,trt,con,CUTreads,adjust_methods){
    xxxx=.calDRUD(mutiraw, trt,con,CUTreads)
    log2col=paste0('indiRUD_',trt,'_',con)
    pcol=paste0('pv_',trt,'_',con)
    tycol=paste0('APAreg_',trt,'_',con)
	xxxx$p_adj	= p.adjust(xxxx[,pcol], method = adjust_methods)
    xxxx=xxxx[,c('gene_symbol',log2col,pcol,'p_adj',tycol)]
    colnames(xxxx)=c('gene_symbol','RED','pvalue','p_adj','APAreg')
    xxxx=xxxx[!duplicated(xxxx),]
    return(xxxx)
}

.trimIPAraw<-function(dfinraw,trt,con,CUTreads){
    col1=paste0(trt,'_IPA_UPRPK')
    col2=paste0(con,'_IPA_UPRPK')
    col9=paste0(trt,'_IPA_DNRPK')
    col10=paste0(con,'_IPA_DNRPK')
    col5=paste0(trt,'_IPA_UPreads')
    col6=paste0(con,'_IPA_UPreads')
    col7=paste0(trt,'_LEreads')
    col8=paste0(con,'_LEreads')
    readCUTcols=c(col5,col6,col7,col8)
    for (i in seq_len(length(readCUTcols))){
        dfinraw=dfinraw[which(dfinraw[,readCUTcols[i]]>CUTreads),]
    }
    dfinraw=dfinraw[which(dfinraw[,col1]>dfinraw[,col9] &
                                dfinraw[,col2]>dfinraw[,col10]),]
    return(dfinraw)

}

.trimIPAraw2<-function(dfinraw,col9S,col10S){
    for (colp in c(col9S,col10S)){
        dfinraw=dfinraw[which(dfinraw[,colp]!=Inf & dfinraw[,colp]!=-Inf),]
    }
    return(dfinraw)
}

.calDRUDIPA<-function(dfinraw,trt,con,CUTreads){
    dfin=.trimIPAraw(dfinraw,trt,con,CUTreads)
    dfin=.calIPADAbn(dfin,trt,con)
    col1=paste0(trt,'_IPA_RE')
    col2=paste0(con,'_IPA_RE')
    col5=paste0(trt,'_IPA_UPreads')
    col6=paste0(con,'_IPA_UPreads')
    col7=paste0(trt,'_LEreads')
    col8=paste0(con,'_LEreads')
    log2col=paste0('indiRUD_',trt,'_',con)
    pcol=paste0('pv_',trt,'_',con)
    tycol=paste0('APAreg_',trt,'_',con)
    dfsub=dfin[,c('gene_symbol','IPAID','PASid',
                    col5,col6,col7,col8,col1,col2,'DAbn')]
    dfsub$pv = apply(dfsub[,c(col5,col6,col7,col8)], 1,
                    function(x) fisher.test(matrix(x,nrow=2))$p.value)
    dfsub$RED=dfsub[,col1]-dfsub[,col2]
    dfsub=.caltype4(dfsub)
    dfsub=dfsub[c('gene_symbol','PASid','RED','pv','type')]
    colnames(dfsub)=c('gene_symbol','PASid',log2col,pcol,tycol)
    dfsub=dfsub[!duplicated(dfsub),]
    dfinraw=merge(dfinraw,dfsub,by=c('gene_symbol','PASid'),all.x=TRUE)
    return(dfinraw)
}


.APAIPA_muti<-function(dfinput, dfpair, CUTreads){
    for (l in seq_len(length(dfpair$trtsamples))){
            trt=as.character(dfpair$trtsamples[l])
            con=as.character(dfpair$consamples[l])
            dfinput=.calDRUDIPA(dfinput, trt,con,CUTreads)
        }
    return(dfinput)
}

.final_tbl_IPAmutil<-function(mutiraw,dfpair,CUTreads,proCUT){
    xxxx=.APAIPA_muti(mutiraw, dfpair, CUTreads)
    RUDcols=c(grep("indiRUD_", colnames(xxxx)))
    pvcols=c(grep("pv_", colnames(xxxx)))
    typecols=c(grep("APAreg_", colnames(xxxx)))
    xxxx$RED=rowMeans(xxxx[,RUDcols], na.rm = TRUE)
    xxxx$Min_pv=apply(xxxx[, pvcols], 1, .absMIN)
    xxxx$coutUP=rowSums(xxxx[, typecols]=='UP')
    xxxx$coutDN=rowSums(xxxx[, typecols]=='DN')
    sizeCUT=.calsizeCUT(dfpair,proCUT)
    xxxx=.finaltype(xxxx,sizeCUT)
    xxxx=xxxx[,c('gene_symbol','PASid','RED','Min_pv','APAreg')]
    xxxx=xxxx[!duplicated(xxxx),]
    return(xxxx)
}

.final_tbl_IPAmutil2<-function(mutiraw,trtsamples,consamples,CUTreads,p_methods){
    col5S=paste0(trtsamples,'_IPA_UPreads')
    col6S=paste0(consamples,'_IPA_UPreads')
    col7S=paste0(trtsamples,'_LEreads')
    col8S=paste0(consamples,'_LEreads')
    col9S=paste0(trtsamples,'_IPA_RE')
    col10S=paste0(consamples,'_IPA_RE')
    dfsubXXX=mutiraw[,c('gene_symbol','PASid',col5S,col6S,col7S,
                            col8S,col9S,col10S)]
    READSCOLS=c(col5S,col6S,col7S,col8S)    
    dfsubXXX=.trimCommon(dfsubXXX,READSCOLS,CUTreads)
    dfsubXXX=.trimIPAraw2(dfsubXXX,col9S,col10S)
    dfsubXXX$trt_RE=rowMeans(dfsubXXX[,col9S], na.rm = TRUE)
    dfsubXXX$con_RE=rowMeans(dfsubXXX[,col10S], na.rm = TRUE)
    dfsubXXX$RED=dfsubXXX$trt_RE-dfsubXXX$con_RE
    dfsubXXX=.calt_p(dfsubXXX,col9S,col10S,p_methods)    
    dfsubXXX=.caltype3(dfsubXXX)
    dfsubXXX=dfsubXXX[c('gene_symbol','PASid','RED','pvalue','p_adj','APAreg')]
    dfsubXXX=dfsubXXX[!duplicated(dfsubXXX),]
    return(dfsubXXX)        
}

.final_tbl_IPAsingl<-function(mutiraw,trt,con,CUTreads,adjust_methods){
    xxxx=.calDRUDIPA(mutiraw, trt,con,CUTreads)
    log2col=paste0('indiRUD_',trt,'_',con)
    pcol=paste0('pv_',trt,'_',con)
    tycol=paste0('APAreg_',trt,'_',con)
	xxxx$p_adj	= p.adjust(xxxx[,pcol], method = adjust_methods)
    xxxx=xxxx[,c('gene_symbol','PASid',log2col,pcol,'p_adj',tycol)]
    colnames(xxxx)=c('gene_symbol','PASid','RED','pvalue','p_adj','APAreg')
    xxxx=xxxx[!duplicated(xxxx),]
    return(xxxx)
}

APAdiff<-function(sampleTable,mutiraw, conKET='NT',trtKEY='KD',
                    PAS='3UTR',CUTreads=0,p_methods="fdr"){
    consamples=.getSPs(sampleTable,conKET)
    trtsamples=.getSPs(sampleTable,trtKEY)
    reptypeRAW=.judge_rep(trtsamples,consamples)
    if(reptypeRAW=='multi' & PAS=='3UTR'){
        APA_diff=.final_tbl_3mutil2(mutiraw,trtsamples,consamples,CUTreads,p_methods)
    }

    if(reptypeRAW=='single' & PAS=='3UTR'){
        APA_diff=.final_tbl_3singl(mutiraw,trtsamples,consamples,CUTreads,p_methods)
    }

    if(reptypeRAW=='multi' & PAS=='IPA'){
        APA_diff=.final_tbl_IPAmutil2(mutiraw,trtsamples,consamples,CUTreads,p_methods)
    }
    if(reptypeRAW=='single' & PAS=='IPA'){
        APA_diff=.final_tbl_IPAsingl(mutiraw,trtsamples,consamples,CUTreads,p_methods)
    }
    if(reptypeRAW=='ERROR'){
    print("Sample matrix, error, please check your sample table")
    }
    APA_diff=APA_diff[!is.na(APA_diff$RED),]
    return(APA_diff)
}
