#### APA Volcano plotting ####
APAVolcano <- function(df, Pcol = "pvalue",PAS='3UTR',
                        top = -1, markergenes = NULL,
                        y_cutoff = 0.05,xlab = "RED", ylab = "-Log10(P-value)",
                        PAScolor = c("gray80", "red", "blue"),
                        alpha = 0.75, plot_title = NULL,
						width = 4, height = 2.5){
	x = "RED"
	y = Pcol
	force = 0.1   
	if(PAS=='3UTR'){
			gg = df[, c("gene_symbol",x, y,'APAreg')]
			gg[, y] = -log10(gg[, y])
			Label = "gene_symbol"
			gg$Label = gg[, Label]
		  }	else if (PAS=='IPA') {
			gg = df[, c("gene_symbol","PASid",x, y,'APAreg')]
			gg[, y] = -log10(gg[, y])		  
			gg$Label=paste(gg$gene_symbol,gg$PASid,sep=":")
		  }
	
    if(!(top==0 & is.null(markergenes))){
      
      gg = gg[order(gg[,y], abs(gg[,x]), decreasing = TRUE), ]
      idx1 = idx2 = c()
      if(top>0){
        idx1 = which(gg$APAreg=="UP")[1:min(top, sum(gg$APAreg=="UP"))]
        idx2 = which(gg$APAreg=="DN")[1:min(top, sum(gg$APAreg=="DN"))]
      }
      idx = unique(c(idx1, idx2, which(gg$Label %in% markergenes)))
      gg$Label = as.character(gg$Label)
      gg$Label[setdiff(1:nrow(gg), idx)] = ""
    }	
	
    gg$color = gg$APAreg
    gg$color[gg$Label!=""] = "black"
    PAScolor = c(PAScolor, "black")
    names(PAScolor) = c("NC", "UP", "DN", "black")
    p = ggplot(gg, aes(x=gg[,x], y=gg[,y], label=Label))
    p = p + geom_point(aes(fill=APAreg), shape = 21, alpha=alpha)
    p = p + geom_point(aes(colour=color), shape = 21, alpha=alpha, show.legend = FALSE)
    if(!(top==0 & is.null(markergenes))){
      p = p + ggrepel::geom_text_repel(aes(gg[,x], gg[,y], color = APAreg, label = Label),
	  force = force, fontface = 'bold', size = 3,
                                       nudge_y = 0.2,
									   nudge_x = 0.2,
										direction    = "x", show.legend = FALSE,
										segment.size = 0.2)											
	}	
    p = p + scale_color_manual(values=PAScolor)
    p = p + scale_fill_manual(values=PAScolor)
    p = p + theme(text = element_text(colour="black",size = 14),
                  plot.title = element_text(hjust = 0.5, size=16),
                  axis.text = element_text(colour="gray10"))
    p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank())
    p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
    p = p + geom_vline(xintercept = 0, linetype = "dotted")
    p = p + labs(x=xlab, y=ylab, title=plot_title)

								   
    p = p + theme(legend.position = "top")
    return(p)
}