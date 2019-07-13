# About APAlyzer
Alternative Polyadenylation(APA) analysis can be challenging due the high difficulty of de-novo identification of polyadenylation site(PAS) using bulk-RNA seq. To overcome this limitation, we developed a toolkit, APAlyzer, which takes the advantage of the PAS database [PolyA_DB3] (http://exon.njms.rutgers.edu/polya_db/v3/) to analysis APA regulation using bulk RNA-seq. It starts with bam files from RNA-seq mapping, extract the reads mapped to corresponding 3UTR APA, intronic PAS and CDS regions, then calculate their relative expression, compare APA patterns between different groups, and label those significantly regulated APA genes.

## Authors
* **Ruijia Wang**
* **Bin Tian**

## License
This project is licensed under the LGPL-3 License.

# Installing APAlyzer
APAlyzer should be installed as follows using BiocManager:

```{r eval=FALSE}
if (!"BiocManager" %in% rownames(installed.packages()))
     install.packages("BiocManager")
BiocManager::install("APAlyzer")
```
Alternatively, it can also be installed through:
```
R CMD INSTALL APAlyzer.zip
```

# obtain sample datasets 

## RNA-seq bam file
The package needs to read the bam file(s) to obtain the expression information around PAS regions. To demonstrate the analysis process in this vignettes, we first determine the paths to the example BAM files in the `r BiocStyle::Biocpkg("TBX20BamSubset")` data package.

We can extrat sample bamfiles (mouse RNA-seq, mapped to mm9) from:
```{r eval=FALSE}
library("TBX20BamSubset")
library("Rsamtools")
flsall <- getBamFileList()
flsall
```

## Genomic reference regions
Since we are extracting the expression information around PAS regions. Genomic reference of PAS (both 3'UTR APA and IPA) are also required by our package. We have prebuilt the corresponding reference in mouse (mm9), which can be loaded from `extdata`:
```{r}
extpath <- system.file("extdata", "mm9_REF.RData", package="APAlyzer")
load(extpath, verbose=TRUE)
```
This `extdata` covers 3'UTR APA region(refUTRraw), IPA region(dfIPA), and last 3'exon region(dfLE). The `refUTRraw` is a data frame containing 6 columns covering genomic information of 3'UTR PASs:
```{r}
head(refUTRraw,2)
```
`dfIPA` is a data frame containing 8 colmuns for Intronic PASs; 'upstreamSS' means closest 5'/3' splice site to IPA, 'downstreamSS' means closest 3' splice site:
```{r}
head(dfIPA,2)
```
`dfLE` is a data frame containing 5 colmuns for 3' least exon; 'LEstart' means the start genomic position of last 3' exon.
```{r}
head(dfLE,2)
```

In additions to mouse mm9, our package has also prebuild the corresponding regions for human hg19 genome:
```{r}
extpath <- system.file("extdata", "hg19_REF.RData", package="APAlyzer")
load(extpath, verbose=TRUE)
```


# Analyze 3'UTR APA

## Build reference aUTR and cUTR regions
To calculate 3'UTR APA relative expression (RE). We first need to define the refence region of aUTR and cUTR using `refUTRraw`.
```{r eval=FALSE}
refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
UTRdbraw=REF3UTR(refUTRraw)
```
The `REF3UTR` function returns a genomic ranges containing aUTR(pPAS to dPAS) and cUTR(cdsend to pPAS) regions for each gene:
```{r, echo = FALSE}
library(APAlyzer)
extpath <- system.file("extdata", "mm9_REF.RData", package="APAlyzer")
load(extpath)
refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
UTRdbraw=REF3UTR(refUTRraw)
```
```{r}
head(UTRdbraw,2)
```

## Extract reads count and calculte relative expression (RE) of aUTR and cUTR 
Once a/cUTR region is defined, the RE of 3'UTR APA of each gene can be quickly calculated by `PASEXP_3UTR`:
```{r, echo = FALSE}
options(warn=-1)
suppressMessages(library(APAlyzer))
suppressMessages(library("TBX20BamSubset"))
suppressMessages(library("Rsamtools"))
extpath <- system.file("extdata", "mm9_REF.RData", package="APAlyzer")
load(extpath)
flsall <- getBamFileList()
refUTRraw=refUTRraw[which(refUTRraw$Chrom=='chr19'),]
UTRdbraw=REF3UTR(refUTRraw)
```

```{r}
DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
```
The `PASEXP_3UTR` requires two part of input, one is the reference aUTR and cUTR regions, and the other one is the path of bamfile(s). In additions to input, you can also define the strandness of the sequencing using `Strandtype`. The detailed the usage can also be checked through `?PASEXP_3UTR`.
The output data frame covers reads count (in a/cUTRs), RPKM(in a/cUTRs) and RE for each gene:
```{r}
head(DFUTRraw,2)
```

# Analyze Intronic PAS (IPA) regulation

## Build reference IPA regions
Analysis of IPA requires two genomic regions: IPA regions and last 3'exon regions. As mentioned in the sample data section, the two regions in mouse and human are prebuilt in this package:

```{r eval=FALSE}
#mouse(mm9):
extpath <- system.file("extdata", "mm9_REF.RData", package="APAlyzer")
load(extpath)
```

```{r eval=FALSE}
#human(hg19):
extpath <- system.file("extdata", "hg19_REF.RData", package="APAlyzer")
load(extpath)
```

## Extract reads count and calculte relative expression (RE) of IPA 
Similar to 3'UTR APA, RE of IPAs can also be simply calculate using `PASEXP_IPA`:
```{r eval=FALSE}
IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=1)
load(extpath)
```
Note as a specific feature in IPA part, you can set more threads using 'nts=' to increase the calculation speed. And the detailed the usage can also be checked through `?PASEXP_IPA`.

The output data frame contains reads count and reads density IPA upstream (a), IPA down stream (b) and 3'exon regions(c). The RE of IPA is calculated as IPA-RE = log2(a - b)/c.

```{r, echo = FALSE}
options(warn=-1)
suppressMessages(library(APAlyzer))
extpath <- system.file("extdata", "mm9_TBX20.APAout.RData", package="APAlyzer")
load(extpath)
```

```{r}
head(IPA_OUTraw,2)
```

# Identify genes and PAS under significant APA regulation

Once we obtained the expression information of PASs in each sample, we can then compare the APA regulation difference between two different groups. In the analysis, we usually we will meet two type of experimental design: the ones with multi-replicates and the ones without replicates. To fit the situation of both cases, let's first generate two sample tables accordingly:

```{r, echo = FALSE}
extpath <- system.file("extdata", "mm9_TBX20.APAout.RData", package="APAlyzer")
load(extpath)
```

```{r}
# Build the sample table with replicates
sampleTable1 = data.frame(samplename = c(names(flsall)),
					 condition = c(rep("NT",3),rep("KD",3)))
sampleTable1
load(extpath)
```

```{r}
# Build the sample table without replicates
sampleTable2 = data.frame(samplename = c("SRR316184","SRR316187"),
					 condition = c("NT","KD")) 
sampleTable2					 
```


## identify genes with signficantly 3'UTR shortening or lengthening

We then start to analysis 3'UTR APA between KD and NT group without replicates using `sampleTable2`. The function we used here to do the analysis is called `APAdiff` (detailed information can be checked using `?APAdiff`). It will fist to go through the sample table to determine whether it is a multi-replicate design or non-replicate design. Then the APA compassion will be performed:
```{r}
# Analysis 3'UTR APA between KD and NT group using non-repilicate design
test_3UTRsing=APAdiff(sampleTable2,DFUTRraw, conKET='NT',trtKEY='KD',PAS='3UTR',CUTreads=0)
```
The `APAdiff` requires two inputs. One is sampleable defining the groups/conditions of the samples, the other one is the expression information of aUTR and cUTR, which can be obtained using `PASEXP_3UTR` in the previous step. The group name of treat and control group can be defined using `trtKEY=` and `conKET=`; the PAS type analyzed should be defined using `PAS=`; and the reads cutoff used for aUTR and cUTR region can alse be defined using `CUTreads=`, the default value is 0. In the non-replicate design, the APA pattern will be compared between two samples and output into a data frame:
```{r}
head(test_3UTRsing,2)
table(test_3UTRsing$APAreg)
```
The output contains 4 columns: 'gene symbol' describes gene information; 'RUD' describes the delta relative expression between two groups; 'pvalue' (Fisher's exact test) describes the statistical significance; and 'APAreg' point the 3'UTR APA regulation patten in the gene. We defined 3 types in 'APAreg', 'UP' means aUTR abundance in treatment group ('KD' in this case) is at least 5% higher than in control ('NT' in this case), and  'pvalue'<0.05; 'DN' means aUTR abundance is 5% in treatment lower than control with p-value<0.05; 'NC' are the left genes. In the 3'UTR APA level, 'UP' stands for 3'UTR shortening, and 'DN' stands for 3'UTR lengthening.

In the multi-replicate design all the possible paired samples are first compared, then the significance and regulation pattern are evaluated among all the comparison:
```{r}
# Analysis 3'UTR APA between KD and NT group using multi-repilicate design
test_3UTRmuti=APAdiff(sampleTable1,DFUTRraw, conKET='NT',trtKEY='KD',PAS='3UTR',CUTreads=0,proCUT=0.5)
head(test_3UTRmuti,2)
table(test_3UTRmuti$APAreg)
```
In the multi-replicate design, 'RUD' is averaged from all the comparison; 'Min_pv' is the minimum p-value from all the comparison. 
In additions to regulation setting, `APAreg` will use consistency as an extra cutoff to determine APA regulation pattern, which can be defined using `proCUT` (default is 0.5). And in this case 'UP' is defined as up-regulation in at least 50% of the comparison and average 'RUD' is greater than 0; while 'DN' is opposite; and 'NC' are the left genes.

## identify IPAs with signficantly activation or suprresiion
IPA comparison is similar to 3'UTR APA using `APAdiff`, except (1) use IPA expression as input, (2) 'PAS=' needs to define to 'IPA', and (3) the analysis is performed on each IPA.

Analysis IPA between KD and NT group without replicates:
```{r} 
test_IPAsing=APAdiff(sampleTable2,IPA_OUTraw, conKET='NT',trtKEY='KD',PAS='IPA',CUTreads=0) 
head(test_IPAsing,2)
```

Analysis IPA between KD and NT group using muti-replicates and relax consistency cutoff:
```{r}
test_IPAmuti=APAdiff(sampleTable1,IPA_OUTraw, conKET='NT',trtKEY='KD',PAS='IPA',proCUT=0.2,CUTreads=0)
head(test_IPAmuti,2)
```  

# Plot the APA results
The APA comparison can be plotted using either violin plot or CDF curves. Let's use the previous 3'UTR APA and IPA comparison output as an example, first build the plotting data frame: 
```{r, echo = FALSE}
options(warn=-1)
extpath <- system.file("extdata", "mm9_TBX20.APAdiff_OUT.RData", package="APAlyzer")
load(extpath)
```
```{r}
test_3UTRmuti$APA="3'UTR"
test_IPAmuti$APA="IPA"
dfplot=rbind(test_3UTRmuti[,c('RUD','APA')],test_IPAmuti[,c('RUD','APA')])
```
We can then plot violin and CDF curves using `r BiocStyle::Biocpkg("ggplot2")`:
```{r}
library(ggplot2)
ggplot(dfplot, aes(x = APA, y = RUD)) + geom_violin(trim = FALSE) + geom_boxplot(width = 0.2)+theme_classic()+ geom_hline(yintercept=0, linetype="dashed", color = "red")
ggplot(dfplot, aes( x = RUD, color = APA)) + stat_ecdf(geom = "step") +theme_classic()+ylab("cumulative fraction")+ geom_vline(xintercept=0, linetype="dashed", color = "gray")+ geom_hline(yintercept=0.5, linetype="dashed", color = "gray")
```

# Analyze expression of coding genes
APA is frequently involved in gene expression regulation. To compare the gene expression versus APA in different samples. Our package also provide a simple function to assess the expression profile using RNA-seq samples.

## Build reference CDS regions using 
```{r eval=FALSE}
library("GenomicFeatures")
library("org.Mm.eg.db")
txdb <- makeTxDbFromUCSC(genome="mm9", tablename="refGene")
seqlevels(txdb) <-"chr19"
IDDB <- org.Mm.eg.db
CDSdbraw=REFCDS(txdb,IDDB)
```

## Extract reads count and calculte TPM for each gene 
```{r eval=FALSE}
DFGENEraw=GENEXP_CDS(CDSdbraw, flsall, Strandtype="forward")
```
