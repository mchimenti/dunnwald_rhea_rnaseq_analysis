## Analysis of Rhea/Dunnwald RNA-Seq Human hip/lip/oral syndromic / non-syndromic
## Date: 6.1.2018
## Author: Michael Chimenti
## Organism: hg38 / human
## Aligners: hisat2 / salmon
## Design: Case/control (syndromic vs. non) single timepoint, age batch, tissue-type batch 
## Reps: 4

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")

#negative binomial GLM
library(DESeq2)
library(calibrate)

#annotation
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library(dplyr)
library(pcaExplorer)
#pathway
library(pathview)
library(gage)
library(gageData)
library(ggplot2)

setwd("~/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/analysis") 

#######################################
## FeatureCounts > DESeq2 > PCAExplorer
#######################################


## Annotation function 
get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

## setup paths and create count matrix and metadata table; create DESeq object

matrix_path <- "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/dunnwall_rhea_annotated_comb.counts"
count_mat <- read.table(matrix_path, header=TRUE, row.names = 'id')

annot <- count_mat[28]
count_mat <- count_mat[,1:27]
count_mat <- count_mat[,c(6,1:5,7:27)]  ## reorder columns to match coldata order 
count_mat <- as.matrix(count_mat)
#count_mat <- count_mat[,c(2:5,7:24)]


colData <- read.csv("/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/sample_metadata.csv")
colData <- arrange(colData, sname)
colData$age <- as.factor(colData$age)

## hip and age are confounded, deseq2 won't let me make a dds with this
## solution: drop age column, and drop hip tissue



# colData <- filter(colData, tissue != 'hip')
# colData <- select(colData, c("sname","tissue","syndromic"))
# 
# dds <- DESeqDataSetFromMatrix(countData = count_mat,
#                               colData = colData,
#                               design = ~ tissue + syndromic)
# 
# ## prefiltering and re-leveling
# 
# dds <- dds[ rowSums(counts(dds)) > 5, ]
# dds$syndromic <- relevel(dds$syndromic, ref='n')
# 
# ## Differential expression analysis with contrasts
# 
# dds <- DESeq(dds)
# 
# ##---------------launch PCA Explorer on dds object 
# anno <- get_annotation(dds, 'hsapiens_gene_ensembl','ensembl_gene_id')
# anno <- na.omit(anno)
# rld <- rlog(dds, blind=FALSE)
# pcaExplorer(dds=dds,annotation=anno,rlt=rld)

## 13A and 13B are outliers, drop and reanalyze

colData <- filter(colData, !(sname %in% c("X13A_20180513000","X13B_20180513000")))
count_mat <- count_mat[,-c(2,3)]

colData <- colData %>% mutate(group = as.factor(paste0(colData$tissue, colData$syndromic)))

dds_drop <- DESeqDataSetFromMatrix(countData = count_mat,
                                   colData = colData,
                                   design = ~ group)

dds_drop <- dds_drop[ rowSums(counts(dds_drop)) > 5, ]
dds_drop$group <- relevel(dds_drop$group, ref='lipn')

dds_drop <- DESeq(dds_drop)

##---------------launch PCA Explorer on dds object 
anno <- get_annotation(dds_drop, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rld_drop <- rlog(dds_drop, blind=FALSE)
pcaExplorer(dds=dds_drop,annotation=anno,rlt=rld_drop)





## results 

res_lip_syn <- results(dds_drop, contrast=c("group", "lipy", "lipn"))
res_oral_syn <- results(dds_drop, contrast=c("group", "oraly", "oraln"))
res_hip_syn <- results(dds_drop, contrast=c("group", "hipy", "hipn"))

sig <- 0.2

## extract syndromic comparisons 
res_lip_syn <- na.omit(res_lip_syn)  #drop NA rows
res_lip_syn_sig <- res_lip_syn[res_lip_syn$padj < sig & res_lip_syn$baseMean > 5.0,]
res_lip_syn_ord <- res_lip_syn_sig[order(res_lip_syn_sig$padj),]

res_oral_syn <- na.omit(res_oral_syn)  #drop NA rows
res_oral_syn_sig <- res_oral_syn[res_oral_syn$padj < sig & res_oral_syn$baseMean > 5.0,]
res_oral_syn_ord <- res_oral_syn_sig[order(res_oral_syn_sig$padj),]

res_hip_syn <- na.omit(res_hip_syn)
res_hip_syn_sig <- res_hip_syn[res_hip_syn$padj < sig & res_hip_syn$baseMean > 5.0,]
res_hip_syn_ord <- res_hip_syn_sig[order(res_hip_syn_sig$padj),]



# ## adaptive shrinkage 
# 
# res_lip_syn_shrink <- lfcShrink(dds_drop, contrast=c("group", "lipy", "lipn"), type = 'ashr')
# res_lip_syn_shrink <- na.omit(res_lip_syn_shrink)  #drop NA rows
# res_lip_syn_shrink_sig <- res_lip_syn_shrink[res_lip_syn_shrink$padj < sig & res_lip_syn_shrink$baseMean > 5.0,]
# res_lip_syn_shrink_ord <- res_lip_syn_shrink_sig[order(res_lip_syn_shrink_sig$padj),]

## Independent hypothesis weighting 

# library("IHW")
# res_lip_syn_IHW <- results(dds_drop, contrast=c("group", "lipy","lipn"), filterFun = ihw)
# res_lip_syn_IHW <- na.omit(res_lip_syn_IHW)  #drop NA rows
# res_lip_syn_IHW_sig <- res_lip_syn_IHW[res_lip_syn_IHW$padj < sig & res_lip_syn_IHW$baseMean > 5.0,]
# res_lip_syn_IHW_ord <- res_lip_syn_IHW_sig[order(res_lip_syn_IHW_sig$padj),]
# 
# res_lip_syn_IHW_ord$gene <- anno[row.names(res_lip_syn_IHW_ord),"gene_name"]
# 
# res_oral_syn_IHW <- results(dds_drop, contrast=c("group", "oraly","oraln"), filterFun = ihw)
# res_oral_syn_IHW <- na.omit(res_oral_syn_IHW)  #drop NA rows
# res_oral_syn_IHW_sig <- res_oral_syn_IHW[res_oral_syn_IHW$padj < sig & res_oral_syn_IHW$baseMean > 5.0,]
# res_oral_syn_IHW_ord <- res_oral_syn_IHW_sig[order(res_oral_syn_IHW_sig$padj),]
# 
# res_oral_syn_IHW_ord$gene <- anno[row.names(res_oral_syn_IHW_ord), "gene_name"]

## extract tissue comparisons

sig = 0.01
res_oral_lip_nosyn <- results(dds_drop, contrast=c("group","oraln","lipn"))
res_oral_lip_nosyn <- na.omit(res_oral_lip_nosyn)
res_oral_lip_nosyn_sig <- res_oral_lip_nosyn[res_oral_lip_nosyn$padj < sig & res_oral_lip_nosyn$baseMean > 5.0,]
res_oral_lip_nosyn_ord <- res_oral_lip_nosyn_sig[order(res_oral_lip_nosyn_sig$padj), ]

res_oral_lip_nosyn_ord$ext_gene <- anno[row.names(res_oral_lip_nosyn_ord), "gene_name"]

write.csv(res_oral_lip_nosyn_ord, 
          file = "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_lip_NonSyn_padj_0p05.csv")

res_oral_lip_syn <- results(dds_drop, contrast=c("group","oraly","lipy"))
res_oral_lip_syn <- na.omit(res_oral_lip_syn)
res_oral_lip_syn_sig <- res_oral_lip_syn[res_oral_lip_syn$padj < sig & res_oral_lip_syn$baseMean > 5.0,]
res_oral_lip_syn_ord <- res_oral_lip_syn_sig[order(res_oral_lip_syn_sig$padj), ]

res_oral_lip_syn_ord$gene <- anno[row.names(res_oral_lip_syn_ord), "gene_name"]

write.csv(res_oral_lip_syn_ord,
          file = "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_lip_Syn_padj_0p05.csv")

## syndromic DE genes plots

#top lip syndromic
plotCounts(dds_drop, gene = "ENSG00000188856", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000262406", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000237506", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000234449", intgroup = 'group')

#top oral syndromic
plotCounts(dds_drop, gene = "ENSG00000171199", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000205592", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000162078", intgroup = 'group')
plotCounts(dds_drop, gene = "ENSG00000188856", intgroup = 'group')


## plotting

plotMA(dds_drop)

plotPCA(rld_drop, intgroup = "group")


## volcano plot (based on Stephen Turner code)
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

## syndromic DE genes in Oral 
res_oral_syn_ord_counts <- merge(as.data.frame(res_oral_syn_ord), as.data.frame(counts(dds_drop, normalized=TRUE)), 
                                 by="row.names", sort=FALSE)
names(res_oral_syn_ord_counts)[1] <- "Gene"
res_oral_syn_ord_counts$ext_gene <- anno[res_oral_syn_ord_counts$Gene, "gene_name"]

png("oral_syn_diffexpr-volcanoplot.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_oral_syn_ord_counts, main = "Volcano Plot: Syndromic DE Genes in Oral Tissue", lfcthresh=1.0, sigthresh=0.1, textcx=.5, xlim=c(-12, 7), ylim = c(3,6))
dev.off()

## Tissue DE genes, non-syndromic 

png("oral_lip_nosyn_diffexpr-volcanoplot.png", 1200, 1500, pointsize=15, res=150)
volcanoplot(res_oral_lip_nosyn_ord, main = "Volcano Plot: DE genes Oral vs. Lip, NonSyndromic", lfcthresh=2.5, sigthresh=0.001, textcx=.25, xlim=c(-12, 12), ylim = c(4,50))
dev.off()


## write gene lists 
write.csv(res_lip_syn_IHW_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_lip_syn_IHW_padj_0p1.csv")
write.csv(res_oral_syn_IHW_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_syn_IHW_padj_0p1.csv")
write.csv(res_oral_lip_syn_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_lip_syn_padj_0p05.csv")
write.csv(res_oral_lip_nosyn_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_lip_nosyn_padj_0p05.csv")
write.csv(res_hip_syn_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_hip_syn_padj_0p05.csv")
write.csv(res_oral_syn_ord, file = 
            "/Users/mchimenti/iihg/RNA_seq/dunnwald_lab/project_rhea_jun2018/res_oral_syn_padj_0p2.csv")
