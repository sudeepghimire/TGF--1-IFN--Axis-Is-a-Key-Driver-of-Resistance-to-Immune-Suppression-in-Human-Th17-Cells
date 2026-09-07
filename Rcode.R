
library(DESeq2)
library(limma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(tibble)
library(tidyverse)
library(pheatmap)

# Loading data and preprocess it 
metadata<-read.csv("metadata3samples.csv")
head(metadata)
metadata
dim(metadata)
counts<-read.csv("counts3samples.csv") 
head(counts)
str(counts)

#combinging genes with same names in the row
dim(counts)
df<-aggregate(.~genes,data=counts,FUN=sum)
df
dim(df)
row.names(df)
row.names(df)<-df$genes
head(df)
data <-df[, -1]
head(data)
str(df)
colData<-colnames(data)
colData
table(colnames(data)==metadata$sample) ##confirming the names of columns and metadata
dim(data)

# 3. Form a Deseq2 object called dds
library(DESeq2)
data
dim(data)
counts <-round(data)
dds <- DESeqDataSetFromMatrix(countData= counts,  
                              colData = metadata,
                              design = ~ com1367vothers)
dds 
colData(dds)
rownames(dds)
head(dds)

colSums(assay(dds)) #total number of reads in each sample

#filtering low expressed genes
is_expressed<-assay(dds) >=5
head(is_expressed)
hist(rowSums(is_expressed),main="Number of samples a gene is expressed in",xlab="Sample Count")


keep<-rowSums(assay(dds)>=5)>=3 ##keeping those with "TRUE" in at least 3 samples
table(keep)
keep

dds<-dds[keep,] ##filtereing some samples 
dds
boxplot(assay(dds)) #visualizing the counts distributions
boxplot(log10(assay(dds)))

# Get log2 counts
vsd <- vst(dds,blind=TRUE) 
# Check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions", col="orange")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd)), col="blue")

library(ggpubr)
plotPCA(vsd,intgroup="com1367vothers")+stat_stars()+stat_ellipse(level = 0.95,linetype=2)+ 
  theme_classic()+geom_point(size=1)+theme(axis.text.x = element_text(size = 12, color = "black"), 
                                           axis.text.y = element_text(size = 12, color="black"), 
                                           text = element_text(size=12, color = "black"))
# Get PCA data from the VST object
pcaData <- plotPCA(vsd, intgroup = "com1367vothers", returnData = TRUE)
metadata <- as.data.frame(colData(vsd))
metadata$sample <- rownames(metadata)

pcaData$sample <- pcaData$name
pcaData <- left_join(pcaData, metadata, by = "sample")
head(pcaData)

vars <- c(
  "group",
  "condition",
  "origin"
)

results <- data.frame()

for (v in vars) {
  
  # Remove NA
  dat <- pcaData[!is.na(pcaData[[v]]), ]
  
  # PC1
  model1 <- aov(PC1 ~ dat[[v]], data = dat)
  anova1 <- summary(model1)[[1]]
  
  ss1 <- anova1[1, "Sum Sq"]
  total_ss1 <- sum(anova1[, "Sum Sq"])
  
  # PC2
  model2 <- aov(PC2 ~ dat[[v]], data = dat)
  anova2 <- summary(model2)[[1]]
  
  ss2 <- anova2[1, "Sum Sq"]
  total_ss2 <- sum(anova2[, "Sum Sq"])
  
  results <- rbind(
    results,
    data.frame(
      variable = v,
      PC1_R2 = ss1 / total_ss1,
      PC1_pvalue = anova1[1, "Pr(>F)"],
      PC2_R2 = ss2 / total_ss2,
      PC2_pvalue = anova2[1, "Pr(>F)"]
    )
  )
}

results

results %>%
  dplyr::arrange(desc(PC1_R2)) %>%
  dplyr::mutate(
    PC1_R2_percent = round(PC1_R2 * 100, 2),
    PC2_R2_percent = round(PC2_R2 * 100, 2)
  )

#PCA plot showing origin of samples 
library(ggpubr)
library(ggrepel)
head(pcaData)

ggplot(pcaData, aes(PC1, PC2, color = com1367vothers.x)) +
  geom_point(size = 2) +
  geom_text_repel(
    aes(label = origin),
    size = 3,
    max.overlaps = Inf
  ) +
  stat_ellipse(level = 0.95, linetype = 2) +
  theme_classic() +
  labs(
    x = paste0("PC1: ", round(attr(pcaData, "percentVar")[1] * 100), "% variance"),
    y = paste0("PC2: ", round(attr(pcaData, "percentVar")[2] * 100), "% variance")
  )

# compute pairwise correlation values
head(assay(vsd))
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
vsd_cor
library(pheatmap)
pheatmap(vsd_cor, labels_row = vsd$condition) #showed clustering of 1367 separately

# 7. normalize read counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds
rownames(dds)
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
normlzd_dds <- counts(dds, normalized=T)
head(normlzd_dds)

##differential gene expression
design(dds) <- ~com1367vothers
rownames(dds)
dds <- DESeq(dds, test = "Wald")
res <-results(dds,contrast=c("com1367vothers","Th17_1_3_6_7",  "Th17_2_4_5_8_9"))
summary(res, alpha=0.05)
res

res[order(res$padj),]
#write.csv(as.data.frame(res[order(res$padj), ]), file="Th17_1_3_6_7 vs Th17_others2.csv")


##making an upregulated downregulated table 
library(tidyverse)
library(ggrepel)
library("kableExtra")
# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}


xyz<-as.data.frame(res)
head(xyz)
yz <- xyz %>% 
  mutate(
    Expression = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(yz) %>% 
  knitr_table()

head(yz)

yz %>% 
  count(Expression) %>% 
  knitr_table()

# Make a basic volcano plot #https://www.r-bloggers.com/2014/05/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/
yz$Gene<-rownames(yz)
head(yz)
with(yz, plot(log2FoldChange, -log10(pvalue), pch=20, main="Th17_1_3_6_7 vs Th17_others"))
# Add colored points: blue if padj<0.05, red if log2FC>2 and padj<0.05)
with(subset(yz, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(yz, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
library("calibrate")
with(subset(yz, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.4))

# 9. setting differentially expressed gene selection criteria
pThr<-0.05   ##for adj. p-value threshold
logFCThr<-3    ##for log2 fold change 
baseMeanThr<-20    ##for average expression level (at least 20 counts)

idx=which(res$padj<=pThr &
            abs(res$log2FoldChange) >=logFCThr &
 res$baseMean >=baseMeanThr)

sigRes = res[idx, ]
dim(res)
dim(sigRes) 
head(sigRes)
dim(sigRes)
rownames(sigRes)
head(sigRes)

#11. Annotating genes names for more information for both differentially expressed as well as all genes
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(dplyr))
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
sigRes

anno_sigRes <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(sigRes),
                                     columns=c("SYMBOL","GENENAME", "ENSEMBL", "ALIAS"),
                                     keytype="ALIAS")
head(anno_sigRes)
rownames(sigRes)
sigRes2<-cbind(SYMBOL=rownames(sigRes), sigRes)
head(sigRes2)
sigRes2
sigtab<-left_join(as.data.frame(sigRes2), anno_sigRes)
head(sigtab)
sigtab
dim(sigtab)##this file contains all signifiantly altering genes with more annotations
head(sigtab)


# 12.plot heatmap of DE genes
mat<-assay(vsd) 
head(mat)
sigRes
sigRes$gene<-row.names(sigRes)
head(sigRes)
idx<-sigRes$gene
head(idx)
DEgenes<-mat[idx, ]
DEgenes
dim(DEgenes)

annotation<-as.data.frame(colData(vsd)[,"com1367vothers"])
rownames(annotation)<-colnames(DEgenes)
pheatmap(DEgenes, scale="row", show_rownames = T, clustering_distance_rows = "correlation",
         annotation_col = annotation, fontsize_col = 6 ,cellwidth = 6,
         border_color = "grey75" )




