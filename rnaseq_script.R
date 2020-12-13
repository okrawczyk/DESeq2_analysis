# Load DESeq2 library
library("DESeq2")
library("ggplot2")

# Prefix for outfile name
outfile = 'HVJFNBGXG'
directory <- "~/HVJFNBGXG/HTSeq/COUNTS"
setwd(directory)

# Tabel -> SampleID, CountsFileName, sampleName, Condition
sampleTable <-read.csv('tabel.csv')


# Running DESeq2
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ Condition)
colData(ddsHTSeq)$Condition <- factor(colData(ddsHTSeq)$Condition, levels =  c('A', 'B'))
dds <- DESeq(ddsHTSeq)
#res <- results(dds)

# Total number of raw counts per sample
rawCountsPerSample <- colSums(counts(dds))

# Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

# Plot dispersion estimates
plotDispEsts(dds)

contrast <- c("condition", "ID", "IB")
res <- results(dds, alpha=0.05, contrast=c("Condition", "A", "B"))
resultsNames(dds)

library("apeglm")
resLFC <- lfcShrink(dds, coef="Condition_A_vs_B", type="apeglm")

# MA plot 
plotMA(res, ylim=c(-2,2), cex=1, colSig='red')
plotMA(resLFC, ylim=c(-2,2), cex=1, colSig='red')

# Identify genes on MA plot by pressing mouse button
idx <- identify(resLFC$baseMean, resLFC$log2FoldChange)
rownames(resLFC)[idx]

## Get differential expression results

# Order by the smallest p-value & output significant results
# threshold = 0.05
table(res$padj < 0.05)
#Order by adjusted p-value
resOrdered <- res[order(res$padj, decreasing = FALSE), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

#Save only names of DE_genes
DE_genes=subset(resdata, resdata$padj < 0.05)
Genes <- DE_genes$Gene
Genes <- as.data.frame(Genes)
write.csv(Genes, file="diffexpr-genes.csv")

# Examine plot of p-values
#hist(resLFC$pvalue, breaks=50, col='grey')

# Gene with the smalles p-value from resdata table
plotCounts <- plotCounts(dds, gene=which.min(resLFC$padj), intgroup="Condition", returnData = TRUE)
ggplot(plotCounts, aes(x=Condition, y=count)) +
       geom_point(position = position_jitter(w=0.1, h=0)) +
       scale_y_log10(breaks = c(25,100,400))


# Transform raw counts into normalized values:
# rlog transformed  
rld <- rlogTransformation(dds, blind=T)

# variance stabilization - for heatmaps, etc.
vsd <- varianceStabilizingTransformation(dds, blind=T)


# Exploratory data analysis - (PCA & Hierarchical Clustering) - identifying sources of variation in the data

# PCA 
library("genefilter")
library("grDevices")
library('RColorBrewer')
library('gplots')
library('pheatmap')

# Sample-to-sample distance 
sampleDists <- dist(t(assay(rld)))
mat <- as.matrix(sampleDists)
rownames(mat) <- colnames(mat) <- with(colData(dds), SampleID)
colors <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
pheatmap(mat, color=colors, clustering_distance_rows=sampleDists, 
         clustering_distance_cols = sampleDists)
#dev.copy(png, paste0(outfile, "-clustering.png"))
#dev.off()

# Hierarchical clustering
condition <- c('A', 'B')
data <- plotPCA(vsd, intgroup = "Condition" , returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
#png(file=paste0(outfile, "-PCA.png"), width=1200,height=1000)
(pcaplot <- ggplot(data, aes(PC1, PC2, color = Condition)) +
    geom_point(size=5) +
    ggtitle("Principal Component Analysis")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")))
#dev.copy(png, paste0(outfile, "-PCA.png"))
#dev.off()



#Heatmap of data

# 1000 top expressed genes with heatmap.2
select <- rownames(resOrdered)[1:1000]
mat1 <- assay(vsd)[select, ]

my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
#png(file=paste0(outfile, "-HEATMAP.png"), width=2000,height=1500)
heatmap.2(mat1, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6,
          main="Top Expressed Genes Heatmap")
#dev.copy(png, paste0(outfile, "-HEATMAP.png"))
#dev.off()

