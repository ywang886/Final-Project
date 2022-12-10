# Milestone 2

# Input the data

## (A) Pre-filtering

This is a second filtering of the data. By removing rows with low gene expression reads, I can save memory size when running `dds`. I will use only the data in a row of more than ten reads.

```
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

## (B) Note on factor level

In my project, there are two groups that I want to compare. The data from the YOUNG group is contrasted with the OLD group's. Therefore, I set two factor levels here to tell DESeq2 which reference level I want to compare against.

`dds$condition <- factor(dds$condition, levels = c("untreated","treated"))`

By entering this code I can remove the levels that do not have samples in the current _DESeqDataSet_.

`dds$condition <- droplevels(dds$condition)`

# Differential expression analysis

Here I created a function `DESeq` and gave an argument `dds` to it (Love, Huber, and Anders 2014). Additionally, by using the function `results`, the result tables can be generated with log2 fold change, p values and adjusted p values.

```
dds <- DESeq(dds)
res <- results(dds)
res
```

![](https://github.com/ywang886/Pictures/blob/main/Result_table.png?raw=true)


Then, set the contrast for the result table.

```
res <- results(dds, name="condition_YOUNG_vs_OLD")
res <- results(dds, contrast=c("condition","YOUNG","OLD"))
```

## (A) Visualization and ranking

To shrink the effect size, I then used the function `resultNames` here.

```
resLFC <- lfcShrink(dds, coef="condition_YOUNG_vs_OLD", type="apeglm")
resLFC
```

![](https://github.com/ywang886/Pictures/blob/main/Log_fold_change_shrinkage.png?raw=true)

Load the BiocParallel package to speed up the analysis.

```
library("BiocParallel")
register(MulticoreParam(4))
```

## (B) Get adjusted p-values

Obtain the adjusted p-value from the result table. Adjusted p-value is a statistical method used to correct for false positives.

```
res05 <- results(dds, alpha=0.05)
summary(res05)
```

![](https://github.com/ywang886/Pictures/blob/main/adjusted_p-value.png?raw=true)

There were 43 p-values less than 0.05.

`sum(res05$padj < 0.05, na.rm=TRUE)`

[1] 43

# Exploring and exporting results

## (A) MA-plot

Make a MA plot with the difference of log fold along y-axis and normalized mean count expression along x-axis to visualize the gene differential expression.

`plotMA(res, ylim=c(-2,2))`

![](https://github.com/ywang886/Pictures/blob/main/MAplot1.png?raw=true)

According to the plot below, the blue points represent the genes that have high differential expression levels, meaning that the result is acceptable.

`plotMA(resLFC, ylim=c(-2,2))`

![](https://github.com/ywang886/Pictures/blob/main/MAplot2.png?raw=true)

Use different estimators to shrink the datasets and visualize the resulting data.

```
resultsNames(dds)

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

![](https://github.com/ywang886/Pictures/blob/main/apeglm_normal_ashr_MAplot.png?raw=true)

## (B) Plot counts

The most and least differentially expressed genes can be obtained from the following codes.

```
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene=which.max(res$padj), intgroup="condition")
```

![](https://github.com/ywang886/Plots/blob/main/plot_counts_min_max.png?raw=true)

The least differentially expressed gene is ENSG00000109321 while the most differentially expressed gene is ENSG00000109321.

Since I filter the gene of TCGA-BRCA at first, making plot counts helps me to examine the counts of the single genes', BRCA1 and BRCA2, reads.

The plots show no significant differences in gene expression between two age groups.

```
d <- plotCounts(dds, "ENSG00000012048.23", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("BRCA1")
```

![](https://github.com/ywang886/Pictures/blob/main/plot_count_BRCA1.png?raw=true)

```
d <- plotCounts(dds, "ENSG00000139618.16", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) +
  ggtitle("BRCA2")
```

![](https://github.com/ywang886/Pictures/blob/main/plot_count_BRCA2.png?raw=true)

MORE INFORMATION ON RESULTS COLUMNS

`mcols(res)$description`

![](https://github.com/ywang886/Pictures/blob/main/more_info.png?raw=true)

# Exporting results to CSV files

Only if the threshold is within a 0.05 value of the p-value results will be written by _write.csv_ function.

```
write.csv(as.data.frame(resOrdered), 
          file="condition_age_results.csv")
resSig_0.05 <- subset(resOrdered, padj < 0.05)
resSig_0.05[which(resSig_0.05$log2FoldChange > 0), "up_down"] <- "up"
resSig_0.05[which(resSig_0.05$log2FoldChange < 0), "up_down"] <- "down"
resSig_up <- subset(resSig_0.05, log2FoldChange > 0)
resSig_down <- subset(resSig_0.05, log2FoldChange < 0)
resSig_0.05
```

![](https://github.com/ywang886/Pictures/blob/main/csv_file_writing.png?raw=true)

# Data transformations and visualization

## (A) Count data transformation

There are two functions, `vst` and `rlog`, that eliminate the variance of the log fold-change when the mean of the count data is very low. Here I used the function `vst` instead of `rlog` to extract transformed values because I had many samples in my dataset. It took less time to run the count transformation using `vst`.

`vsd <- vst(dds, blind=FALSE)`

![](https://github.com/ywang886/Pictures/blob/main/extract_values.png?raw=true)

## (B) Effects of transformations on the variance

By typing the code, I can make the plot of standard deviation of the transformed data. As can be seen from the plot, the variance has been eliminated where the mean is low, making the dataset more concentrated.

```
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```

![](https://github.com/ywang886/Pictures/blob/main/transformed_plot_ntd.png?raw=true)

`meanSdPlot(assay(vsd))`

![](https://github.com/ywang886/Pictures/blob/main/transformed_plot_vsd.png?raw=true)

# Heat Map

The heat map helps to compare the X-axis and Y-axis. As shown, the blue region in the middle has different gene expression in the two age groups. This gene name can be found in the "ntd" table and can be known to have lower expression in younger patients.

```
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

![](https://github.com/ywang886/Pictures/blob/main/heatmap_ntd.png?raw=true)

```
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

![](https://github.com/ywang886/Pictures/blob/main/heatmap_vsd.png?raw=true)

# Heatmap of the sample-to-sample distances

This heat map plot illustrates the correlation of gene expression between the two groups, which shows some differences in gene expression between younger and older patients.

```
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](https://github.com/ywang886/Pictures/blob/main/heatmap_sampledistance.png?raw=true)

# Principal component plot of the samples

`plotPCA(vsd, intgroup="condition")`

![](https://github.com/ywang886/Pictures/blob/main/plotPCA.png?raw=true)

# Deliverable

### R MarkDown

# References

1. Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinfo

2. Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. 10.1093/biostatistics/kxw041