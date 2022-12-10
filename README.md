# TRGN510 Final Project

# Title

## Differential Gene Expression in Breast Invasive Ductal Carcinoma by Age Using DeSEQ2.


# Author

## Yu Hsuan (Annie) Wang


# Overview of project

### I will identify differentially expressed genes between older and younger breast cancer patients. This analysis will utilize the package DeSEQ2 (http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and complete the entire vignette. For this analysis, I’ll utilize the TCGA cohort, and have identified 385 STAR counts files for tumors that fit within my cohort with people under 50 years of age and those over 50 years of age, according to the CDC website.

# Data

### I will use the data from [https://portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository). Examining clinical data, there are 695 breast ductal carcinoma samples, and 496 are individuals greater than or equal to 50 years of age, and 198 are identified as patients under 50. The specific files are available are [here](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22set_id%3AegaBY4QBHCyNF9U_zGgg%22%5D%7D%2C%22op%22%3A%22IN%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22transcriptome%20profiling%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22tsv%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=cases).


# Milestone 1

# Screening files

The topic I would like to study is the differential expression of BRCA gene in infiltrating ductal breast cancer at different age levels. The related data will be downloaded from the [GDC Data Portal](https://portal.gdc.cancer.gov/). First, I go to the "Explore" page and set filters in the "Cases" and "Clinical" sections as shown below.

1. Primary Site: breast

2. Program: TCGA

3. Project: TCGA-BRCA

4. Sample type: Primary tumor

5. Gender: female

6. Race: white

7. Ethnicity: not hispanic or latino

8. Vital Status: alive

9. Primary Diagnosis: infiltrating duct carcinoma, nos

10. Age at Diagnosis: (1) 30y-39y (2) 70y-79y

According to the [American Cancer Society](https://www.cancer.org/cancer/breast-cancer/about/how-common-is-breast-cancer.html#:~:text=Breast%20cancer%20mainly%20occurs%20in,cancer%20are%20younger%20than%2045.), the median age of breast cancer diagnosis in women is 62 years old. Because I am interested in knowing if age has an effect on gene expression, I decide to compare two groups that are below and above this age.

![](https://github.com/ywang886/Plots/blob/main/exploration.png?raw=true)

After setting these filters in Exploration, I then click the bottom "View Files in Repository" to set the data type I will need to use.

1. Data Category: transcriptome profiling

2. Data type: Gene expression quantification

3. Experimental Strategy: RNA-Seq

4. Workflow Type: STAR-counts

5. Data format: tsv

6. Access: open

![](https://github.com/ywang886/Plots/blob/main/repository.png?raw=true)

The total files I have gotten are (1) 25 files / 22 cases, and (2) 36 files / 33 cases.

# Data Ddownloads

The data I will need to use are _manifest.txt_ files and _clinical tsv_ files.

# Organizing counts files

## (A) Installation of GDC Data Transfer Tool

Since the genome data stored at GDC Data Portal might be large sizes to arrange, a high-performance data download tool is useful in this case. I will use GDC Data Transfer Tool Client to execute data downloads and submissions.

1. First, go to [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) download and install the GDC-client tool in MAC OSX x64 version.

2.  Type the command in the terminal `export PATH="directory path:$PATH"`

3. Check if the tool is installed succefully by typing `gdc-client -version` 

## (B) Download the Counts files

1. Manifest files were downloaded from GDC Data Protal. I create two directories, BCyoung_data and BCold_data to save each manifest file.

2. By using GDC-client tool, I can read the manifest.txt file through typing `gdc-client download -m gdc_manifest.young.txt` in the terminal.

3. Each manifest.txt file can generate two kinds of files at the same time, the target tsv file and the log parcel file.

## (C) Organize the Counts files

Set the working directory to the current directory.

`setwd("/Users/Directory path")`

Because the downloaded file contains a log file within the folder, the target file should be extracted to the "BCyoung_data" and "BDold_data" folders that I create.

List all the files in the working directory.

```
i <- list.dirs()

i
```

Move the target files into the folder.

```
m = i[2:51]
for(n in m){
  x.path=paste(n,list.files(n),sep='/')
  file.copy(x.path,'./data',recursive = T)}
```

Now I have 22 tsv files in BCyoung_data tsv folder and 33 files in BCold folder, these data need to be merged in one table. First, I change the file names to their case ID names using clinical information. Second, I extract the gene ID and gene data type, HTseq-fpkm, from the _rna_seq.augmented_star_gene_counts.tsv_ files and manually put them into a table as _gene_counts2.csv_. The csv file shows as below.

![](https://github.com/ywang886/Plots/blob/main/gene_counts.png?raw=true)

Since my pasilla package cannot be installed successfully, I created a _sample_annotation.csv_ file containing the patients' conditions by myself. The csv file shows as below.

![](https://github.com/ywang886/Plots/blob/main/sample_annotation.png?raw=true)

# Inputting data to the vignette

## (A) Install packages

1. Install DESeq2 package in [Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

2. Install ggplot2 package by typing `intall.packages(ggplot2)`

3. Install all the packages below and prepare to plot the results.

	3.1 apeglm
	
	3.2 BiocParallel
	
	3.3 ashr
	
	3.4 vsn
	
	3.5 pheatmap
	
	3.6 BioManager

## (B) Input and screen the data

Input the data to R studio.

`cts <- read.delim("gene_counts2.csv", header = TRUE, row.names = 1, sep = ",")`

There are 60660 rows in cts. Now I have to do is to remove the empty gene expression data in the table.

`cts <- counts[which(rowSums(counts) > 0),]`

After removing the empty rows, the total row is 53906.

Then inputing the annotation file to set the condition for each case.

`coldata <- read.delim("sample_annotation2.csv", header = TRUE, row.names = 1, sep = ",")
`

Converting the data into integer mode because _DESeqDataSetFromMatrix_ needs integered data to work.

`cts <- round(cts)`

![](https://github.com/ywang886/Plots/blob/main/cts.png?raw=true)

Setting the condition for the matrix.

```
condition <- factor(c("OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD","OLD",
"YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG","YOUNG"))
```

Now I can run _DESeqDataSetFromMatrix_ to organize my data by making a matrix for plots.

`dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)`


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

Since I filtered the gene of TCGA-BRCA at first, making plot counts helps me to examine the counts of the single genes', BRCA1 and BRCA2, reads.

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

# Known issues

Problem: An error was encountered while installing the pasilla package as well as the PNG package was not working in R studio.

Solution: Follow the instruction in the vignette and learn by searching for similar projects online to create the same files manually.

```
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
```


# Conclusions

I downloaded 25 files for breast cancer cases aged 30-39 years and 36 files for cases aged 70-79 years. These files were then organized into tsv and csv tables and analyzed using DESeqDataSetFromMatrix. After a second screening for low gene expression and the shrinkage of the dataset of adjusted p-value < 0.05, the data were ready for plotting. Different plots require different levels of calculation and illustrate some significance of differential gene expressions. Gene expression in younger and older patients did not show significant differences in single plot counts and transformed data plots. However, MA plots and heat maps did show some differences between young and old in breast cancer.

In the process of gene expression analysis, I learned to become more proficient in applying R studio and improved my problem-solving skills. In my project, although age has some effect on gene expression differences, I suppose there are other gene contrasts that can show more genetic differences in breast cancer patients, such as ER, PR ,and HER2. Moreover, I have learned different perspectives on analyzing gene expression based on the production of plots. By plotting, the results can be visualized more effectively and I will continue to learn to improve my skills.