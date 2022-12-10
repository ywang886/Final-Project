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

