# Milestone 1

## (A) Screening files

First, I went to "Exploration" page and set the filters in "Cases" and "Clinical" sections as below.

1. Primary Site: breast

2. Program: TCGA

3. Gender: female

4. Race: white

5. Ethnicity: not hispanic or latino

6. Vital Status: alive

7. Primary Diagnosis: infiltrating duct carcinoma, nos

8. Age at Diagnosis: (1) 26y-50y (2) 51y-90y 

![](https://github.com/ywang886/Pictures/blob/main/Filter%202.png?raw=true)

After setting these filters in Exploration, I then clicked the bottom "View Files in Repository" to set the data type I would need to use.

1. Data Category: transcriptome profiling

2. Experimental Strategy: RNA-Seq

3. Workflow Type: STAR-counts

4. Access: open

![](https://github.com/ywang886/Pictures/blob/main/Filter%201.png?raw=true)

The total files I have gotten are (1) 110 files / 98 cases, and (2) 232 files / 213 cases. I would download and arrange the files for data analysis later.

## (B) Organizing counts files

### Installation of GDC Data Transfer Tool

Since the genome data stored at GDC Data Portal might be large sizes to arrange, a high-performance data download tool is useful in this case. I would use GDC Data Transfer Tool Client to execute data downloads and submissions.

1. First, go to [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) download and install the GDC-client tool in MAC OSX x64 version.

2.  Type the command in the terminal `export PATH="directory path:$PATH"`

3. Check if the tool is installed succefully by typing `gdc-client -version` 

### Download the Counts files

1. The Manifest files were downloaded from GDC Data Protal. I created two directories, BC_younger and BC_older to save each manifest file.

2. By using GDC-client tool, I need to download the Manifest file through typing `gdc-client download -m gdc_manifest.2022-11-15.txt` in the terminal.

3. Then I organized the downloaded data into two directories, gdc_younger_data and gdc_older_data.

### Organize the Counts files

The total counts files I downloaded are 342. Each folder contains a .tsv count file and a logs folder. Since the .tsv file is needed only, I use R studio to collect all the count files.

Here, I will take 3 files from each group 26y-50y and group 51y-90y to test running the vignette.

I will select cases in groups 20y-30y and 80y-90y to observe the differential gene expression under the extreme situation.

## (C) Inputting data to the vignette

### Install DEseq2 from Bioconductor

Install the package in [Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

### Input the data

I confronted error when typing `dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)`. The error is supposed to be the different type of count files. The STAR count files may not be analyzed here since the command is for HTseq count files.

![](https://github.com/ywang886/Pictures/blob/main/R%20TEST.png?raw=true)