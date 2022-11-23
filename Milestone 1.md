# TRGN510 Final Project

# Title

## Differential Gene Expression in Breast Invasive Ductal Carcinoma by Age Using DeSEQ2.


# Author

## Yu Hsuan (Annie) Wang


## Milestone 1

## (A) Screening files

### I will obtain my data from [CDC Data Protal](https://portal.gdc.cancer.gov/repository) by setting several filters, such as cancer primary diagnosis, gender, race, data category, workflow type, etc.

First, I went to "Exploration" page and set the filters in "Cases" and "Clinical" sections as below.

1. Primary Site: breast

2. Program: TCGA

3. Gender: female

4. Race: white

5. Ethnicity: not hispanic or latino

6. Vital Status: alive

7. Primary Diagnosis: infiltrating duct carcinoma, nos

![](https://github.com/ywang886/Pictures/blob/main/Filter%202.png?raw=true)

After setting these filters in Exploration, I then clicked the bottom "View Files in Repository" to set the data type I would need to use.

1. Data Category: transcriptome profiling

2. Experimental Strategy: RNA-Seq

3. Workflow Type: STAR-counts

4. Access: open

![](https://github.com/ywang886/Pictures/blob/main/Filter%201.png?raw=true)

The total files I have gotten is 385. I would download and arrange the JSON file for data analysis later.

## (B) Organizing counts files

### Since the file downloaded from the CDC data portal is a STAR-count file, I will need to convert it to a HTseq-count file to facilitate the subsequent use of DEseq to run the data.

## (C) Inputting data to the vignette

### I will input the HTseq-count file and contruct the *DESeqDataSet*.

# Milestone 2

## (A) Differential expression analysis

### I will generate results table by using data frames for operating differential expression results.

## (B) Exploring and exporting results

### MA-plot

### Plot counts

### More information on results columns

### Exporting results to CSV files

## (C) Data transformations and visualization

### I will build discrete distributions by using the raw counts to test for differential expression.


# Deliverable

### R MarkDown
