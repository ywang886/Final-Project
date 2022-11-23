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

## (C) Input data to the vignette

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
