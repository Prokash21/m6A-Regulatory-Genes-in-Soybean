#
#
#
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
getwd()
#set the working directory
setwd("D:/DWCT/DOT_apu_project/aBiotic/GSE57252")

#load the count data

count_data <- read.csv("GSE57252_DESeq_normalized_counts.csv", header=TRUE,row.names = 1)

colnames(count_data)
head(count_data)

#load the sample info
#sample_info1 <- read.csv("meta_keratinocyte.csv")

sample_info <- read.csv("meta_data1.csv", header = TRUE,row.names = 1)

#colData1 <- read.csv("meta_keratinocyte.csv", header = T, sep = '\t', 
#                     stringsAsFactors = TRUE)
#colData <- read.csv("meta_keratinocyte.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$treatment <- factor(sample_info$treatment)
#sample_info$Time_point <- factor(sample_info$Time_point)
#sample_info$tissue <- factor(sample_info$tissue)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)


unique(sample_info$treatment)
#unique(sample_info$Time_point)
#unique(sample_info$tissue)


dim(count_data)
dim(new_count_data)
dim(sample_info)
#design <- ~ agent + Time_point + agent:Time_point

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ treatment )

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
#dds$treatment <- factor(dds$treatment, levels = c ("Control","Water deficit","Heat stress","Combined water deficit and heat stress")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$treatment <- relevel(dds$treatment , ref = "Co")
dds$treatment

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"Salt stress_vs_control.csv")



#Identify available coefficient names
coeff_names <- resultsNames(dds)

#Print the coefficient names
print(coeff_names)

#[1] "treatment_De_vs_Co" "treatment_Na_vs_Co"                                                   
#[2] "treatment_Combined.water.deficit.and.heat.stress._vs_Control."
#[3] "treatment_Heat.stress._vs_Control."                           
#[4] "treatment_Water.deficit._vs_Control."                         

######################################## "Time_point_12h_vs_0h"  ################################################

resLFC <- lfcShrink(dds, coef ="treatment_De_vs_Co", type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
head(resLFC)

XX<-read.csv("Dot_aPU_45_gene.csv")
top <- XX$geneid
top <- as.character(top)

Y<- resLFC[top, ]

write.csv(Y, file = "treatment_De_vs_Co.csv")


######################################## "Time_point_24h_vs_0h" ################################################
resLFC <- lfcShrink(dds, coef ="treatment_Na_vs_Co"  , type = "apeglm")
resLFC<- as.data.frame(resLFC)
#change resLFC to a dataframe
resLFC1 <- as.data.frame(resLFC)

XX<-read.csv("Dot_aPU_45_gene.csv")
top <- XX$geneid
print(top)
top <- as.character(top)
Y<- resLFC1[top, ]

write.csv(Y, file = "treatment_Na_vs_Co.csv")


######################################## "Time_point_48h_vs_0h" ################################################
# Load necessary libraries if not already loaded
library(DESeq2)

# Perform shrinkage on log2 fold changes
resLFC <- lfcShrink(dds, coef = "treatment_Heat.stress_vs_Control", type = "apeglm")

# Convert resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

# Load gene IDs to select top genes
XX <- read.csv("Dot_aPU_45_gene.csv")
top <- as.character(XX$Gene_id)

# Filter resLFC for top genes
filtered_resLFC <- resLFC[top, ]

# Write filtered results to CSV
write.csv(filtered_resLFC, file = "treatment_Heat.stress_vs_Control.csv", row.names = FALSE)
