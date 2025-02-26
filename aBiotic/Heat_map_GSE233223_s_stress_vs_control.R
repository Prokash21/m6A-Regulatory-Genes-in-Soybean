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
setwd("D:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL")

#load the count data

count_data <- read.csv("Count_Data_final.csv", header=TRUE,row.names = 1)

colnames(count_data)
head(count_data)

#load the sample info
#sample_info1 <- read.csv("meta_keratinocyte.csv")

sample_info <- read.csv("Meta_Data_final.csv", header = TRUE,row.names = 1)

#colData1 <- read.csv("meta_keratinocyte.csv", header = T, sep = '\t', 
#                     stringsAsFactors = TRUE)
#colData <- read.csv("meta_keratinocyte.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$agent <- factor(sample_info$agent)
sample_info$Time_point <- factor(sample_info$Time_point)
sample_info$tissue <- factor(sample_info$tissue)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)


unique(sample_info$agent)
unique(sample_info$Time_point)
unique(sample_info$tissue)


dim(count_data)
dim(new_count_data)
dim(sample_info)
design <- ~ agent + Time_point + agent:Time_point

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ agent )

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
dds$Treatment <- factor(dds$agent, levels = c ("control","Salt stress")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$agent <- relevel(dds$agent , ref = "control")
dds$agent

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

#[1] "Intercept"                    "agent_Salt.stress_vs_control"

resLFC <- lfcShrink(dds, coef ="agent_Salt.stress_vs_control"  , type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

##################18th oct


resLFC_p_cut <- as.data.frame(resLFC[ resLFC$padj<0.05,])
write.csv(resLFC_p_cut, file = "resLFC_p_cut.csv")
#################
# Create a Volcano Plot
ggplot(resLFC,aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels and add the header
  labs(
    title = "HeatMap of agent_Salt.stress_vs_control", # Add the header here
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

########################################## padj #############################################

Upregulated <- resLFC[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05, ]

top_UP_hits1 <- Upregulated[order(Upregulated$padj ),]
top_UP_hits <- Upregulated[order(Upregulated$padj),][1:20,]
top_UP_hits1 <-row.names(top_UP_hits1)
top_UP_hits <-row.names(top_UP_hits)
top_UP_hits

write.csv(top_UP_hits, file = "UP_top_20_agent_Salt.stress_vs_control.csv", row.names = FALSE)

Downregulated <- resLFC[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05, ]
top_hits1 <- Downregulated[order(Downregulated$padj),]
top_hits <- Downregulated[order(Downregulated$padj),][1:20,]
top_hits1 <-row.names(top_hits1)
top_hits <-row.names(top_hits)
top_hits

write.csv(top_hits, file = "DOWN_top_20_agent_Salt.stress_vs_control.csv", row.names = FALSE)




#Hratmap


# Calculate Z-scores for the top 20 genes
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

zscore_all <- t(apply(normalize_counts, 1, cal_z_score))
XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
top <- as.character(top)
zscore_subset <- zscore_all[top, ]

# Create a heatmap of Z-scores with gene names as row names
rownames(zscore_subset) <- top
colnames(zscore_subset) <- sample_info$Time_point
pheatmap(zscore_subset,
         main = "Heatmap with gene Ids")

write.csv(zscore_subset,file = "Z_Score_Afser_vai_heatmap.csv")



