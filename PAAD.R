library("TCGAbiolinks")
library(SummarizedExperiment)
library("limma")
library("edgeR")
library("EDASeq")
library(gplots)
library(biomaRt)



# Project information
getProjectSummary("TCGA-PAAD")

# Download and preprocess data
paadQ <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification")
GDCdownload(paadQ)
paad.data <- GDCprepare(paadQ)
head(paad.data)

# Explore some metadata information
paad.data$tissue_type
paad.data$tumor_descriptor
paad.data$barcode

# Create a simple metadata for our use case
simpleMeta <- data.frame("barcode" = paad.data$barcode, 
                         "tissue_type" = paad.data$tissue_type, 
                         "tumor_type" = paad.data$tumor_descriptor)

# Select unstranded dataset
paad.raw.data <- assays(paad.data) # assays function extracts the expression data
dim(paad.raw.data$unstranded)

# Let's downsize our data to 10 primary and 10 normal tissue samples
selectedBarcodes<- c(subset(simpleMeta, tissue_type== "Tumor")$barcode[c(1:10)],
                     subset(simpleMeta, tissue_type=="Normal")$barcode[c(1:10)])

# Filter the selected barcodes to include only those present in paad.raw.data$unstranded
validBarcodes <- colnames(paad.raw.data$unstranded)
selectedBarcodes <- selectedBarcodes[selectedBarcodes %in% validBarcodes]

# Select the data for the valid barcodes
selectData <- paad.raw.data$unstranded[, selectedBarcodes]
dim(selectData)

# Data normalization and filtering
NormData <- TCGAanalyze_Normalization(tabDF = selectData, geneInfo = geneInfoHT, method = "geneLength")

# Filtering
filtData <- TCGAanalyze_Filtering(tabDF = NormData, method = "quantile", qnt.cut = 0.25)

# Display dimensions of filtered data
dim(filtData)

# Differential expressed analysis (DEA)
selectResults <- TCGAanalyze_DEA(mat1 = filtData[, selectedBarcodes[1:10]], 
                                 mat2 = filtData[, selectedBarcodes[11:14]], 
                                 Cond1type = "Normal", 
                                 Cond2type = "Tumor",
                                 pipeline = "edgeR", 
                                 fdr.cut = 0.01, 
                                 logFC.cut = 2)

# Plot the results
plot(selectResults$logFC, -log10(selectResults$FDR), 
     xlab = "log Fold Change", 
     ylab = "-log10 FDR", 
     main = "Volcano Plot of Differential Expression")

# Differential expression analysis with treatment levels
selectResults.Level <- TCGAanalyze_LevelTab(selectResults, 
                                            "Normal", 
                                            "Tumor",
                                            filtData[, selectedBarcodes[1:10]],
                                            filtData[, selectedBarcodes[11:14]])

head(selectResults.Level)

# Now we can visualize with a heatmap
head(selectResults.Level)
dim(selectResults.Level)

heat.data <- filtData[rownames(selectResults.Level),]

# Color the plot by the kind of tumor
cancer.type <- c(rep("Normal", 10), rep("Tumor", 10))

ccodes <- c()

for (i in cancer.type) {
  if (i == "Normal") {
    ccodes <- c(ccodes, "red")
  } else {
    ccodes <- c(ccodes, "blue")
  }
}

heatmap.2(x = as.matrix(heat.data),
          col = hcl.colors(10, palette = 'Blue-Red 2'), # search hcl colors in r
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = 'none',
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Normal vs Tumor",
          na.color = 'black',
          ColSideColors = ccodes)
