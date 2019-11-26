library(dplyr)
library(Seurat)

# load data from Seurat package
data(cc.genes)

# segregate this list into markers of G2/M phase and markers of S phase
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes

# Load the dataset (header = TRUE, as.is = TRUE, row.names = 1)
feature_counts.data <- read.table(file = "Desktop/Seurat/S/Data/feature_counts.table.tsv", 
                                  sep = "\t",
                                  header = TRUE,
                                  quote = "",
                                  stringsAsFactors = FALSE)

dim(x=feature_counts.data)
#head(x = rownames(x = feature_counts.data))
#head(x = colnames(x = feature_counts.data))
#names(x = feature_counts.data)

# Convert gene names
geneInfoFile <- read.table(file = "Desktop/Seurat/S/Data/geneInfo.QC.method.seurat.SCORPIUS.tsv", 
                           sep = "\t",
                           header = TRUE,
                           quote = "",
                           stringsAsFactors = FALSE)

dim(x=geneInfoFile)
#head(x = rownames(x = geneInfoFile))
#head(x = colnames(x = geneInfoFile))

geneInfoNames = geneInfoFile %>%  select(Geneid,external_gene_name)
#head(geneInfoNames)
dim(x=geneInfoNames)
#head(x = rownames(x = geneInfoNames))
#head(x = colnames(x = geneInfoNames))
geneInfo = inner_join(geneInfoNames, feature_counts.data, by = c("Geneid" = "Geneid"))

dim(x=geneInfo)
#head(x = rownames(x = geneInfo))
head(x = colnames(x = geneInfo))
#head(geneInfo)
geneInfo = geneInfo[,-1]
dim(x=geneInfo)
head(geneInfo)
head(x = colnames(x = geneInfo))

#remove data that is not sc
removeIDs <- function(geneInfo) {
  column_indices_ID <- c()
  
  for(i in 1:ncol(geneInfo)) {
    if(grepl(".dil.", colnames(geneInfo[i]))) {
      column_indices_ID <- c(column_indices_ID, i)
    }
    if(grepl(".FACS.", colnames(geneInfo[i]))) {
      column_indices_ID <- c(column_indices_ID, i)
    }
  }
  feature_counts_removedColumns <- c()
  feature_counts_removedColumns <- geneInfo[,-column_indices_ID]
  return(feature_counts_removedColumns)
}

geneInfo_SC <- removeIDs(geneInfo)
dim(geneInfo_SC)
head(x = colnames(x = geneInfo_SC))
head(x = rownames(x = geneInfo_SC))
rownames(geneInfo_SC) =  geneInfo_SC$external_gene_name

dim(x=geneInfo)
#geneInfo
feature_counts.data = geneInfo[,-1:-7]
ncol(feature_counts.data)
nrow(feature_counts.data)

#remove genes that are not expressed
rowSums(feature_counts_SC)
keep <- rowSums(feature_counts_SC>1)
feature_counts_SC <- feature_counts_SC[keep,]
dim(feature_counts_SC)

#Create a Seurat object out of the table
feature_counts_seurat <- CreateSeuratObject(counts = feature_counts_SC)
feature_counts_seurat
#summary(feature_counts.data)
feature_counts_seurat <- NormalizeData(object = feature_counts_seurat)
feature_counts_seurat <- FindVariableFeatures(object = feature_counts_seurat, selection.method = "vst" )
feature_counts_seurat
feature_counts_seurat <- ScaleData(feature_counts_seurat, features = rownames(feature_counts_seurat))
feature_counts_seurat <- RunPCA(object = feature_counts_seurat)
feature_counts_seurat

DimHeatmap(feature_counts_seurat)


#Assign cell-cycle scores
feature_counts_seurat <- CellCycleScoring(feature_counts_seurat, s.features = s.genes, g2m.features = g2m.genes)

# view cell cycle scores and phase assignments
head(feature_counts_seurat[[]])