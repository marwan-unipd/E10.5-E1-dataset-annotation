---
title: "R Notebook"
output: html_notebook
---

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")
```
```{r}
install.packages('Seurat')
```

```{r}
#very basic Seurat preprocessing
library(dplyr)
library(Seurat)
 Raw_data <- Read10X(data.dir ='.//new Dataset/M1RE105/M1RE105/embryo1/filtered_feature_bc_matrix/' )
 seuset_data <- CreateSeuratObject(counts = Raw_data, min.cells = 3, min.features = 200)
 seuset_data[["percent.mt"]] <- PercentageFeatureSet(seuset_data, pattern = "mt-")
 lb <- quantile(seuset_data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.01)
 ub <- quantile(seuset_data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.99)
 seuset_data <- seuset_data[, seuset_data[["nFeature_RNA"]] > lb & seuset_data[["nFeature_RNA"]] < ub & seuset_data[["percent.mt"]] < 25] 
 seuset_data <- NormalizeData(object = seuset_data, verbose = FALSE)
 seuset_data <- FindVariableFeatures(object = seuset_data, nfeatures = 3000, verbose = FALSE, selection.method = 'vst')
 seuset_data <- ScaleData(seuset_data, verbose = FALSE)
 seuset_data <- RunPCA(seuset_data, npcs = 30, verbose = FALSE)
 seuset_data <- FindNeighbors(seuset_data, dims = 1:30)
 seuset_data <- FindClusters(seuset_data, resolution = 0.3)
 seuset_data <- RunUMAP(seuset_data, reduction = "pca", dims = 1:30)





data <- seuset_data


DimPlot(data, reduction = "umap", label=TRU




library(SingleR)



ref <- celldex::MouseRNAseqData()



results <- SingleR(test = as.SingleCellExperiment(data), ref = ref, labels = ref$label.main)



data$singlr_labels <- results$labels



DimPlot(data, reduction = 'umap', group.by = 'singlr_labels', label = TRUE)



FeaturePlot(data, features = c("Ptprc", "Cd3e"))



if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRNAseq")


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scuttle")


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TabulaMurisData")


library(ExperimentHub)


eh <- ExperimentHub()


query(eh, "TabulaMurisData")


eh[['EH1617']]



lung_ref <- eh[['EH1617']]

lung_ref <- lung_ref[,!is.na(lung_ref$cell_ontology_class)]
```

```{r}
lung_ref
```

```{r}
library(scuttle)
```

```{r}
lung_ref <- logNormCounts(lung_ref)
```


```{r}
results <- SingleR(test = as.SingleCellExperiment(data), ref = lung_ref, labels = lung_ref$cell_ontology_class)
```


```{r}
data$singlr_label <- results$labels
```

```{r}
 DimPlot(data, reduction = "umap", group.by = 'singlr_label', label = FALSE)


















