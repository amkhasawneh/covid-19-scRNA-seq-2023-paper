#################################Loading########################################
#Loading libraries:
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)


#Preparing the environment:
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()

#################################Split-Scaling & Azimuth########################

covid <- readRDS("../data/byproducts/01-covid-filtered.rds")

#Splitting Seurat object by sample to perform scaling on all samples:
#First splitting into samples:
split.covid <- SplitObject(covid, split.by = "sample")

#Scaling:
for (i in 1:length(split.covid)) {
  split.covid[[i]] <- NormalizeData(split.covid[[i]], verbose = T)
  split.covid[[i]] <- FindVariableFeatures(split.covid[[i]], verbose = T)
  split.covid[[i]] <- ScaleData(split.covid[[i]], verbose = T)
}

#Saving separate RDS files for each sample to upload to Azimuth:
for (i in names(split.covid)) {
  saveRDS(split.covid[[i]], file = paste0("../data/byproducts/", names(split.covid[i]), ".rds"))
  names(split.covid[i])
}

#After uploading the separate RDS files to Azimuth, incorporating the results
#by importing Azimuth's results for each sample:
for (i in names(split.covid)) {
  split.covid[[i]]$azimuthNames <- read.table(paste0("../data/byproducts/", i, "_azimuth_pred.tsv"), sep = "\t", header = T)$predicted.celltype.l2
}

gc()

#Saving the split object:
saveRDS(split.covid, "../data/byproducts/02-covid-split-scaled.rds")
remove(critical119_Patient5, critical120_Patient6, critical238_Patient4, critical293_Patient1,
       critical213_Patient3, critical308_Patient2,
       mild227_Patient4, mild186_Patient3,
       moderate124_Patient6, moderate138_Patient5, moderate272_Patient1, moderate303_Patient2,
       healthy1_control1, healthy2_control2, healthy3_control3)

gc()

#################################Integration####################################

#Importing the split file:
split.covid <- readRDS("../data/byproducts/02-split-scaled.rds")

#Selecting the most variable features for integration:
integ.feat <- SelectIntegrationFeatures(object.list = split.covid, nfeatures = 5000)

#Performing canonical correlation analysis (CCA), and identifying and filtering anchors:
anchors <- FindIntegrationAnchors(object.list = split.covid, anchor.features = integ.feat,
                                  verbose = T)                  #
#These two steps had to be done
saveRDS(anchors, "../data/byproducts/anchors.rds")                                 #on the Linux machine,
#as an RScript,
#due to memory restrictions.
remove(covid, split.covid, integ.feat)                          ###(They need a ton of RAM!)
gc()                                                            
#
#
anchors <- readRDS("../data/byproducts/anchors.rds")                               #
#  
#Integrating across conditions:                             #    
covid <- IntegrateData(anchorset = anchors, verbose = T)###             

remove(anchors)
gc()

saveRDS(covid, "../data/byproducts/03-covid-integrated.rds")
covid <- readRDS("../data/byproducts/03-covid-integrated.rds")

DefaultAssay(covid) <- "integrated"

#Scaling the data:
covid <- ScaleData(covid, features = rownames(covid))

#Linear dimensional reduction:
covid <- RunPCA(covid, features = VariableFeatures(object = covid), verbose = T)

#Examining PCA results:
print(covid[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(covid, dims = 1:2, reduction = "pca")
DimPlot(covid, reduction = "pca", split.by = "severity")
DimPlot(covid, reduction = "pca")
DimHeatmap(covid, dims = 1:15, cells = 500, balanced = TRUE)

#Determining the dimensionality of the data:
ElbowPlot(covid, ndims = 40) #This implies going for around 30 might be OK.

saveRDS(covid, "../data/byproducts/03-covid-integrated.rds")


#################################UMAP Clustering################################

covid <- readRDS("../data/byproducts/03-covid-integrated.rds")
gc()

#Graph-based clustering:
covid <- FindNeighbors(covid, dims = 1:40, verbose = T)
covid <- FindClusters(covid, verbose = T,
                      resolution = 0.8)

Idents(covid) <- "integrated_snn_res.0.8" #This one seems to make the most sense.

#Non-linear dimensional reduction:
covid <- RunUMAP(covid, dims = 1:40, return.model = T)
Idents(covid) <- "integrated_snn_res.0.8"
DimPlot(covid, reduction = "umap", label = T, repel = T, raster = T, group.by = "azimuthNames") + NoLegend()
DimPlot(covid, reduction = "umap", split.by = "severity") + NoLegend()
DimPlot(covid, reduction = "umap", split.by = "sample") + NoLegend()
gc()

saveRDS(covid, "../data/byproducts/04-covid-clustered.rds")
