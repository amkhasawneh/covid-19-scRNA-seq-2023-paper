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
  saveRDS(split.covid[[i]], file = paste0(names(split.covid[i]), ".rds"))
  names(split.covid[i])
}

#After uploading the separate RDS files to Azimuth, incorporating the results
#by importing Azimuth's results for each sample:
for (i in names(split.covid)) {
  split.covid[[i]]$azimuthNames <- read.table(paste0(i, "_azimuth_pred.tsv"), sep = "\t", header = T)$predicted.celltype.l2
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

