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

#Loading data:
mild227 <- Read10X("../data/from_cellranger/mild227/sample_feature_bc_matrix")
moderate124 <- Read10X("../data/from_cellranger/moderate124/sample_feature_bc_matrix")
mild186 <- Read10X("../data/from_cellranger/mild186/sample_feature_bc_matrix")
critical213<- Read10X("../data/from_cellranger/critical213/sample_feature_bc_matrix")
moderate303 <- Read10X("../data/from_cellranger/moderate303/sample_feature_bc_matrix")
critical308 <- Read10X("../data/from_cellranger/critical308/sample_feature_bc_matrix")
moderate138 <- Read10X("../data/from_cellranger/moderate138/sample_feature_bc_matrix")
moderate272 <- Read10X("../data/from_cellranger/moderate272/sample_feature_bc_matrix")
critical119 <- Read10X("../data/from_cellranger/critical119/sample_feature_bc_matrix")
critical120 <- Read10X("../data/from_cellranger/critical120/sample_feature_bc_matrix")
critical238 <- Read10X("../data/from_cellranger/critical238/sample_feature_bc_matrix")
critical293 <- Read10X("../data/from_cellranger/critical293/sample_feature_bc_matrix")
healthy1 <- Read10X("../data/from_cellranger/healthy1/sample_feature_bc_matrix")
healthy2 <- Read10X("../data/from_cellranger/healthy2/sample_feature_bc_matrix")
healthy3 <- Read10X("../data/from_cellranger/healthy3/sample_feature_bc_matrix")

#Creating Seurat objects:
mild227 <- CreateSeuratObject(counts = mild227, min.cells = 3, min.features = 100)
mild186 <- CreateSeuratObject(counts = mild186, min.cells = 3, min.features = 100)
moderate124 <- CreateSeuratObject(counts = moderate124, min.cells = 3, min.features = 100)
moderate138 <- CreateSeuratObject(counts = moderate138, min.cells = 3, min.features = 100)
moderate272 <- CreateSeuratObject(counts = moderate272, min.cells = 3, min.features = 100)
moderate303 <- CreateSeuratObject(counts = moderate303, min.cells = 3, min.features = 100)
critical119 <- CreateSeuratObject(counts = critical119, min.cells = 3, min.features = 100)
critical120 <- CreateSeuratObject(counts = critical120, min.cells = 3, min.features = 100)
critical213 <- CreateSeuratObject(counts = critical213, min.cells = 3, min.features = 100)
critical238 <- CreateSeuratObject(counts = critical238, min.cells = 3, min.features = 100)
critical293 <- CreateSeuratObject(counts = critical293, min.cells = 3, min.features = 100)
critical308 <- CreateSeuratObject(counts = critical308, min.cells = 3, min.features = 100)
healthy1 <- CreateSeuratObject(counts = healthy1, min.cells = 3, min.features = 100)
healthy2 <- CreateSeuratObject(counts = healthy2, min.cells = 3, min.features = 100)
healthy3 <- CreateSeuratObject(counts = healthy3, min.cells = 3, min.features = 100)

covid <- merge(x = critical119, y = c(critical120, critical213, critical238,
                                      critical293, critical308,
                                      mild186, mild227, 
                                      moderate124, moderate138, moderate272, moderate303,
                                      healthy1, healthy2, healthy3), 
               add.cell.ids = c("critical119_Patient5", "critical120_Patient6", "critical213_Patient3",
                                "critical238_Patient4", "critical293_Patient1", "critical308_Patient2", 
                                "mild186_Patient3", "mild227_Patient4", 
                                "moderate124_Patient6", "moderate138_Patient5", "moderate272_Patient1", "moderate303_Patient2",
                                "healthy1_control1", "healthy2_control2", "healthy3_control3"),
               project = "covid-19")

head(covid@meta.data)

#Saving file for upload to Azimuth:
saveRDS(covid, "00-covid-raw.rds")

remove(critical119, critical120, critical238, critical213, critical293, critical308,
      mild186, mild227, moderate124, moderate138, moderate272, moderate303,
      healthy1, healthy2, healthy3)
gc()

#################################Environment####################################
#Preparing the environment:
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
gc()

