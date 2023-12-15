################################Loading#########################################
library(scRepertoire) 
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(stringr)
library(circlize)
library(scales)
library(dplyr)
library(pals)
library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(cowplot)
library(vegan)

#Loading the data:
critical119 <- read.csv("../data/from_cellranger/critical119/vdj_b/filtered_contig_annotations.csv")
critical120 <- read.csv("../data/from_cellranger/critical120/vdj_b/filtered_contig_annotations.csv")
critical213 <- read.csv("../data/from_cellranger/critical213/vdj_b/filtered_contig_annotations.csv")
critical238 <- read.csv("../data/from_cellranger/critical238/vdj_b/filtered_contig_annotations.csv")
critical293 <- read.csv("../data/from_cellranger/critical293/vdj_b/filtered_contig_annotations.csv")
critical308 <- read.csv("../data/from_cellranger/critical308/vdj_b/filtered_contig_annotations.csv")
mild186 <- read.csv("../data/from_cellranger/mild186/vdj_b/filtered_contig_annotations.csv")
mild227 <- read.csv("../data/from_cellranger/mild227/vdj_b/filtered_contig_annotations.csv")
moderate124 <- read.csv("../data/from_cellranger/moderate124/vdj_b/filtered_contig_annotations.csv")
moderate138 <- read.csv("../data/from_cellranger/moderate138/vdj_b/filtered_contig_annotations.csv")
moderate272 <- read.csv("../data/from_cellranger/moderate272/vdj_b/filtered_contig_annotations.csv")
moderate303 <- read.csv("../data/from_cellranger/moderate303/vdj_b/filtered_contig_annotations.csv")
healthy1 <- read.csv("../data/from_cellranger/healthy1/vdj_b/filtered_contig_annotations.csv")
healthy2 <- read.csv("../data/from_cellranger/healthy2/vdj_b/filtered_contig_annotations.csv")
healthy3 <- read.csv("../data/from_cellranger/healthy3/vdj_b/filtered_contig_annotations.csv")


#Combinging the data frames:
contig.list <- list(critical119, critical120, critical213,
                    critical238, critical293, critical308,
                    mild186, mild227, 
                    moderate124, moderate138, moderate272, moderate303,
                    healthy1, healthy2, healthy3)

B.combined <- combineBCR(contig.list, 
                         samples = c("critical119", "critical120", "critical213",
                                     "critical238", "critical293", "critical308",
                                     "mild186","mild227",
                                     "moderate124", "moderate138", "moderate272", "moderate303",
                                     "healthy1", "healthy2", "healthy3"),
                         ID = c("Patient5", "Patient6", "Patient3",
                                "Patient4", "Patient1", "Patient2",
                                "Patient3", "Patient4", 
                                "Patient6", "Patient5", "Patient1", "Patient2",
                                "control1", "control2", "control3"))

B.combined <- addVariable(B.combined, name = "severity", 
                          variables = c("critical", "critical", "critical",
                                        "critical", "critical", "critical",
                                        "mild", "mild",
                                        "moderate", "moderate", "moderate", "moderate",
                                        "healthy", "healthy", "healthy"))
head(B.combined)

#Saving the contigs:
saveRDS(B.combined, "../data/byproducts/BCR-combined-contigs.rds")

remove(critical119, critical120, critical213,
       critical238, critical293, critical308,
       mild186, mild227, 
       moderate124, moderate138, moderate272, moderate303,
       healthy1, healthy2, healthy3, contig.list)
gc()
