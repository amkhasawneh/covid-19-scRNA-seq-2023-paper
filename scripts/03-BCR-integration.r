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

################################Integration with Seurat object##################

B.combined <- readRDS("../data/byproducts/BCR-combined-contigs.rds")

#Loading clustered Seurat object:
covid <- readRDS("../data/byproducts/04-covid-clustered.rds")

#Changing the current Idents:
Idents(covid) <- "azimuthNames"

#Subsetting B-cells:
BCR <- covid[,grep("(B)|(Plasmablast)", covid$azimuthNames)]

#Combining Seurat object with BCR table:
BCR <- combineExpression(B.combined, BCR, cloneCall="gene", filterNA = T)
head(BCR@meta.data)
gc()

#Adjusting the sample variable's order:
BCR$sample <- factor(BCR$sample, levels = c("healthy1_control1", "healthy2_control2", "healthy3_control3",
                                            "moderate272_Patient1", "critical293_Patient1", 
                                            "moderate303_Patient2", "critical308_Patient2",
                                            "mild186_Patient3", "critical213_Patient3",
                                            "mild227_Patient4", "critical238_Patient4",
                                            "critical119_Patient5", "moderate138_Patient5",
                                            "critical120_Patient6", "moderate124_Patient6"))

#Modifying cell type order:
BCR$azimuthNames <- factor(BCR$azimuthNames, levels = c("B intermediate", "B memory", "B naive", "Plasmablast"))

#Adding V and J gene usage:
BCR$HChain <- vapply(strsplit(BCR$CTgene, "[_]"), "[", "", 1)
BCR$HChain <- sub("(NA)", NA, BCR$HChain)
BCR$v_gene <- vapply(strsplit(BCR$HChain, "[.]"), "[", "", 1)
BCR$j_gene <- vapply(strsplit(BCR$HChain, "[.]"), "[", "", 2)
BCR$c_gene <- vapply(strsplit(BCR$HChain, "[.]"), "[", "", 4)

BCR$LChain <- vapply(strsplit(BCR$CTgene, "[_]"), "[", "", 2)
BCR$LChain <- sub("(NA)", NA, BCR$LChain)
BCR$lkv_gene <- vapply(strsplit(BCR$LChain, "[.]"), "[", "", 1)
BCR$lkj_gene <- vapply(strsplit(BCR$LChain, "[.]"), "[", "", 2)
BCR$lkc_gene <- vapply(strsplit(BCR$LChain, "[.]"), "[", "", 3)

saveRDS(BCR, "../data/byproducts/05-BCR-combined.rds")

BCR <- readRDS("../data/byproducts/05-BCR-combined.rds")

remove(covid, B.combined, abundace, top.clonotypes, longest)
gc()

#Saving some tables:
write.table(BCR@meta.data[order(BCR@meta.data[["Frequency"]],decreasing=TRUE),],
            file = "../results/tables/BCR-frequency.tsv", sep="\t", append = FALSE, quote=FALSE, 
            row.names = FALSE, col.names = TRUE)

#Tables with cell numbers for each sample:
cells <- BCR@meta.data %>%
  group_by(sample, azimuthNames) %>% count() %>% 
  spread(key = sample, value = n) %>% as.data.frame()
rownames(cells) <- cells$azimuthNames
cells$azimuthNames <- NULL
cells <- as.matrix(cells)
cells <- proportions(as.matrix(cells), margin = 2) * 100
write.table(cells, file = "../results/tables/cell-proportions-samples.tsv", sep = "\t", col.names = NA)

