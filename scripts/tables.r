

#Table S10: cell numbers for each sample:
cells <- BCR@meta.data %>%
  group_by(sample, azimuthNames) %>% dplyr::count() %>% 
  spread(key = sample, value = n) %>% as.data.frame()
rownames(cells) <- cells$azimuthNames
cells$azimuthNames <- NULL
cells <- as.matrix(cells)
write.table(cells, file = "../results/tables/cell-numbers-samples.tsv", sep = "\t", col.names = NA)
cells <- proportions(as.matrix(cells), margin = 2) * 100
write.table(cells, file = "../results/tables/cell-proportions-samples.tsv", sep = "\t", col.names = NA)

