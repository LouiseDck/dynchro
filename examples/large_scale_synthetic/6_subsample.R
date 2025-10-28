library(tidyverse)
library(dyno)
library(dyngen)

modelb3 <- readRDS(
  "examples/large_scale_synthetic/branching_data_vsc42417/model_wt_exp_b3.rds"
)
modelb2 <- readRDS(
  "examples/large_scale_synthetic/branching_data_vsc42417/model_wt_exp_b2.rds"
)

datasetb3 <- readRDS(
  "examples/large_scale_synthetic/branching_data_vsc42417/dataset_b3.rds"
)
datasetb2 <- readRDS(
  "examples/large_scale_synthetic/branching_data_vsc42417/dataset_b2.rds"
)


library(anndata)
ad_b3 <- as_anndata(modelb3)
ad_b2 <- as_anndata(modelb2)

ad_b3$obs$from <- modelb3$experiment$cell_info$from
ad_b3$obsm$lmds <- datasetb3$dimred

ad_b2$obs$from <- modelb2$experiment$cell_info$from
ad_b2$obsm$lmds <- datasetb2$dimred

write_h5ad(
  ad_b3,
  "examples/large_scale_synthetic/branching_data_vsc42417/ad_b3.h5ad"
)
write_h5ad(
  ad_b2,
  "examples/large_scale_synthetic/branching_data_vsc42417/ad_b2.h5ad"
)

modelb2$gold_standard
