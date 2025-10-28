library(tidyverse)
library(dyngen)
library(dyno)
library(wrapr)
library(anndata)
source("examples/large_scale_synthetic/0_diverging_kinetics.R")

DEBUG <- FALSE

backbone <- backbone_linear()
config <-
  initialise_model(
    backbone = backbone,
    num_cells = 10000,
    num_tfs = nrow(backbone$module_info),
    num_targets = 2000,
    num_hks = 5000,
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 100)
    )
  )

if (DEBUG) {
  config <-
    initialise_model(
      backbone = backbone,
      num_cells = 1000,
      num_tfs = nrow(backbone$module_info),
      num_targets = 50,
      num_hks = 50,
      verbose = interactive(),
      download_cache_dir = tools::R_user_dir("dyngen", "data"),
      simulation_params = simulation_default(
        census_interval = 5,
        ssa_algorithm = ssa_etl(tau = .01),
        experiment_params = simulation_type_wild_type(num_simulations = 10)
      )
    )
}

num_cells <- 10000
num_features <- 2000
num_tfs <- 100
num_targets <- round((num_features - num_tfs) / 2)
num_hks <- num_features - num_targets - num_tfs

out <-
  initialise_model(
    backbone = backbone,
    num_tfs = num_tfs,
    num_targets = num_targets,
    num_hks = num_hks,
    num_cells = num_cells,
    gold_standard_params = gold_standard_default(
      census_interval = 1,
      tau = 100 / 3600
    ),
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(
        num_simulations = num_cells / 10
      )
    ),
    verbose = TRUE
  ) %>%
  generate_dataset(make_plots = TRUE)

write_rds(out, "dataset.rds")


#############################
# This generates 2 datasets with a 99 procent batch effect, gradient?
# Generates paired datasets
#############################

generate_batch_effect <- function(index, config) {
  print(index)
  model_common <- config %>%
    generate_tf_network() %>%
    generate_feature_network()

  model_a <- model_common %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells()

  model_b <- model_common %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells()

  model_a <- model_a %>% generate_experiment()
  dataset_a <- model_a %>% as_dyno()

  # a1 <- dataset_a$milestone_percentages %>% filter(milestone_id == "sA")
  # a1 <- a1[orderv(a1[, "percentage"]), ][1, 1]$cell_id
  # adata_ds <- as_anndata(model_a)
  # adata_ds$obs[["model"]] <- dataset_a$cell_info[["model"]]
  # adata_ds$obs[["milestones"]] <- dataset_a$progressions$to
  # adata_ds$uns[["iroot"]] <- a1
  # anndata::write_h5ad(
  #     adata_ds,
  #     paste0("examples/large_scale_synthetic/data/batcheffect_dataseta", index, ".h5ad"))

  saveRDS(
    dataset_a,
    paste0(
      "examples/large_scale_synthetic/data/batcheffect_dataseta",
      index,
      ".rds"
    )
  )

  model_between2 <- generate_diverging_kinetics(model_a, model_b, 0.9)

  dataset_between <- model_between2 %>% as_dyno()

  # b1 <- dataset_between$milestone_percentages %>% filter(milestone_id == "sA")
  # b1 <- b1[orderv(b1[, "percentage"]), ][1, 1]$cell_id

  # adata_btwn <- as_anndata(model_between2)
  # adata_btwn$obs[["model"]] <- dataset_between$cell_info[["model"]]
  # adata_btwn$obs[["milestones"]] <- dataset_between$progressions$to
  # adata_btwn$uns[["iroot"]] <- b1
  # anndata::write_h5ad(
  #     adata_btwn,
  #     paste0("examples/large_scale_synthetic/data/batcheffect_datasetb", index, ".h5ad"))

  saveRDS(
    dataset_between,
    paste0(
      "examples/large_scale_synthetic/data/batcheffect_datasetb",
      index,
      ".rds"
    )
  )
}
set.seed(9356)

lapply(1:10, function(x) generate_batch_effect(x, config))
lapply(11:20, function(x) generate_batch_effect(x, config))
lapply(21:30, function(x) generate_batch_effect(x, config))
lapply(31:40, function(x) generate_batch_effect(x, config))
lapply(41:50, function(x) generate_batch_effect(x, config))


generate_batch_effect(0, config)
generate_batch_effect(1, config)
generate_batch_effect(2, config)
