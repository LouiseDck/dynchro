library(tidyverse)
library(dyngen)
library(dyno)
library(assertthat)
library(tibble)

set.seed(87961)
DEBUG <- TRUE

path <- Sys.getenv("VSC_DATA")

backbone <- backbone_bifurcating()

num_cells <- 100000
num_features <- 2000
num_tfs <- 100
num_targets <- round((num_features - num_tfs) / 2)
num_hks <- num_features - num_targets - num_tfs

if (DEBUG) {
  num_cells <- 100
  num_features <- 200
  num_tfs <- 10
  num_targets <- round((num_features - num_tfs) / 2)
  num_hks <- num_features - num_targets - num_tfs
}

model <- initialise_model(
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
)

#########################
# Generate common model #
#########################

model_common_tf <- generate_tf_network(model)
saveRDS(model_common_tf, paste0(path, "/model_common_tf.rds"))
model_common_ft <- generate_feature_network(model_common_tf)
saveRDS(model_common_ft, paste0(path, "/model_common_ft.rds"))
model_common_kn <- generate_kinetics(model_common_ft)
saveRDS(model_common_kn, paste0(path, "/model_common_kn.rds"))
model_common_gs <- generate_gold_standard(model_common_kn)
# p_common_gs <- plot_gold_mappings(model_common_gs, do_facet = FALSE)
ggsave(paste0(path, "/p_common_gs.png"), p_common_gs)
saveRDS(model_common_gs, paste0(path, "/model_common_gs.rds"))

#####################
# Generate WT model #
#####################

simulation_type_knockdown_ <- function(
  num_simulations,
  timepoint = runif(num_simulations),
  genes = "*",
  num_genes = sample(1:5, num_simulations, replace = TRUE, prob = 0.25^(1:5)),
  multiplier = runif(num_simulations, 0, 1),
  seed = sample.int(10 * num_simulations, num_simulations)
) {
  if (num_simulations == 0) {
    NULL
  } else {
    assert_that(
      (is.list(genes) && length(genes) == num_simulations) ||
        is.character(genes) ||
        is.factor(genes)
    )
    if (!is.list(genes)) genes <- as.list(rep(list(genes), num_simulations))

    tibble(
      type = "knockdown",
      timepoint,
      genes,
      num_genes,
      multiplier,
      seed
    )
  }
}

knockout <- function(model, module) {
  genes <- model$feature_info %>%
    filter(module_id == module) %>%
    pull(feature_id)

  nr_sim <- 100L
  genes_list <- rep(list(genes), nr_sim)

  model$simulation_params$experiment_params <- simulation_type_knockdown_(
    num_simulations = nr_sim,
    timepoint = 0,
    genes = genes_list,
    num_genes = length(genes),
    multiplier = 0
  )
  return(model)
}

model_b3 <- model_common_gs %>%
  knockout("B3")
saveRDS(model_b3, paste0(path, "/model_wt_b3.rds"))
model_b3_cells <- model_b3 %>% generate_cells()
saveRDS(model_b3_cells, paste0(path, "/model_wt_cells_b3.rds"))
model_b3_exp <- model_b3_cells %>% generate_experiment()
saveRDS(model_b3_exp, paste0(path, "/model_wt_exp_b3.rds"))
ds_b3 <- as_dyno(model_b3_exp)
saveRDS(ds_b3, paste0(path, "/dataset_b3.rds"))

model_b2 <- model_common_gs %>%
  knockout("B2")
saveRDS(model_b2, paste0(path, "/model_wt_b2.rds"))
model_b2_cells <- model_b2 %>% generate_cells()
saveRDS(model_b2_cells, paste0(path, "/model_wt_cells_b2.rds"))
model_b2_exp <- model_b2_cells %>% generate_experiment()
saveRDS(model_b2_exp, paste0(path, "/model_wt_exp_b2.rds"))
ds_b2 <- as_dyno(model_b2_exp)
saveRDS(ds_b2, paste0(path, "/dataset_b2.rds"))
