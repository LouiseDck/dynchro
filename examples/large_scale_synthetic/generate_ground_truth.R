library(tidyverse)
library(dyngen)
library(dyno)
library(ggplot2)

hpc <- FALSE

if (hpc) {
  backbone <- backbone_branching()
  num_cells <- 10000
  num_features <- 2000
  num_tfs <- 100
} else {
  backbone <- backbone_branching()
  num_cells <- 1000
  num_features <- 200
  num_tfs <- 10
}
num_targets <- round((num_features - num_tfs) / 2)
num_hks <- num_features - num_targets - num_tfs

model <- initialise_model(backbone = backbone,
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

p0 <- plot_backbone_modulenet(model)
ggsave("p0.png", p0)
p0
model_tf <- generate_tf_network(model)
p1 <- plot_feature_network(model_tf)
ggsave("p1.png", p1)
saveRDS(model_tf, "model_tf.rds")

model_ft <- generate_feature_network(model_tf)
p2 <- plot_feature_network(model_ft)
ggsave("p_feature.png", p2)
saveRDS(model_ft, "model_ft.rds")

model_kn <- generate_kinetics(model_ft)
p3 <- plot_feature_network(model_kn)
ggsave("p_kinetics.png", p3)
saveRDS(model_kn, "model_kn.rds")

model_gs <- generate_gold_standard(model_kn)
p4 <- plot_gold_simulations(model_gs)
ggsave("p_gold.png", p4)
saveRDS(model_gs, "model_gs.rds")

model_cl <- generate_cells(model_gs)
p5 <- plot_simulations(model_cl)
ggsave("p_sim.png", p5)
saveRDS(model_cl, "model_cl.rds")

model_xp <- generate_experiment(model_cl)
saveRDS(model_xp, "model_xp.rds")

dataset <- as_dyno(model_xp)
saveRDS(dataset, "dataset.rds")

library(dynplot)
p6 <- plot_dimred(dataset)
ggsave("p_dataset.png", p6)
