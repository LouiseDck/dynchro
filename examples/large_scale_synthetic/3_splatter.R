library(dyngen)
library(dyno)
library(splatter)
library(scater)
library(SingleCellExperiment)
library(dynplot)

# params <- newSplatParams()
# params <- setParam(params, "nGenes", 2000)
# params <- setParam(params, "nCells", 10000)

ds <- readRDS("examples/large_scale_synthetic/vsc42417/dataset.rds")
model <- readRDS("examples/large_scale_synthetic/vsc42417/model_xp.rds")

plot_dimred(ds)

counts <- t(as.matrix(ds$counts))
scec <- SingleCellExperiment(assays = list(counts = counts))
scec <- logNormCounts(scec)
gt <- scater::runPCA(scec)
plotPCA(gt)

params_estimated <- splatEstimate(counts)
saveRDS(params_estimated, "examples/large_scale_synthetic/params_estimated.rds")

sim <- splatSimulate(
  params_estimated,
  batchCells = c(1000, 1000),
  group.prob = c(0.65, 0.1, 0.25),
  path.from = c(0, 1, 2),
  method = "paths"
)
sim <- logNormCounts(sim)
plotPCA(sim, colour_by = "Batch")
saveRDS(sim, "examples/large_scale_synthetic/sim.rds")
sim <- logNormCounts(sim)
saveRDS(sim, "examples/large_scale_synthetic/sim_lognorm.rds")

sim <- scater::runPCA(sim)
plotPCA(sim, colour_by = "Batch")

sim2 <- splatSimulate(
  params_estimated,
  batchCells = c(1000, 1000, 1000, 1000, 1000, 250, 1250, 3000),
  method = "paths"
)
sim2 <- logNormCounts(sim2)
saveRDS(sim2, "examples/large_scale_synthetic/sim2_lognorm.rds")
sim2 <- scater::runPCA(sim2)
plotPCA(sim2, colour_by = "Batch")
