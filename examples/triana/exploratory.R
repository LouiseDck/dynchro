library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Load the data
wd <- "Work/dynchro/examples/triana"
aml <- readRDS("Work/dynchro/examples/triana/AML.rds")
healthy <- readRDS("Work/dynchro/examples/triana/Healthy.rds")
wta_projected <- readRDS("Work/dynchro/examples/triana/WTA_projected.rds")
ab200_projected <- readRDS("Work/dynchro/examples/triana/200AB_projected.rds")

DimPlot(healthy)

DimPlot(wta_projected)
DimPlot(ab200_projected)

DimPlot(healthy, reduction = "MOFAUMAP")

healthy_wta <- merge(healthy, wta_projected, add.cell.ids = c("Healthy", "WTA"))
healthy_wta_umap <-
  DimPlot(healthy_wta)


healthy_umap <- healthy@reductions$MOFAUMAP@cell.embeddings
wta_projected_umap <- wta_projected@reductions$ProjectedUMAP@cell.embeddings
ab200_projected_umap <- ab200_projected@reductions$ProjectedUMAP@cell.embeddings

test1 <- rbind(healthy_umap, wta_projected_umap)
test <- rbind(test1, ab200_projected_umap)
test <- as.data.frame(test)

sample <- c(
  rep("Healthy", nrow(healthy_umap)),
  rep("WTA", nrow(wta_projected_umap)),
  rep("200AB", nrow(ab200_projected_umap))
)

p <- ggplot(test, aes(x = MOFAUMAP_1, y = MOFAUMAP_2, color = sample)) +
  geom_point()
p

p2 <- ggplot(healthy_umap, aes(x = MOFAUMAP_1, y = MOFAUMAP_2)) + geom_point()
p2
