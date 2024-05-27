library("MouseGastrulationData")
library("ExperimentHub")
##########################
# READ IN DATA
##########################
sce <- EmbryoAtlasData(type = "processed")
sce_raw <- EmbryoAtlasData(type = "raw")
wt <- WTChimeraData()
tal1 <- Tal1ChimeraData()
tch <- TChimeraData()

##########################
# RECREATE FIG1
##########################

singlets <- which(!(colData(sce)$doublet | colData(sce)$stripped))
par(mar=c(4, 4, 4, 15), xpd=TRUE)
plot(
    x = reducedDim(sce, "umap")[singlets, 1],
    y = reducedDim(sce, "umap")[singlets, 2],
    col = EmbryoCelltypeColours[colData(sce)$celltype[singlets]],
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)
coord <- par("usr")
color_factors <- as.factor(colData(sce)$celltype[singlets])
legend_colors <- as.vector(EmbryoCelltypeColours[levels(color_factors)])
legend(x = coord[2] * 1.05, y = coord[4], legend = levels(color_factors), fill = EmbryoCelltypeColours, xpd = TRUE)

##########################
# SELECT BLOOD LINEAGE
##########################

# Ery BP Haem EC Mes Mk BP4 My Haem4 EC7
# Blood progenitors 1
# Blood progenitors 2
# Endothelium
# Erythroid1
# Erythroid2
# Erythroid3
# Haematoendothelial progeniros
# Mixed mesoderm
bloodlabels <- c(
    "Blood progenitors 1",
    "Blood progenitors 2",
    "Endothelium",
    "Erythroid1",
    "Erythroid2",
    "Erythroid3",
    "Haematoendothelial progenitors",
    "Mixed mesoderm"
)
bloodcells <- which(colData(sce)$celltype %in% bloodlabels)

plot(
    x = reducedDim(sce, "umap")[bloodcells, 1],
    y = reducedDim(sce, "umap")[bloodcells, 2],
    col = EmbryoCelltypeColours[colData(sce)$celltype[bloodcells]],
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)

##########################
# FORCE DIRECTED LAYOUT
##########################

library("zellkonverter")
writeH5AD(sce, "sce.h5ad")
adata <- SCE2AnnData(sce)

plot(
    x = reducedDim(sce, "pca.corrected")[,1],
    y = reducedDim(sce, "pca.corrected")[,2],
    col = EmbryoCelltypeColours[colData(sce)$celltype],
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)

plot(
    x = reducedDim(tal1, "pca.corrected")[,1],
    y = reducedDim(tal1, "pca.corrected")[,2],
    col = "black",
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)

plot(
    x = reducedDim(tch, "pca.corrected.E7.5")[,1],
    y = reducedDim(tch, "pca.corrected.E7.5")[,2],
    col = "black",
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)

tal1_singlets <- which(!(colData(tal1)$celltype.mapped == "Doublet"))
xs1 <- reducedDim(sce, "pca.corrected")[singlets, 1]
ys1 <- reducedDim(sce, "pca.corrected")[singlets, 2]
xs2 <- reducedDim(tal1, "pca.corrected")[tal1_singlets, 1]
ys2 <- reducedDim(tal1, "pca.corrected")[tal1_singlets, 2]
xs <- c(xs1, xs2)
ys <- c(ys1, ys2)

colors <- c(
     1 * nrow(sce),
    2 * nrow(tal1)
)

plot(
    x = xs,
    y = ys,
    col = colors,
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
)


colors_celltype <- c(
    EmbryoCelltypeColours[colData(sce)$celltype[singlets]],
    EmbryoCelltypeColours[colData(tal1)$celltype.mapped[tal1_singlets]]
)
color_factors = as.factor(c(
    colData(sce)$celltype[singlets],
    colData(tal1)$celltype.mapped[tal1_singlets]
))

par(mar=c(15, 4, 4, 15), xpd=TRUE)
plot(
    x = xs,
    y = ys,
    col = colors_celltype,
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "PCA1", ylab = "PCA2"
)
legend_colors = as.vector(EmbryoCelltypeColours[levels(color_factors)])
legend("topright", legend = levels(color_factors), fill = legend_colors)
coord <- par("usr")
legend(x = coord[2] * 1.05, y = coord[4],legend = color_factors, fill = legend_colors, xpd = TRUE )

EmbryoCelltypeColours
scater::runUMAP(tch)
