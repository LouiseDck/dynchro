library(condimentsPaper)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(slingshot)
library(condiments)
library(tradeSeq)

# remove scores, slingshot, tradeSeq, crv

write_sce_parts <- function(sce, basefile, colData = TRUE){

    counts <- as(counts(sce), "sparseMatrix")
    counts_triples <- summary(counts)

    # logcounts <- as(logcounts(sce), "sparseMatrix")
    # logcounts_triples <- summary(logcounts)

    # rownames (genes)
    rownames <- rownames(counts)

    # colnames (cells, barcodes)
    colnames <- colnames(counts)

    if(colData){
        # metadata
        # colData names
        colDataNames <- colnames(colData(sce))
        colDataSce <- colData(sce)
        colDataSce_rm <- within(colDataSce, rm("scores", "slingshot", "tradeSeq", "crv"))
        write.csv(colDataSce_rm, file = glue("{basefile}_colData.csv"))
    }

    write.csv(rownames, file = glue("{basefile}_rownames.csv"))
    write.csv(colnames, file = glue("{basefile}_colnames.csv"))

    fwrite(list(counts_triples$i), sep = ",", file = glue("{basefile}_counts_triples_i.csv"))
    fwrite(list(counts_triples$j), sep = ",", file = glue("{basefile}_counts_triples_j.csv"))
    fwrite(list(counts_triples$x), sep = ",", file = glue("{basefile}_counts_triples_x.csv"))

    # fwrite(list(logcounts_triples$i), sep = ",", file = glue("{basefile}_logcounts_triples_i.csv"))
    # fwrite(list(logcounts_triples$j), sep = ",", file = glue("{basefile}_logcounts_triples_j.csv"))
    # fwrite(list(logcounts_triples$x), sep = ",", file = glue("{basefile}_logcounts_triples_x.csv"))

}

write_sls <- function(object, file){
  o2 <- as.data.frame(object)
  write.csv(o2, file = file)
}

write_csv(pathStats(sds_1)$pseudotime, "sds_1.csv")

write_sce_parts(kras, "kras")

theme_set(theme_classic())

Sys.setenv(VROOM_CONNECTION_SIZE=500072)
kras <- condimentsPaper::import_KRAS()
data("kras", package = "condimentsPaper")

rm1 <- readLines("examples/kras/kras1_remove.csv")
rm2 <- readLines("examples/kras/kras2_remove.csv")
kras <- kras[, !(colnames(kras) %in% rm1)]
kras <- kras[, !(colnames(kras) %in% rm2)]

cols <- c(brewer.pal(5, "Blues")[2:5], brewer.pal(3, "Greens")[2:3], brewer.pal(3, "Reds")[2:3],
          brewer.pal(3, "Oranges")[2], "Grey")
names(cols) <- c(3, 5, 4, 1, 8, 2, 9, 10, 6, 7)
kras$Cluster <- as.character(kras$Cluster)
reducedDim(kras, "TSNE") <- colData(kras)[, c("tSNE1", "tSNE2")] %>% as.matrix()
df <- colData(kras)[, 1:98] %>% as.data.frame() %>%
  sample_frac(1)
p1 <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = Batch)) +
  geom_point(size = .7) +
  scale_color_brewer(palette = "Accent") +
  labs(col = "Type")
p1

p2 <- ggplot(df, aes(x = tSNE1, y = tSNE2, fill = Cluster)) +
  geom_point(size = 1, alpha = .65, col = "grey70", shape = 21) +
  scale_fill_manual(values = cols) +
  labs(fill = "Cell Clusters")
p2

kras <- imbalance_score(Object = kras, conditions = "Batch", dimred = "TSNE")
df$scores <- kras[, df$X1]$scores$scaled_scores
p3 <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Score")
p3


kras <- slingshot(kras, reducedDim = "TSNE", clusterLabels = kras$Cluster,
                 start.clus = 7, extend = "n", reweight = FALSE, reassign = FALSE)
topologyTest(kras, conditions = "Batch", rep = 50)


p4 <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = Batch)) +
  geom_point(size = 1.5) +
  scale_color_brewer(palette = "Accent") +
  labs(col = "Pseudotime") +
  geom_path(data = slingCurves(kras, as.df = TRUE),
            col = "black", size = 1.5)
p4

sdss <- slingshot_conditions(kras, kras$Batch, approx_points = FALSE,
                             extend = "n", reweight = FALSE, reassign = FALSE)
sdss$condition_id <- names(sdss)
sdss$mapping <- matrix(rep(1:3, each = 3), nrow = 3, ncol = 3, byrow = TRUE)
sds <- do.call(merge_sds, sdss)

sds

df <- full_join(
  df %>% select(X1, tSNE1, tSNE2, Cluster, Batch) %>%
    dplyr::rename("cells" = "X1"),
  slingPseudotime(sds) %>% 
    as.data.frame() %>%
    mutate(cells = rownames(.))
) %>%
  pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "pst")

p4 <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = Batch)) +
  geom_point(size = .7, alpha = .1) +
  scale_color_brewer(palette = "Accent")
for (batch in unique(kras$Batch)) {
  sds_cond <- sdss[[batch]]
  for (i in 1:3) {
    p4 <- p4 +  
      geom_path(data = slingCurves(sds_cond)[[i]]$s[slingCurves(sds_cond)[[i]]$ord, ] %>%
                  as.data.frame() %>%
                  mutate(Batch = batch), 
                size = 1.5)   
  }
}
position <- data.frame(
  "tSNE1" = c(40, -30, 45),
  "tSNE2" = c(50, -50, -50),
  "Batch" = "H2122A",
  "text" = paste0("Lineage ", 1:3)
)
p4 <- p4 + 
  geom_text(data = position, col = "black", aes(label = text), size = 5)
p4

p <- ggplot(df, aes(x = tSNE1, y = tSNE2, col = Batch)) +
  geom_point(size = .7, alpha = .1) +
  scale_color_brewer(palette = "Accent")
cls <- df %>% group_by(Cluster, Batch) %>%
  dplyr::summarise(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2), .groups = NULL)

p <- p +
  geom_point(data = cls, size = 3)
edges <- lapply(slingLineages(sds), function(lin) {
  from <- lin[1:(length(lin) - 1)]
  to <- lin[2:length(lin)]
  return(data.frame("from" = from, "to" = to))
}) %>% bind_rows()
for (batch in unique(kras$Batch)) {
  cl_batch <- cls %>% filter(Batch == batch)
  edges_batch <- left_join(edges, cl_batch, by = c("from" = "Cluster")) %>%
    left_join(cl_batch %>% dplyr::rename("tSNE1_end" = "tSNE1",
                                  "tSNE2_end" = "tSNE2") %>%
                select(-Batch), by = c("to" = "Cluster") )
    p <- p +
     geom_segment(data = edges_batch, aes(xend = tSNE1_end, yend = tSNE2_end),
                  size = 2)
}
p

progressionTest(sds, conditions = kras$Batch, lineages = TRUE)

p5 <- ggplot(df, aes(x = pst)) +
  geom_density(alpha = .4, aes(fill = Batch), col = "transparent") +
  geom_density(aes(col = Batch), fill = "transparent", size = 1.5) +
  guides(col = "none") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  labs(x = "Pseudotime", fill = "Type") +
  facet_wrap(~ Curve, scales = "free_x")
p5

fateSelectionTest(sds, conditions = kras$Batch, pairwise = TRUE)

weights <- condiments:::.sling_reassign(sds)
df <- df %>%
  full_join(weights %>% 
              as.data.frame() %>%
              mutate(cells = rownames(.)) %>%
              dplyr::rename("Lineage1" = V1, "Lineage2" = V2, "Lineage3" = V3) %>%
              pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "weights")
      )

df_w <- df %>%
  select(-pst) %>%
  group_by(cells) %>%
  mutate(weights = weights / sum(weights)) %>%
  pivot_wider(names_from = "Curve", values_from = "weights")
p <- ggplot(df_w, aes(x = Lineage1, y = Lineage3)) +
  geom_hex() +
  scale_fill_viridis_c(direction = -1) +
  facet_wrap(~Batch, scales = "free") +
  geom_abline(slope = -1, intercept = 1, linetype = "dotted") +
  geom_abline(slope = -1, intercept = 2/3, linetype = "dotted") +
  geom_abline(slope = -1, intercept = 1/3, linetype = "dotted") +
  annotate("text", x = .53, y = .53, label = "w3 = 0", angle = -52) +
  annotate("text", x = .62, y = .1, label = "w3 = 1/3", angle = -52) +
  annotate("text", x = .14, y = .14, label = "w3 = 2/3", angle = -52) +
  theme(legend.position = "bottom") +
  labs(x = "Weights for Lineage 1 (w1)", y = "Weights for Lineage 2 (w2)",
       fill = "counts per hexagon")
p

df_w <- df %>%
  select(-pst) %>%
  group_by(cells) %>%
  mutate(weights = weights / sum(weights)) %>%
  ungroup() %>%
  group_by(Batch, Curve) %>%
  summarise(weights = mean(weights), .groups = NULL)

p2 <- ggplot(df_w, aes(x = Curve, fill = Batch, y = weights)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Accent") +
  theme(legend.position = c(.7, .7)) +
  labs(x = "", y = "Mean weight")
p2

set.seed(3)
filter <- apply(counts(kras), 1, function(g) {
    sum(g >= 5) >= 10
})
kras <- kras[filter, ]

set.seed(3)
library(BiocParallel)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4
icMat <- evaluateK(counts = as.matrix(assays(kras)$counts),
                   pseudotime = slingPseudotime(sds, na = FALSE),
                   cellWeights = weights,
                   conditions = factor(colData(kras)$Batch),
                   nGenes = 300,
                   k = 3:7,
                   parallel = TRUE,
                   BPPARAM = BPPARAM)

set.seed(3)
kras@int_metadata$slingshot <- sds
kras <- fitGAM(counts = kras,
               conditions = factor(colData(kras)$Batch),
               parallel = TRUE,
               BPPARAM = BPPARAM,
               nknots = 6)

condRes <- conditionTest(kras, l2fc = log2(2), lineages = TRUE)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

# plot genes
scales <- c("#7FC97F", "#57A059", "#2E7935",
            "#BEAED4", "#9687AC", "#706285",
            "#FDC086", "#CC945C", "#9D6A34")
scales <- scales[c(1, 4, 7, 2, 5, 8, 3, 6, 9)]
oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
p6 <- plotSmoothers(kras, counts(kras),
                    gene = rownames(assays(kras)$counts)[oo[1]],
                    alpha = 1, curvesCols = scales, sample = .3) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(kras)$counts)[oo[1]])

# Second most significant gene
p7 <- plotSmoothers(kras, assays(kras)$counts,
              gene = rownames(assays(kras)$counts)[oo[2]],
              alpha = 1, curvesCols = scales, sample = .3) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(kras)$counts)[oo[2]])

# least significant gene
p8 <- plotSmoothers(kras, assays(kras)$counts,
              gene = rownames(assays(kras)$counts)[oo[nrow(kras)]],
              alpha = 1, curvesCols = scales, sample = .3) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(kras)$counts)[oo[nrow(kras)]])

p6
p7
p8

### based on mean smoother
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
conditionGenes_lineage1 <- rownames(condRes)[condRes$padj_lineage1 <= 0.05]
conditionGenes_lineage1 <- conditionGenes_lineage1[!is.na(conditionGenes_lineage1)]

yhatSmooth <- predictSmooth(kras, gene = conditionGenes_lineage1, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth[, c(1:50, 151:200, 301:350)], 1, scales::rescale))
heatSmooth_H2122A <- pheatmap(yhatSmoothScaled[, 1:50],
  cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H2122A", legend = FALSE,
  silent = TRUE
)

matchingHeatmap_H358A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 51:100],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H358A",
  legend = FALSE, silent = TRUE
)

matchingHeatmap_SW1573A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 101:150],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "SW1573A",
  legend = FALSE, silent = TRUE
)

p9 <- plot_grid(heatSmooth_H2122A[[4]], matchingHeatmap_H358A[[4]], matchingHeatmap_SW1573A[[4]],
                NULL, NULL, NULL, ncol = 3, rel_widths = c(1.4, 1, 1), rel_heights = c(10, 1)) +
  draw_text("Lineage 1", x = .5, y = .05)
p9

condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
conditionGenes_lineage2 <- rownames(condRes)[condRes$padj_lineage2 <= 0.05]
conditionGenes_lineage2 <- conditionGenes_lineage1[!is.na(conditionGenes_lineage2)]

yhatSmooth <- predictSmooth(kras, gene = conditionGenes_lineage2, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth[, c(51:100, 201:250, 351:400)],1, scales::rescale))

heatSmooth_H2122A <- pheatmap(yhatSmoothScaled[, 1:50],
  cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H2122A", legend = FALSE,
  silent = TRUE
)

matchingHeatmap_H358A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 51:100],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H358A",
  legend = FALSE, silent = TRUE
)

matchingHeatmap_SW1573A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 101:150],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "SW1573A",
  legend = FALSE, silent = TRUE
)

p10 <- plot_grid(heatSmooth_H2122A[[4]], matchingHeatmap_H358A[[4]], matchingHeatmap_SW1573A[[4]],
                NULL, NULL, NULL, ncol = 3, rel_widths = c(1.4, 1, 1), rel_heights = c(10, 1)) +
  draw_text("Lineage 2", x = .5, y = .05)
p10

condRes$padj_lineage3 <- p.adjust(condRes$pvalue_lineage3, "fdr")
conditionGenes_lineage3 <- rownames(condRes)[condRes$padj_lineage3 <= 0.05]
conditionGenes_lineage3 <- conditionGenes_lineage3[!is.na(conditionGenes_lineage3)]

yhatSmooth <- predictSmooth(kras, gene = conditionGenes_lineage3, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth[, c(101:150, 251:300, 401:450)],1, scales::rescale))
heatSmooth_H2122A <- pheatmap(yhatSmoothScaled[, 1:50],
  cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H2122A", legend = FALSE,
  silent = TRUE
)

matchingHeatmap_H358A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 51:100],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "H358A",
  legend = FALSE, silent = TRUE
)

matchingHeatmap_SW1573A <- pheatmap(yhatSmoothScaled[heatSmooth_H2122A$tree_row$order, 101:150],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "SW1573A",
  legend = FALSE, silent = TRUE
)

p11 <- plot_grid(heatSmooth_H2122A[[4]], matchingHeatmap_H358A[[4]], matchingHeatmap_SW1573A[[4]],
                 NULL, NULL, NULL, ncol = 3, rel_widths = c(1.4, 1, 1), rel_heights = c(10, 1)) +
  draw_text("Lineage 3", x = .5, y = .05)
p11
