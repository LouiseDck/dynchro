require(ggplot2)
require(Seurat)
require(reshape2)
require(plyr)
require(DESeq2)
require(ggrepel)
require(SingleCellExperiment)
require(scmap)
require(parallel)
require(randomForest)
aml <- readRDS("Work/dynchro/examples/triana/AML.rds")
healthy <- readRDS("Work/dynchro/examples/triana/Healthy.rds")
Healthy <- healthy
reproduce <- TRUE

AMLs <- aml

#replace this line by a command that loads your seurart object. To make sure you use the same normalization that we use , run Seurat normalization. In this case, the antibody counts are CLR normalizaed and the RNA counts were Log Normalized.
#The object needs to have a metadata column Batch that encodes sample identity
if (!reproduce) {
  merged <- merge(Healthy, y = AMLs, add.cell.ids = c("reference", "query"), project = "abAML")
} else {
  merged <- AMLs
}
#by

normalized_counts <- GetAssayData(merged, assay = "BOTH", slot="data")
counts <- GetAssayData(merged, assay = "BOTH", slot="counts") 

#replace this line by a command that loads you patient-level meta data. metadata.tsv file is contained in this git directory.
#(row names should correspond to the Batch IDs used in the cell-level meta data, and there should be a column called Status that determines if the individuals are "healthy" or "diseased")
patient.metadata <- read.table("Work/dynchro/examples/triana/metadata.tsv", sep="\t",header=T, row.names = 1)

merged$status <- patient.metadata[as.character(merged$Batch),"Status"]

counts.projection <- normalized_counts[intersect(rownames(normalized_counts), rownames(Healthy)),]

#set up query data
sce_Culture <- SingleCellExperiment(assays = list(normcounts =  as.matrix(counts.projection)))
logcounts(sce_Culture) <- normcounts(sce_Culture)
rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)

sce_All <- SingleCellExperiment(assays = list(normcounts = as.matrix(Healthy@assays$BOTH@data[rownames(counts.projection),])))

logcounts(sce_All) <- normcounts(sce_All)
# use gene names as feature symbols
rowData(sce_All)$feature_symbol <- rownames(sce_All)
# remove features with duplicated names
sce_All <- sce_All[!duplicated(rownames(sce_All)), ]

sce_Culture<-setFeatures(sce_Culture,features =  rownames(sce_Culture))
sce_Culture <- indexCell(sce_Culture)

sce_All<-setFeatures(sce_All,features =  rownames(sce_All))
sce_All <- indexCell(sce_All)


Culture_Map <- scmapCell(
  projection = sce_Culture,
  index_list = list(
    sce_All = metadata(sce_All)$scmap_cell_index
  ),
  w = 5)

#it makes sense to save the result


Calc<-function(id,cult){
  u <- cult[,id]
  xcoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,1][u]
  ycoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,2][u]
  x=mean(xcoords)
  y=mean(ycoords)
  ct.t=table(Idents(Healthy)[u])
  ct.t<- ct.t[order(ct.t,decreasing = T)]
  nearest <- Healthy@assays$BOTH@data[rownames(sce_Culture),u]
  query <- normcounts(sce_Culture)[,id]
  cornn<- apply(nearest,2,function(x) cor(x, query))
  data.frame(row.names = id, x = x, y =y, ct = names(ct.t)[1], score = mean(cornn))
  
}


library(parallel)
mapped <- mclapply(colnames(Culture_Map$sce_All[[1]]), Calc, cult = Culture_Map$sce_All[[1]], mc.cores = 6)
mapped_ <- mapped
mapped <- do.call(rbind,mapped_)

Idents(merged) <- mapped$ct
mapped$sample <- merged$Batch

# mapped_amls <- subset(mapped, sample %in% rownames(patient.metadata)[patient.metadata$Status == "diseased"])
options <- rownames(patient.metadata)[patient.metadata$Status == "diseased"] %>% as.vector()
# options <- factor(options, levels = options)
# other <- mapped$sample %>% as.vector() %>% unique()
# mapped$sample <- factor(mapped$sample, levels = other)
mapped_amls <- mapped[mapped$sample %in% options,]

forcluster <- table(as.character(mapped_amls$ct), as.character(mapped_amls$sample))

#the following line needs to be modified to include only 
forcluster <- forcluster[apply(forcluster,1,sum)> 100 & !grepl("T cell|NK cell|B cell", rownames(forcluster)),]
forcluster <- t(t(forcluster)/apply(forcluster,2,sum))

cl <- pheatmap::pheatmap(forcluster, clustering_method = "ward.D2")

bckgr <- data.frame(x = Embeddings(Healthy, reduction = "MOFAUMAP")[,1],
                    y = Embeddings(Healthy, reduction = "MOFAUMAP")[,2])

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

mapped <- split(mapped, mapped$sample)
mapped <- lapply(mapped, function(x) {
  x$d <- get_density(x$x, x$y)
  x$dnorm <- x$d / max(x$d)
  x
})
mapped <- do.call(rbind, mapped)

qplot(x = x, y =y, data= bckgr, color = I("#BBBBBB"),size=I(0.1)) + geom_point(aes(color = dnorm),data =mapped,size=I(0.1)) + facet_wrap(~ sample, nrow=3) + scale_color_gradientn(colours = c("black","blue","red","orange")) +
  theme_bw() + theme(panel.grid = element_blank(), axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank()) 



if (reproduce) {
  Healthy$status <- "healthy"
  merged <- merge(Healthy, y = merged, add.cell.ids = c("reference", "query"), project = "abAML")
}

cat("Samples included into DE comparison (min 20 cells per sample):\n")

getTest <- function(id, assay = "AB") {
  ab <- GetAssayData(merged, assay = assay, slot = "counts")
  
  summed <-sapply(unique(merged$Batch), function(btc) {
    use <- merged$Batch == btc & Idents(merged)  == id
    if (sum(use)==1) ab[,merged$Batch  == btc &Idents(merged)  == id] else if (sum(use)==0) rep(0, nrow(ab)) else apply(ab[,merged$Batch  == btc &Idents(merged)  == id],1,sum)
  }) 
  ncells <- sapply(unique(merged$Batch), function(btc) sum(merged$Batch == btc & Idents(merged) == id))
  rownames(summed) <- rownames(ab)
  sums <- apply(summed,2,sum)
  summed <- summed[,ncells > 20]
  cat(id, ":", paste(colnames(summed),collapse = ","),"\n")
  
  #adjust to include your cavariates
  dds <- data.frame(row.names = colnames(summed), 
                    individual =colnames(summed),
                    status = factor(patient.metadata[colnames(summed),"Status"]),
                    age = scale(patient.metadata[colnames(summed),"Age"]),
                    gender = factor(patient.metadata[colnames(summed),"Gender"]))
  
  statuscounts <- table(dds$status)
  runaml <- statuscounts["diseased"] > 2 & statuscounts["healthy"] > 2
  
  tryCatch(
    {
      if (runaml) {
        #adjust to include your covariates
        deseq <- DESeqDataSetFromMatrix(summed, dds, ~  status + age + gender )
        deseq <- estimateSizeFactors(deseq)
        deseq <- estimateDispersions(deseq)
        deseq <- nbinomWaldTest(deseq)
        results <- as.data.frame(results(deseq, contrast = c("status","healthy","diseased")))
        results$gene <- rownames(results)
        results$ct <- id
        results$nAML <- statuscounts["diseased"]
        results$nHealthy <- statuscounts["healthy"] 
        aml <- T
        
        return(results)
      } else return(NULL)
      
      
      
      
    }, error = function(e) {warning("In ", id, e); NULL}
  )
  
}

#if you want to look at RNA differential expression, replace assay by RNA
all.tests.covariate <- lapply(unique(as.character(Idents(merged) )), getTest, assay="AB")
all.tests.covariate <- do.call(rbind, all.tests.covariate)

head(all.tests.covariate)


AMLs <- subset(merged, status == "diseased")
Healthy <- subset(merged, status != "diseased")

#train RF 
getRF <- function(ct) {
  AMLs.myeloid <- subset(AMLs, idents = ct)
  ab <- GetAssayData(AMLs.myeloid, assay = "AB", slot = "counts")
  ncells <- sapply(unique(AMLs.myeloid$Batch), function(btc) sum(AMLs.myeloid$Batch== btc & Idents(AMLs.myeloid) == ct))
  
  usepat <- unique(AMLs.myeloid$Batch)[ncells>20] #same fitlering criterion as above
  if (length(usepat) <= 2) return(NULL)
  AMLs.myeloid <- subset(AMLs.myeloid, Batch %in% usepat)
  pat <- factor(AMLs.myeloid$Batch)
  pred <- t(GetAssayData(AMLs.myeloid, "AB", slot = "data"))
  
  if (length(unique(pat)) > 1) {
    rf <- randomForest(as.matrix(pred), y = pat, importance = T)
    data.frame(row.names = rownames(rf$importance), ct = ct, gene = rownames(rf$importance), 
               MeanDecreaseAccuracy = rf$importance[,"MeanDecreaseAccuracy"],
               MeanDecreaseGini = rf$importance[,"MeanDecreaseGini"],
               MeanDecreaseAccuracySD = rf$importanceSD[,"MeanDecreaseAccuracy"],
               norm =  rf$importance[,"MeanDecreaseAccuracy"] / max(rf$importance[,"MeanDecreaseAccuracy"]))  
  } else NULL
  
  
}

ct <- Idents(AMLs)
usect <- table(ct)
usect <- names(usect)[usect > 10]

rfs <- lapply(usect, getRF)
rfs <- do.call(rbind, rfs)

showcelltypes <- c("HSCs & MPPs", "Myelocytes")
usecounts.AMLs <- GetAssayData(AMLs, assay = "AB", slot="counts")
usecounts.Healthy <- GetAssayData(Healthy, assay = "AB", slot="counts")


means <-lapply(showcelltypes, function(x) {
  data.frame(ct = x,
             aml = apply(usecounts.AMLs[,Idents(AMLs)  == x], 1, mean),
             healthy = apply(usecounts.Healthy[,Idents(Healthy)  == x ], 1, mean),
             gene  = rownames(usecounts.Healthy))
  
})

means <- do.call(rbind, means)

sig_means <- merge(means, subset(all.tests.covariate, padj < 0.1 & ct %in% showcelltypes))
means <- merge(means, rfs)
complete <- merge(means, sig_means, all=T)

qplot(x = healthy, y = aml, size=-log10(padj), data =complete, log="xy" , color = norm)+facet_wrap(~ ct) + 
  geom_text_repel(size = 3, aes(label = gsub("-AB","",gene)),data = subset(sig_means, padj < 0.0001), color ="black")+ theme_bw() + theme(panel.grid = element_blank()) + 
  scale_color_gradientn(colors = c("grey","black","blue","red","orange"),name = "Inter-patient\nvariability") + scale_size_area(name = "-log10 p\nDE", max_size = 3,na.value = 0.5) + scale_shape_discrete(name = "LSC marker\n(Hanekamp 2017)") +
  xlab("Counts in healthy individuals") + ylab("Counts in leukemic individuals")

