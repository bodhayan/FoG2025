# Code for Festival of Genomics event
setwd("C:/Users/bodha/Desktop/FoG WS/")
rm(list = ls())

library(Seurat)
library(tidyverse)
library(dplyr)

# Plot 
library(ggplot2)
library(gridExtra)
library(ggrepel)

# Dataset: https://www.10xgenomics.com/datasets/10k-bone-marrow-mononuclear-cells-bmmncs-5-v2-0-without-intronic-reads-2-standard
expr_mat <- Read10X_h5('10k_BMMNC_5pv2_nextgem_10k_BMMNC_5pv2_nextgem_count_sample_feature_bc_matrix.h5')

seu_obj = CreateSeuratObject(counts = expr_mat,
                             min.cells = 1,
                             project = "bmmnc")

seu_obj # 22715 features across 8288 samples



# Data1: https://www.10xgenomics.com/datasets/frozen-bmm-cs-healthy-control-1-1-standard-1-1-0
# Data2: https://www.10xgenomics.com/datasets/frozen-bmm-cs-healthy-control-2-1-standard-1-1-0


donor1_dir <- 'frozen_bmmc_healthy_donor1_filtered_gene_bc_matrices/filtered_matrices_mex/hg19'
donor2_dir <- 'frozen_bmmc_healthy_donor2_filtered_gene_bc_matrices/filtered_matrices_mex/hg19'

list.files(donor1_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
list.files(donor2_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx

donor1_expression_matrix <- Read10X(data.dir = donor1_dir)
donor2_expression_matrix <- Read10X(data.dir = donor2_dir)

donor1_seu_obj = CreateSeuratObject(counts = donor1_expression_matrix,
                                    project = "donor1")
donor2_seu_obj = CreateSeuratObject(counts = donor2_expression_matrix,
                                    project = "donor2")

donor1_seu_obj # 32738 features across 1985 samples
donor2_seu_obj # 32738 features across 2472 samples


# FILT ------------
donor1_seu_obj = CreateSeuratObject(counts = donor1_expression_matrix,
                                    min.features = 200,
                                    min.cells = 1,
                                    project = "donor1")
donor2_seu_obj = CreateSeuratObject(counts = donor2_expression_matrix,
                                    min.features = 200,
                                    min.cells = 1,
                                    project = "donor2")

donor1_seu_obj # 16259 features across 1985 samples
donor2_seu_obj # 16248 features across 2471 samples



# merge datasets
merged_seurat <- merge(donor1_seu_obj,
                       y = donor2_seu_obj,
                       add.cell.ids = c('donor1','donor2'),
                       project = 'proj')

View(merged_seurat@meta.data)
table(merged_seurat$orig.ident)



# % MT reads
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
View(merged_seurat@meta.data)

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# Filtering -----------------
merged_seurat # 17737 features across 4456 samples
merged_seurat <- subset(merged_seurat, subset = percent.mt < 25) 
merged_seurat # 17737 features across 4451 samples
table(merged_seurat$orig.ident)

# Tranform and PCA pipeline (needed for integrated join) -----------------------
#merged_seurat <- NormalizeData(object = merged_seurat) %>%
#                 FindVariableFeatures() %>%
#                 ScaleData() %>%
#                 RunPCA()

merged_seurat <- SCTransform(merged_seurat) %>%
  RunPCA()

# determine dimensionality of the data
ElbowPlot(merged_seurat)

# Merge with integration
merged_obj <- IntegrateLayers(object = merged_seurat,
                              method = HarmonyIntegration, # CCAIntegration, RPCAIntegration
                              assay = "SCT",
                              orig.reduction = "pca",
                              new.reduction = "harmony",
                              verbose = FALSE)


# determine dimensionality of the data
ElbowPlot(merged_obj)

# Clustering and dim redn
merged_obj <- FindNeighbors(merged_obj,
                            reduction = "harmony",
                            dims = 1:30) %>%
  FindClusters(resolution = c(0.4,0.6,0.8,1.0,1.2,1.4)) %>% # parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
  RunUMAP(reduction = "harmony",
          dims = 1:30)

# plot
p1 <- DimPlot(merged_obj,
              #reduction = 'umap',
              group.by = 'orig.ident',
              cols = c('indianred','lavender','lightgreen','dodgerblue')) + ggtitle("Integrated")
p1

# Clustering and dim redn - in original dataset
merged_seurat <- FindNeighbors(merged_seurat,
                               reduction = "pca",
                               dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca",
          dims = 1:30)

p2 <- DimPlot(merged_seurat,
              #reduction = 'umap',
              group.by = 'orig.ident',
              cols = c('indianred','lavender','lightgreen','dodgerblue')) + ggtitle("As it is")
p2

grid.arrange(p2, p1, ncol = 2, nrow = 1)

# Macrophage detection through AUCell ##########################################

FeaturePlot(merged_obj, 'CSF1R') +
FeaturePlot(merged_obj, 'ITGAM') # Cd11b

genes <- c("ITGAM", "CSF1R")
VlnPlot(merged_obj, features = genes)+ DimPlot(merged_obj, label = T)

# PagloDB
mph_mk <- read.csv2('PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t')
unique(mph_mk$species)
mm <- c('Mm Hs', 'Hs')
mm_mk <- filter(mph_mk, species %in% mm)
mm_mph <- mm_mk$official.gene.symbol

# CellMarker
# https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
#ACT http://biocc.hrbmu.edu.cn/ACT/download.jsp
mm_bm <- read.csv2('Human_cell_markers.txt', sep = '\t')
View(mm_bm)
mph_mk <- filter(mm_bm, cellName == 'Macrophage' & tissueType == 'Bone marrow' & cancerType == 'Normal')[1,]$cellMarker
mph_mk
mph_mk<-str_split(mph_mk, ", ")[[1]]
mph_mk

# AUCell on CellMarker
library(AUCell)
exprMatrix <- LayerData(merged_obj, assay = "SCT", layer = "counts")
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(mph_mk, cells_rankings)
set.seed(333)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
cells_assignment$geneSet$aucThr$thresholds
th <- cells_assignment$geneSet$aucThr$thresholds[1,1]
th
#th <- 0.5
tot_cells <- names(which(getAUC(cells_AUC)["geneSet",]>=0))
length(tot_cells) # 58857
new_cells <- names(which(getAUC(cells_AUC)["geneSet",]>th))
length(new_cells) # 30655
merged_obj$is_mph <- ifelse(colnames(merged_obj) %in% new_cells, "mph", "non_mph")
View(merged_obj@meta.data)
DimPlot(merged_obj, group.by = "is_mph")

merged_obj <- SetIdent(merged_obj, value = 'is_mph')
VlnPlot(merged_obj, features = genes, idents = 'is_mph')
VlnPlot(merged_obj, genes)

table(merged_obj$is_mph, merged_obj$seurat_clusters)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.0.4)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.0.6)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.0.8)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.1)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.1.2)
table(merged_obj$is_mph, merged_obj$SCT_snn_res.1.4)

a <- as.data.frame(table(merged_obj$seurat_clusters, merged_obj$is_mph))
b <- reshape(a, idvar = "Var1", timevar = "Var2", direction = "wide")
b$tot <- b$Freq.mph + b$Freq.non_mph
b$mph <- b$Freq.mph/b$tot

# CellType Identification ###################################################
library(SingleR)
# get reference data -----------
# https://stackoverflow.com/questions/77370659/error-failed-to-collect-lazy-table-caused-by-error-in-db-collect-using
# devtools::install_version("dbplyr", version = "2.3.4") # downgrade dbplyr
ref <- celldex::MonacoImmuneData() # Use
#refm <- scRNAseq::GrunHSCData() 
View(as.data.frame(colData(ref)))
#View(as.data.frame(colData(refm)))

pred <- SingleR(test = exprMatrix,
                ref = ref,
                labels = ref$label.main)
pred
merged_obj$label <- pred$labels[match(rownames(merged_obj@meta.data), rownames(pred))]
DimPlot(merged_obj,
        group.by = 'label',
        label = T)

# Annotation diagnostics ----------
# ...Based on the scores within cells -----------
pred
pred$scores
plotScoreHeatmap(pred)
plotScoreHeatmap(pred, cluster_rows = F, cluster_cols = F)
ct <- sort(unique(pred$labels))
plotScoreHeatmap(pred, rows.order = ct)
# ...Based on deltas across cells ----------
plotDeltaDistribution(pred)

# Choosing the right clustering
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.0.4)
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.0.6)
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.0.8)
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.1)
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.1.2)
janitor::tabyl(merged_obj@meta.data,
               is_mph,
               SCT_snn_res.1.4)

p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.0.4, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r, labels = p0.4$SCT_snn_res.0.4)
p1<-plot(p0.4$SCT_snn_res.0.4, p0.4$r)
p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.0.6, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r)
p2<-plot(p0.4$SCT_snn_res.0.6, p0.4$r)
p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.0.8, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r)
p3<-plot(p0.4$SCT_snn_res.0.8, p0.4$r)
p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.1, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r)
p4<-plot(p0.4$SCT_snn_res.1, p0.4$r)
p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.1.2, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r)
p5<-plot(p0.4$SCT_snn_res.1.2, p0.4$r)
p0.4 <- janitor::tabyl(merged_obj@meta.data, SCT_snn_res.1.4, is_mph)
p0.4$r <- 100*p0.4$mph/(p0.4$mph+p0.4$non_mph)
#beeswarm::beeswarm(p0.4$r)
p6<-plot(p0.4$SCT_snn_res.1.4, p0.4$r)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
p1+p2+p3+p4+p5+p6

janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.0.4)
janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.0.6)
janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.0.8)
janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.1)
janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.1.2)
janitor::tabyl(merged_obj@meta.data,
               label,
               SCT_snn_res.1.4)

library(pheatmap)
a<-table(merged_obj$label,merged_obj$SCT_snn_res.0.4)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))
a<-table(merged_obj$label,merged_obj$SCT_snn_res.0.6)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))
a<-table(merged_obj$label,merged_obj$SCT_snn_res.0.8)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))
a<-table(merged_obj$label,merged_obj$SCT_snn_res.1)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))
a<-table(merged_obj$label,merged_obj$SCT_snn_res.1.2)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))
a<-table(merged_obj$label,merged_obj$SCT_snn_res.1.4)
pheatmap(log10(a+10), color = colorRampPalette(c('white','dodgerblue'))(10))

p <- janitor::tabyl(merged_obj@meta.data,
                    label,
                    is_mph)
p$r <- 100*p$mph/(p$mph+p$non_mph)
View(p)

janitor::tabyl(merged_obj@meta.data,
               SCT_snn_res.1.4,
               SCT_snn_res.0.4)


# ...Comparing to unsupervised clustering ------------
tab <- table(Assigned=pred$labels, Clusters=merged_obj$seurat_clusters)
library(pheatmap)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','dodgerblue'))(10))


# Save an object to a file #####################################################
saveRDS(merged_obj, file = "mer_obj.rds")

rm(list=ls())
################################################################################
# Restore the object
obj <- readRDS(file = "mer_obj.rds")
DimPlot(obj)
DimPlot(obj, split.by = 'orig.ident')

p1 <- DimPlot(obj,
              group.by = 'SCT_snn_res.0.4',
              label = T)

p2 <- DimPlot(obj,
              group.by = 'label',
              label = T) 
p3 <- DimPlot(obj,
              group.by = 'is_mph',
              label = T) 
p4 <- DimPlot(obj,
              group.by = 'orig.ident',
              label = T)
grid.arrange(p1, p2, p3,p4, ncol = 2, nrow = 2)

View(obj@meta.data)

obj[["RNA"]] <- JoinLayers(obj) #merging layer

# Observing distribution
# 0.4
dict_0p4 <- as.data.frame(table(obj$SCT_snn_res.0.4, obj$is_mph))
dict_0p4 <- reshape(dict_0p4, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_0p4$tot <- dict_0p4$Freq.mph+dict_0p4$Freq.non_mph
dict_0p4$mpc_pc <- 100*dict_0p4$Freq.mph/dict_0p4$tot
plot(dict_0p4$Var1, dict_0p4$mpc_pc)
dict_0p4$mph <- ifelse(dict_0p4$mpc_pc > 50, 1,0)
group_by(dict_0p4, mph) %>% summarise(dist = sum(tot))
# 0.6
dict_0p6 <- as.data.frame(table(obj$SCT_snn_res.0.6, obj$is_mph))
dict_0p6 <- reshape(dict_0p6, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_0p6$tot <- dict_0p6$Freq.mph+dict_0p6$Freq.non_mph
dict_0p6$mpc_pc <- 100*dict_0p6$Freq.mph/dict_0p6$tot
plot(dict_0p6$Var1, dict_0p6$mpc_pc)
dict_0p6$mph <- ifelse(dict_0p6$mpc_pc > 50, 1,0)
group_by(dict_0p6, mph) %>% summarise(dist = sum(tot))
# 0.8
dict_0p8 <- as.data.frame(table(obj$SCT_snn_res.0.8, obj$is_mph))
dict_0p8 <- reshape(dict_0p8, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_0p8$tot <- dict_0p8$Freq.mph+dict_0p8$Freq.non_mph
dict_0p8$mpc_pc <- 100*dict_0p8$Freq.mph/dict_0p8$tot
plot(dict_0p8$Var1, dict_0p8$mpc_pc)
dict_0p8$mph <- ifelse(dict_0p8$mpc_pc > 50, 1,0)
group_by(dict_0p8, mph) %>% summarise(dist = sum(tot))
# 1.0
dict_1p0 <- as.data.frame(table(obj$SCT_snn_res.1, obj$is_mph))
dict_1p0 <- reshape(dict_1p0, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_1p0$tot <- dict_1p0$Freq.mph+dict_1p0$Freq.non_mph
dict_1p0$mpc_pc <- 100*dict_1p0$Freq.mph/dict_1p0$tot
plot(dict_1p0$Var1, dict_1p0$mpc_pc)
dict_1p0$mph <- ifelse(dict_1p0$mpc_pc > 50, 1,0)
group_by(dict_1p0, mph) %>% summarise(dist = sum(tot))
# 1.2
dict_1p2 <- as.data.frame(table(obj$SCT_snn_res.1.2, obj$is_mph))
dict_1p2 <- reshape(dict_1p2, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_1p2$tot <- dict_1p2$Freq.mph+dict_1p2$Freq.non_mph
dict_1p2$mpc_pc <- 100*dict_1p2$Freq.mph/dict_1p2$tot
plot(dict_1p2$Var1, dict_1p2$mpc_pc)
dict_1p2$mph <- ifelse(dict_1p2$mpc_pc > 50, 1,0)
group_by(dict_1p2, mph) %>% summarise(dist = sum(tot))
# 1.4
dict_1p4 <- as.data.frame(table(obj$SCT_snn_res.1.4, obj$is_mph))
dict_1p4 <- reshape(dict_1p4, idvar = "Var1", timevar = "Var2", direction = "wide")
dict_1p4$tot <- dict_1p4$Freq.mph+dict_1p4$Freq.non_mph
dict_1p4$mpc_pc <- 100*dict_1p4$Freq.mph/dict_1p4$tot
plot(dict_1p4$Var1, dict_1p4$mpc_pc)
dict_1p4$mph <- ifelse(dict_1p4$mpc_pc > 40, 1,0)
group_by(dict_1p4, mph) %>% summarise(dist = sum(tot))

library(ggbeeswarm)
v1 <- ggplot(dict_0p4, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v2 <- ggplot(dict_0p6, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v3 <- ggplot(dict_0p8, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v4 <- ggplot(dict_1p0, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v5 <- ggplot(dict_1p2, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v6 <- ggplot(dict_1p4, aes(x=1, y=mpc_pc)) + geom_violin() + geom_beeswarm(groupOnX = FALSE) + geom_hline(yintercept = c(40, 50))
v1+v2+v3+v4+v5+v6

obj$res.0.4 <- obj$SCT_snn_res.0.4 # Saving cluster

# Create cluster-wise MPH label
mph_clust <- dict_0p4[dict_0p4$mph == 1,]$Var1
obj$res.0.4_mph <- ifelse(obj$SCT_snn_res.0.4 %in% mph_clust, "mph", "others")

DimPlot(obj, group.by = "res.0.4_mph") + DimPlot(obj, group.by = "is_mph") +
  DimPlot(obj, group.by = "res.0.4_mph")+ DimPlot(obj, group.by = "SCT_snn_res.0.4", label = T) 

# Extract mph
DimPlot(obj, group.by = "orig.ident")
mph <- subset(obj, subset = res.0.4_mph == "mph")
mph <- SCTransform(mph) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunTSNE(dims = 1:30)

sort(unique(mph$res.0.4))
View(mph@meta.data)
DimPlot(mph, group.by = "res.0.4", reduction = 'tsne', label = T, pt.size = 1) +
  DimPlot(mph, group.by = "orig.ident", reduction = 'tsne', label = T, pt.size = 1)
DimPlot(mph, reduction = 'tsne', label = T)


table(obj$orig.ident, obj$SCT_snn_res.0.4)
DimPlot(obj, group.by = "orig.ident") + DimPlot(obj, group.by = "SCT_snn_res.0.4")



# DE genes - muskat ##########################################################################
# https://www.10xgenomics.com/analysis-guides/differential-gene-expression-analysis-in-scrna-seq-data-between-conditions-with-biological-replicates
# https://www.bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(SingleCellExperiment)

# check
library(ExperimentHub)
eh <- ExperimentHub()
query(eh, "Kang")
(sce <- eh[["EH2259"]])

# obj.sce <- as.SingleCellExperiment(obj) # Not working
object <-obj
object[["RNA"]] <- as(object[["RNA"]], Class="Assay")
sce <- as.SingleCellExperiment(object)

# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

library(sctransform)
assays(sce)$vstresiduals <- vst(counts(sce), verbosity = FALSE)$y

sce$id <- colnames(sce)
(sce <- prepSCE(sce, 
                kid = "label", # subpopulation assignments
                gid = "orig.ident",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

# compute UMAP using 1st 20 PCs
sce <- runTSNE(sce, pca = 20)
sce <- runUMAP(sce, pca = 20)

# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (dr in c("TSNE", "UMAP"))
  for (col in c("cluster_id", "group_id"))
    .plot_dr(sce[, cs100], dr, col)

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

(pb_mds <- pbMDS(pb))

# use very distinctive shaping of groups & change cluster colors
pb_mds <- pb_mds + 
  scale_shape_manual(values = c(17, 4)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# change point size & alpha
pb_mds$layers[[1]]$aes_params$size <- 5
pb_mds$layers[[1]]$aes_params$alpha <- 0.6
pb_mds

# run DS analysis
res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("stim-ctrl", levels = mm)

# run DS analysis
pbDS(pb, design = mm, contrast = contrast)

# 1st approach
mm <- mmDS(sce, method = "dream",
           n_cells = 10, n_samples = 2,
           min_counts = 1, min_cells = 20)

# 2nd & 3rd approach
mm <- mmDS(sce, method = "vst", vst = "sctransform")
mm <- mmDS(sce, method = "nbinom")

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# view top 2 hits in each cluster
top2 <- bind_rows(lapply(tbl_fil, top_n, 2, p_adj.loc))
format(top2[, -ncol(top2)], digits = 2)

frq <- calcExprFreqs(sce, assay = "counts", th = 0)
# one sheet per cluster
assayNames(frq)

# expression frequencies in each
# sample & group; 1st cluster
t(head(assay(frq), 5))

gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(frq)), 
                function(u) apply(u[, gids] > 0.1, 1, any), 
                logical(nrow(sce)))
t(head(frq10))

tbl_fil2 <- lapply(kids, function(k)
  dplyr::filter(tbl_fil[[k]], 
                gene %in% names(which(frq10[, k]))))

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil2, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# tidy format; attach pre-computed expression frequencies
resDS(sce, res, bind = "row", frq = frq)

# big-table (wide) format; attach CPMs
resDS(sce, res, bind = "col", cpm = TRUE)

# compute expression frequencies on the fly
resDS(sce, res, frq = TRUE)

library(UpSetR)
de_gs_by_k <- map(tbl_fil, "gene")
upset(fromList(de_gs_by_k))

# pull top-8 DS genes across all clusters
top8 <- bind_rows(tbl_fil) %>% 
  slice_min(p_adj.loc, n = 8, 
            with_ties = FALSE) %>% 
  pull("gene")

# for ea. gene in 'top8', plot t-SNE colored by its expression 
ps <- lapply(top8, function(g)
  .plot_dr(sce[, cs100], "TSNE", g) + 
    ggtitle(g) + theme(legend.position = "none"))

# arrange plots
plot_grid(plotlist = ps, ncol = 4, align = "vh")

plotExpression(sce[, sce$cluster_id == "B cells"],
               features = tbl_fil$`B cells`$gene[seq_len(6)],
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# top-5 DS genes per cluster
pbHeatmap(sce, res, top_n = 5)

# top-20 DS genes for single cluster
pbHeatmap(sce, res, k = "B cells")

# single gene across all clusters
pbHeatmap(sce, res, g = "ISG20")

# Amy's analyses ###############################################################################

cts=as.matrix(read.csv("C:/Users/bodha/OneDrive - University of Glasgow/Tommy/Joana_Amy/05_Analysis/noNeut/counts_noNeut_noRedBlood_2020-11-20.csv", check.names = F, row.names = 1))

cellsOfInterest=data.frame("cell"=colnames(cts))
coldata = read.csv("C:/Users/bodha/OneDrive - University of Glasgow/Tommy/Joana_Amy/05_Analysis/coldata_noNeut2020-01-21.csv") %>% dplyr::filter(cell %in% cellsOfInterest$cell )
coldata=as.matrix(coldata)

print("Are the rownames of coldata the same as the colnames in the count matrix?")
all(rownames(coldata) == colnames(cts))

library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(counts=cts), 
  colData = coldata
)
sce

#first check how many genes would remain
ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)
sce <- sce[keep,]

is.spike <- grepl("^ERCC-", rownames(sce))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
altExpNames(sce)

library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

mito <- which(rowData(sce)$CHR=="chrM")

par(mfrow=c(1,2))
hist(counts(sce), breaks=40, main="", col="grey80", ylim = c(0,500), xlab="counts")
hist(log((counts(sce))), breaks=40, main="", col="grey80", ylim = c(0,100000), xlab="logcounts")

df <- perCellQCMetrics(sce, subsets=list(Mito=mito))
colnames(df)

#adding cell names to this vector, so I can individually ID cells
rownames(df)=as.character(sce@colData$cell)

cellsDet=data.frame("cell"=rownames(df),
                    "features"=df$detected,
                    "ERCC"=df$altexps_ERCC_sum,
                    "mitoc_genes"=df$subsets_Mito_sum)

par(mfrow=c(2,2))
hist(df$subsets_Mito_percent, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df$altexps_ERCC_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df$sum/1e6, xlab="sum of counts (^1e6)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(df$detected, xlab="features detected", main="",
     breaks=20, col="grey80", ylab="Number of cells")

par(mfrow=c(1,2))
plot(df$altexps_ERCC_percent,df$subsets_Mito_percent, ylab="Mitochondrial proportion (%)",
     xlab="ERCC proportion (%)",main="", col="grey50")
plot(df$detected,df$subsets_Mito_percent, ylab="Mitochondrial proportion (%)",
     xlab="features detected",main="", col="grey50")
qc.param <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent",
                                                 "altexps_ERCC_percent"), nmads=4)
colSums(as.matrix(qc.param))

colData(sce) <- cbind(colData(sce), df)
sce$plate <- factor(sce$plate)
sce$condition <- ifelse(grepl("control", sce$condition),
                        "control", "CML")
sce$discard <- qc.param$discard
gridExtra::grid.arrange(
  plotColData(sce, x="plate", y="sum", colour_by="discard",
              other_fields="condition", shape_by = "condition") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="plate", y="detected", colour_by="discard", 
              other_fields="condition", shape_by = "condition") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="plate", y="subsets_Mito_percent", 
              colour_by="discard", other_fields="condition", shape_by = "condition") + 
    ggtitle("Mito percent"),
  plotColData(sce, x="plate", y="altexps_ERCC_percent", 
              colour_by="discard", other_fields="condition", shape_by = "condition") + 
    ggtitle("ERCC percent"),
  ncol=2
)

plotColData(sce, x="sum", y="subsets_Mito_percent", 
            colour_by="discard", other_fields=c("plate", "condition"), shape_by = "plate") +
  geom_text(aes(label=ifelse(df$subsets_Mito_percent>20,as.character(rownames(df)),'')),hjust=0,vjust=0, size=2) +
  facet_grid(~condition) +
  theme(panel.border = element_rect(color = "grey"))

plotColData(sce, x="altexps_ERCC_percent", y="subsets_Mito_percent",
            colour_by="discard", other_fields=c("plate", "condition"), shape_by = "plate") + 
  geom_text(aes(label=ifelse(df$subsets_Mito_percent>20,as.character(rownames(df)),'')),hjust=0,vjust=0, size=2) +
  facet_grid(~condition) + 
  theme(panel.border = element_rect(color = "grey"))

filt_sce <- sce[,!qc.param$discard]

print(paste(as.character(length(sce@colData$cell)), as.character(length(filt_sce@colData$cell))))

df_filtered <- perCellQCMetrics(filt_sce, subsets=list(Mito=mito))
#adding cell names to this vector, so I can individually ID cells
rownames(df_filtered)=as.character(filt_sce@colData$cell)
cellsDetFilt=data.frame("cell"=rownames(df_filtered),
                        "features"=df_filtered$detected,
                        "ERCC"=df_filtered$altexps_ERCC_sum,
                        "mitoc_genes"=df_filtered$subsets_Mito_sum)

par(mfrow=c(2,2))
hist(df_filtered$subsets_Mito_percent, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df_filtered$altexps_ERCC_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df_filtered$sum/1e6, xlab="sum of counts (^1e6)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(df_filtered$detected, xlab="features detected", main="",
     breaks=20, col="grey80", ylab="Number of cells")

logcounts(filt_sce) <- log(counts(filt_sce)+1)

filt_sce <- scater::runPCA(filt_sce)

set.seed(2302)
filt_sce <- scater::runTSNE(filt_sce, perplexity=50)

reducedDims(filt_sce)

filt.counts=counts(filt_sce)

library(scran)
dec <- modelGeneVar(filt_sce)#, block=assignments$phases)
par(mfrow=c(1,1))
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")

output <- getClusteredPCs(reducedDim(filt_sce))
#npcs <- metadata(output)$chosen
npcs <- output@metadata$chosen #Bodhayan by commenting previous
reducedDim(filt_sce, "PCAsub") <- reducedDim(filt_sce, "PCA")[,1:npcs,drop=FALSE]
npcs

g <- buildSNNGraph(filt_sce, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership
filt_sce$cluster <- factor(cluster)
table(filt_sce$cluster)

set.seed(1)
filt_sce <- runTSNE(filt_sce, dimred="PCAsub")
plotTSNE(filt_sce, colour_by="cluster", shape_by="condition")

# Change ofcolor
new_pallete <- c("skyblue3", "red3", "salmon", "lavender", "skyblue", "indianred3")
p <- plotTSNE(filt_sce, colour_by="cluster", shape_by="condition", point_size=3)
p  + scale_color_manual(values = new_pallete) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))

new_pallete <- c("skyblue3", "red3", "salmon", "lavender", "skyblue", "indianred3")
p1 <- plotTSNE(filt_sce, colour_by="Cd36", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p2 <- plotTSNE(filt_sce, colour_by="Lgals1", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p3 <- plotTSNE(filt_sce, colour_by="Cd14", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p1+p2+p3

set.seed(100)
filt_sce <- runUMAP(filt_sce, dimred="PCAsub")
plotUMAP(filt_sce, colour_by="cluster", shape_by="condition")

plotTSNE(filt_sce, colour_by="cluster", shape_by="plate")

plotTSNE(filt_sce, colour_by="cluster", shape_by="plate", point_size=3) + scale_color_manual(values = new_pallete) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))

set.seed(2302)
plotTSNE(filt_sce, colour_by="Cd36", shape_by="condition")

ratio <- bluster::pairwiseModularity(g, cluster, as.ratio=TRUE)
library(pheatmap)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))

ass.prob <- bluster::bootstrapStability(filt_sce, FUN=function(x) {
  g <- buildSNNGraph(x, use.dimred="PCAsub")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce$cluster)
pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=colorRampPalette(c("white", "blue"))(100))

markers <- findMarkers(filt_sce, filt_sce$cluster)
markers.table<-markers[[1]]

top.markers <- rownames(markers.table)[markers.table$Top <= 10]
plotHeatmap(filt_sce, features=top.markers, columns=order(filt_sce$cluster),
            colour_columns_by=c("cluster", "condition"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5), 
            show_colnames = FALSE) 

library(gplots)
ann_colors = list(cluster = c("1"="skyblue3","2"="red3","3"="salmon","4"="lavender","5"="skyblue", "6"="indianred3"),
                  condition = c("control" = "lightblue3", "CML" = "pink"))
plotHeatmap(filt_sce, features=top.markers, color_columns_by = c("cluster","condition"), column_annotation_colors = ann_colors, 
            color = bluered(256), center = T, cex=0.9, 
            columns = order(filt_sce$cluster), cluster_cols = F, zlim=c(-5, 5))

# Saving Amy obj locally
saveRDS(filt_sce, file = "amy_sce_obj.rds")
amy_obj = readRDS("amy_sce_obj.rds")
plotTSNE(amy_obj) # Check

# Analyses in Seurat #########################################################################
assayNames(filt_sce)
amy_meta <- as.matrix(colData(filt_sce))
amy_seu <- CreateSeuratObject(counts = counts(filt_sce))
View(amy_seu@meta.data)
all(rownames(amy_meta) == rownames(amy_seu@meta.data))
amy_seu$cell = rownames(amy_seu@meta.data)
amy_seu@meta.data <- right_join(amy_seu@meta.data, amy_meta, by="cell", copy = T)
rownames(amy_seu@meta.data) <- amy_seu$cell

# Basic pipeline
amy_seu <- SCTransform(amy_seu) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:7, method ="igraph") %>%
  FindClusters(resolution = 1.5) %>%
  RunTSNE(dims = 1:7)
DimPlot(amy_seu, label = T, pt.size = 3) + DimPlot(amy_seu, label = T, pt.size = 3, group.by = "cluster")

DimPlot(amy_seu, label = T, pt.size = 3) +
  DimPlot(amy_seu, label = T, pt.size = 3, group.by = "cluster") +
  DimPlot(amy_seu, label = T, pt.size = 3, group.by = "condition") +
  DimPlot(amy_seu, label = T, pt.size = 3, group.by = "plate")

table(amy_seu$seurat_clusters, amy_seu$cluster)

# Intergrating for plate
amy_seu[["RNA"]] <- split(amy_seu[["RNA"]], f = amy_seu$plate)
amy_seu
amy_seu <- IntegrateLayers(object = amy_seu,
                           method = HarmonyIntegration,
                           orig.reduction = "pca",
                           new.reduction = "harmony",
                           verbose = FALSE)                # Not working?

# Amy pipeline on Zsombor ########################################################

mph <- JoinLayers(mph)


GetAssayData(mph, ass)
cts=LayerData(mph, assay = "RNA", layer = "counts")

cellsOfInterest=data.frame("cell"=colnames(cts))
coldata = read.csv("C:/Users/bodha/OneDrive - University of Glasgow/Tommy/Joana_Amy/05_Analysis/coldata_noNeut2020-01-21.csv") %>% dplyr::filter(cell %in% cellsOfInterest$cell )
coldata=as.matrix(coldata)

print("Are the rownames of coldata the same as the colnames in the count matrix?")
all(rownames(coldata) == colnames(cts))

library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(counts=cts), 
  colData = coldata
)
sce


sce <- as.SingleCellExperiment(mph)
#first check how many genes would remain
ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)
sce <- sce[keep,]

is.spike <- grepl("^ERCC-", rownames(sce))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
altExpNames(sce)

library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

mito <- which(rowData(sce)$CHR=="chrM")

par(mfrow=c(1,2))
hist(counts(sce), breaks=40, main="", col="grey80", ylim = c(0,500), xlab="counts")
hist(log((counts(sce))), breaks=40, main="", col="grey80", ylim = c(0,100000), xlab="logcounts")

df <- perCellQCMetrics(sce, subsets=list(Mito=mito))
colnames(df)

#adding cell names to this vector, so I can individually ID cells
rownames(df)=as.character(sce@colData$cell)

cellsDet=data.frame("cell"=rownames(df),
                    "features"=df$detected,
                    "ERCC"=df$altexps_ERCC_sum,
                    "mitoc_genes"=df$subsets_Mito_sum)

par(mfrow=c(2,2))
hist(df$subsets_Mito_percent, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df$altexps_ERCC_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df$sum/1e6, xlab="sum of counts (^1e6)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(df$detected, xlab="features detected", main="",
     breaks=20, col="grey80", ylab="Number of cells")

par(mfrow=c(1,2))
plot(df$altexps_ERCC_percent,df$subsets_Mito_percent, ylab="Mitochondrial proportion (%)",
     xlab="ERCC proportion (%)",main="", col="grey50")
plot(df$detected,df$subsets_Mito_percent, ylab="Mitochondrial proportion (%)",
     xlab="features detected",main="", col="grey50")
qc.param <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent",
                                                 "altexps_ERCC_percent"), nmads=4)
colSums(as.matrix(qc.param))

colData(sce) <- cbind(colData(sce), df)
sce$plate <- factor(sce$plate)
sce$condition <- ifelse(grepl("control", sce$condition),
                        "control", "CML")
sce$discard <- qc.param$discard
gridExtra::grid.arrange(
  plotColData(sce, x="plate", y="sum", colour_by="discard",
              other_fields="condition", shape_by = "condition") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="plate", y="detected", colour_by="discard", 
              other_fields="condition", shape_by = "condition") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="plate", y="subsets_Mito_percent", 
              colour_by="discard", other_fields="condition", shape_by = "condition") + 
    ggtitle("Mito percent"),
  plotColData(sce, x="plate", y="altexps_ERCC_percent", 
              colour_by="discard", other_fields="condition", shape_by = "condition") + 
    ggtitle("ERCC percent"),
  ncol=2
)

plotColData(sce, x="sum", y="subsets_Mito_percent", 
            colour_by="discard", other_fields=c("plate", "condition"), shape_by = "plate") +
  geom_text(aes(label=ifelse(df$subsets_Mito_percent>20,as.character(rownames(df)),'')),hjust=0,vjust=0, size=2) +
  facet_grid(~condition) +
  theme(panel.border = element_rect(color = "grey"))

plotColData(sce, x="altexps_ERCC_percent", y="subsets_Mito_percent",
            colour_by="discard", other_fields=c("plate", "condition"), shape_by = "plate") + 
  geom_text(aes(label=ifelse(df$subsets_Mito_percent>20,as.character(rownames(df)),'')),hjust=0,vjust=0, size=2) +
  facet_grid(~condition) + 
  theme(panel.border = element_rect(color = "grey"))

filt_sce <- sce[,!qc.param$discard]

print(paste(as.character(length(sce@colData$cell)), as.character(length(filt_sce@colData$cell))))

df_filtered <- perCellQCMetrics(filt_sce, subsets=list(Mito=mito))
#adding cell names to this vector, so I can individually ID cells
rownames(df_filtered)=as.character(filt_sce@colData$cell)
cellsDetFilt=data.frame("cell"=rownames(df_filtered),
                        "features"=df_filtered$detected,
                        "ERCC"=df_filtered$altexps_ERCC_sum,
                        "mitoc_genes"=df_filtered$subsets_Mito_sum)

par(mfrow=c(2,2))
hist(df_filtered$subsets_Mito_percent, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df_filtered$altexps_ERCC_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(df_filtered$sum/1e6, xlab="sum of counts (^1e6)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(df_filtered$detected, xlab="features detected", main="",
     breaks=20, col="grey80", ylab="Number of cells")

logcounts(filt_sce) <- log(counts(filt_sce)+1)

filt_sce <- scater::runPCA(filt_sce)

set.seed(2302)
filt_sce <- scater::runTSNE(filt_sce, perplexity=50)

reducedDims(filt_sce)

filt.counts=counts(filt_sce)

library(scran)
dec <- modelGeneVar(filt_sce)#, block=assignments$phases)
par(mfrow=c(1,1))
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")

output <- getClusteredPCs(reducedDim(filt_sce))
#npcs <- metadata(output)$chosen
npcs <- output@metadata$chosen #Bodhayan by commenting previous
reducedDim(filt_sce, "PCAsub") <- reducedDim(filt_sce, "PCA")[,1:npcs,drop=FALSE]
npcs

g <- buildSNNGraph(filt_sce, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership
filt_sce$cluster <- factor(cluster)
table(filt_sce$cluster)

set.seed(1)
filt_sce <- runTSNE(filt_sce, dimred="PCAsub")
plotTSNE(filt_sce, colour_by="cluster", shape_by="condition")

# Change ofcolor
new_pallete <- c("skyblue3", "red3", "salmon", "lavender", "skyblue", "indianred3")
p <- plotTSNE(filt_sce, colour_by="cluster", shape_by="condition", point_size=3)
p  + scale_color_manual(values = new_pallete) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))

new_pallete <- c("skyblue3", "red3", "salmon", "lavender", "skyblue", "indianred3")
p1 <- plotTSNE(filt_sce, colour_by="Cd36", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p2 <- plotTSNE(filt_sce, colour_by="Lgals1", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p3 <- plotTSNE(filt_sce, colour_by="Cd14", shape_by="condition", point_size=3) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))
p1+p2+p3

set.seed(100)
filt_sce <- runUMAP(filt_sce, dimred="PCAsub")
plotUMAP(filt_sce, colour_by="cluster", shape_by="condition")

plotTSNE(filt_sce, colour_by="cluster", shape_by="plate")

plotTSNE(filt_sce, colour_by="cluster", shape_by="plate", point_size=3) + scale_color_manual(values = new_pallete) + theme(axis.text = element_text(size=10), axis.title = element_text(size=14, face="bold"))

set.seed(2302)
plotTSNE(filt_sce, colour_by="Cd36", shape_by="condition")

ratio <- bluster::pairwiseModularity(g, cluster, as.ratio=TRUE)
library(pheatmap)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))

ass.prob <- bluster::bootstrapStability(filt_sce, FUN=function(x) {
  g <- buildSNNGraph(x, use.dimred="PCAsub")
  igraph::cluster_walktrap(g)$membership
}, clusters=sce$cluster)
pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE,
         col=colorRampPalette(c("white", "blue"))(100))

markers <- findMarkers(filt_sce, filt_sce$cluster)
markers.table<-markers[[1]]

top.markers <- rownames(markers.table)[markers.table$Top <= 10]
plotHeatmap(filt_sce, features=top.markers, columns=order(filt_sce$cluster),
            colour_columns_by=c("cluster", "condition"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5), 
            show_colnames = FALSE) 

library(gplots)
ann_colors = list(cluster = c("1"="skyblue3","2"="red3","3"="salmon","4"="lavender","5"="skyblue", "6"="indianred3"),
                  condition = c("control" = "lightblue3", "CML" = "pink"))
plotHeatmap(filt_sce, features=top.markers, color_columns_by = c("cluster","condition"), column_annotation_colors = ann_colors, 
            color = bluered(256), center = T, cex=0.9, 
            columns = order(filt_sce$cluster), cluster_cols = F, zlim=c(-5, 5))

# Saving Amy obj locally
saveRDS(filt_sce, file = "amy_sce_obj.rds")
amy_obj = readRDS("amy_sce_obj.rds")
plotTSNE(amy_obj) # Check





# Differential analyses ##################################################################
# Ima vs Veh
obj <- PrepSCTFindMarkers(obj)
obj <- SetIdent(obj, value = 'orig.ident')

FindConservedMarkers(obj,
                     ident.1 = 1,
                     grouping.var = 'orig.ident')
mk <- FindMarkers(obj,
                  ident.1 = 'ima',
                  ident.2 = 'veh')

iv_mk <- rownames(mk[mk$p_val_adj < 0.05 & (mk$avg_log2FC >= 1 | mk$avg_log2FC <= -1),])

mk <- FindMarkers(obj,
                  ident.1 = 'dtg',
                  ident.2 = 'stg')
ds_mk <- rownames(mk[mk$p_val_adj < 0.05 & (mk$avg_log2FC >= 1 | mk$avg_log2FC <= -1),])

mk <- FindMarkers(obj,
                  ident.1 = 'veh',
                  ident.2 = 'stg')
vs_mk <- rownames(mk[mk$p_val_adj < 0.05 & (mk$avg_log2FC >= 1 | mk$avg_log2FC <= -1),])

head(mk, n = 5)

iv_mk
ds_mk
vs_mk

length(iv_mk)
length(ds_mk)
length(vs_mk)

s2_int <- intersect(ds_mk, iv_mk)
length(s2_int)
s2_rest <- setdiff(ds_mk, iv_mk)
length(s2_rest)
length(intersect(s2_rest, vs_mk))
s2_qc <- setdiff(s2_rest, vs_mk)
length(s2_qc)

DotPlot(obj, features = s2_qc, split.by = "label", cols = 1:15) + RotatedAxis()

DotPlot(obj, features = s2_qc, split.by = "orig.ident", cols = 1:15) #+ RotatedAxis()

#Volcano
mk$diffexpressed <- "NO"
mk$diffexpressed[mk$avg_log2FC > 1 & mk$p_val_adj < 0.1] <- "UP"
mk$diffexpressed[mk$avg_log2FC < -1 & mk$p_val_adj < 0.1] <- "DOWN"

# Higlight all diff ex genes
mk$gene_symbol <- row.names(mk)
mk$delabel <- NA
mk$delabel[mk$diffexpressed != "NO"] <- mk$gene_symbol[mk$diffexpressed != "NO"]

# OR

#Highlight specific genes
sp <- c('SHMT1','SHMT2')
mk$gene_symbol <- row.names(mk)
mk$delabel <- NA
mk$delabel[mk$gene_symbol %in% sp] <- mk$gene_symbol[mk$gene_symbol %in% sp]


ggplot(data=mk, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values=c("dodgerblue", "black", "indianred")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red")





# GSEA https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
library(clusterProfiler)
library(enrichplot)

# Interested in Phagocytosis, Efferocytosis

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

# Prepare Input
df <- mk[mk$p_val_adj < 0.05,]
# df <- mk[rownames(mk) %in% s2_qc,]
original_gene_list <- df$avg_log2FC
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# GO
gse <- gseGO(geneList=gene_list, 
             ont ="BP", #ALL, BP, CC, MF
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Encrichment Map
gse1 <- pairwise_termsim(gse) # Not present
emapplot(gse1, showCategory = 10)

# Category Netplot
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3) # categorySize can be either 'pvalue' or 'geneNum'

# Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

# GSEA plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

View(gse@result)

# PubMed trend of enriched terms
library(europepmc)
terms <- gse$Description[1:10]
pmcplot(terms, 2010:2024, proportion=FALSE)

# KEGG Gene Set Enrichment Analysis ####################################################
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df2 = df[rownames(df) %in% dedup_ids$SYMBOL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$avg_log2FC
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "mmu" #hsa
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

# Dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# Encrichment map
kk21 <- pairwise_termsim(kk2)
emapplot(kk21)

#Category Netplot:
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list) # categorySize can be either 'pvalue' or 'geneNum'

# Ridgeplot
ridgeplot(kk2) + labs(x = "enrichment distribution")

# GSEA Plot
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

# Pathview ###############################################################################
library(pathview)
mmu <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04820", species = kegg_organism)
mmu <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05167", species = kegg_organism)
mmu <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04820", species = kegg_organism, kegg.native = F)

knitr::include_graphics("mmu04130.pathview.png")

# Over-Representation Analysis #########################################################
#prepare input
genes <- names(gene_list)[abs(gene_list) > 2] # filter on min log2fold change (log2FoldChange > 2)

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
# Upset Plot
library(enrichplot)
library(ggupset)
upsetplot(go_enrich)

# Wordcloud
library(wordcloud)
wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)

# Barplot
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

# Dotplot
dotplot(go_enrich)

# Encrichment map
go_enrich1 <- pairwise_termsim(go_enrich)
emapplot(go_enrich1)

# Enriched GO induced graph:
goplot(go_enrich, showCategory = 10)

# Category Netplot
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list) # categorySize can be either 'pvalue' or 'geneNum'

# KEGG Pathway Enrichment #######################################################
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df2 = df[rownames(df) %in% dedup_ids$SYMBOL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$avg_log2FC
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_sig_genes_df = subset(df2, p_val_adj < 0.05)
kegg_genes <- kegg_sig_genes_df$avg_log2FC
names(kegg_genes) <- kegg_sig_genes_df$Y
kegg_genes <- na.omit(kegg_genes)
kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]

# Create enrichKEGG object
kegg_organism = "mmu"
kk <- enrichKEGG(gene=kegg_genes,
                 universe=names(kegg_gene_list),
                 organism=kegg_organism,
                 pvalueCutoff = 0.05,
                 keyType = "ncbi-geneid")

# Barplot
barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

# Dotplot
dotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

# Category Netplot
cnetplot(kk, categorySize="pvalue", foldChange=gene_list)

# Pathview
ids<-bitr(names(gene_list), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb=organism) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
names(gene_list) <- ids$ENSEMBL
mmu <- pathview(gene.data=gene_list, # requires ENSEMBL
                pathway.id="mmu04080",
                species = kegg_organism,
                gene.idtype=gene.idtype.list[3])
dme <- pathview(gene.data=gene_list, pathway.id="mmu04080", species = kegg_organism, gene.idtype=gene.idtype.list[3], kegg.native = F)

knitr::include_graphics("mmu04080.pathview.png")




#Saving h5ad - Not working
library(SeuratDisk)
SaveH5Seurat(obj,
             filename = "zk_int.h5Seurat",
             overwrite = TRUE)
Convert("zk_int.h5Seurat", dest = "h5ad")
