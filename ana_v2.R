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
