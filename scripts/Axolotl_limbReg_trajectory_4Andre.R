##########################################################################
##########################################################################
# Project: Andre's cultured cell trajectory 
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Mar 31 11:49:38 2026
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
#library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)
#library("viridis")

version.analysis = '_axolotl20260331'
resDir = paste0("../results/scRNAseq_axolotle_limbReg_Diego", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/Andre/'

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')


########################################################
########################################################
# Section I: Import the new scRNA-seq data for axolotl limb regeneration
# with dpa0, dpa5, dpa11, dpa18
########################################################
########################################################

processing_annot = FALSE
if(processing_annot){
  annot = rtracklayer::import(paste0("/groups/tanaka/People/current/Diego/Projects/",
                                     "8_Transxolotl/5_AxolotlT2T/data/00_assemblies/1_hifiasm/",
                                     "231106_22xFlowcells+UL/assembly/primary_curated/250329/",
                                     "annotation/geneAnnotation.with-chrM.gtf"))
  annot = data.frame(annot$gene_id, annot$name)
  colnames(annot) = c('gene_id', 'gene')
  annot = annot[!is.na(annot$gene), ]
  
  annot$name = paste0(annot$gene, '_', annot$gene_id)
  
  names_uniq = unique(annot$name)
  
  annot = annot[match(names_uniq, annot$name), ]
  
  jj = which(annot$gene_id != annot$gene)
  
  xx = annot 
  genes_uniq = unique(xx$gene)
  for(n in 1:length(genes_uniq))
  {
    kk = which(xx$gene == genes_uniq[n])
    
    if(length(kk) > 1) {
      cat(n, '--', genes_uniq[n], '--', xx$name[kk], '\n')
      xx$gene[kk] = xx$name[kk]
    }
  }
  
  annot = xx
  rm(xx)
  
  saveRDS(annot, file = paste0(RdataDir, 'primary_curated_250329_geneAnnotation.withChrM.rds'))
  
}

annot = readRDS(file = paste0(RdataDir, 'primary_curated_250329_geneAnnotation.withChrM.rds'))



aa = readRDS(paste0(dataDir, '250701_LimbRegeneration.RDS'))
metadata = aa@meta.data
counts = aa@assays$GenesExons$counts

mm = match(rownames(counts), annot$gene_id)

rownames(counts) = annot$gene[mm]

xx <- CreateSeuratObject(counts = counts, project = "axoLimbReg", meta.data = metadata)

xx[['integrated']] = aa[['integrated']]
xx[['tsne']] = aa[['tsne']]
xx[['umap3d']] = aa[['umap3d']]
xx[['tsne3d']] = aa[['tsne3d']]
xx[['umap_old']] = aa[['umap']]

rm(aa)
aa = xx
rm(xx)

aa$condition = aa$orig.ident

aa$condition[which(aa$condition == 'ForelimbRegeneration_UpperArmMature_1')] = 'dpa0'
aa$condition[which(aa$condition == 'ForelimbRegeneration_UpperArm5dpa_1')] = 'dpa5'
aa$condition[which(aa$condition == 'ForelimbRegeneration_UpperArm11dpa_1')] = 'dpa11'
aa$condition[which(aa$condition == 'ForelimbRegeneration_UpperArm18dpa_1')] = 'dpa18'

aa$condition = factor(aa$condition, levels = c('dpa0', 'dpa5', 'dpa11', 'dpa18'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)


p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'cell_type', reduction = 'umap', label = TRUE, repel = TRUE)

p1 /p2


mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4", "ATP8", "MT-CO1", "COI")
mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))

mtgenes = rownames(aa)[!is.na(match(rownames(aa), mtgenes))]

xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "RNA", features = mtgenes)
aa[['percent.mt']] = xx$percent.mt
rm(xx)

# Visualize QC metrics as a violin plot
VlnPlot(aa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'condition', ncol = 1)
VlnPlot(aa, features = 'nFeature_RNA', y.max = 7000, group.by = 'condition') +
  geom_hline(yintercept = c(500,1000))

VlnPlot(aa, features = 'nCount_RNA', y.max = 50000, group.by = 'condition') +
  geom_hline(yintercept = c(500,1000))

VlnPlot(aa, features = 'percent.mt', y.max = 20, group.by = 'condition') +
  geom_hline(yintercept = c(5,10))


## second time cell filtering 
aa <- subset(aa, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)

#aa = subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

saveRDS(aa, file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_QCsfiltered.rds'))

##########################################
# discard the doublets 
##########################################
library(DoubletFinder)

aa = readRDS(file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_QCsfiltered.rds'))

aa$DF_out = NA

#aa$condition = factor(aa$condition)
Idents(aa) = aa$condition

cc = unique(aa$condition)

for(n in 2:length(cc))
{
  # n = 1
  cat(n, '-----', as.character(cc[n]), '\n')
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, verbose = TRUE)
  subs <- FindNeighbors(subs, dims = 1:30)
  subs <- FindClusters(subs, resolution = 0.5)
  
  subs <- RunUMAP(subs, dims = 1:30)
  
  sweep.res.list_nsclc <- paramSweep(subs, PCs = 1:30)
  sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
  
  pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- subs@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  
  
  nExp_poi <- round(0.076*nrow(subs@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  subs <- DoubletFinder::doubletFinder(subs, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                                          reuse.pANN = FALSE , sct = FALSE)
  
  df_out = subs@meta.data
  subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
  
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
          raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], version.analysis, '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], version.analysis,  
                              '.rds'))
  
}

DimPlot(aa, group.by = 'DF_out')

saveRDS(aa, file = paste0(RdataDir, 
                          '/ForeLimb_upperArmRegeneration_QCsfiltered_DFout.rds'))

##########################################
# redo processing after DF filtering 
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           '/ForeLimb_upperArmRegeneration_QCsfiltered_DFout.rds'))


DimPlot(aa, group.by = 'DF_out')

ggsave(filename = paste0(resDir, '/all_doubletFinder_out.pdf'), 
       width = 12, height = 8)


aa = subset(aa, cells = colnames(aa)[which(aa$DF_out == 'Singlet')])

aa$condition = factor(aa$condition, levels = c('dpa0', 'dpa5', 'dpa11', 'dpa18'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)


p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'cell_type', reduction = 'umap', label = TRUE, repel = TRUE)

p1 /p2

ggsave(filename = paste0(resDir, '/UMAP_condition_cellTypes_filtered.pdf'), 
       width = 12, height = 16)

##########################################
# subset CT  
##########################################
aa = subset(aa, cells = colnames(aa)[which(aa$cell_type == 'Mesenchyme')])

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa, features = rownames(aa))
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)


p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'cell_type', reduction = 'umap', label = TRUE, repel = TRUE)

p1 /p2


aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'curated_clusters', reduction = 'umap', label = TRUE, repel = TRUE) + 
  NoLegend()

p1 + p2

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE) 

p1 + p2

ggsave(filename = paste0(resDir, '/UMAP_CT_newClusters.pdf'), width = 16, height = 8)



markers =  toupper(c('Procr', 'Dpt', 'Pi16', 'Col1a2', 'Acta2', 'Lum', 'Col3a1', 'Col1a1', 'Mmp2', 
                     'Pdgfra',  'Twist1', 'Vim', 'Mmp9', 'PRRX1'))
mm = match(markers, rownames(aa))
markers = markers[!is.na(mm)]

FeaturePlot(aa, features = markers)

ggsave(filename = paste0(resDir, '/UMAP_CT_featurePlots.pdf'), 
       width = 16, height = 12)

oupMarker <- FindAllMarkers(aa, logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE)
oupMarker = oupMarker[grep('^AME', oupMarker$gene, invert = TRUE), ]

oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, "heatmap_CT_newClusters_topMarkers.pdf"),
       width = 12, height = 20)


FeaturePlot(aa, features = c('OTOS', 'SOX9', "ACAN",'COMP'))

aa$cell_type[which(aa$seurat_clusters == "10")] = 'CT_cartilage'
DimPlot(aa, group.by = 'cell_type', reduction = 'umap', label = TRUE, repel = TRUE)

FeaturePlot(aa, features = c('CDH1', "KRT8", 'KRT7', 'CLDN1'), reduction = 'umap')

aa$cell_type[which(aa$seurat_clusters == "11")] = 'CT_epithelia'

markers =  toupper(c('Dpt', 'Col1a2',  'Lum', 'Col3a1', 
                     'Pdgfra',  'Twist1',  'PRRX1'))
mm = match(markers, rownames(aa))
markers = markers[!is.na(mm)]

FeaturePlot(aa, features = markers)

FeaturePlot(aa, features = c('nCount_RNA', 'nFeature_RNA'))

VlnPlot(aa, features = c('nCount_RNA', 'nFeature_RNA'), group.by = 'seurat_clusters')

ggsave(filename = paste0(resDir, '/UMAP_CT_featurePlots_lowCoverage.pdf'), 
       width = 16, height = 8)

aa = subset(aa, cells = colnames(aa)[which(aa$seurat_clusters != '4' & aa$seurat_clusters != '11' &
                                             aa$seurat_clusters != '10')])


aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa, features = rownames(aa))
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE)

p1  + p2

oupMarker <- FindAllMarkers(aa, logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE)
oupMarker = oupMarker[grep('^AME', oupMarker$gene, invert = TRUE), ]

oupMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, "heatmap_CT_newClusters_topMarkers_v2.pdf"),
       width = 12, height = 20)


aa = subset(aa, cells = colnames(aa)[which(aa$seurat_clusters != '9')])

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000) # find subset-specific HVGs

aa <- ScaleData(aa, features = rownames(aa))
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE)

p1  + p2


saveRDS(aa, file = paste0(RdataDir, 
                          '/ForeLimb_upperArmRegeneration_QCsfiltered_DFout_CTfiltered.rds'))


########################################################
########################################################
# Section II: Trajectory analysis 
# 
########################################################
########################################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

aa = readRDS(file = paste0(RdataDir, 
                           '/ForeLimb_upperArmRegeneration_QCsfiltered_DFout_CTfiltered.rds'))

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE)

p1  + p2

Idents(aa) = aa$condition

sce = as.SingleCellExperiment(aa)
#rm(aa)
dec <- modelGeneVar(sce)

nb_features = 5000; n_neighbors = 200;

top.hvgs <- getTopHVGs(dec, n=nb_features)

sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
# reducedDimNames(sce)
ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]


tic()
dm <- DiffusionMap(ll.pca, sigma = 'local', k = n_neighbors, n_eigs = 20, distance = 'cosine')

toc()

cells = names(dm$DC1)
metadata = aa@meta.data
dcs = data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, 
                 DC4 = dm$DC4, DC5 = dm$DC5, stringsAsFactors = FALSE)
dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)
aa[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(aa))

rm(metadata)

p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2), group.by = 'condition')
p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), group.by = 'condition')
#p3 = DimPlot(aa, reduction = 'DC', dims = c(2, 3), group.by = 'condition')
p1 + p2


aa <- FindNeighbors(aa, reduction = 'DC', dims = 1:5)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.3)


p1 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'DC', label = TRUE, repel = TRUE)

p1  + p2

aa = subset(aa, cells = colnames(aa)[which(aa$seurat_clusters != '10')])

saveRDS(aa, file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_QCsfiltered_DFout_CTfiltered_DM.rds'))

##########################################
# run slingshot using the seurat clusters 
##########################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

aa = readRDS(file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_QCsfiltered_DFout_CTfiltered_DM.rds'))

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa, features = rownames(aa))
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, repel = TRUE)

p1  + p2

p1 = DimPlot(aa, group.by = 'condition', reduction = 'pca', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'pca', label = TRUE, repel = TRUE)

p1  + p2

#aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'pca', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'pca', label = TRUE, repel = TRUE)

p1  + p2


p1 = DimPlot(aa, group.by = 'condition', reduction = 'DC', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', reduction = 'DC', label = TRUE, repel = TRUE)

p1  + p2

## start the slingshot 
## example code found https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
sce = as.SingleCellExperiment(aa)

# dec <- modelGeneVar(sce)
# 
# nb_features = 3000
# top.hvgs <- getTopHVGs(dec, n=nb_features)
# 
# sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 50, scale = FALSE)

# reducedDimNames(sce)
rd1= reducedDim(sce, 'PCA')[, c(1:2)]
rd2 = Embeddings(aa, reduction = "pca")[, c(1:2)]

# FQnorm <- function(counts){
#   rk <- apply(counts,2,rank,ties.method='min')
#   counts.sort <- apply(counts,2,sort)
#   refdist <- apply(counts.sort,1,median)
#   norm <- apply(rk,2,function(r){ refdist[r] })
#   rownames(norm) <- rownames(counts)
#   return(norm)
# }
# 
# 
# assays(sce)$norm <- FQnorm(assays(sce)$counts)
# 
# pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
# rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

cl1 <- colData(sce)$condition
#colData(sce)$GMM <- cl1

#library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

plot(rd2, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'condition', reducedDim = 'PCA', start.clus = 'dpa0',
                 end.clus = 'dpa18')

pt <- slingPseudotime(sce)
aa$pseudotime <- pt[, 1]

saveRDS(sce, file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_slingshotOut.rds'))

library(grDevices)

pdf(paste0(resDir, "Slingshot_getLineage_withStarting.cluster.pdf"),
    height = 6, width =10, useDingbats = FALSE)

p1 = DimPlot(aa, group.by = 'condition', reduction = 'pca', label = TRUE, repel = TRUE, pt.size = 1.5,
             label.size = 6)
plot(p1)

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$condition], pch=16)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

plot(reducedDims(sce)$PCA, col = plotcol, pch=16)
lines(SlingshotDataSet(sce), lwd=3, col='black')

p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1.,
             label.size = 6)
p2 = FeaturePlot(aa, features = 'pseudotime') +  
  scale_color_viridis_c(option = "magma") +
  #scico::scale_color_scico(palette = "vik") + 
  ggtitle("Slingshot pseudotime")

plot(p1 + p2)

dev.off()


##########################################
# test driven genes for trajectory using code from
# https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html
##########################################
#library(tradeSeq)
#library(TSCAN)
library(gam)

sce = readRDS(file = paste0(RdataDir, 'ForeLimb_upperArmRegeneration_slingshotOut.rds'))

# Only look at the 3,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(sce) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:3000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
pt <- slingPseudotime(sce)
sce$pseudotime <- pt[, 1]
t <- sce$pseudotime

gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

drivers = data.frame(gene = names(gam.pval), gam.pval, stringsAsFactors = FALSE)
drivers = drivers[order(drivers$gam.pval), ]
drivers$gam.fdr = p.adjust(drivers$gam.pval)

write.csv2(drivers, file = paste0(resDir, 'drivers_trajectory_gam.csv'), quote = FALSE, 
           row.names = FALSE)
# Identify genes with the most significant time-dependent model fit.


topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:300]
# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
require(clusterExperiment)
heatdata = aa@assays$RNA$data
heatdata <- as.matrix(heatdata[rownames(heatdata) %in% topgenes, order(t, na.last = NA)])

heatclus <- aa$condition[order(t, na.last = NA)]

png(paste0(resDir, "heatmap_time_genes.png"), width=10, height=10, units = "in", res=200)
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", 
                               visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)
dev.off()


########################################################
########################################################
# Section III: test Andre's vitro data 
# 
########################################################
########################################################
aa = readRDS(file = paste0("/groups/tanaka/People/current/Andre/260314_FibroblastActivation.RDS"))

DimPlot(aa)

DefaultAssay(aa) = 'GenesExons'


aa$condition = colnames(aa)
aa$condition = sapply(aa$condition, 
                      function(x){xx = unlist(strsplit(x, '_')); paste0(xx[1:(length(xx)-3)], '_')})


aa$condition = gsub('FibroblastActivation_PrimaryLimbCells_', '', aa$condition)



DimPlot(aa, reduction = 'umap', group.by = 'condition')

DimPlot(aa, features = '')


