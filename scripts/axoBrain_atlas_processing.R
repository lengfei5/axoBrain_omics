##########################################################################
##########################################################################
# Project: axolotl brain omics data processing 
# Script purpose: processing axoBrain_atlas data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Apr  2 16:32:43 2026
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

library(zellkonverter)
library(SingleCellExperiment)

inputDir = '/groups/tanaka/Collaborations/axolotlBrain_Omics/data_shared_Mateja/'
outDir = '/groups/tanaka/Collaborations/axolotlBrain_Omics/Rdata/'

if(!dir.exists(outDir)) dir.create(outDir)

resDir = paste0("../results/axoBrain_Atlas_analysis/")
RdataDir = paste0('../results/Rdata/')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

##########################################
# convert the h5ad file to seurat object 
##########################################
#library(SingleCellExperiment)
#library(SeuratDisk)

# Convert h5ad file to h5seurat format
#h5ad_file = paste0(inputDir, 'Retina.h5ad')
#sce = readH5AD(h5ad_file, use_hdf5 = TRUE)

#Convert(h5ad_file, dest = "h5seurat")

# load Seurat file
#aa <- LoadH5Seurat(paste0(inputDir, 'Retina.h5seurat'))

##########################################
# creat Seurat objec for Retina 
##########################################
counts = read.csv(paste0(outDir, 'Retina_adata_layer_counts.csv'), header = TRUE, sep = ',', row.names = c(1))
metadata = read.csv(paste0(outDir, 'Retina_adata_metacell.csv'), header = TRUE, sep = ',', row.names = c(1))
pca = read.csv(paste0(outDir, 'Retina_adata_rd_pca.csv'), header = TRUE, sep = ',', row.names = c(1))
pca_harmony = read.csv(paste0(outDir, 'Retina_adata_rd_pca_harmony.csv'), header = TRUE, sep = ',', 
                       row.names = c(1))

umap = read.csv(paste0(outDir, 'Retina_adata_rd_umap.csv'), header = TRUE, sep = ',', 
                       row.names = c(1))

aa <- CreateSeuratObject(counts = t(as.matrix(counts)), project = "axoBrain", assay = "RNA",
                         min.cells = 3, min.features = 50, meta.data = metadata)

colnames(umap) <- paste0("UMAP_", 1:ncol(umap))
aa[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(aa))

colnames(pca) = paste0("PCA_", 1:ncol(pca))
aa[['pca_orig']] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PCAorig_", assay = DefaultAssay(aa))

colnames(pca_harmony) = paste0("PCA_harmony", 1:ncol(pca_harmony))
aa[['pca_harmony']] <- CreateDimReducObject(embeddings = as.matrix(pca_harmony), 
                                            key = "PCAharmony_", assay = DefaultAssay(aa))

DimPlot(aa, group.by = 'anno', reduction = 'umap', label = TRUE, repel = TRUE, label.size = 4)

ggsave(filename = paste0(resDir, '/Retina_umap_annotation.pdf'), 
       width = 12, height = 8)

saveRDS(aa, file = paste0(outDir, 'SeuratObject_v5_Retina.rds'))

##########################################
# creat Seurat objec for Hindbrain
##########################################
brainSection = 'Midbrain'
counts = read.csv(paste0(outDir, brainSection, '_adata_layer_counts.csv'), header = TRUE, 
                  sep = ',', row.names = c(1))
metadata = read.csv(paste0(outDir, brainSection, '_adata_metacell.csv'), header = TRUE, 
                    sep = ',', row.names = c(1))
pca = read.csv(paste0(outDir, brainSection, '_adata_rd_pca.csv'), header = TRUE, 
               sep = ',', row.names = c(1))
pca_harmony = read.csv(paste0(outDir, brainSection, '_adata_rd_pca_harmony.csv'), header = TRUE, sep = ',', 
                       row.names = c(1))
umap = read.csv(paste0(outDir, brainSection, '_adata_rd_umap.csv'), header = TRUE, sep = ',', 
                row.names = c(1))

aa <- CreateSeuratObject(counts = t(as.matrix(counts)), project = "axoBrain", assay = "RNA",
                         min.cells = 3, min.features = 50, meta.data = metadata)

colnames(umap) <- paste0("UMAP_", 1:ncol(umap))
aa[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(aa))

colnames(pca) = paste0("PCA_", 1:ncol(pca))
aa[['pca_orig']] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PCAorig_", assay = DefaultAssay(aa))

colnames(pca_harmony) = paste0("PCA_harmony", 1:ncol(pca_harmony))
aa[['pca_harmony']] <- CreateDimReducObject(embeddings = as.matrix(pca_harmony), 
                                            key = "PCAharmony_", assay = DefaultAssay(aa))

DimPlot(aa, group.by = 'batch', reduction = 'umap', label = TRUE, repel = TRUE, label.size = 4)

ggsave(filename = paste0(resDir, '/', brainSection, '_umap_batchInfo.pdf'), 
       width = 12, height = 8)

saveRDS(aa, file = paste0(outDir, 'SeuratObject_v5_', brainSection, '.rds'))


VlnPlot(aa, features = c('nFeature_RNA'), group.by = 'batch') +
  geom_hline(yintercept = c(500,1000))



