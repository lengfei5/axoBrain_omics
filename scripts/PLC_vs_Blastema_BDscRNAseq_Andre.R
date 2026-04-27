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
#resDir = paste0("../results/scRNAseq_axolotle_limbReg_Diego", version.analysis, '/')
RdataDir = paste0('../results/Rdata/')
resDir = paste0("../results/scRNAseq_axolotle_PrimaryLimbCells", version.analysis, '/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/Andre/'

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

annot = readRDS(file = paste0(RdataDir, 'primary_curated_260226_geneAnnotation.withChrM.rds'))

levels = c("MatLimb_0dpa_1", "Blastema_11dpa_1", 
           "PLC_1dpd_1", "PLC_3dpd_1", "PLC_5dpd_0","PLC_5dpd_1",
           "PLC_8dpd_1", "PLC_10dpd_1", "PLC_12dpd_1", "PLC_15dpd_1",
           "PLC_8dpd_IL11treated_1", "PLC_8dpd_Dexatreated_1")
cols = c('#054674', '#801517',
         '#31C53F', '#2FF18B', '#28CECA','#4B4BF7',
         '#25aff5', '#D4D915','#ff9a36','#B95FBB',
         '#CCB1F1', '#AC8F14')

#names(cols) = levels
names(cols) = levels


c('#F68282','#31C53F','#1FA195','#B95FBB','#D4D915',
  '#28CECA', '#ff9a36', '#2FF18B',
  "#054674", '#25aff5', "#4d7ea9", '#D4D915','#ff9a36','#B95FBB')
c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
  '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
  '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
  '10'='#E6C122')



########################################################
########################################################
# Section III: test Andre's vitro data 
# 
########################################################
########################################################

resDir = paste0("../results/scRNAseq_axolotle_PrimaryLimbCells", version.analysis, '/')
if(!dir.exists(resDir)) dir.create(resDir)
annot = readRDS(file = paste0(RdataDir, 'primary_curated_260226_geneAnnotation.withChrM.rds'))

##########################################
# process the matched gene annotation 
##########################################
processing_annot = FALSE
if(processing_annot){
  annot_file = paste0("/groups/tanaka/People/current/Diego/Projects/",
                      "8_Transxolotl/5_AxolotlT2T/data/00_assemblies/1_hifiasm/",
                      "231106_22xFlowcells+UL/assembly/primary_curated/260226/",
                      "annotation/geneAnnotation.with-chrM.gtf")
  annot = rtracklayer::import(annot_file)
  
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
  
  saveRDS(annot, file = paste0(RdataDir, 'primary_curated_260226_geneAnnotation.withChrM.rds'))
  
}


##########################################
# import and overview of Andre's data
##########################################
aa = readRDS(file = paste0("/groups/tanaka/People/current/Andre/260314_FibroblastActivation.RDS"))


DefaultAssay(aa) = 'GenesExons'

aa$condition = aa$orig.ident

aa$condition = gsub('FibroblastActivation_', '', aa$condition)
aa$condition = gsub('PrimaryLimbCells_', 'PLC_', aa$condition)
aa$condition = gsub('Tissue_Blastema', 'Blastema_', aa$condition)
aa$condition = gsub('Tissue_MatureLimb', 'MatLimb_0dpa', aa$condition)
aa$condition = gsub('DexamethasoneRuxolitinibtreated', 'Dexatreated', aa$condition)

DimPlot(aa, reduction = 'umap', group.by = 'condition', label = TRUE, repel = TRUE, raster=FALSE)

ggsave(filename = paste0(resDir, "Andre_PLC_scRNAseq_UMAP.pdf"),
       width = 12, height = 20)


metadata = aa@meta.data
counts = aa@assays$GenesExons$counts

mm = match(rownames(counts), annot$gene_id)

rownames(counts) = annot$gene[mm]

xx <- CreateSeuratObject(counts = counts, project = "axoPrimaryLimbCells", meta.data = metadata)

xx[['integrated']] = aa[['integrated']]
xx[['tsne']] = aa[['tsne']]
xx[['umap3d']] = aa[['umap3d']]
xx[['tsne3d']] = aa[['tsne3d']]
xx[['umap_old']] = aa[['umap']]

rm(aa)
aa = xx
rm(xx)


aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)


p1 = DimPlot(aa, group.by = 'condition', reduction = 'umap', label = TRUE, repel = TRUE)

p1

ggsave(filename = paste0(resDir, "Andre_PLCscRNAseq_UMAP_firstCheck.pdf"),
       width = 12, height = 8)

mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4", "ATP8", "MT-CO1", "COI")
mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
mtgenes = rownames(aa)[!is.na(match(rownames(aa), mtgenes))]

xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "RNA", features = mtgenes)
aa[['percent.mt']] = xx$percent.mt

rm(xx)

# Visualize QC metrics as a violin plot
#VlnPlot(aa, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), group.by = 'condition', ncol = 1)
VlnPlot(aa, features = 'nFeature_RNA', y.max = 7000, group.by = 'condition', pt.size = 0.0) +
  geom_hline(yintercept = c(500, 800, 1000))

VlnPlot(aa, features = 'nCount_RNA', y.max = 50000, group.by = 'condition', pt.size = 0.0) +
  geom_hline(yintercept = c(500,1000))

VlnPlot(aa, features = 'percent.mt', y.max = 20, group.by = 'condition', pt.size = 0.0) +
  geom_hline(yintercept = c(5,10))


## second time cell filtering 
aa <- subset(aa, subset = nFeature_RNA > 500 & percent.mt < 5)

#aa = subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)
saveRDS(aa, file = paste0(RdataDir, 'Andre_PLCscRNAseq_QCsfiltered.rds'))

##########################################
# identify doublet
##########################################
library(DoubletFinder)

aa = readRDS(file = paste0(RdataDir, 'Andre_PLCscRNAseq_QCsfiltered.rds'))

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
  
  sweep.res.list_nsclc <- DoubletFinder::paramSweep(subs, PCs = 1:30)
  sweep.stats_nsclc <- DoubletFinder::summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- DoubletFinder::find.pK(sweep.stats_nsclc)
  
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
                          '/Andre_PLCscRNAseq_QCsfiltered_DFout.rds'))


##########################################
# clean the doublet
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           '/Andre_PLCscRNAseq_QCsfiltered_DFout.rds'))

aa = subset(aa, subset = DF_out == "Singlet")

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:50, n.neighbors = 50,  min.dist = 0.3)




aa$condition = factor(aa$condition, levels = levels)

DimPlot(aa, group.by = 'condition', reduction = 'umap', cols = cols, label = TRUE, repel = TRUE)


ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered_umap.pdf"),
       width = 12, height = 8)

FeaturePlot(aa, features = 'nFeature_RNA', reduction = 'umap', label = TRUE, repel = TRUE)

saveRDS(aa, file = paste0(RdataDir, 
                          '/Andre_PLCscRNAseq_QCsfiltered_rmDFout.rds'))

########################################################
########################################################
# Section IV: manually annot clusters and select CT cells
# 
########################################################
########################################################

##########################################
# clusters and marker genes 
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           '/Andre_PLCscRNAseq_QCsfiltered_rmDFout.rds'))

ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

DimPlot(aa,  reduction = 'umap',  label = TRUE, repel = TRUE)

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered_clusters_30PCs.res0.5.pdf"),
       width = 12, height = 8)

saveRDS(aa, file = paste0(RdataDir,
                          '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_30PCs.res0.5Clusters.rds'))


markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

saveRDS(markers, file = paste0(RdataDir, 
                               'Andre_PLCscRNAseq_QCsfiltered_rmDFout_markers_30PCs.res0.5.rds'))


## reload the calculated clusters and markers
aa = readRDS(file = paste0(RdataDir,
                           '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_30PCs.res0.5Clusters.rds'))

markers = readRDS(file = paste0(RdataDir, 
                                'Andre_PLCscRNAseq_QCsfiltered_rmDFout_markers_30PCs.res0.5.rds'))


cat(length(unique(markers$gene)), ' markers found in ax\n')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#ggs = top10$gene[which(top10$cluster == 'immune_cells')]
#ggs = ggs[grep('^AME|^LOC', ggs, invert = TRUE)]
#ggs = ggs[order(ggs)]
xx = subset(aa, downsample = 200)

DoHeatmap(xx, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, 
                         '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_markers_30PCs.res0.5_heatmap.pdf'), 
       width = 30, height = 40)


##########################################
# just select the CT cells
##########################################
features = c('Prrx1', 'Mfap5',  'Fbn1',   'Pdgfra', ## Li et al., 2020
             'Procr', 'Dpt', 'Pi16', 'Col1a2', 'Acta2', 'Lum', 'Col3a1', 'Col1a1', 'Mmp2', 'Pdgfra', ## Nastya
             'Col1a2', 'Vim', 'Fstl1', 'DDR2', 'Acta2', 'Postn', 'Tcf21', 'Pdgfra', 'Col3a1', 'Col1a1', 'Gsn', 
             'Fbln2', 'Sparc', 'MMP2', 'Rspo1', 'Lum', 'Col8a1', # Elad
             'Lum', 'Dpt', 'Postn', # whited paper
             'Dpt', 'Pi16' # Sabine
)
features = unique(toupper(features))
features = features[!is.na(match(features, rownames(aa)))]

features = c("PRRX1", "PDGFRA", "DPT", "ACTA2", "LUM", "COL1A2", "MMP2", 'POSTN', "COL3A1")

FeaturePlot(aa, features = features)

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered", 
                         "_clusters.30PCs.res0.5_CTfeatures.pdf"),
       width = 16, height = 12)


VlnPlot(aa, features = "PRRX1")
VlnPlot(aa, features = "DPT")
VlnPlot(aa, features = "LUM")

## manually select CT clusters and subset it 
aa$celltype = NA

aa$celltype[!is.na(match(aa$seurat_clusters, 
                         c('0', '1', '5', '10', '13', '16',
                           '4', '7', '6', '9', '12', '25', 
                           '3', '8')))] = 'CT'

DimPlot(aa, group.by = 'celltype')

aa = subset(aa, subset = celltype == 'CT')

p1 = DimPlot(aa, group.by = 'seurat_clusters')
p2 = DimPlot(aa, group.by = 'condition', reduction = 'umap', cols = cols, 
             label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered", 
                         "_1roundCTselect.pdf"),
       width = 18, height = 6)


aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50,  min.dist = 0.3)

DimPlot(aa, group.by = 'condition', reduction = 'umap', cols = cols, label = TRUE, repel = TRUE)

saveRDS(aa, file = paste0(RdataDir,
                          '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.rds'))


##########################################
# downsample 3k cells for each condition and 2rd round of CT cleaning 
##########################################
aa = readRDS(file = paste0(RdataDir,
                           '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.rds'))

#aa = subset(aa, cells = colnames(aa)[which(aa$condition != "PLC_8dpd_IL11treated_1" & 
#                                             aa$condition != "PLC_8dpd_Dexatreated_1")])

#aa$condition = droplevels(aa$condition)
Idents(aa) = aa$condition
aa = subset(x = aa, downsample = 3000)

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50,  min.dist = 0.3)

DimPlot(aa,  reduction = 'umap',  group.by = 'condition', label = TRUE, repel = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:30)

aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

DimPlot(aa,  reduction = 'umap',  label = TRUE, repel = TRUE)



p1 = DimPlot(aa,  reduction = 'umap', cols = cols,
             group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa,  reduction = 'umap',  label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_rmDFout_CTselected.downsampled_umap_clusters", 
                         ".pdf"),
       width = 18, height = 6)

saveRDS(aa, file = paste0(RdataDir,
                          '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.downsampled.rds'))


markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

saveRDS(markers, file = paste0(RdataDir, 
                               'Andre_PLCscRNAseq_QCsfiltered_CTselected.downsampled.rds'))

## reload the calculated clusters and markers
cat(length(unique(markers$gene)), ' markers found in ax\n')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.0) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#ggs = top10$gene[which(top10$cluster == 'immune_cells')]
#ggs = ggs[grep('^AME|^LOC', ggs, invert = TRUE)]
#ggs = ggs[order(ggs)]
xx = subset(aa, downsample = 200)

DoHeatmap(xx, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, 
                         '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected_heatmap.pdf'), 
       width = 20, height = 30)

## double check the CT markers
features = c("PRRX1", "PDGFRA", "DPT", "ACTA2", "LUM", "COL1A2", "MMP2", 'POSTN', "COL3A1")
FeaturePlot(aa, features = features)

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered", 
                         "_CTselected_CTfeatures.pdf"),
       width = 16, height = 12)

FeaturePlot(aa, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered", 
                         "_CTselected_lowSequencingClusters.pdf"),
       width = 12, height = 8)

VlnPlot(aa, features = 'nFeature_RNA', group.by = 'condition', pt.size = 0)

ggsave(filename = paste0(resDir, "AndrePLCscRNAseq_QCs_doubletFiltered", 
                         "_CTselected_nCount_RNA.pdf"),
       width = 12, height = 8)



aa = subset(aa, cells = colnames(aa)[which(aa$seurat_clusters != '17' &
                                             aa$seurat_clusters != '18')])

#aa = subset(aa, cells = colnames(aa)[which(aa$condition != "PLC_5dpd_0" &
#                                             aa$condition != "PLC_10dpd_1" &
#                                             aa$condition != "PLC_12dpd_1")])

#aa$condition = droplevels(aa$condition)

aa = NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50,  min.dist = 0.3)

DimPlot(aa,  reduction = 'umap',  group.by = 'condition', label = TRUE, repel = TRUE, cols = cols)

ggsave(paste0(resDir, 'Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.downsampled_clean.pdf'), 
       width = 12, height = 8)

saveRDS(aa, file = paste0(RdataDir,
                          '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.downsampled_clean.rds'))


##########################################
# test batch correction
##########################################
aa = readRDS(paste0(RdataDir,
                    '/Andre_PLCscRNAseq_QCsfiltered_rmDFout_CTselected.downsampled_clean.rds'))

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source(paste0(functionDir, '/functions_dataIntegration.R'))

aa$batch = 'PLC'
aa$batch[aa$condition == 'MatLimb_0dpa_1' | aa$condition == "Blastema_11dpa_1"] = 'blastema' 

aa$batch = factor(aa$batch, levels = c('blastema', 'PLC'))

xx = IntegrateData_Seurat_RPCA(seuratObj = aa, group.by = 'batch', nfeatures = 3000)

# Visualization
ElbowPlot(xx, ndims = 50)
xx = RunUMAP(xx, reduction = "pca", dims = 1:30, n.neighbors = 50, 
             min.dist = 0.3) 
DimPlot(xx, reduction = "umap", group.by = "condition", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols, label.size = 5) 

ggsave(paste0(resDir, '/Integration_PLC_Blastema_SeuratRPCA_3000HVGs.pdf'), 
       width = 12, height = 8)


xx = IntegrateData_Seurat_CCA(seuratObj = aa, group.by = 'batch', nfeatures = 3000)

# Visualization
ElbowPlot(xx, ndims = 50)
xx = RunUMAP(xx, reduction = "pca", dims = 1:30, n.neighbors = 50, 
             min.dist = 0.3) 
DimPlot(xx, reduction = "umap", group.by = "condition", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols, label.size = 5) 

ggsave(paste0(resDir, '/Integration_PLC_Blastema_SeuratRPCA_3000HVGs.pdf'), 
       width = 12, height = 8)


xx = IntegrateData_Seurat_CCA(seuratObj = aa, group.by = 'batch', nfeatures = 3000)

# Visualization
ElbowPlot(xx, ndims = 50)
xx = RunUMAP(xx, reduction = "pca", dims = 1:30, n.neighbors = 30, 
             min.dist = 0.3) 
DimPlot(xx, reduction = "umap", group.by = "condition", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols, label.size = 5) 

ggsave(paste0(resDir, '/Integration_PLC_Blastema_SeuratCCA_3000HVGs.pdf'), 
       width = 12, height = 8)


#Run harmony to correct for chemistry
library(harmony)
tic()
xx<- RunHarmony(aa, group.by.vars = 'batch', 
                
                #nclust = 10, 
                max.iter.harmony = 10, 
                #epsilon.harmony = -Inf,
                epsilon.harmony = 0.0001,
                verbose = TRUE,
                #sreference_values = 'blastema',
                plot_convergence = TRUE)
toc()

#Identify the highest contributing PCs
ElbowPlot(xx, ndims = 50)
npcs = 30
xx <- RunUMAP(xx, reduction = 'harmony', dims = 1:npcs, 
              reduction.name = 'umap_harmony')

DimPlot(xx, reduction = "umap_harmony", label = TRUE, group.by = 'condition',
        cols = cols)

ggsave(paste0(resDir, '/Integration_PLC_Blastema_RunHarmony_3000HVGs.pdf'), 
       width = 12, height = 8)







