##########################################################################
##########################################################################
# Project: axolotl WE (wound epidermis) trajectory 
# Script purpose: process and analyze trajectory of Leigh et al., 2018
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri May 19 11:01:01 2023
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library("viridis")

version.analysis = '_axolotl_Leigh2018_20230519'
resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

dataDir = '/groups/tanaka/People/current/jiwang/projects/limbRegeneration_scRNA/raw_NGS/axolotl/Leigh_2018'
metaDir = paste0(dataDir, '/indrops/Leigh_et_al_2018_Supplementary_R_code')

########################################################
########################################################
# Section I: import processed Seurat object from Steven Blair (Whited lab)
# 
########################################################
########################################################
load(paste0(dataDir, '/processed_data/nLeigh_seuratConverted_natCom2018.rdata'))

## mature sample cell types
meta = read.table(file = paste0(metaDir, '/homeostasis.meta.data.txt'), header = TRUE)
aa = seu_intact
aa$celltypes = NA
aa$celltypes = meta$Cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)

seu_intact = aa

rm(aa)

# wound healing
meta = read.table(file = paste0(metaDir, '/wound_healing.meta.data.txt'), header = TRUE)
aa = seu_3dpa
aa$celltypes = NA
aa$celltypes = meta$cell_type[match(colnames(aa), meta$Cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_3dpa = aa

rm(aa)

# 14dpa
meta = read.table(file = paste0(metaDir, '/early_bud_blastema.meta.data.txt'), header = TRUE)
aa = seu_14dpa
aa$celltypes = NA
aa$celltypes = meta$cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_14dpa = aa

rm(aa)

# 23dpa
meta = read.table(file = paste0(metaDir, '/medium_bud_blastema.meta.data.txt'), header = TRUE)
aa = seu_23dpa
aa$celltypes = NA
aa$celltypes = meta$Cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_23dpa = aa

rm(aa)
rm(meta)

saveRDS(seu_intact, file = paste0(RdataDir, '/seuratObj_processedWhited_0dpa.rds'))
saveRDS(seu_3dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_3dpa.rds'))
saveRDS(seu_14dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_14dpa.rds'))
saveRDS(seu_23dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_23dpa.rds'))


##########################################
# change the gene names 
##########################################
#seu_dpa0 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_0dpa.rds'))
#aa = seu_dpa0
#counts = aa@assays$RNA@counts
#metadata = aa@meta.data
#seu_dpa3 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_3dpa.rds'))
#seu_dpa14 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_14dpa.rds'))
#seu_dpa23 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_23dpa.rds'))

########################################################
########################################################
# Section II: check marker genes 
# 
########################################################
########################################################
seu_dpa0 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_0dpa.rds'))
#aa = seu_dpa0
#counts = aa@assays$RNA@counts
#metadata = aa@meta.data
seu_dpa3 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_3dpa.rds'))
seu_dpa14 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_14dpa.rds'))
seu_dpa23 = readRDS(file = paste0(RdataDir, '/seuratObj_processedWhited_23dpa.rds'))


p1 = DimPlot(seu_dpa0, group.by = 'celltypes', label = TRUE, repel = TRUE) + ggtitle('dpa0') +
  NoLegend()
p2 = DimPlot(seu_dpa3, group.by = 'celltypes', label = TRUE, repel = TRUE) + ggtitle('dpa3') +
  NoLegend()
p3 = DimPlot(seu_dpa14, group.by = 'celltypes', label = TRUE, repel = TRUE) + ggtitle('dpa14') +
  NoLegend()
p4 = DimPlot(seu_dpa23, group.by = 'celltypes', label = TRUE, repel = TRUE) + ggtitle('dpa23') +
  NoLegend()

(p1 + p2)/(p3 + p4) 

ggsave(filename = paste0(resDir, '/tSNE_Whited_celltypes.pdf'), 
       width = 16, height = 12)

### dpa0 Wnt1, Wnt2, Wnt8a and Wnt8b
aa = seu_dpa0
features = c(rownames(aa)[grep('WNT3A', rownames(aa))],
             rownames(aa)[grep('WNT5A', rownames(aa))],
             rownames(aa)[grep('WNT1_', rownames(aa))],
             rownames(aa)[grep('WNT2_', rownames(aa))],
             rownames(aa)[grep('WNT8_', rownames(aa))],
             rownames(aa)[grep('WNT8A_', rownames(aa))],
             rownames(aa)[grep('WNT8B_', rownames(aa))],
             rownames(aa)[grep('TGFB', rownames(aa))],
             rownames(aa)[grep('INHA', rownames(aa))])

ggs = sapply(features, function(x){x = unlist(strsplit(as.character(x), '[|]'))[3]; 
                                    x = gsub("\\^.*",'', x); 
                                    as.character(x)})
#ggs[length(ggs)] = "INHA"

p = FeaturePlot(aa, 
            features = features,
            ncol = 3)

# Loop over the individual plots (assuming p[[1]] to p[[N]])
for(i in  seq_along(p)){
  p[[i]] <- p[[i]] + labs(title = ggs[i])
}

print(p)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa0_featureplot.pdf'), 
       width = 12, height = 10)

p1 = VlnPlot(aa, features = features, group.by = 'celltypes', ncol = 1)
for(i in  seq_along(p1)){
  p1[[i]] <- p1[[i]] + labs(title = ggs[i])
}

print(p1)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa0_VlnPlot.pdf'), 
       width = 8, height = 28)

write.csv2(features, file = paste0(resDir, '/Whiteddataset_dpa0_featuresList.csv'), 
           row.names = FALSE)


### dpa3
aa = seu_dpa3
features = c(rownames(aa)[grep('WNT3A', rownames(aa))],
             rownames(aa)[grep('WNT5A', rownames(aa))],
             rownames(aa)[grep('WNT1_', rownames(aa))],
             rownames(aa)[grep('WNT2_', rownames(aa))],
             rownames(aa)[grep('WNT8_', rownames(aa))],
             rownames(aa)[grep('WNT8A_', rownames(aa))],
             rownames(aa)[grep('WNT8B_', rownames(aa))],
             rownames(aa)[grep('TGFB', rownames(aa))],
             rownames(aa)[grep('INHA', rownames(aa))])

ggs = sapply(features, function(x){x = unlist(strsplit(as.character(x), '[|]'))[3]; 
x = gsub("\\^.*",'', x); 
as.character(x)})
#ggs[length(ggs)] = "INHA"

p = FeaturePlot(aa, 
                features = features,
                ncol = 3)

# Loop over the individual plots (assuming p[[1]] to p[[N]])
for(i in  seq_along(p)){
  p[[i]] <- p[[i]] + labs(title = ggs[i])
}

p1 = VlnPlot(aa, features = features, group.by = 'celltypes', ncol = 1)
for(i in  seq_along(p1)){
  p1[[i]] <- p1[[i]] + labs(title = ggs[i])
}

print(p)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa3_featureplot.pdf'), 
       width = 12, height = 10)

print(p1)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa3_VlnPlot.pdf'), 
       width = 8, height = 28)

write.csv2(features, file = paste0(resDir, '/Whiteddataset_dpa3_featuresList.csv'), 
           row.names = FALSE)

### dpa14
aa = seu_dpa14
features = c(rownames(aa)[grep('WNT3A', rownames(aa))],
             rownames(aa)[grep('WNT5A', rownames(aa))],
             rownames(aa)[grep('WNT1_', rownames(aa))],
             rownames(aa)[grep('WNT2_', rownames(aa))],
             rownames(aa)[grep('WNT8_', rownames(aa))],
             rownames(aa)[grep('WNT8A_', rownames(aa))],
             rownames(aa)[grep('WNT8B_', rownames(aa))],
             rownames(aa)[grep('TGFB', rownames(aa))],
             rownames(aa)[grep('INHA', rownames(aa))])

ggs = sapply(features, function(x){x = unlist(strsplit(as.character(x), '[|]'))[3]; 
x = gsub("\\^.*",'', x); 
as.character(x)})
#ggs[length(ggs)] = "INHA"

p = FeaturePlot(aa, 
                features = features,
                ncol = 3)

# Loop over the individual plots (assuming p[[1]] to p[[N]])
for(i in  seq_along(p)){
  p[[i]] <- p[[i]] + labs(title = ggs[i])
}

p1 = VlnPlot(aa, features = features, group.by = 'celltypes', ncol = 1)
for(i in  seq_along(p1)){
  p1[[i]] <- p1[[i]] + labs(title = ggs[i])
}

print(p)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa14_featureplot.pdf'), 
       width = 12, height = 10)


plot(p1)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa14_VlnPlot.pdf'), 
       width = 8, height = 28)

write.csv2(features, file = paste0(resDir, '/Whiteddataset_dpa14_featuresList.csv'), 
           row.names = FALSE)


### dpa23
aa = seu_dpa23
features = c(rownames(aa)[grep('WNT3A', rownames(aa))],
             rownames(aa)[grep('WNT5A', rownames(aa))],
             rownames(aa)[grep('WNT1_', rownames(aa))],
             rownames(aa)[grep('WNT2_', rownames(aa))],
             rownames(aa)[grep('WNT8_', rownames(aa))],
             rownames(aa)[grep('WNT8A_', rownames(aa))],
             rownames(aa)[grep('WNT8B_', rownames(aa))],
             rownames(aa)[grep('TGFB', rownames(aa))],
             rownames(aa)[grep('INHA', rownames(aa))])

ggs = sapply(features, function(x){x = unlist(strsplit(as.character(x), '[|]'))[3]; 
x = gsub("\\^.*",'', x); 
as.character(x)})

ggs[length(ggs)] = "INHA"


p = FeaturePlot(aa, 
                features = features,
                ncol = 3)

# Loop over the individual plots (assuming p[[1]] to p[[N]])
for(i in  seq_along(p)){
  p[[i]] <- p[[i]] + labs(title = ggs[i])
}

p1 = VlnPlot(aa, features = features, group.by = 'celltypes', ncol = 1)
for(i in  seq_along(p1)){
  p1[[i]] <- p1[[i]] + labs(title = ggs[i])
}

print(p)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa23_featureplot.pdf'), 
       width = 12, height = 14)

plot(p1)

ggsave(filename = paste0(resDir, '/Whiteddataset_tSNE_dpa23_VlnPlot.pdf'), 
       width = 8, height = 28)

write.csv2(features, file = paste0(resDir, '/Whiteddataset_dpa23_featuresList.csv'), 
           row.names = FALSE)


