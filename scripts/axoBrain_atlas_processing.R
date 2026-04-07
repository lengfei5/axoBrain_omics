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

inputDir = '/groups/tanaka/Collaborations/axolotlBrain_Omics/data_shared_Mateja/'
outDir = '/groups/tanaka/Collaborations/axolotlBrain_Omics/Rdata/'

if(!dir.exists(outDir)) dir.create(outDir)

##########################################
# convert the h5ad file to seurat object 
##########################################
library(SingleCellExperiment)
library(SeuratDisk)

# Convert h5ad file to h5seurat format
h5ad_file = paste0(inputDir, 'Retina.h5ad')
Convert(h5ad_file, dest = "h5seurat")

# load Seurat file
aa <- LoadH5Seurat(paste0(inputDir, 'Retina.h5seurat'))






