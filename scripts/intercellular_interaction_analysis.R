rm(list = ls())
load('./data/scRNA_harmony_annotation.Rdata')
library(SingleCellExperiment)
scRNA_harmony.sce = as.SingleCellExperiment(scRNA_harmony)

library(Seurat)
library(SeuratDisk)

####pre-processing####
# extract information from seurat object with annotation
# convert Seurat object to 'h5ad' format as counts file for cellphonedb
SaveH5Seurat(scRNA_harmony, filename = './data/scRNA_harmony.h5seurat')
Convert('./data/scRNA_harmony.h5seurat', dest = 'h5ad')

# generate annotation file for cellphonedb
anno_file = cbind(rownames(scRNA_harmony@meta.data), 
                  scRNA_harmony@meta.data[, 'singleR', drop = F])

write.table(as.matrix(anno_file), './data/scRNA_harmony_anno.txt', sep = '\t',
            quote = F, row.names = F)


####visulization####
library(CellChat)
library(tidyr)
library(ktplots)

pvals <- read.delim("./results/cellphonedb_output/pvalues.txt", check.names = FALSE)
means <- read.delim("./results/cellphonedb_output/means.txt", check.names = FALSE)
decon <- read.delim("./results/cellphonedb_output/deconvoluted.txt", check.names = FALSE)
sig_means <- read.delim("./results/cellphonedb_output/significant_means.txt", check.names = FALSE)

heatmap_cpdb = plot_cpdb_heatmap(pvals, cellheight = 30, cellwidth = 30, symmetrical = F)
ggsave("./results/cpdb_heatmap.pdf", heatmap_cpdb, width = 10, height = 5)


chord_cpdb <- plot_cpdb3(cell_type1 = 'Endothelial_cells', cell_type2 = 'Macrophage|T_cells|B_cell',
                scdata = scRNA_harmony.sce,
                celltype_key = 'singleR', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = decon, # new options from here on specific to plot_cpdb3
                keep_significant_only = TRUE,
                standard_scale = TRUE,
                remove_self = TRUE)


