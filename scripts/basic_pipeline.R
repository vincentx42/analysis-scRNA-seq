rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)

setwd("/Users/vincent/bioinfomatics_projects/scRNA_analysis_reproduction")

####load expression matrices as Seurat objects####
source_dir = "./data/GSE148071/"
files = list.files(path = source_dir, pattern = "*_exp.txt", all.files = F)

n = 1 
for (i in files) {
  counts = read.table(paste0(source_dir, i))
  tmp <- CreateSeuratObject(counts,
                            min.cells = 3,
                            project = 'scRNA_lung_cancer',
                            min.features = 300)
  assign(paste0("scRNA_P", n), tmp)
  n = n + 1
}

# combine all objects together
mergelist = grep("scRNA_P*", ls(), value =T)
mergelist = mergelist[mergelist != "scRNA_P1"]
seuratlist = list()
for (scname in mergelist) {
  obj = get(scname)
  seuratlist = append(seuratlist, obj)
}

# merge
scRNA_comb = merge(x=scRNA_P1, y=seuratlist, 
                   add.cell.ids = seq(1:42), 
                   project = "scRNA_lung_cancer")

# view cell counts
table(scRNA_comb@meta.data$orig.ident)

####Quality control and Normalization####
# calculate proportion of mitochondria
scRNA_comb[['percent.mt']] <- PercentageFeatureSet(scRNA_comb, 
                                                   pattern = '^MT-')

# calculate proportion of RBCs and remove related genes
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_comb@assays$RNA)) 
HB.genes <- rownames(scRNA_comb@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 

scRNA_comb[["percent.HB"]] <- PercentageFeatureSet(scRNA_comb, 
                                                   features = HB.genes) 

# cell composition visualization
# scatter plot
plot1 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot1 + plot2 + plot3

# apply filter
scRNA_afterqc <- subset(scRNA_comb, subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 5000 &
                          percent.mt < 30 & 
                          percent.HB < 3 & 
                          nCount_RNA < 30000)

# normalization
scRNA_normalized <- NormalizeData(scRNA_afterqc, 
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)
# save normalized data
save(scRNA_normalized, file = './data/scRNA_normalized.Rdata')

####Feature selection####
# load(file = './data/scRNA_normalized.RData')
scRNA <- FindVariableFeatures(scRNA_normalized, 
                              selection.method = "vst", 
                              nfeatures = 600) 

# view top10 features
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 

plot3 = plot1 + plot2
ggsave("./results/high_variable_genes.pdf", plot3, width = 10, height = 5)

####Scaling the data####
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

####PCA####
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 

# choose number of PCs
ElbowPlot(scRNA, ndims=50, reduction="pca") 

####Eliminate batch effect with harmony####
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")

# use harmony data to find clusters
scRNA_harmony <- FindNeighbors(scRNA_harmony, 
                               reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1.0)

table(scRNA_harmony@meta.data$seurat_clusters)
####Dimension reduction with UMAP####
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

plot1 = DimPlot(scRNA_harmony, reduction = "umap", label = T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
plotc <- plot1 + plot2
plotc

# save harmony data
save(scRNA_harmony, file = './data/scRNA_harmony.Rdata')

####find markers####
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers = markers %>% 
  dplyr::select(gene, everything()) %>% 
  subset(p_val<0.05)

top10 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

scRNA_harmony@meta.data$seurat_clusters

####singleR cell types annotation####
library(SingleR)

# load reference dataset
refdata = celldex::HumanPrimaryCellAtlasData()

# extract assay data and cluster info
counts_sc <- GetAssayData(scRNA_harmony, slot="data")
clusters <- scRNA_harmony@meta.data$seurat_clusters
# invoke singler
cell_anno <- SingleR(test = counts_sc, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

# export annotation info
celltype = data.frame(ClusterID = rownames(cell_anno), 
                      celltype = cell_anno$labels, 
                      stringsAsFactors = FALSE)

write.csv(celltype, 
          "./results/anno_singler.csv",
          row.names = FALSE)

# add annotation to seurat object
scRNA_harmony@meta.data$singleR = celltype[match(clusters, celltype$ClusterID), 'celltype']
# view cell types after annotation
p1 = DimPlot(scRNA_harmony, group.by="singleR", label=T, label.size=5)

p2 = DimPlot(scRNA_harmony, group.by="seurat_clusters", label=T, label.size=5, reduction='umap')

ggsave("./results/singler_umap_celltype.pdf", p1, width=10, height=5)
ggsave("./results/find_cluster_umap_celltype.pdf", p2, width=10, height=5)

# view marker genes
genelist = top10[top10$p_val == 0, ]
genelist = arrange(genelist, -genelist$avg_log2FC)
genelist = genelist[1:20, ]

p1 = DotPlot(scRNA_harmony, features = genelist$gene, group.by = "singleR") +
  RotatedAxis() +
  scale_x_discrete("") +
  scale_y_discrete("")
ggsave("./results/marker_genes_20.pdf", p1, width=10, height=5)


save(scRNA_harmony, file = './data/scRNA_harmony_annotation.Rdata')



