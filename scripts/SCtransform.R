rm(list = ls())

# load packages
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)
library(readxl)

setwd("/Users/vincent/bioinfomatics_projects/scRNA_analysis_reproduction")

####Preparation####
# sample list
samples <- read_excel("./data/sc_LUAD/patients_metadata.xlsx", 
                      range = cell_cols("A:A")) %>% .$sample_id

### import cellranger files from different data sets
for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("./data/sc_LUAD/", samples[i], "/filtered_feature_bc_matrix")))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6, seu_obj7, seu_obj8, seu_obj9, seu_obj10, seu_obj11, seu_obj12, seu_obj13, seu_obj14, seu_obj15, seu_obj16, seu_obj17, seu_obj18, seu_obj19, seu_obj20), add.cell.ids = samples, project = "lung")

# save(seu_obj, file = "./data/sc_LUAD/scRNA.Rdata")

#load("./data/sc_LUAD/scRNA.Rdata") 
#scRNA_comb = seu_obj

####Quality control and Normalization####
# calculate proportion of mitochondria / RBCs / ribosome
scRNA_comb[['percent.mt']] <- PercentageFeatureSet(scRNA_comb, 
                                                   pattern = '^MT-')
scRNA_comb[['percent.HB']] <- PercentageFeatureSet(scRNA_comb, 
                                                   pattern = "^HBA|^HBB")
scRNA_comb[['percent.RP']] <- PercentageFeatureSet(scRNA_comb, 
                                                   pattern = "^RPS|^RPL")

# cell composition visualization
# scatter plot
plot1 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 = FeatureScatter(scRNA_comb, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot1 + plot2 + plot3

# violin plot
VlnPlot(scRNA_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
ggsave2("violin_percent_mt.pdf", path = "./results2/", width = 10, height = 5)

# apply filter
scRNA_afterqc <- subset(scRNA_comb, subset = nFeature_RNA > 500 & 
                          nFeature_RNA < 10000 &
                          percent.mt < 30 & 
                          percent.HB < 5 & 
                          nCount_RNA > 1000 &
                          nCount_RNA < 100000)

# use SCtransform to normalize data
scRNA <- SCTransform(scRNA_afterqc, verbose = T, 
                     vars.to.regress = c("nCount_RNA", "percent.mt"), 
                     conserve.memory = T)

save(scRNA, file = "./data/sc_LUAD/scRNA_SCtransform.Rdata")


####Dimensionality reduction####
scRNA <- RunPCA(scRNA)
ElbowPlot(scRNA, ndims = 50)
scRNA <- RunUMAP(scRNA, dims = 1:25, verbose = FALSE)

# find clusters
scRNA <- FindNeighbors(scRNA, dims = 1:25, verbose = FALSE)
scRNA <- FindClusters(scRNA, resolution = 0.5, verbose = FALSE)
plot1 = DimPlot(scRNA, label = TRUE) + NoLegend()
ggsave("./results2/umap.pdf", plot1, width = 10, height = 5)

####Find markers####
DefaultAssay(scRNA) = "SCT"
scRNA = PrepSCTFindMarkers(scRNA)

markers <- FindAllMarkers(object = scRNA, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers = markers %>% 
  dplyr::select(gene, everything()) %>% 
  subset(p_val<0.05 & abs(markers$avg_log2FC) > 1)

top10 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# visualization
p1 = DoHeatmap(scRNA, features = top10$gene)
ggsave("./results2/heatmap_genes.pdf", p1, width = 30, height = 15)

####singleR cell types annotation####
library(SingleR)

# load reference dataset
refdata = celldex::HumanPrimaryCellAtlasData()

# extract assay data and cluster info
counts_sc <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
# invoke singler
cell_anno <- SingleR(test = counts_sc, ref = refdata, labels = refdata$label.main, 
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")

# export annotation info
celltype = data.frame(ClusterID = rownames(cell_anno), 
                      celltype = cell_anno$labels, 
                      stringsAsFactors = FALSE)

# add annotation to seurat object
scRNA@meta.data$singleR = celltype[match(clusters, celltype$ClusterID), 'celltype']
# view cell types after annotation
p1 = DimPlot(scRNA, group.by="singleR", label=T, label.size=5, reduction = 'umap')

p2 = DimPlot(scRNA, group.by="seurat_clusters", label=T, label.size=5, reduction='umap')

ggsave("./results2/singler_umap_celltype.pdf", p1, width=10, height=5)
ggsave("./results2/find_cluster_umap_celltype.pdf", p2, width=10, height=5)

# view marker genes
genelist = unique(top10$gene)
p1 = DotPlot(scRNA, features = genelist, group.by = "singleR") +
  RotatedAxis() +
  scale_x_discrete("") +
  scale_y_discrete("")

ggsave("./results2/dot_genes_annotation.pdf", p1, width = 30, height = 10)

p2 = DoHeatmap(scRNA, features = genelist, group.by = "singleR")

ggsave("./results2/heatmap_genes_annotation.pdf", p2, width = 30, height = 15)

# save(scRNA, file = "./data/sc_LUAD/scRNA_final.Rdata")

