rm(list = ls())
load('./data/scRNA_harmony_annotation.Rdata')

####trajectory analysis with monocle3#####
library(monocle3)
expression_matrix = GetAssayData(scRNA_harmony, assay = 'RNA', slot = 'counts')
cell_metadata = scRNA_harmony@meta.data
gene_annotation = data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) = gene_annotation$gene_short_name

# construct CDS
cds = new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)
save(cds, file = './results/initial_cds.Rdata')

# pre-processing
load('./results/initial_cds.Rdata')
cds = preprocess_cds(cds, num_dim = 50)

# umap dimension reduction
cds = reduce_dimension(cds, preprocess_method = 'PCA')

p1 = plot_cells(cds, reduction_method = 'UMAP', color_cells_by = "singleR") +
  ggtitle('cds-umap')

ggsave('./results/UMAP_monocle3.pdf', p1, width = 10, height = 5)

# cluster the cells with unsupervised method
cds = cluster_cells(cds)
p2 = plot_cells(cds, color_cells_by = "partition")
p2

# trajectory analysis
cds = learn_graph(cds)
p3 = plot_cells(cds, color_cells_by = "singleR", label_groups_by_cluster = F,
                label_leaves = F, label_branch_points = F)

ggsave('./results/trajectory.pdf', p3, width = 10, height = 5)

# pseudotime analysis
cds = order_cells(cds)
p4 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
                label_leaves = F, label_branch_points = F)

ggsave('./results/pseudotime.pdf', p4, width = 10, height = 5)