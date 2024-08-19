library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

x<-readRDS("rsi_3-split-250.rds")

#subset x
Idents(x)<-x$limb
f<-subset(x, idents="F")
Idents(f)<-f$cell.type3
f2<-f
f2<-subset(f, idents=c("EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast"))
cds.rf <- as.cell_data_set(f2)
fData(cds.rf)$gene_short_name <- rownames(fData(cds.rf))
head(fData(cds.rf))

h<-subset(x, idents="H")
Idents(h)<-h$cell.type3
h2<-h
h2<-subset(h, idents=c("EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast"))
DefaultAssay(f6)<-"RNA"


##retrieve umap coordinate from seurat
recreate.partitions <- c(rep(1, length(cds.rf@colData@rownames)))
names(recreate.partitions) <- cds.rf@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds.rf@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- f2@active.ident

cds.rf@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds.rf@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- f2@reductions$umap@cell.embeddings
plot_cells(cds.rf, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")

cds.rf <- learn_graph(cds.rf, use_partition = F, close_loop = F)
plot_cells(cds.rf, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)+

plot_cells(cds.rf, color_cells_by = "stage", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
cds.rf <- order_cells(cds.rf, reduction_method = "UMAP")
plot_cells(cds.rf, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = T,
           label_roots = F, 
           group_label_size = 6, 
           label_cell_groups=T,
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 0)

cds.rf$monocle3_pseudotime <- pseudotime(cds.rf)
f.pseudo <- as.data.frame(colData(cds.rf))
FeaturePlot(f2, features = "pseudotime")

##hind
#h2
DefaultAssay(h2)<-"RNA"
cds.rh <- as.cell_data_set(h2)
fData(cds.rh)$gene_short_name <- rownames(fData(cds.rh))
head(fData(cds.rh))

##retrieve umap coordinate from seurat
recreate.partitions <- c(rep(1, length(cds.rh@colData@rownames)))
names(recreate.partitions) <- cds.rh@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds.rh@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- h2@active.ident
cds.rh@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds.rh@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- h2@reductions$umap@cell.embeddings
plot_cells(cds.rh, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")#+facet_wrap(~cell.type2)

cds.rh <- learn_graph(cds.rh, use_partition = F)
cds.rh <- learn_graph(cds.rh, close_loop = F, use_partition = F)
plot_cells(cds.rh, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)+

plot_cells(cds.rh, color_cells_by = "stage", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
cds.rh <- order_cells(cds.rh, reduction_method = "UMAP")
plot_cells(cds.rh, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = T,
           label_roots = F, 
           group_label_size = 6, 
           label_cell_groups=T,
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 0)

cds.rh$monocle3_pseudotime <- pseudotime(cds.rh)
h.pseudo <- as.data.frame(colData(cds.rh))
FeaturePlot(h2, features = "pseudotime")


