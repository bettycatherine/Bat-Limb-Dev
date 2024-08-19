library(cellAlign)
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
h<-subset(x, idents="H")
Idents(h)<-h$cell.type3
f.b<-subset(f, idents=c("EMP", "OP", "Osteoblast"))
##or chondrogenesis as: f.b<-subset(f, idents=c("EMP", "CP", "Chondrocyte"))
h.b<-subset(h, idents=c("EMP", "OP", "Osteoblast"))

DefaultAssay(f.b)<-"RNA"
cds.rf <- as.cell_data_set(f.b)
fData(cds.rf)$gene_short_name <- rownames(fData(cds.rf))
head(fData(cds.rf))
cds.rf <- preprocess_cds(cds.rf, num_dim = 50)
cds.rf <- align_cds(cds.rf, alignment_group = "method")
cds.rf <- reduce_dimension(cds.rf, preprocess_method = "Aligned")
cds.rf <- cluster_cells(cds.rf, reduction_method = "Aligned")
cds.rf <- cluster_cells(cds.rf, reduction_method = "UMAP")
cds.rf <- learn_graph(cds.rf, learn_graph_control=list(ncenter=2), use_partition = F)
plot_cells(cds.rf,color_cells_by = "stage")+
plot_cells(cds.rf,label_principal_points = T)
cds.rf <- order_cells(cds.rf, root_pr_nodes = 'Y_1')
plot_cells(cds.rf,color_cells_by = "pseudotime")
gene_test.rf <- graph_test(cds.rf, neighbor_graph="principal_graph")
whichgene.rf <- row.names(subset(gene_test.rf, q_value < 0.05))
traj.f<- cds.rf@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

DefaultAssay(h.b)<-"RNA"
cds.rh <- as.cell_data_set(h.b)
fData(cds.rh)$gene_short_name <- rownames(fData(cds.rh))
head(fData(cds.rh))
cds.rh <- preprocess_cds(cds.rh, num_dim = 50)
cds.rh <- align_cds(cds.rh, alignment_group = "method")
cds.rh <- reduce_dimension(cds.rh, preprocess_method = "Aligned")
cds.rh <- cluster_cells(cds.rh, reduction_method = "Aligned")
cds.rh <- cluster_cells(cds.rh, reduction_method = "UMAP")
cds.rh <- learn_graph(cds.rh, learn_graph_control=list(ncenter=2), use_partition = F)
plot_cells(cds.rh,color_cells_by = "stage")+
  plot_cells(cds.rh,label_principal_points = T)
cds.rh <- order_cells(cds.rh, root_pr_nodes = 'Y_1')
plot_cells(cds.rh,color_cells_by = "pseudotime")
gene_test.rh <- graph_test(cds.rh, neighbor_graph="principal_graph", cores=4)
whichgene.rh <- row.names(subset(gene_test.rh, q_value < 0.05))
#whichgene.rh <- row.names(gene_test.rh)[order(gene_test.rh$qval)][1:500]
traj.h<- cds.rh@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

DefaultAssay(f.b)<-"RNA"
f.b <- FindVariableFeatures(f.b, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(f.b)
f.b <- ScaleData(f.b, features = all.genes)
expr.f <- GetAssayData(f.b, assay='RNA', slot='scale.data')
expr.f<-as.data.frame(expr.f)
#expr.f[1:4,1:5]

DefaultAssay(h.b)<-"RNA"
h.b <- FindVariableFeatures(h.b, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(h.b)
h.b <- ScaleData(h.b, features = all.genes)
expr.h <- GetAssayData(h.b, assay='RNA', slot='scale.data')
expr.h<-as.data.frame(expr.h)

whichgene<-intersect(whichgene.rf,whichgene.rh)
whichgene<-unlist(whichgene)
expr.f<-expr.f[whichgene,]
expr.h<-expr.h[whichgene,]

expr.f<-expr.f[whichgene.rh,]
expr.h<-expr.h[whichgene.rh,]

cluster<-as.data.frame(cluster)
whichgene.local<-row.names(cluster[which(cluster$cluster==4),])
whichgene.local<-unlist(whichgene.local)
expr.f.local<-expr.f[whichgene.local,]
expr.h.local<-expr.h[whichgene.local,]

traj.f.norm<-(traj.f-min(traj.f))/(max(traj.f)-min(traj.f))
traj.h.norm<-(traj.h-min(traj.h))/(max(traj.h)-min(traj.h))

traj.f.norm <- as.numeric(t(traj.f.norm))
names(traj.f.norm) <- colnames(expr.f)
traj.h.norm <- as.numeric(t(traj.h.norm))
names(traj.h.norm) <- colnames(expr.h)

traj.f.norm.local <- as.numeric(t(traj.f.norm))
names(traj.f.norm.local) <- colnames(expr.f.local)
traj.h.norm.local <- as.numeric(t(traj.h.norm))
names(traj.h.norm.local) <- colnames(expr.h.local)

##interpolation
numPts <- 200
inter.f <- cellAlign::interWeights(expDataBatch = expr.f, trajCond = traj.f.norm,
                                         winSz = 0.1, numPts = numPts)
inter.h <- cellAlign::interWeights(expDataBatch = expr.h, trajCond = traj.h.norm,
                                         winSz = 0.1, numPts = numPts)
#scale the interpolated data (Recommended):
interScaled.f = cellAlign::scaleInterpolate(inter.f)
interScaled.h = cellAlign::scaleInterpolate(inter.h)
alignment = globalAlign(interScaled.h$scaledData,interScaled.f$scaledData, 
                        scores = list(query = interScaled.h$traj, 
                                      ref = interScaled.f$traj),
                        sigCalc = F, numPerm = 20
)
plotAlign(alignment)




