library(Seurat)

x<-readRDS("rsi_3-split-250.rds")
Idents(x)<-x$limb
f<-subset(x, idents = "F")
h<-subset(x, idents = "H")
Idents(f)<-f$cell.type3
f
Idents(h)<-h$cell.type3
h

f2<-subset(f, idents=c("EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast"))
h2<-subset(h, idents=c("EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast"))

DefaultAssay(f2)<-"integrated"
f2 <- RunPCA(object = f2, verbose = FALSE)
f2 <- FindNeighbors(object = f2, dims = 1:30, verbose = FALSE)
f2 <- RunUMAP(object = f2, dims = 1:30, verbose = FALSE)
DimPlot(f2)
DimPlot(f2, split.by = "cell.type3")

DefaultAssay(h2)<-"integrated"
h2 <- RunPCA(object = h2, verbose = FALSE)
h2 <- FindNeighbors(object = h2, dims = 1:30, verbose = FALSE)
h2 <- RunUMAP(object = h2, dims = 1:30, verbose = FALSE)
DimPlot(h2)

f2$barcode <- colnames(f2)
f2$UMAP_1 <- f2@reductions$umap@cell.embeddings[,1]
f2$UMAP_2 <- f2@reductions$umap@cell.embeddings[,2]
write.csv(f2@meta.data, file='./velo/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(f2, assay='RNA', slot='counts')
writeMM(counts_matrix, file='./velo/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(f2@reductions$pca@cell.embeddings,file='./velo/pdmp/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./velo/gene_names.csv',
  quote=F,row.names=F,col.names=F
)

h2$barcode <- colnames(h2)
h2$UMAP_1 <- h2@reductions$umap@cell.embeddings[,1]
h2$UMAP_2 <- h2@reductions$umap@cell.embeddings[,2]
write.csv(h2@meta.data, file='./velo/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(h2, assay='RNA', slot='counts')
writeMM(counts_matrix, file='./velo/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(h2@reductions$pca@cell.embeddings,file='./velo/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./velo/gene_names.csv',
  quote=F,row.names=F,col.names=F
)


