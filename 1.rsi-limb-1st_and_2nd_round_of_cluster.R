###f+h
library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(Rtsne)
library(sctransform)
library(ggplot2)

setwd("D:/whatever/working/project/analysis/rsi-limb/sc")
options(future.globals.maxSize=7000*1024^2)
hex.data<-read.table("rsi-hex-rsi_3_rescued_stages.txt", header = T, row.names = 1)
dt.data<-read.table("rsi-dt-rsi_3_rescued_stages.txt", header = T, row.names = 1)
hex<-CreateSeuratObject(counts = hex.data)
dt<-CreateSeuratObject(counts = dt.data)
dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^mt-")
hex[["percent.mt"]] <- PercentageFeatureSet(hex, pattern = "^mt-")
VlnPlot(hex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(dt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hex <- subset(hex, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5)
dt <- subset(dt, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 5)

rsi.merged <- merge(dt, y =hex, add.cell.ids = c("dt","hex"))
aaa <- stringr::str_split(colnames(rsi.merged),'_', simplify = T)[,1] 
head(aaa)
names(aaa) <- colnames(rsi.merged) 
table(aaa)
rsi.merged <- AddMetaData(object = rsi.merged, metadata = aaa, col.name = 'primer')

rsi.list <- SplitObject(rsi.merged, split.by = 'primer')
for (i in 1:length(rsi.list)) {
  rsi.list[[i]] <- SCTransform(rsi.list[[i]], verbose = FALSE)
}
rsi.features <- SelectIntegrationFeatures(object.list = rsi.list, nfeatures = 3000)
rsi.list <- PrepSCTIntegration(object.list = rsi.list, anchor.features = rsi.features, verbose = FALSE)
rsi.anchors <- FindIntegrationAnchors(object.list = rsi.list, normalization.method = "SCT", 
                                      anchor.features = rsi.features, verbose = FALSE)
rsi.combined <- IntegrateData(anchorset = rsi.anchors, normalization.method = "SCT", verbose = FALSE)

#saveRDS(rsi.combined, "rsi_3-split.rds")

aaa <- stringr::str_split(colnames(rsi.combined),'[.]', simplify = T)[,2]
head(aaa)
bbb<-stringr::str_split(aaa,'_', simplify = T)[,1]
head(bbb)
names(bbb) <- colnames(rsi.combined) 
table(bbb)
rsi.combined <- AddMetaData(object = rsi.combined, metadata = bbb, col.name = 'limb')

aaa <- stringr::str_split(colnames(rsi.combined),'[.]', simplify = T)[,1]
head(aaa)
table(aaa)
bbb<-stringr::str_split_fixed(aaa,'_', 2)[,-1]
names(bbb) <- colnames(rsi.combined) 
table(bbb)
rsi.combined <- AddMetaData(object = rsi.combined, metadata = bbb, col.name = 'stage')

aaa <- stringr::str_split(colnames(rsi.combined),'CS',simplify = T)[,2]
head(aaa)
table(aaa)
bbb<-stringr::str_split(aaa,'_', 2,simplify = T)[,1]
names(bbb) <- colnames(rsi.combined) 
table(bbb)
rsi.combined <- AddMetaData(object = rsi.combined, metadata = bbb, col.name = 'stage.2')

Idents(rsi.combined)<-"stage"
rsi.no15<-subset(rsi.combined, idents = c("CS16","CS18","CS20"))

rsi.no15 <- RunPCA(object = rsi.no15, verbose = FALSE)
rsi.no15 <- RunUMAP(object = rsi.no15, dims = 1:20, verbose = FALSE)
rsi.no15 <- FindNeighbors(object = rsi.no15, dims = 1:20, verbose = FALSE)
rsi.no15 <- RunTSNE(rsi.no15, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)

DefaultAssay(rsi.no15)<-"integrated"
rsi.no15 <- FindClusters(object = rsi.no15, verbose = FALSE, resolution = 0.05)
DimPlot(rsi.no15)

cell.type <- c("MP", "chondrocyte","epithelia","muscle","endothelia","neuron","WBC")
names(cell.type) <- levels(rsi.no15)
rsi.no15 <- RenameIdents(rsi.no15, cell.type)
#DimPlot(f, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
cell.type<-rsi.no15@active.ident
names(cell.type) <- colnames(rsi.no15)
rsi.no15 <- AddMetaData(object = rsi.no15, metadata = cell.type, col.name = 'cell.type')
DimPlot(rsi.no15)

sub_cluster <- c("MP_2","MP_10","MP_13", "MP_8","MP_9","MP_3","MP_1","epithelia", 
                 "MP_12", "chondrocyte",
"neuron","muscle", "MP_6","MP_5","MP_4","endothelia","MP_7","MP_0","chondrocyte", 
"chondrocyte","MP_11","WBC","chondrocyte")
names(sub_cluster) <- levels(x)
x <- RenameIdents(x, sub_cluster)
#DimPlot(f, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
sub_cluster<-x@active.ident
names(sub_cluster) <- colnames(x)
x <- AddMetaData(object = x, metadata = sub_cluster, col.name = 'sub_cluster')
DimPlot(x)

DefaultAssay(rsi.no15)<-"RNA"
rsi.no15 <- NormalizeData(rsi.no15)
rsi.no15 <- FindVariableFeatures(rsi.no15, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rsi.no15)
rsi.no15 <- ScaleData(rsi.no15, features = all.genes)

markers.to.plot <- c("EBF3", "COL3A1", "DKK2", "ATRNL1", "COL12A1", "GPC5","COL2A1", "COL11A1", "COL9A1", 
                     "CXCL14", "DSP", "COL17A1", "G2E3", "DLGAP5", "NDC80","LDB2", "ARHGAP26", "EMCN",
                     "SNTB1", "ADAMTS20", "CDH6","TTN", "NEB", "ACTC1",      
                     "DAB2", "RBPJ", "MRC1")
DotPlot(rsi.no15, features = markers.to.plot, dot.scale = 8,cols = "RdYlBu") + RotatedAxis()

markers <- FindAllMarkers(h, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers, "rsi_3-split-h-markers.txt")
top10 <- markers %>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC)
DoHeatmap(rsi.no15, features = top10$gene) + NoLegend()
DotPlot(rsi.no15, features = unique(top10$gene))+  RotatedAxis()

saveRDS(x, "rsi_3-split-250.rds")

#all
markers.to.plot <- c("COL3A1", "DKK2", "ATRNL1", "COL2A1", "COL11A1", "COL9A1", 
                     "CXCL14", "DSP", "GRHL2", "TTN", "NEB", "ACTC1", "LDB2", 
                     "FLT1", "EMCN", "SNTB1", "ADAMTS20", "CDH6", "DAB2", "RBPJ", "MRC1")
DotPlot(x, features = markers.to.plot, dot.scale = 8,cols = "RdYlBu") + RotatedAxis()


levels(x)<-c("MP_0", "MP_1", "MP_2", "MP_3", "MP_4", "MP_5",  "MP_6", "MP_7", 
             "MP_8", "MP_9", "MP_10",  "MP_11", "MP_12", "MP_13", "chondrocyte",
             "epithelia","muscle", "endothelia","neuron", "WBC")
#MP
markers.to.plot <- c("PNISR", "AKAP9", "USP47",
                     "GNAS", "ZNF804A", "SASH1",
                     "EBF2",  "EPHA4", "HPSE2",
                     "ARHGAP24", "PLXDC2", "ISM1",
                     "CPA6","BMPER","RNF220",
                     "SEMA3E","COL12A1","PDGFD",
                     "DIAPH3","NDC80","CENPE",
                     "HS3ST3A1","SEMA3D","PRDM16",
                     "MEIS2","HMGA2","FREM1",
                     "KCND2","NLGN1","AC097634.4",
                     "COL1A1","COL1A2","ABI3BP",
                     "RUNX2","SNED1","BMP5",
                     "HTR1F","ZFHX3","PHLDB2",
                     "DCC","TOX","SPATS2L",
                     "COL2A1", "COL11A1", "COL9A1", 
                     "CXCL14", "DSP", "GRHL2", "TTN", "NEB", "ACTC1", "LDB2", 
                     "FLT1", "EMCN", "SNTB1", "ADAMTS20", "CDH6", "DAB2", "RBPJ", "MRC1")
DotPlot(x, features = markers.to.plot, dot.scale = 8,cols = "RdYlBu") + RotatedAxis()

##Acording to correlation results combining MPs

Idents(x)<-x$cell.type3
c<-subset(x, idents = "MP")
DefaultAssay(c)<-"integrated"
c <- RunPCA(object = c, verbose = FALSE)
c <- FindNeighbors(object = c, dims = 1:30, verbose = FALSE)
c <- RunUMAP(object = c, dims = 1:30, verbose = FALSE)
#c <- RunTSNE(c, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
c <- FindClusters(object = c, verbose = FALSE, resolution = 0.5)

DefaultAssay(c)<-"RNA"
c <- NormalizeData(c)
c <- FindVariableFeatures(c, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(c)
c <- ScaleData(c, features = all.genes)
markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%group_by(cluster) %>%top_n(n = 5, wt = avg_log2FC)
DoHeatmap(x, features = top5$gene) + NoLegend()
DotPlot(x, features = unique(top5$gene),dot.scale = 8,cols = "RdYlBu")+  RotatedAxis()

cell.type2 <- c("cop","osteoblast","MP_12", "MP_8","cop","MP_3","MP_0","epithelia", 
                "MP_12", "chondrocyte",
                "neuron","muscle", "pe","MP_5","cop","endothelia","cp",
                "MP_0","op","WBC")
names(cell.type2) <- levels(x)
x <- RenameIdents(x, cell.type2)
#DimPlot(f, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
cell.type2<-x@active.ident
names(cell.type2) <- colnames(x)
x <- AddMetaData(object = x, metadata = cell.type2, col.name = 'cell.type3')
DimPlot(x)

levels(x)<-c("MP_0", "cop","MP_3", "MP_5",  "MP_8", "MP_12", "cp","chondrocyte",
             "op",  "osteoblast", "pe", "epithelia", "muscle", "endothelia", "neuron", "WBC" )

x@reductions$umap1<-x@reductions$umap
x <- RunUMAP(object = x, dims = 1:30, verbose = FALSE)
x@reductions$umap2<-x@reductions$umap
x@reductions$tsne1<-x@reductions$tsne
x <- RunTSNE(x, dims = 1:30, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)

DimPlot(x)
#figure size:6*7

x<-RenameIdents(x, 'cop'='EMP','osteoblast'='Osteoblast',  "CTP_12"='ZMP', "CTP_8"="MMP",
                "CTP_3"="AMP", "CTP_0"='PNMP', 'epithelia'= "Epithelia", "chondrocyte"='Chondrocyte',
                "neuron"='Neuron', "muscle"='Muscle', "pe"='PE', "CTP_5"= 'PDMP', 
                "endothelia"='Endothelia', "cp"='CP', "op"='OP', "WBC"='BM' )
levels(x)<-c("PNMP", "EMP", "AMP", "PDMP", "MMP", "ZMP", "CP", "Chondrocyte", "OP", "Osteoblast",
             "PE", "Epithelia", "Muscle", "Endothelia", "Neuron", "BM")

write.csv(FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, 
                         logfc.threshold = 0.25),"rsi_3-markers-combined-2.csv")

levels(x)<-c("BM", "Neuron", "Endothelia","Muscle","Epithelia","PE","Osteoblast", "OP", 
              "Chondrocyte",  "CP","ZMP", "MMP",  
             "PDMP", "AMP", "EMP", "PNMP"  )
markers.to.plot <- c("PNISR", "AKAP9", "USP47",
  "EBF2","HPSE2","EPHA4",
  "ARHGAP24", "PLXDC2", "ISM1",
  "SEMA3E","COL12A1","PDGFD",
  "MEIS2","HMGA2","FREM1",
  "HTR1F","ZFHX3","DCC",
  "SOX5","TRPS1","FOXP2",
  "COL2A1", "COL11A1", "COL9A1",
  "RUNX2","SNED1","CPED1",
  "COL1A1","COL1A2","ABI3BP",
  "DIAPH3","NDC80","CENPE",
  "CXCL14", "DSP", "GRHL2", 
  "TTN", "NEB", "ACTC1", 
  "LDB2", "FLT1", "EMCN", 
  "SNTB1", "ADAMTS20", "CDH6", 
  "DAB2", "RBPJ", "MRC1")
DotPlot(x, features = markers.to.plot, dot.scale = 8,cols = "RdYlBu") + RotatedAxis()

DotPlot(x, features = markers.to.plot) + 
  scale_colour_gradient2(low = "navy", high = "darkred", na.value = NA)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



