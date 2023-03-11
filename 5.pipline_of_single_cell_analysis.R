

if(T){
  library(Seurat)
  library(ChIPseeker)
  library(clusterProfiler)
  library(GenomicRanges)
  library(GenomicRanges)
  library(tidyverse)
  library(reshape2)
  library(stringr)
  library(clusterProfiler)
  library(ReactomePA)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  library(org.Mm.eg.db)
  
}


path <- ""
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   


data <- lapply(filePath, function(x){
  data <- Read10X(x,gene.column = 1)
  data <- CreateSeuratObject(data,project = "Axin2", min.cells = 3, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  data <- NormalizeData(data,verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data, features = rownames(data), verbose = FALSE)
  data <- RunPCA(data, features = VariableFeatures(object =data),dims = 1:30)
  data <- FindNeighbors(data, dims = 1:30)
  data <- FindClusters(data, resolution = 0.5,random.seed = 12580)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30,n.neighbors = 10, n.components = 2, min.dist=0.3,seed.use = 12580)
  return(data)
})  


if(F){
  data <- lapply(data, function(x){
    x <- subset(x,subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 7 & nCount_RNA < 307000)
  })
}



data<- lapply(X = data, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list =data)
data <- lapply(X = data, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


single_cell_BGI_pipline <- function(data){
  data <- NormalizeData(data,verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data, features = rownames(data), verbose = FALSE)
  data <- RunPCA(data, features = VariableFeatures(object =data),dims = 1:30)
  data <- FindNeighbors(data, dims = 1:30)
  data <- FindClusters(data, resolution = 0.5,random.seed = 12580)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30,n.neighbors = 10, n.components = 2, min.dist=0.3,seed.use = 12580)
  return(data)
}

single_cell_integrate_pipline <- function(data){
  data <- ScaleData(data, features = rownames(data), verbose = FALSE)
  data <- RunPCA(data, features = features,dims = 1:30)
  data <- FindNeighbors(data, dims = 1:30)
  data <- FindClusters(data, resolution = 0.5,random.seed = 12580)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30,n.neighbors = 10, n.components = 2, min.dist=0.3,seed.use = 12580)
  return(data)
}


mouse.anchors <- FindIntegrationAnchors(object.list = data, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
Axin2<- IntegrateData(anchorset = mouse.anchors)
DefaultAssay(Axin2) <- "integrated"

Axin2 <- single_cell_integrate_pipline(Axin2)
saveRDS(Axin2,file = "./Axin2_integrate.rds")

DefaultAssay(Axin2) <- "integrated"
Axin2 <- FindClusters(Axin2, resolution = 0.1,random.seed = 12580)
DimPlot(Axin2,label = T) 
VlnPlot(Axin2,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0,ncol = 1)


######################################### Subset Sox2 NPC cell for analysis
Idents(Axin2) <- Axin2$integrated_snn_res.0.1
DimPlot(Axin2,label = T)
NPC <- subset(Axin2, idents= c("1","2","3"))
saveRDS(NPC,file = "./NPC.rds")


NPC <- single_cell_BGI_pipline(NPC)
DimPlot(NPC_wt,label = T)
NPC <- FindClusters(NPC_wt, resolution = 0.1,random.seed = 12580)
FeaturePlot(NPC,features = c("Axin2","Pax6","Sox2","Emx1"))
DimPlot(NPC,label = T)

NPC_clean <- subset(NPC_clea,subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 7 & nCount_RNA < 30000)
NPC_clean <- single_cell_BGI_pipline(NPC_clean)
NPC_clean <- FindClusters(NPC_clean, resolution = 0.1,random.seed = 12580)
Idents(NPC_clea) <- NPC_clea$RNA_snn_res.0.1
DimPlot(NPC_clean,label=T)
FeaturePlot(NPC_clean,features = c("Axin2","Pax6","Sox2","Emx1"))
VlnPlot(NPC_clean,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0)

saveRDS(NPC_clean,file = "./NPC_clean.rds")




