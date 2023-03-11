
#### PACKAGES
if(T){
  library(ggplot2)
  library(reshape2)
  library(DESeq2)
  library(FactoMineR)
  library(factoextra)
  library(Rtsne)
  library(ggpubr)
  library(pheatmap)
  library(tidyverse)
  library(ggrepel)
}


#####  DATA
if(T){
  rm(list = ls())
  setwd("path") 
  load("data")
  GR <- c(G,R)
  
  cellID_STAR <- cellID_STAR[cellID_STAR$number %in% GR,]
  TPM_STAR <- TPM_STAR[,colnames(TPM_STAR) %in% GR]
  geneid
}

#### qc
if(T){
  QC_data <- TPM_STAR
 
  QC_data=QC_data[apply(QC_data,1, function(x) sum(x>1) > 5),]
  head(QC_data)
  dim(QC_data)
  cellID_STAR$GR <- factor(cellID_STAR$GR)
  boxplot(log2(QC_data+1),outline=FALSE, notch=T,col=as.factor(cellID_STAR$GR),las=2,labels= cellID_STAR$GR)+
    legend(list(x= 0,y = 13),inset = 0.1,c("G","R"),pch = c(15,15),col = c("black","red"))
  cell_sample = apply(QC_data,2,function(x) sum(x>1))
  hist(cell_sample,nclass = 30,labels = T,xlab = "Cell gene expression total number")
  TPM_qc <- QC_data
  TPM_qc <- cbind(TPM_qc[,colnames(TPM_qc) %in% G],TPM_qc[,colnames(TPM_qc) %in% R])
  cellID_STAR <- rbind(cellID_STAR[cellID_STAR$number %in% G,],cellID_STAR[cellID_STAR$number %in% R,])
  
  
  
}


####  PCA 
if(T){
    #
    dat <- TPM_qc 
    dat=as.data.frame(t(dat)) 
    dat=cbind(dat,cellID_STAR) #
   
    library("FactoMineR")
    library("factoextra") 
   
    dat.pca <- PCA(dat[,-(12976:12980)], graph = FALSE)
    P <- fviz_pca_ind(dat.pca,repel =T,
                      geom.ind = c("point","text"), 
                      col.ind = cellID_STAR$GR, 
                      addEllipses = T, 
                      pointsize = 3,
                      legend.title = "Groups",
                      title = "PCA - Axin2 single cell")
    P
    P+theme_bw()+ 
      theme(panel.grid.major = element_blank(),   
            panel.grid.minor = element_blank(),
            panel.border = element_blank())
    P
    

}





#### Hierarchical clustering
if(T){
  hc_data <- TPM_qc
  hc=hclust(dist(t(hc_data))) 
  plot(hc,labels =colnames(hc_data))
  if(F){
    clus = cutree(hc, 2) 
    table(clus) 
    cellID$hc_cluster = as.factor(clus) #
  }
}

##### DEG
if(T){
  library(DESeq2)
  coldata <- as.data.frame(cellID_STAR$GR)
  rownames(coldata) <- rownames(cellID_STAR)
  colnames(coldata) <- "celltype"
  coldata$celltype <- as.factor(coldata$celltype)
  condition <- as.factor(coldata$celltype)
  dds <- DESeqDataSetFromMatrix(countData = ceiling(TPM_qc), colData = DataFrame(condition), design = ~condition) 
  
  dds1 <- DESeq(dds, minReplicatesForReplace = 7, parallel = FALSE) 
  res <- results(dds1, contrast = c('condition', 'G', 'R')) 
  res<-res[order(res$padj),] 
  sum(res$padj < 0.05, na.rm=TRUE) 
  summary(res)
  res1 <- as.data.frame(res)
  res1$Ensemble <- rownames(res1)
  res1 <- merge(res1,geneid,by = "Ensemble")
 
  plotMA(res, ylim=c(-10,10),main="TC_NG", alpha = 0.05)
 
  find.significant.genes <- function(de.result,alpha=0.05) {
   
    filtered <- de.result[(de.result$padj < alpha) & !is.infinite(de.result$log2FoldChange) 
                          & !is.nan(de.result$log2FoldChange) & !is.na(de.result$padj),]
   
    sorted <- filtered[order(filtered$padj),c(1,2,6)]
  }
  de2.genes <- find.significant.genes(res) 
  
 
  DEG <- res1
  DEG$Significant <- ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) >= 0.5, 
                            ifelse(DEG$log2FoldChange > 0.5, "Up", "Down"), "Stable")
  
  DEG$Significant <- as.factor(DEG$Significant)
  
  DEG <- na.omit(DEG)
  
  
}


if(T){
  res_table_thres <- as.data.frame(res1)
  res_table_thres <- na.omit(res_table_thres)
 
  res_table_thres$Significant <- ifelse(res_table_thres$pvalue < 0.05 & abs(res_table_thres$log2FoldChange) >= 0.4, 
                                        ifelse(res_table_thres$log2FoldChange > 0.4, "Up", "Down"), "Stable")
  
  res_table_thres$Significant <- as.factor(res_table_thres$Significant)

  p <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = Significant )) +
    scale_color_manual(values = c("#619CFF","#00BA38","#F8766D"))+
    ggtitle("Volcano plot of G relative to R") +
    xlab("log2 fold change") + 
    ylab("-log10 p-value") +
    scale_y_continuous(limits = c(0,8)) + 
    scale_x_continuous(limits = c(-3,3)) +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25)))+
    theme_classic()+
    geom_text_repel(data = subset(res_table_thres, (padj < 0.05 & abs(log2FoldChange) >= 0.8)),
                    aes(log2FoldChange, -log10(pvalue), label = Symble),
                    size = 4, box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE)
  #pdf("./volcano_plot2.pdf")
  print(p)
  # dev.off()
  
}




















