#########  enrichment ##############
{
  if(T){
    library(clusterProfiler)
    library(DOSE)
    library(org.Mm.eg.db)
    library(ggplot2)
    library(enrichplot)
    library(stringr)
    library(topGO)
    library(ggupset)
    library("pathview")
    library(enrichMap)
  }
}


### id transfer
{
  funtion_analysis <- as.data.frame(res1)

  gene.df <- bitr(funtion_analysis$Ensemble, fromType = "ENSEMBL" , 
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  
  funtion_analysis <- merge(funtion_analysis,gene.df,by.x="Ensemble" ,by.y="ENSEMBL") 
  RG <- funtion_analysis[(funtion_analysis$pvalue < 0.05 & abs(funtion_analysis$log2FoldChange) > 0.5),]

    R_up <- funtion_analysis[funtion_analysis$pvalue < 0.05 & funtion_analysis$log2FoldChange <0,] 
    G_up <- funtion_analysis[funtion_analysis$pvalue < 0.05 & funtion_analysis$log2FoldChange >0,] 
  }
  
} 


{
  
  ego_R<-enrichGO(gene       = R_up$Ensemble,
                  OrgDb      = org.Mm.eg.db,
                  keyType    = 'ENSEMBL',
                  ont        = "ALL", #
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)
  
  ego_G<-enrichGO(gene       = G_up$Ensemble,
                  OrgDb      = org.Mm.eg.db,
                  keyType    = 'ENSEMBL',
                  ont        = "ALL", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)
  
  ego_GR<-enrichGO(gene       = RG$Ensemble,
                   OrgDb      = org.Mm.eg.db,
                   keyType    = 'ENSEMBL',
                   ont        = "ALL", 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 0.5)
  
  
  dotplot(ego_G, showCategory=30)
  barplot(ego_G,showCategory = 30,title="The GO enrichment analysis of G up ")+ 
    scale_size(range=c(2, 12))+
    scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 70))
  
  barplot(ego_R,showCategory = 30,title="The GO enrichment analysis of R up ")+ 
    scale_size(range=c(2, 12))+
    scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 70))
  
  # Gene-Concept Network 
  edox <- setReadable(ego_G, 'org.Mm.eg.db', 'ENSEMBL') ## convert gene ID to Symbol
  cnetplot(edox, foldChange=G_up$log2FoldChange,showCategory = 10)
  cnetplot(edox, foldChange=G_up$log2FoldChange,circular = TRUE,colorEdge = TRUE)## Gene-Concept Network
  emapplot(ego_G)
  
  edox2 <- setReadable(ego_R, 'org.Mm.eg.db', 'ENSEMBL') ## convert gene ID to Symbol
  cnetplot(edox2, foldChange=R_up$log2FoldChange,showCategory = 10)
  cnetplot(edox2, foldChange=R_up$log2FoldChange,circular = TRUE,colorEdge = TRUE)## Gene-Concept Network
  emapplot(ego_R)
  
  ## kegg 分析
  kk_G <- enrichKEGG(gene     = G_up$ENTREZID, 
                     organism     = 'mmu',
                     pAdjustMethod = 'BH', 
                     #universe     = gene_all,
                     pvalueCutoff = 0.5,
                     qvalueCutoff =0.5)
  kk_R <- enrichKEGG(gene     = R_up$ENTREZID, 
                     organism     = 'mmu',
                     pAdjustMethod = 'BH', 
                     #universe     = gene_all,
                     pvalueCutoff = 0.5,
                     qvalueCutoff =0.5)
  kk_GR <- enrichKEGG(gene     = RG$ENTREZID, 
                      organism     = 'mmu',
                      pAdjustMethod = 'BH', 
                      #universe     = gene_all,
                      pvalueCutoff = 0.5,
                      qvalueCutoff =0.5)
  
  
  barplot(kk_G, showCategory=30)
  dotplot(kk.ego,showCategory = 25, title="The KEGG enrichment analysis of Male Neurons")+
    scale_size(range=c(2, 12))+
    scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))
  barplot(kk_R,showCategory = 30, title="The KEGG enrichment analysis of R ")+
    scale_size(range=c(2, 12))+
    scale_x_discrete(labels=function(kk) str_wrap(kk,width = 50))
  barplot(kk_G,showCategory = 30, title="The KEGG enrichment analysis of G ")+
    scale_size(range=c(2, 12))+
    scale_x_discrete(labels=function(kk) str_wrap(kk,width = 50))
  
  browseKEGG(kk_G, 'mmu04310')
  
  
}
{ 
  genelist <- G_up$log2FoldChange
  names(genelist) <- G_up$ENTREZID
  genelist2 <- R_up$log2FoldChange
  names(genelist2) <- R_up$ENTREZID
  genelist <- sort(genelist, decreasing = TRUE)
  genelist2 <- sort(genelist2, decreasing = TRUE)
  gsemf <- gseGO(genelist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "ENTREZID",
                 ont="MF",
                 pvalueCutoff = 0.05)
  gsemf2 <- gseGO(genelist2,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont="CC",
                  pvalueCutoff = 0.5) 
}
head(summary(gsecc))



