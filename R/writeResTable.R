writeResTable <- function(res, fileN = paste0("../results/dds_all_paired_aster_", K,".tex")){
  df <- data.frame(Topic = rownames(res), 
                   lfc = res$log2FoldChange, 
                   lfcSE= res$lfcSE, 
                   WTS = res$stat, 
                   pvalue = res$pvalue, 
                   p.adj = res$padj)
  df$lfc <- round(df$lfc, digits = 2)
  df$lfcSE <- round(df$lfcSE, digits = 2)
  df$WTS <- round(df$WTS, digits = 2)
  df$p.adj <- round(df$p.adj, digits = 4)
  df$p.adj[which(df$p.adj == 0)] <- "<.0001"
  print(xtable(df, 
               type = "latex", 
               digits = c(0,0,2, 2,2,2,4)), 
        file = fileN)
}