phyloseqTransformation <- function(ps, 
                                   c = .4){
  dds <- phyloseq_to_deseq2(ps, design = ~ 1)
  dds <-  estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds, fitType = "local")
  abund <- assay(dds)
  abund_temp <- matrix(nrow = nrow(abund), ncol = ncol(abund))
  
  # Column co-factor
  
  for(i in 1:dim(abund)[1]){
    for(j in 1:dim(abund)[2]){
      abund_temp[i,j] <- abund[i,j]/sizeFactors(dds)[j]
    }
  }
  
  
  # Row co-factor
  
  # k > 2
  ind <- which(dispersions(dds) > 2)
  # K <= 2
  ind2 <- which(dispersions(dds) <= 2)
  
  for(i in ind){
    for(j in 1:dim(abund_temp)[2]){
      if(abund_temp[i,j] > 0){
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt((abund_temp[i, j]+ c)/(dispersions(dds)[i]-2*c)))
      }else{
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt(abund_temp[i, j]))
      }
      
    }
  }
  
  
  if(length(ind2) > 0){
    for(i in ind2){
      for(j in 1:dim(abund_temp)[2]){
        if(abund_temp[i,j] > 0){
          abund_temp[i,j] <- log2(abund_temp[i, j] + .5*dispersions(dds)[i])
        }
        
      }
    }
  }
  
  
  rownames(abund_temp) <- rownames(abund)
  colnames(abund_temp) <- colnames(abund)
  
  
  ps <- phyloseq(otu_table(abund_temp, taxa_are_rows = TRUE),
                 sample_data(ps), 
                 tax_table(ps), 
                 phy_tree(ps))
  
  return(ps)
  
}






# phyloseqTransformation <- function(ps, 
#                                   c = .4){
#   dds <- phyloseq_to_deseq2(ps, design = ~ 1)
#   dds = estimateSizeFactors(dds, type = "poscounts")
#   dds = estimateDispersions(dds, fitType = "local")
#   abund = assay(dds)
#   abund_temp = sweep(abund, MARGIN = 2, sizeFactors(dds), "/")
#   
#   abund[abund > 0] = abund + c
#   den = dispersions(dds)-2*c
#   abund_temp = sweep(abund_temp, MARGIN = 1, den, "/")
#   # k > 2
#   ind = which(dispersions(dds) > 2)
#   abund_temp[ind, ] = sqrt(dispersions(dds)[ind]- .5)*asinh(sqrt(abund_temp[ind, ]))
#   #abund_temp[ind, ] = asinh(sqrt(abund_temp[ind, ]))
#   # K <= 2
#   ind2 = which(dispersions(dds) <= 2)
#   if(length(ind2) > 0){
#     abund_temp[ind2, ][abund_temp[ind2, ] > 0] = log(abund_temp[ind2, ] + .5*dispersions(dds)[ind2])
#   }
#   
#   ps <- phyloseq(otu_table(abund_temp, taxa_are_rows = TRUE),
#                  sample_data(ps), 
#                  tax_table(ps), 
#                  phy_tree(ps))
#   
#   return(ps)
#   
# }
# 
# 
