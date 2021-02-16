betaAligned <- function(beta, 
                         K, 
                         aligned, 
                         iter = 2000, 
                         chain = 4){
  # align topics between chains
  # no need to switch array dimnension
  beta_chain <- list()
  beta_chain[[1]] <- beta[1:(iter/2), ,]
  beta_chain[[1]] <- beta_chain[[1]][, aligned[,1], ]
  
  for(ch in 2:chain){
    beta_chain[[ch]] <- beta[((ch-1)*(iter/2)+1):(ch*(iter/2)),,]
    beta_chain[[ch]] <- beta_chain[[ch]][,aligned[,ch],]
  }
  
  beta_aligned <- abind(beta_chain[[1]], beta_chain[[2]], along = 1)
  for(ch in 3:chain){
    beta_aligned <- abind(beta_aligned, beta_chain[[ch]], along = 1)
  }
  
  return(beta_aligned)
}