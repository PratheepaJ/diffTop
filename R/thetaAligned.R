thetaAligned <- function(theta, 
                         K, 
                         aligned, 
                         iter = 2000, 
                         chain = 4){
  # align topics between chains
  # switch samples and Topic dimension in array
  theta <- aperm(theta, c(1,3,2))
  
  theta_chain <- list()
  theta_chain[[1]] <- theta[1:(iter/2), ,]
  theta_chain[[1]] <- theta_chain[[1]][, aligned[,1],]
  
  for(ch in 2:chain){
    theta_chain[[ch]] <- theta[((ch-1)*(iter/2)+1):(ch*(iter/2)),,]
    theta_chain[[ch]] <- theta_chain[[ch]][,aligned[,ch],]
  }
  
  theta_aligned <- abind(theta_chain[[1]], theta_chain[[2]], along = 1)
  for(ch in 3:chain){
    theta_aligned <- abind(theta_aligned, theta_chain[[ch]], along = 1)
  }
  # switch back samples and Topic dimension in array
  theta_aligned <- aperm(theta_aligned, c(1,3,2))
  
  
  return(theta_aligned)
}