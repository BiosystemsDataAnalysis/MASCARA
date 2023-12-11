MASCARA4_test_1 <- function(resids, baits, spikes){
  #MIXOMICS
  spls_res <- spls(resids[,-which(colnames(resids) %in% baits)], 
                   resids[,which(colnames(resids) %in% baits)], all.outputs = TRUE)
  
  #####
  y <- colMeans(spls_res$loadings$Y)
  U <- spls_res$loadings$X
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = TRUE)])
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_bait_cands,UTP, spls_res))
}



MASCARA4_test_2 <- function(resids, baits, spikes){
  #target projection
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                          resids[,which(colnames(resids) %in% baits)], ncomp = 2)
  
  #####
  y <- colMeans(spls_res$Yloadings)
  U <- spls_res$projection
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = TRUE)])
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_bait_cands,UTP,spls_res,y))
}




MASCARA4_test_8 <- function(resids, baits, spikes){
  #vip
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                         resids[,which(colnames(resids) %in% baits)], ncomp = 2)
  
  
  #####
  y <- colMeans(spls_res$Yloadings)
  U <- spls_res$projection
  
  #CHECK spls_res$residuals for variance explained 
  
  
  x <- Yloadings(spls_res)
  ve_y <- colSums(x^2)/nrow(x)
  


  sPLS_cands <- cbind.data.frame(abs(spls_res$projection) %*% ve_y, spls_res$projection)
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]

  RP <- prod(which(rownames(sPLS_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_cands))
}





MASCARA4_test_3 <- function(resids, baits, spikes){
  #Target projection PLS1 model on u1 of bait resids
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                         prcomp(resids[,which(colnames(resids) %in% baits)])$x[,1], ncomp = 2)
  
  #####
  y <- spls_res$Yloadings[]
  U <- spls_res$projection
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = TRUE)])
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_bait_cands,UTP))
}


#selectivity ratio (SR)
MASCARA4_test_7 <- function(resids, baits, spikes, ncomp = 2){
  
  X <- resids[,-which(colnames(resids) %in% baits)]
  Y <- resids[,which(colnames(resids) %in% baits)]
  spls_res <- simpls.fit(X, Y, ncomp = ncomp)
  
  #####
  
  B <- spls_res$projection %*% t(spls_res$Yscores)
  y_hat <- X %*% B
  
  Y_u <- prcomp(y_hat)$x[,1:ncomp]   
  
  X_hat <- (Y_u %*% solve(t(Y_u) %*% Y_u) %*% t(Y_u)) %*% X
  
  X_res <- X - X_hat
  
  SR <- colSums(X_hat[]^2)/colSums(X_res[]^2)
  
  SR <- cbind.data.frame(colnames(X),"Rank" = SR)
  SR <- SR[order(SR$Rank, decreasing = TRUE),]
  
  sPLS_bait_cands <- SR
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_bait_cands))
}


