library(devtools)
library(ggplot2)
library(MetStaT)
library(gASCA)
library(patchwork)
library(readxl)
library(ggfortify)
library(ggrepel)
library(grid)
library(gridExtra)
library(tidyr)
library(tidyverse)
library(reshape2)
library(scales)
library(matrixStats)
library(DESeq2)
library(MASS)
library(glmnet)
library(foreach)
library(doParallel)    
library(MUVR) 
library(pheatmap)
library(caret)
library(mixOmics)
library(pracma)
library(WGCNA)
library(parallel)


# source("../DATA_SIM_FUNCS.R")

source("/zfs/omics/personal/fwhite/MASCARA/DATA_SIM_FUNCS.R")


# Create_Core_DEV <- Create_Core_DEV_2

F1_plot <- F1_plot_no_plot
# nreps <- 3
alphaN <- 2
betaN <- 4
ncands <- 500


splsR <- function(resids, baits, spikes){
  #vip
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                         resids[,which(colnames(resids) %in% baits)], ncomp = 2)
  
  
  #####
  y <- colMeans(spls_res$Yloadings)
  U <- spls_res$projection
  
  x <- Yloadings(spls_res)
  ve_y <- colSums(x^2)/nrow(x)
  
  sPLS_cands <- cbind.data.frame(abs(spls_res$projection) %*% ve_y, spls_res$projection)
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
  
  RP <- prod(which(rownames(sPLS_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP, sPLS_cands))
}




MASCARA4_test <- function(resids, baits, spikes){
  #target projection
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                         resids[,which(colnames(resids) %in% baits)], ncomp = 2)
  
  #####
  y <- colMeans(spls_res$Yloadings)
  U <- spls_res$projection
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = TRUE)])
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  
  return(list(RP, sPLS_bait_cands,UTP))
}

WGCNA4_test <- function(data, baits, spikes){
  
  softPower <- 6
  adjacency <- adjacency(data, power = softPower)
  
  TOM <- TOMsimilarity(adjacency)
  TOM.dissimilarity <- 1-TOM
  geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
  Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 12)
  ModuleColors <- labels2colors(Modules)
  MElist <- moduleEigengenes(data, colors = ModuleColors)
  MEs <- MElist$eigengenes
  ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete")
  
  merge <- mergeCloseModules(data, ModuleColors, cutHeight = .25)
  # mergedColors = merge$colors
  
  mergedMEs = merge$newMEs
  
  
  sKME <- signedKME(data, mergedMEs)
  
  #feature selection
  q_bar <- colMeans(sKME[baits,])
  R <- sKME[!rownames(sKME) %in% baits,]
  R_TP <- (dot(t(R),q_bar)/dot(q_bar,q_bar))
  
  Candidates <- as.data.frame(R_TP[order(R_TP, decreasing = TRUE)])
  
  RP <- prod(which(rownames(Candidates) %in% spikes))^(1/length(spikes))
  
  return(list(RP, Candidates, mergedMEs, sKME))
  
  
}


baits <- paste0("X_",c(1985:1988))   #501:510
spikes <- paste0("X_",c(1989:2000))


Baits <- c(baits,spikes)


nreps <- 3
meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = nreps))
colnames(meta) <- c("growth_condition","time")
meta$ID <- paste(meta$growth_condition, meta$time, sep = "_")



sim_data <- Create_Core_DEV(nreps = 3, meta = meta, plot = F, 
                            EffectSize = c(X_a_ab = 1, time = 1, E = 0.5, mu = 1, Struc_E = 1.5), 
                            struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                            b_sigma = c(1,0.9), irr_spikes = TRUE, SCORE_SEED = 1027)







nreps <- 15

X_30 <- list()

big_meta <- cbind.data.frame(rep(c(1,-1), each = betaN*nreps), rep(c(1:betaN), each = nreps))

colnames(big_meta) <- c("growth_condition","time")
big_meta$ID <- paste(big_meta$growth_condition, big_meta$time, sep = "_")


big_meta$cond <- paste0(big_meta[,1],"_",big_meta[,2])
conds <- unique(big_meta$cond)
big_meta$reps <- paste0(big_meta$cond,"_",c(1:nreps))


#sample the X_30 dataset here with reps from x
X_funced <- list()
meta <- list()

norms <- list()

l <- 1


# y <- rep(c(0.0004,0.002, 0.01),3)      
# x <- rep(30,9)  #100
# nse <- rep(c(0.1,0.3,0.5), each = 3)
# 

# c(0.25,0.5,0.75)/30

# y <- rep(c(0.008333333, 0.016666667, 0.025000000),3)      
# x <- rep(30,9)  #100
# nse <- rep(c(0.2,0.4,0.6), each = 3)
# 
# 
# 
# main_lev <- x*y
# ERs <- x
# 
# #######
# 
# idlist <- list()
# 
# a <- NULL
# 
# # set.seed(40001)
# 
# for(a in 1:length(ERs)){
#   
#   i <- NULL
#   
#   for(i in 1:20){
#     
#     
#     X_30[[i]] <- Create_Core_DEV(nreps = nreps, meta = big_meta, plot = FALSE, 
#                     EffectSize = c(X_a_ab = main_lev[a], time = 0.2, E = nse[a], mu = 0.2, Struc_E = 0.076), 
#                     struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
#                     b_sigma = c(1,0), irr_spikes = FALSE, SCORE_SEED = (250021 + i), Experiment_responders = ERs[a], ts = 3 + i) 
#     
#     
#     
#     
#     k <- NULL
#     for(k in 3:nreps){
#       
#       IDS <- NULL
#       
#       j <- NULL
#       for(j in 1:length(conds)){
#         
#         set.seed(40001 + l*2)
#         ids <- sample(big_meta$reps[which(big_meta$cond == conds[j])], k)
#         IDS <- c(IDS,ids)
#         
#         idlist[[l]] <- IDS
#         
#         
#         
#       }
#       
#       X_funced[[l]] <- X_30[[i]][[1]][which(big_meta$reps %in% IDS),]
#       meta[[l]] <- big_meta[which(big_meta$reps %in% IDS),]
#       
#       norms[[l]] <- X_30[[i]][[6]]
#       
#       l <- l + 1
#       
#       
#       
#     }
#     
#     
#     
#   }
#   
#   
#   
#   
# }

# y <- rep(c(0.008333333, 0.016666667, 0.025000000),3)
# x <- rep(30,9)  #100
# nse <- rep(c(0.2,0.4,0.6), each = 3)


m_l <- c(1:5)/10




x <- rep(30,25)  #100  number of Experiment responders (CE) genes  URP

y <- m_l/x[1]     #rep(c(0.008333333, 0.016666667, 0.025000000))
y <- rep(y,length(y))



nse <- rep(m_l, each = length(m_l))  #rep(c(0.1,0.2, 0.3, 0.4, 0.5,0.6), each = length(y))


main_lev <- x*y
ERs <- x

library(parallel)

# Initialize global lists
X_30 <- list()
X_funced <- list()
meta <- list()
norms <- list()
idlist <- list()

# Define the number of cores
num_cores <- 16#detectCores() - 1

# Initialize the cluster
cl <- makeCluster(num_cores)

# Export necessary variables and functions to the cluster
# clusterExport(cl, varlist = c("Create_Core_DEV", "nreps", "big_meta", "conds", "main_lev", "ERs", "nse"))
# source("DATA_SIM_FUNCS.R")

# If Create_Core_DEV is defined in a separate script, source it on the cluster nodes
# clusterEvalQ(cl, )
clusterEvalQ(cl, {
  library(devtools)
  library(ggplot2)
  library(MetStaT)
  library(gASCA)
  library(patchwork)
  library(readxl)
  library(ggfortify)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(tidyr)
  library(tidyverse)
  library(reshape2)
  library(scales)
  library(matrixStats)
  library(DESeq2)
  library(MASS)
  library(glmnet)
  library(foreach)
  library(doParallel)    
  library(MUVR) 
  library(pheatmap)
  library(caret)
  library(mixOmics)
  library(pracma)
  library(WGCNA)
  source("/zfs/omics/personal/fwhite/MASCARA/DATA_SIM_FUNCS.R")
})

# Create_Core_DEV <- Create_Core_DEV_2


# Loop over 'a'
for (a in 1:length(ERs)) {
  
  # Prepare a list of indices for 'i'
  i_indices <- 1:20
  
  clusterExport(cl, varlist = c("Create_Core_DEV_2", "nreps", "big_meta", "conds", "main_lev", "ERs", "nse","a"))
  
  # Use parLapply to parallelize over 'i'
  results <- parLapply(cl, i_indices, function(i) {
    
    # Ensure necessary variables are in scope
    nreps_local <- nreps
    big_meta_local <- big_meta
    conds_local <- conds
    main_lev_a <- main_lev[a]
    nse_a <- nse[a]
    ERs_a <- ERs[a]
    
    # print(Create_Core_DEV)
    
    # Call Create_Core_DEV
    X_30_i <- Create_Core_DEV_2(
      nreps = nreps_local, meta = big_meta_local, plot = FALSE, 
      EffectSize = c(X_a_ab = main_lev_a, time = 0.2, E = nse_a, mu = 0.2, Struc_E = 0.076), 
      struc_resid = TRUE, e_sigma = c(1.5, 0.8, 0.0), a_sigma = c(1.5, 0.7, 0.6, 0.1),
      b_sigma = c(1, 0), irr_spikes = FALSE, SCORE_SEED = (250021 + i), 
      Experiment_responders = ERs_a, ts = 3 + i
    )
    
    # Initialize local lists
    X_funced_i <- list()
    meta_i <- list()
    norms_i <- list()
    idlist_i <- list()
    l_local <- 1  # Local counter
    
    # Loop over 'k'
    for (k in 3:nreps_local) {
      
      IDS <- c()
      
      # Loop over 'j'
      for (j in 1:length(conds_local)) {
        
        set.seed(40001 + a * 1000 + i * 10 + l_local * 2)
        ids <- sample(big_meta_local$reps[which(big_meta_local$cond == conds_local[j])], k)
        IDS <- c(IDS, ids)
      }
      
      idlist_i[[l_local]] <- IDS
      X_funced_i[[l_local]] <- X_30_i[[1]][which(big_meta_local$reps %in% IDS), ]
      meta_i[[l_local]] <- big_meta_local[which(big_meta_local$reps %in% IDS), ]
      norms_i[[l_local]] <- X_30_i[[6]]
      
      l_local <- l_local + 1
    }
    
    # Return results from this iteration
    list(
      X_30_i = X_30_i,
      X_funced_i = X_funced_i,
      meta_i = meta_i,
      norms_i = norms_i,
      idlist_i = idlist_i
    )
    
  })
  
  # Combine results from all 'i' iterations
  for (res in results) {
    X_30 <- c(X_30, list(res$X_30_i))
    X_funced <- c(X_funced, res$X_funced_i)
    meta <- c(meta, res$meta_i)
    norms <- c(norms, res$norms_i)
    idlist <- c(idlist, res$idlist_i)
  }
}

# Stop the cluster after computation
stopCluster(cl)



# Determine the number of cores to use
numCores <- 16
cl <- makeCluster(numCores)



clusterEvalQ(cl, {
  library(devtools)
  library(MetStaT)
  library(gASCA)
  library(reshape2)
  library(scales)
  library(matrixStats)
  library(MASS)
  library(glmnet)
  library(pracma)
  library(WGCNA)
  source("/zfs/omics/personal/fwhite/MASCARA/DATA_SIM_FUNCS.R")
})
Baits <- c(baits, spikes)

# Initialize result data frames as before
ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
COEXP_HIGH_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
PLS_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
MASCARA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
WGCNA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

ASCA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
PLS_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
WGCNA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))

EI <- c()  # To store error indices

lip <- c(1:length(X_funced))

clusterExport(cl, varlist = c("X_funced","meta","baits","spikes","ncands"))


results <- parLapply(cl, lip, function(i) {

  MASCARA4_test <- function(resids, baits, spikes){
    #target projection
    spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                           resids[,which(colnames(resids) %in% baits)], ncomp = 2)
    
    #####
    y <- colMeans(spls_res$Yloadings)
    U <- spls_res$projection
    UTP <- (dot(t(U),y)/dot(y,y))
    
    sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = TRUE)])
    RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
    
    
    return(list(RP, sPLS_bait_cands,UTP))
  }
  
  
  
  F1_plot <- F1_plot_no_plot
  # nreps <- 3
  alphaN <- 2
  betaN <- 4
  ncands <- 500
  
  
  splsR <- function(resids, baits, spikes){
    #vip
    spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                           resids[,which(colnames(resids) %in% baits)], ncomp = 2)
    
    
    #####
    y <- colMeans(spls_res$Yloadings)
    U <- spls_res$projection
    
    x <- Yloadings(spls_res)
    ve_y <- colSums(x^2)/nrow(x)
    
    sPLS_cands <- cbind.data.frame(abs(spls_res$projection) %*% ve_y, spls_res$projection)
    sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
    
    RP <- prod(which(rownames(sPLS_cands) %in% spikes))^(1/length(spikes))
    
    return(list(RP, sPLS_cands))
  }
  
  
  
  WGCNA4_test <- function(data, baits, spikes){
    
    softPower <- 6
    adjacency <- adjacency(data, power = softPower)
    
    TOM <- TOMsimilarity(adjacency)
    TOM.dissimilarity <- 1-TOM
    geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
    Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 12)
    ModuleColors <- labels2colors(Modules)
    MElist <- moduleEigengenes(data, colors = ModuleColors)
    MEs <- MElist$eigengenes
    ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete")
    
    merge <- mergeCloseModules(data, ModuleColors, cutHeight = .25)
    # mergedColors = merge$colors
    
    mergedMEs = merge$newMEs
    
    
    sKME <- signedKME(data, mergedMEs)
    
    #feature selection
    q_bar <- colMeans(sKME[baits,])
    R <- sKME[!rownames(sKME) %in% baits,]
    R_TP <- (dot(t(R),q_bar)/dot(q_bar,q_bar))
    
    Candidates <- as.data.frame(R_TP[order(R_TP, decreasing = TRUE)])
    
    RP <- prod(which(rownames(Candidates) %in% spikes))^(1/length(spikes))
    
    return(list(RP, Candidates, mergedMEs, sKME))
    
    
  }
  
  Ranked_Coexp_HIGH <- function(X, baits, spikes, ncands){
    
    ranked <- ranked_coexp(baits, X)
    RP <- prod(which(rownames(ranked) %in% spikes))^(1/length(spikes))
    
    return(list(RP, ranked))
    
  }
  
  
  ar <- ASCA(X_funced[[i]], m = meta[[i]], Baits = baits, spikes = spikes,
             ncands = ncands, distance_calc = FALSE, Return_Model = FALSE)
  ASCA_RES_i <- ar[[1]]
  ASCA_CANDS_i <- rownames(ar[[2]])
  
  hcr <- Ranked_Coexp_HIGH(X_funced[[i]], baits, spikes, ncands)
  COEXP_HIGH_RES_i <- hcr[[1]]
  COEXP_HIGH_CANDS_i <- rownames(hcr[[2]])
  
  hcr <- splsR(X_funced[[i]], baits, spikes)  
  PLS_RES_i <- hcr[[1]]
  PLS_CANDS_i <- rownames(hcr[[2]])
  
  MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
  MASCARA_RES_i <- MASCARAh[[1]]
  MASCARA_CANDS_i <- rownames(MASCARAh[[2]])
  
  WGCNAh <- WGCNA4_test(X_funced[[i]], baits, spikes = spikes)
  WGCNA_RES_i <- WGCNAh[[1]]
  WGCNA_CANDS_i <- rownames(WGCNAh[[2]])
  
  # Return the results as a list
  list(
    ASCA_RES_i = ASCA_RES_i,
    ASCA_CANDS_i = ASCA_CANDS_i,
    COEXP_HIGH_RES_i = COEXP_HIGH_RES_i,
    COEXP_HIGH_CANDS_i = COEXP_HIGH_CANDS_i,
    PLS_RES_i = PLS_RES_i,
    PLS_CANDS_i = PLS_CANDS_i,
    MASCARA_RES_i = MASCARA_RES_i,
    MASCARA_CANDS_i = MASCARA_CANDS_i,
    WGCNA_RES_i = WGCNA_RES_i,
    WGCNA_CANDS_i = WGCNA_CANDS_i
  )
  

})

# Stop the cluster after computation
stopCluster(cl)

# Assemble the results
for (i in 1:length(results)) {
  res <- results[[i]]
  if (!is.null(res$error_index)) {
    EI <- c(EI, res$error_index)
    next
  }
  if (!is.null(res$ASCA_RES_i)) {
    ASCA_RES[, i] <- res$ASCA_RES_i
    ASCA_CANDS[, i] <- res$ASCA_CANDS_i
  }
  if (!is.null(res$COEXP_HIGH_RES_i)) {
    COEXP_HIGH_RES[, i] <- res$COEXP_HIGH_RES_i
    COEXP_HIGH_CANDS[, i] <- res$COEXP_HIGH_CANDS_i
  }
  if (!is.null(res$PLS_RES_i)) {
    PLS_RES[, i] <- res$PLS_RES_i
    PLS_CANDS[, i] <- res$PLS_CANDS_i
  }
  if (!is.null(res$MASCARA_RES_i)) {
    MASCARA_RES[, i] <- res$MASCARA_RES_i
    MASCARA_CANDS[, i] <- res$MASCARA_CANDS_i
  }
  if (!is.null(res$WGCNA_RES_i)) {
    WGCNA_RES[, i] <- res$WGCNA_RES_i
    WGCNA_CANDS[, i] <- res$WGCNA_CANDS_i
  }
  
  # Optionally, print progress
  print(i)
}

# Proceed as before
Noise_sim <- list(ASCA_RES, PLS_RES, COEXP_HIGH_RES, MASCARA_RES, WGCNA_RES)





# 
# 
# 
# 
# 
# ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
# 
# COEXP_HIGH_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
# 
# PLS_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
# 
# MASCARA_RES  <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
# 
# WGCNA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))
# 
# 
# 
# 
# ASCA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
# 
# PLS_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
# 
# COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
# 
# MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
# 
# WGCNA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))
# 
# 
# Baits <- c(baits,spikes)
# 
# 
# EI <- c()
# i <- 1
# while(i < length(X_funced) + 1){
#   
#   tryCatch({
#     
#     ar <- ASCA(X_funced[[i]], m = meta[[i]], Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...
#     ASCA_RES[,i] <- ar[[1]]
#     ASCA_CANDS[,i] <- rownames(ar[[2]])
#     
#     
#     hcr <- Ranked_Coexp_HIGH(X_funced[[i]], baits, spikes, ncands)
#     COEXP_HIGH_RES[,i] <- hcr[[1]]
#     COEXP_HIGH_CANDS[,i] <- rownames(hcr[[2]])
#     
#     hcr <- sPLSr(X_funced[[i]], meta[[i]], baits, spikes, ncands)  
#     PLS_RES[,i] <- hcr[[1]]
#     PLS_CANDS[,i] <- rownames(hcr[[2]])
#     
#     
#     
#     MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
#     MASCARA_RES[,i] <- MASCARAh[[1]]
#     MASCARA_CANDS[,i] <- rownames(MASCARAh[[2]])
#     
#     WGCNAh <- WGCNA4_test(X_funced[[i]], baits, spikes = spikes)
#     WGCNA_RES[,i] <- WGCNAh[[1]]
#     WGCNA_CANDS[,i] <- rownames(WGCNAh[[2]])
#     
#     
#     
#     i <- i + 1
#     
#   },
#   error=function(e){
#     print("Nog eentje")
#     #log which method and which feature set(?) makes it fail?
#     error_index <- i
#     EI <- c(EI,error_index)
#   },
#   finally = {})
#   
#   print(i)
#   
# }
# 
# 
# 
# 
# Noise_sim <- list(ASCA_RES, PLS_RES, COEXP_HIGH_RES, MASCARA_RES, WGCNA_RES)

###################################################



# saveRDS(Noise_sim, "~/personal/23_03_08_Tests/Noise_Sim_24_01_22_Replicate_Tests.RDS")  #Noise_Sim_23_06_27_Replicate_Tests.RDS

saveRDS(Noise_sim, "~/personal/23_03_08_Tests/Noise_Sim_24_12_08_Replicate_Tests_3.RDS")  #Noise_Sim_23_06_27_Replicate_Tests.RDS
#Noise_Sim_24_05_24_Replicate_Tests.RDS
###################################################

