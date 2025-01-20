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


# source("../DATA_SIM_FUNCS.R")
# setwd("C:/Users/fwhite/Documents/GitHub/MASCARA/R")

setwd("/zfs/omics/personal/fwhite/MASCARA/R")

source("DATA_SIM_FUNCS.R")


Create_Core_DEV <- Create_Core_DEV_2

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
  browser()
  Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 4)
  
  # if(length(Modules > 15))
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



# # spikes <- paste0("X_",c(1977:1988))
# baits <- paste0("X_",c(1985:1988))   #501:510
# spikes <- paste0("X_",c(1989:2000))
# 

meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = 3), rep(c(-1:-4,1:4), each = 3))
colnames(meta) <- c("growth_condition","time","interaction")


sim_data <- Create_Core_DEV(nreps = 3, meta = meta, plot = F, 
                            EffectSize = c(X_a_ab = 1, time = 1, E = 0.5, mu = 1, Struc_E = 1.5), 
                            struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                            b_sigma = c(1,0), irr_spikes = TRUE, SCORE_SEED = 1027)


###################
source("DATA_SIM_FUNCS.R")



# baits 
# spikes


# y <- rep(c(0.008333333, 0.016666667, 0.025000000),3)
# x <- rep(30,9)  #100
# nse <- rep(c(0.2,0.4,0.6), each = 3)
# 
# 
# 
# 
# y <- rep((c(0:10)/50))
# x <-  c(0:10) * 10     #5
# 
# 
# x <- rep(x,each = length(y))
# y <- rep(y, length(unique(x)))
# 
# 
# z <- x*y
# 
# 
# main_lev  <-  z


bx <- c(2:15)
sx <-c(2:50)


POI <- expand_grid(bx,sx)


# E_lev <- c(0:10)/20
# ERs <- x

X_funced <- list()
names <- c()
i <- NULL
l <- 1
for(i in 1:nrow(POI)){
  
  
  X_funced[[l]] <- Create_Core_DEV_3(nreps = 3, meta = meta,
                                   EffectSize = c(X_a_ab = 0.5, time = 1, E = 0.4, mu = 1, Struc_E = 0.05), 
                                   struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                                   b_sigma = c(1,0), irr_spikes = FALSE, SCORE_SEED = 1000,  #_8 = 1027
                                   plot = FALSE, Experiment_responders = 1000, ts = 1234, baits = c(as.matrix(POI[l,1])), spikes = c(as.matrix(POI[l,2])))
  
  names[l] <- paste(baits = POI[l,1], spikes = POI[l,2], sep = "_")
  l <- l + 1
  
  
  
  
  
}




Ranked_Coexp_HIGH <- function(X, baits, spikes, ncands){
  # baits_small <- baits[c(9:12)]
  
  ranked <- ranked_coexp(baits, X)
  
  # F1_scores <- F1_plot(ranked, spikes, ncands, TITLE = ": Ranked Correlation")
  
  RP <- prod(which(rownames(ranked) %in% spikes))^(1/length(spikes))
  
  
  return(list(RP, ranked))
  
}




ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

COEXP_HIGH_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

PLS_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

MASCARA_RES  <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

WGCNA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))




ASCA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 2, ncol = length(X_funced)))

PLS_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 2, ncol = length(X_funced)))

COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 2, ncol = length(X_funced)))

MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 2, ncol = length(X_funced)))

WGCNA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 2, ncol = length(X_funced)))



EI <- c()
i <- 1
while(i < length(X_funced) + 1){
  
  tryCatch({
    
    
    p1 <- sum(POI[i,])
    
    b1 <- c(data.matrix(POI[i,1]))
    s1 <- c(data.matrix(POI[i,2]))
    
    baits <- paste0("X_",c((ncol(X_funced[[i]][[1]])-p1+1): (ncol(X_funced[[i]][[1]])-p1+b1) ))
    
    spikes <- paste0("X_",c((ncol(X_funced[[i]][[1]])-s1+1): (ncol(X_funced[[i]][[1]] )) ))
    

                                    
    ar <- ASCA(X_funced[[i]][[1]], m = meta, Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...
    
    ind <- nrow(ASCA_CANDS)

    ASCA_RES[,i] <- ar[[1]]
    ASCA_CANDS[,i] <- c(rownames(ar[[2]]),rep(0,(ind - length(rownames(ar[[2]])))))
    
    
    hcr <- Ranked_Coexp_HIGH(X_funced[[i]][[1]], baits, spikes, ncands)
    COEXP_HIGH_RES[,i] <- hcr[[1]]
    COEXP_HIGH_CANDS[,i] <- c(rownames(hcr[[2]]),rep(0,(ind - length(rownames(hcr[[2]])))))
    
    hcr <- sPLSr(X_funced[[i]][[1]], meta, baits, spikes, ncands)  
    PLS_RES[,i] <- hcr[[1]]
    PLS_CANDS[,i] <- c(rownames(hcr[[2]]),rep(0,(ind - length(rownames(hcr[[2]])))))
    
    MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
    MASCARA_RES[,i] <- MASCARAh[[1]]
    MASCARA_CANDS[,i] <- c(rownames(MASCARAh[[2]]),rep(0,(ind - length(rownames(MASCARAh[[2]])))))
    
    # WGCNAh <- WGCNA4_test(X_funced[[i]][[1]], baits, spikes = spikes)
    # WGCNA_RES[,i] <- WGCNAh[[1]]
    # WGCNA_CANDS[,i] <- rownames(WGCNAh[[2]])
    # 
    
    
    i <- i + 1
    
  },
  error=function(e){
    print("Nog eentje")
    #log which method and which feature set(?) makes it fail?
    error_index <- i
    EI <- c(EI,error_index)
  },
  finally = {})
  
  print(i)
  
}




Noise_sim <- list(ASCA_RES, PLS_RES, COEXP_HIGH_RES, MASCARA_RES, WGCNA_RES)


# saveRDS(Noise_sim, "Noise_Sim_Coexp_24_01_22.RDS")

saveRDS(Noise_sim, "Noise_Sim_b_s_24_12_04.RDS")   #Noise_Sim_Coexp_24_05_24.RDS






#divide by RP by  (prod(c(1:r))^(1/r))/r  #r is number of spikes (true positives)
