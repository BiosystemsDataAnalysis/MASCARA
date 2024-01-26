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
library(pls)

# source("../DATA_SIM_FUNCS.R")

source("DATA_SIM_FUNCS_DEV.R")

F1_plot <- F1_plot_no_plot
nreps <- 3
ncands <- 500

Create_Core_DEV <- Create_Core_DEV_2


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




# spikes <- paste0("X_",c(1977:1988))
baits <- paste0("X_",c(1985:1988))   #501:510
spikes <- paste0("X_",c(1989:2000))


meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = 3), rep(c(-1:-4,1:4), each = 3))
colnames(meta) <- c("growth_condition","time","interaction")


sim_data <- Create_Core_DEV(nreps = 3, meta = meta, plot = F, 
                            EffectSize = c(X_a_ab = 1, time = 1, E = 0.5, mu = 1, Struc_E = 1.5), 
                            struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                            b_sigma = c(1,0), irr_spikes = TRUE, SCORE_SEED = 1027)



y <- rep((c(0:10)/50))
x <-  c(0:10) * 10     #5


x <- rep(x,each = length(y))
y <- rep(y, length(unique(x)))


z <- x*y


main_lev  <-  z

# E_lev <- c(0:10)/20
ERs <- x

X_funced <- list()
names <- c()
i <- NULL
l <- 1
for(i in 1:length(ERs)){

    
    X_funced[[l]] <- Create_Core_DEV(nreps = 3, meta = meta,
                                     EffectSize = c(X_a_ab = main_lev[i], time = 1, E = 0.3, mu = 1, Struc_E = 0.5), 
                                     struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                                     b_sigma = c(1,0), irr_spikes = FALSE, SCORE_SEED = 1000,  #_8 = 1027
                                     plot = FALSE, Experiment_responders = ERs[i], ts = 1234)
    
    names[l] <- paste(ERs[i],main_lev[i], sep = "_")
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



ASCA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

PLS_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

Baits <- c(baits,spikes)

EI <- c()
i <- 1
while(i < length(X_funced) + 1){
  
  tryCatch({
    
    ar <- ASCA(X_funced[[i]][[1]], m = meta, Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...
    ASCA_RES[,i] <- ar[[1]]
    ASCA_CANDS[,i] <- rownames(ar[[2]])
    
    
    hcr <- Ranked_Coexp_HIGH(X_funced[[i]][[1]], baits, spikes, ncands)
    COEXP_HIGH_RES[,i] <- hcr[[1]]
    COEXP_HIGH_CANDS[,i] <- rownames(hcr[[2]])
    
    hcr <- sPLSr(X_funced[[i]][[1]], meta, baits, spikes, ncands)  
    PLS_RES[,i] <- hcr[[1]]
    PLS_CANDS[,i] <- rownames(hcr[[2]])
    
    MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
    MASCARA_RES[,i] <- MASCARAh[[1]]
    MASCARA_CANDS[,i] <- rownames(MASCARAh[[2]])
    


    
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


# ns <- list(ASCA_DIST_LOW, ASCA_DIST_HIGH, PCA_DIST_LOW, PCA_DIST_HIGH)
# saveRDS(ns, "Noise_Sim_220926_2_redo_dist.RDS")



Noise_sim <- list(ASCA_RES, PLS_RES, COEXP_HIGH_RES, MASCARA_RES)

## # saveRDS(Noise_sim, "Noise_Sim_220926.RDS")

## saveRDS(Noise_sim, "Noise_Sim_Coexp_23_06_28.RDS")

# saveRDS(Noise_sim, "Noise_Sim_Coexp_23_07_26_11.RDS")


saveRDS(Noise_sim, "Noise_Sim_Coexp_24_01_22.RDS")


