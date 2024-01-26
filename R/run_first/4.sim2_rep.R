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

# source("../DATA_SIM_FUNCS.R")

source("DATA_SIM_FUNCS_DEV.R")


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
# names <- c()




#######
# y <- c(0.01, 0.05, 0.1)
# x <- c(100,100,100)

y <- rep(c(0.0004,0.002, 0.01),3)      

x <- rep(30,9)  #100

nse <- rep(c(0.1,0.3,0.5), each = 3)

# y <- c(0.1, 2, 4)
# x <- c(25,125,250) 



main_lev <- x*y
ERs <- x

#######

idlist <- list()

a <- NULL

# set.seed(40001)

for(a in 1:length(ERs)){
  
  i <- NULL
  
  for(i in 1:20){
    
    
    X_30[[i]] <- Create_Core_DEV(nreps = nreps, meta = big_meta, plot = FALSE, 
                    EffectSize = c(X_a_ab = main_lev[a], time = 0.2, E = nse[a], mu = 0.2, Struc_E = 0.05), 
                    struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                    b_sigma = c(1,0), irr_spikes = FALSE, SCORE_SEED = (250021 + i), Experiment_responders = ERs[a], ts = 3 + i) 
    
    
    
    
    k <- NULL
    for(k in 3:nreps){
      
      IDS <- NULL
      
      j <- NULL
      for(j in 1:length(conds)){
        
        set.seed(40001 + l*2)
        ids <- sample(big_meta$reps[which(big_meta$cond == conds[j])], k)
        IDS <- c(IDS,ids)
        
        idlist[[l]] <- IDS
        
        
        
      }
      
      X_funced[[l]] <- X_30[[i]][[1]][which(big_meta$reps %in% IDS),]
      meta[[l]] <- big_meta[which(big_meta$reps %in% IDS),]
      
      norms[[l]] <- X_30[[i]][[6]]
      
      l <- l + 1
      
      
      
    }
    
    
    
  }
  
  
  
  
}





Ranked_Coexp_HIGH <- Ranked_Coexp_HIGH.d <- function(X, baits, spikes, ncands){
  
  ranked <- ranked_coexp(baits, X)
  RP <- prod(which(rownames(ranked) %in% spikes))^(1/length(spikes))
  
  return(list(RP, ranked))
  
}





ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

COEXP_HIGH_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

PLS_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

MASCARA_RES  <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))



ASCA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))

PLS_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))

COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))

MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]]) - 4, ncol = length(X_funced)))

Baits <- c(baits,spikes)


EI <- c()
i <- 1
while(i < length(X_funced) + 1){
  
  tryCatch({
    
    ar <- ASCA(X_funced[[i]], m = meta[[i]], Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...
    ASCA_RES[,i] <- ar[[1]]
    ASCA_CANDS[,i] <- rownames(ar[[2]])
    
    
    hcr <- Ranked_Coexp_HIGH(X_funced[[i]], baits, spikes, ncands)
    COEXP_HIGH_RES[,i] <- hcr[[1]]
    COEXP_HIGH_CANDS[,i] <- rownames(hcr[[2]])
    
    hcr <- sPLSr(X_funced[[i]], meta[[i]], baits, spikes, ncands)  
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




Noise_sim <- list(ASCA_RES, PLS_RES, COEXP_HIGH_RES, MASCARA_RES)

###################################################


# saveRDS(Noise_sim, "Noise_Sim_220926.RDS")

# saveRDS(Noise_sim, "~/personal/23_03_08_Tests/Noise_Sim_23_07_30_Replicate_Tests.RDS")  #Noise_Sim_23_06_27_Replicate_Tests.RDS

saveRDS(Noise_sim, "~/personal/23_03_08_Tests/Noise_Sim_24_01_22_Replicate_Tests.RDS")  #Noise_Sim_23_06_27_Replicate_Tests.RDS


# 
# saveRDS(Noise_cands, "Noise_Sim_221004_Replicate_Cands_8.RDS")
# 
# saveRDS(X_funced, "Datasets_Noise_Sim_221004_Replicate_Tests_8.RDS")

###################################################

