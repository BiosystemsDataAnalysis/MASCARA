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

source("DATA_SIM_FUNCS_DEV.R")
F1_plot <- F1_plot_no_plot
nreps <- 3
ncands <- 500

Create_Core_DEV <- Create_Core_DEV_2

sPLSr <- function(X, meta, baits, spikes, ncands, TITLE = NULL){
  
  
  splsR <- spls(X[,-which(colnames(X) %in% baits)],X[,which(colnames(X) %in% baits)], ncomp = 2)
  
  # splsR <- simpls.fit(X[,-which(colnames(X) %in% baits)],X[,which(colnames(X) %in% baits)], ncomp = 2)
  
  sPLS_cands <- cbind.data.frame(abs(splsR$loadings$X) %*% splsR$prop_expl_var$X, splsR$loadings$X)
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = T),]
  
  F1_scores <- F1_plot(sPLS_cands, spikes, ncands, TITLE = paste0(": ",TITLE))
  
  return(list(F1_scores,sPLS_cands))
  
}



MASCARA4_test <- function(resids,ref, baits, spikes, ncands){
  
  Fref <- ref
  
  spls_res <- spls(resids[,-which(colnames(resids) %in% baits)], 
                   resids[,which(colnames(resids) %in% baits)], all.outputs = T)
  
  # spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)], 
  #                  resids[,which(colnames(resids) %in% baits)], all.outputs = T)
  
  sPLS_cands <- cbind.data.frame(abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X, spls_res$loadings$X) #  * spls_res$prop_expl_var$X
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = T),]
  
  importance <- paste0(round(spls_res$prop_expl_var$X * 100, digits = 2), "%")
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_cands),Fref2$Feature),]
  
  #####

  y <- colMeans(spls_res$loadings$Y)

  U <- spls_res$loadings$X
  
  # UTP <- ((U %*% y)/norm(y, "2")) * y
  # UTP <- ((U %*% y)/dot(y,y)) * y   #_9
  
  # UTP <- (dot(t(U),y)/dot(y,y)) * y  #_10
  
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(UTP, decreasing = T)])
  
  # sPLS_bait_cands <- as.data.frame(UTP[order(UTP[,1], decreasing = T),])  #_9 and below
  
  sPLS_maxima_cands <- sPLS_bait_cands
  
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_maxima_cands),Fref2$Feature),]
  
 
  RP <- prod(which(rownames(sPLS_maxima_cands) %in% spikes))^(1/length(spikes))
  
  

  return(list(RP, sPLS_maxima_cands,UTP))
}



# baits <- paste0("X_",c(1989:1996))   #501:510,
# baits <- paste0("X_",c(1985:1996))   #501:510

# baits <- paste0("X_",c(1989:1992))
# spikes <- paste0("X_",c(1977:1988))
baits <- paste0("X_",c(1985:1988))   #501:510
spikes <- paste0("X_",c(1989:2000))


meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = 3), rep(c(-1:-4,1:4), each = 3))
colnames(meta) <- c("growth_condition","time","interaction")
# meta <- as.data.frame(apply(meta,2,factor))

# meta <- as.data.frame(apply(meta,2,factor))


sim_data <- Create_Core_DEV(nreps = 3, meta = meta, plot = F, 
                            EffectSize = c(X_a_ab = 1, time = 1, E = 0.5, mu = 1, Struc_E = 1.5), 
                            struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                            b_sigma = c(1,0), irr_spikes = TRUE, SCORE_SEED = 1027)

ref <- cbind.data.frame(Feature = colnames(sim_data[[1]]), sim_data[[2]], sim_data[[3]])
cols <- colnames(ref[,-1])
ref$Effect <- do.call(paste, c(ref[cols], sep = "_"))
ref$Effect[which(ref$Effect == "0.5_0_0_0.1")] <- "Low PC1"
ref$Effect[which(ref$Effect == "1_0_0_0.1")] <- "Medium PC1"
ref$Effect[which(ref$Effect == "1.5_0_0_0.1")] <- "High PC1"
ref$Effect[which(ref$Effect == "0_0.5_0.1_0")] <- "PC2"
ref$Effect[grep("^0", ref$Effect)] <- "None"
ref$Effect <- factor(ref$Effect, levels = c("Low PC1", "Medium PC1", "High PC1", "PC2", "None"))
colnames(ref)[8] <- "Baits"

ref$Description <- ref$Baits



ref$Baits[1985:1988] <- "Baits"
ref$Baits[1989:2000] <- "Candidates"

ref$Baits[(2000 - (16 + Experiment_responders)):1984] <- "ab1+ve"
# ref$Baits[1969:1976] <- "ab1+ve"
ref$Baits[501:515] <- "ab2+ve"
ref$Baits[486:500] <- "ab2-ve"
# ref$Baits[1969:1976] <- "Expt_Independent"

ref$Baits[1:24] <- "ab1-ve"



ref$Description <- ref$Baits


# a_ab_p3 <- c(rep(0,185),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1685))
#   a_ab_p4 <- c(rep(0,85),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1885))
# 

# baits <- paste0("X_",c(1977:1992))   #501:510

# E_lev <-  c(0:20)/10
# time_lev <- c(0:20)/10




# y <- rep((c(0:20)/33.3333333333333333333333333333333),40)
# x <-  c(1:10) * 10
# x <- rep(x,each = 84)
# 
# z <- x*y


# y <- rep((c(0:10)/33.3333333333333333333333333333333)*2)
# x <-  c(0:10) * 5
# 
# 
# x <- rep(x,each = length(y))
# y <- rep(y, length(unique(x)))
# 
# 
# z <- x*y

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


# main_lev  <-  c(0:20)/4
# 
# # E_lev <- c(0:10)/20
# ERs <- c(1:10) * 10
# 
# X_funced <- list()
# names <- c()
# i <- NULL
# l <- 1
# for(i in 1:length(ERs)){
#   j <- NULL
#   for(j in 1:length(main_lev)){
#     
#     
#     X_funced[[l]] <- Create_Core_DEV(nreps = 3, meta = meta,
#                                      EffectSize = c(X_a_ab = main_lev[j], time = 1, E = 0.3, mu = 1, Struc_E = 0.5), 
#                                      struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
#                                      b_sigma = c(1,0), irr_spikes = FALSE, SCORE_SEED = 2541, 
#                                      plot = FALSE, Experiment_responders = ERs[i], ts = 3)
#     
#     names[l] <- paste(ERs[i],main_lev[j], sep = "_")
#     l <- l + 1
#     
#     
#   }
#   
#   
# }


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
    
    
    
    MASCARAh <- MASCARA4_test(ar[[3]], ref, baits, spikes = spikes, ncands)
    MASCARA_RES[,i] <- MASCARAh[[1]]
    MASCARA_CANDS[,i] <- rownames(MASCARAh[[2]])
    
    
    ##add MASCA here and below for save RDS 
    
    
    
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

# # saveRDS(Noise_sim, "Noise_Sim_220926.RDS")

# saveRDS(Noise_sim, "Noise_Sim_Coexp_23_06_28.RDS")

saveRDS(Noise_sim, "Noise_Sim_Coexp_23_07_26_11.RDS")


###################################################




# 
# Noise_cands <- list(ASCA_CANDS, PLS_CANDS, COEXP_HIGH_CANDS,MASCARA_CANDS)
# 
# # saveRDS(Noise_sim, "Noise_Sim_220926.RDS")
# 
# saveRDS(Noise_cands, "Noise_Sim_Coexp_23_02_19s_Cands_6.RDS")
# # saveRDS(X_funced, "Datasets_Noise_Sim_220926_Noise_Sim_Tests_6.RDS")
# 
