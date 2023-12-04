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

source("DATA_SIM_FUNCS_DEV.R")
source("MASCARA_Test_Funcs.R")

Create_Core_DEV <- Create_Core_DEV_2

F1_plot <- F1_plot_no_plot
# nreps <- 3
alphaN <- 2
betaN <- 4
ncands <- 500
# baits <- paste0("X_",c(1989:1996))   #501:510,
# baits <- paste0("X_",c(1985:1996))   #501:510





# baits <- paste0("X_",c(1997:2000))
# spikes <- paste0("X_",c(1985:1996))

baits <- paste0("X_",c(1985:1988))   #501:510
spikes <- paste0("X_",c(1989:2000))


Baits <- c(baits,spikes)

Experiment_responders <- 12

nreps <- 3
meta <- cbind.data.frame(rep(c(1,-1), each = 12), rep(c(1:4), each = nreps))
colnames(meta) <- c("growth_condition","time")
meta$ID <- paste(meta$growth_condition, meta$time, sep = "_")

# sim_data <- Create_Core(nreps = 3, meta, irr_spikes = TRUE, struc_resid = FALSE, 
#                         a_sigma = c(1,1), b_sigma = c(1.5,1.5), e_sigma = c(1,0.3),
#                         noise_sd = 0.75, EffectSize = c(X_a_ab = 1, time = 1, E = 1, mu = 1), 
#                         plot = FALSE)

sim_data <- Create_Core_DEV(nreps = 3, meta = meta, plot = F, 
                            EffectSize = c(X_a_ab = 1, time = 1, E = 0.5, mu = 1, Struc_E = 1.5), 
                            struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1.5,0.7,0.6,0.1),
                            b_sigma = c(1,0.9), irr_spikes = TRUE, SCORE_SEED = 1027)

ref <- cbind.data.frame(Feature = colnames(sim_data[[1]]), sim_data[[2]], sim_data[[3]])
cols <- colnames(ref[,-1])
ref$Effect <- do.call(paste, c(ref[cols], sep = "_"))
ref$Effect[which(ref$Effect == "0.5_0_0_0.1")] <- "Low PC1"
ref$Effect[which(ref$Effect == "1_0_0_0.1")] <- "Medium PC1"
ref$Effect[which(ref$Effect == "1.5_0_0_0.1")] <- "High PC1"
ref$Effect[which(ref$Effect == "0_0.5_0.1_0")] <- "PC2"
ref$Effect[grep("^0", ref$Effect)] <- "None"
ref$Effect <- factor(ref$Effect, levels = c("Low PC1", "Medium PC1", "High PC1", "PC2", "None"))
colnames(ref)[6] <- "Baits"

ref$Description <- ref$Baits

ref$Baits <- as.character(ref$Baits)


ref$Baits[1985:1988] <- "Baits"
ref$Baits[1989:2000] <- "Candidates"

ref$Baits[(2000 - (16 + Experiment_responders)):1984] <- "ab1+ve"
# ref$Baits[1969:1976] <- "ab1+ve"
ref$Baits[501:515] <- "ab2+ve"
ref$Baits[486:500] <- "ab2-ve"
# ref$Baits[1969:1976] <- "Expt_Independent"

ref$Baits[1:24] <- "ab1-ve"

ref$Description <- ref$Baits






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
# meta <- list()

norms <- list()

l <- 1
# names <- c()



y <- rep((c(0:10)/50))
x <-  c(0:10) * 10     #5


x <- rep(x,each = length(y))
y <- rep(y, length(unique(x)))


z <- x*y


main_lev  <-  z

# E_lev <- c(0:10)/20
ERs <- x

E_in <- c(5,10,15)


X_funced <- list()
names <- c()
i <- NULL
l <- 1

k <- NULL

for(i in 1:length(E_in)){
  
  for(k in 1:30){


      X_funced[[l]] <- Create_Core_DEV_2(nreps = 3, meta = meta, plot = FALSE,
                        EffectSize = c(X_a_ab = 1.8, time = 1, E = E_in[i], mu = 1, Struc_E = 1),
                        struc_resid = T, e_sigma = c(1.5,0.8,0.0), a_sigma = c(1,0.7,0.6,0.1),
                        b_sigma = c(1,0.9), irr_spikes = TRUE, SCORE_SEED = 12, Experiment_responders = Experiment_responders,
                        struc_seed = (127+k))
  
      
      l <- l + 1
      
      
      # names[l] <- paste(ERs[i],main_lev[i], sep = "_")
      
  }
}






ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))


TP_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

SR_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

VIP_RES <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

SVD1_RES  <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))

# TPSR_RES  <- data.frame(matrix(nrow = ncands, ncol = length(X_funced)))



TP_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

SR_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

VIP_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

SVD1_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))

# TPSR_CANDS <- data.frame(matrix(nrow = ncol(X_funced[[1]][[1]]) - 4, ncol = length(X_funced)))


Baits <- c(baits,spikes)


EI <- c()
i <- 1
while(i < length(X_funced) + 1){
  
  tryCatch({
    
    ar <- ASCA(X_funced[[i]][[1]], m = meta, Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...

    

    
    MASCARAh <- MASCARA4_test_2(ar[[3]], baits, spikes = spikes)
    TP_RES[,i] <- MASCARAh[[1]]
    TP_CANDS[,i] <- rownames(MASCARAh[[2]])
    
    
    MASCARAh <- MASCARA4_test_7(ar[[3]], baits, spikes = spikes)
    SR_RES[,i] <- MASCARAh[[1]]
    SR_CANDS[,i] <- rownames(MASCARAh[[2]])
    
    MASCARAh <- MASCARA4_test_8(ar[[3]], baits, spikes = spikes)
    VIP_RES[,i] <- MASCARAh[[1]]
    VIP_CANDS[,i] <- rownames(MASCARAh[[2]])
    
    MASCARAh <- MASCARA4_test_3(ar[[3]], baits, spikes = spikes)
    SVD1_RES[,i] <- MASCARAh[[1]]
    SVD1_CANDS[,i] <- rownames(MASCARAh[[2]])
    
    # MASCARAh <- MASCARA4_test(ar[[3]], ref, baits, spikes = spikes, ncands)
    # TPSR_RES[,i] <- MASCARAh[[1]]
    # TPSR_CANDS[,i] <- rownames(MASCARAh[[2]])
    # 
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




Noise_sim <- list(TP_RES, SR_RES, VIP_RES, SVD1_RES)  #, TPSR_RES

###################################################


# saveRDS(Noise_sim, "Noise_Sim_220926.RDS")

saveRDS(Noise_sim, "~/personal/23_03_08_Tests/Noise_Sim_23_07_30_FS_Tests.RDS")  #Noise_Sim_23_06_27_Replicate_Tests.RDS

# 
# saveRDS(Noise_cands, "Noise_Sim_221004_Replicate_Cands_8.RDS")
# 
# saveRDS(X_funced, "Datasets_Noise_Sim_221004_Replicate_Tests_8.RDS")

###################################################

