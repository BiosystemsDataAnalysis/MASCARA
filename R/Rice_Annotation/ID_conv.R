library(riceidconverter)
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
library(WGCNA)
library(DESeq2)



setwd("C:/Users/fwhite/Documents/Github/MASCARA/R")  # - WARNING -


source("DATA_SIM_FUNCS.R")


# ref <- read.table("R/Rice_Annotation/Oryza_sativa.IRGSP-1.0.50.gtf", sep = "\t")
# ref <- ref[which(ref$V3 == "transcript"),]
# ref$ID <- gsub(".*gene_id ","",ref$V9)
# ref$ID <- gsub(";.*","",ref$ID)
# 
# ref$transcript <- gsub(".*transcript_id ","",ref$V9)
# ref$transcript <- gsub(";.*","",ref$transcript)
# 
# ref$biotype <- gsub(".*transcript_biotype ","",ref$V9)
# ref$biotype <- gsub(";.*","",ref$biotype)
#


X <- readRDS("Data_Haider/Rice_Counts4ASCA.RDS")
meta <- readRDS("Data_Haider/Rice_Meta4ASCA.RDS")
colnames(meta)[3] <- c("time")

meta$time <- as.character(meta$time)

meta$time[which(meta$time == "8")] <- "1"
meta$time[which(meta$time == "10")] <- "3"
meta$time[which(meta$time == "14")] <- "7"
meta$time[which(meta$time == "15")] <- "8"


X <- t(assay(X))

meta <- meta[order(meta$time),]
meta <- meta[order(meta$growth_condition),]

X <- X[match(meta$ID, rownames(X)),]



#X <- X[,which(colnames(X) %in% rownames(sPLS_bait_cands))]

symbol <- RiceIDConvert(unique(colnames(X)),'TRANSCRIPTID',toType = 'SYMBOL')
symbol$SYMBOL <- gsub("LOC","osa:",symbol$SYMBOL)
symbol <- symbol[grep("^osa",symbol$SYMBOL),]
# KOIDs <- read.table("https://rest.kegg.jp/link/osa/ko")
path <- read.table("https://rest.kegg.jp/link/osa/pathway")
# module <- read.table("https://rest.kegg.jp/link/osa/module")

# react <- read.table("https://rest.kegg.jp/link/osa/reaction")


colnames(path) <- c("PATH","SYMBOL")
merged <- merge(symbol, path, by = "SYMBOL")


# pnum <- as.data.frame(table(merged$PATH))
# 
# # same_path <- readRDS("path_example.RDS")
# # pnum <- pnum[which(pnum[,1] %in% same_path[,1]),]
# 
# pnum <- pnum[-which(pnum$Freq < 4),]
# 
# merged <- merged[which(merged$PATH %in% pnum$Var1),]
# 
# 
# 
# 
# 
# if (!requireNamespace("KEGGREST", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install("KEGGREST")
# }
# 
# library(KEGGREST)
# 
# 
# kegg_ids <- unique(gsub("path:","",merged$PATH))
# 
# # - STUPID FUNCTION BELOW -  ??
# 
# # Function to get pathway description for each ID
# get_kegg_description <- function(kegg_id) {
#   pathway_info <- keggGet(kegg_id)  # Query KEGG
#   if (!is.null(pathway_info[[1]])) {
#     return(pathway_info[[1]]$NAME)  # Extract the name/description
#   } else {
#     return(NA)  # Return NA if the ID is invalid
#   }
# }
# 
# # Apply the function to all KEGG IDs
# kegg_descriptions <- sapply(kegg_ids, get_kegg_description)
# 
# # Combine IDs and descriptions into a data frame
# kegg_mapping <- data.frame(
#   KEGG_ID = kegg_ids,
#   Description = kegg_descriptions,
#   stringsAsFactors = FALSE
# )
# 
# # Print the mapping
# print(kegg_mapping)
# 
# 
# kegg_DESC <- t(do.call(cbind,kegg_descriptions))



#
#   # WIP for reaction ID update
# library(KEGGREST)
# rxn2enz <- keggLink("enzyme", "reaction")
# enz2osa <- keggLink("osa", "enzyme")
#
# rxnrnz <- data.frame(rxn=names(rxn2enz), enz=rxn2enz)
# rxn2enz <- data.frame(rxn=names(rxn2enz), enz=rxn2enz)
# enz2osa <- data.frame(enz=names(enz2osa), osa=enz2osa)
#
# merged2 <- merge(enz2osa, rxn2enz, by = "enz")
#
#

# WIP for reaction ID update
library(KEGGREST)
rxn2enz <- keggLink("enzyme", "reaction")
enz2osa <- keggLink("osa", "enzyme")

#also need KO to osa 

# enz2rxn <- keggLink("osa", "reaction")

rxn2enz <- data.frame(rxn=names(rxn2enz), enz=rxn2enz)
enz2osa <- data.frame(enz=names(enz2osa), osa=enz2osa)

merged2 <- merge(enz2osa, rxn2enz, by = "enz")


merged <- merge(merged2,symbol, by.x = "osa", by.y = "SYMBOL")


# add PATH_ID column here to filter paths better
merged$PATH_RXN <- paste(merged$rxn,merged$TRANSCRIPTID, sep = "_")
merged <- merged[!duplicated(merged$PATH_RXN),]

pnum <- as.data.frame(table(merged$rxn))

# same_path <- readRDS("path_example.RDS")
# pnum <- pnum[which(pnum[,1] %in% same_path[,1]),]



pnum <- pnum[-which(pnum$Freq < 4),]

merged <- merged[which(merged$rxn %in% pnum$Var1),]


# merged <- merged2    # - WARNING -

colnames(merged)[c(3)] <- c("PATH")

X <- X[,which(colnames(X) %in% merged$TRANSCRIPTID)]       #from ~25k to 1.5k IDs




# X <- X[,which(colnames(X) %in% merged$TRANSCRIPTID)]       #from ~25k to 4k IDs



F1_plot <- F1_plot_no_plot
# nreps <- 3
ncands <- 500




# Create_Core_DEV <- Create_Core_DEV_2
ranked_coexp <- function(baits, data){
  
  adj_matrix <- abs(cor(data))
  melt_adj <- corclean(adj_matrix)
  
  ###
  rm(adj_matrix)
  
  i <- NULL
  RES <- list()
  for(i in 1:length(baits)){
    
    RES[[i]] <- baiter(baits[i], melt_adj)
    # print(i)
    
  }
  
  rm(melt_adj)
  names(RES) <- baits
  RES <- RES[c(which(do.call(rbind, lapply(RES,dim))[,1] > 1))]
  l <- lapply(RES, "[[",2)
  l <- lapply(l, as.character)
  
  common <- Reduce(intersect, l)
  rres <- lapply(RES, com_filt, common = common)
  
  
  ####cbind(?) ranks together 
  i <- NULL
  ranked <- data.frame()
  for(i in 1:length(rres)){
    
    if(i == 1){
      
      ranked <- rres[[i]][,c(2,4)]
      
    }else{
      
      ranked <- cbind(ranked,rres[[i]][,4])
      
    }
    
    
  }
  
  colnames(ranked)[-1] <- baits[c(which(do.call(rbind, lapply(RES,dim))[,1] > 1))]
  
  ranked$comb_rank <- rowSums(ranked[,-1])/(dim(ranked)[2] - 1)
  ranked <- ranked[order(ranked$comb_rank),]
  
  rownames(ranked) <- ranked$Feature
  ranked <- ranked[,-1]
  
  return(ranked)
  
}

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



softPower <- 6
adjacency <- adjacency(X, power = softPower)

TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 12)
ModuleColors <- labels2colors(Modules)
MElist <- moduleEigengenes(X, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete")

merge <- mergeCloseModules(X, ModuleColors, cutHeight = .25)
# mergedColors = merge$colors

mergedMEs = merge$newMEs


sKME <- signedKME(X, mergedMEs)


WGCNA4_test <- function(sKME, baits, spikes){
  

  #feature selection
  q_bar <- colMeans(sKME[baits,])
  R <- sKME[!rownames(sKME) %in% baits,]
  R_TP <- (dot(t(R),q_bar)/dot(q_bar,q_bar))
  
  Candidates <- as.data.frame(R_TP[order(abs(R_TP), decreasing = TRUE)])
  
  RP <- prod(which(rownames(Candidates) %in% spikes))^(1/length(spikes))
  
  return(list(RP, Candidates, mergedMEs, sKME))
  
  
}



MASCARA4_test <- function(resids, baits, spikes){
  #target projection
  spls_res <- simpls.fit(resids[,-which(colnames(resids) %in% baits)],
                         resids[,which(colnames(resids) %in% baits)], ncomp = 2)
  
  #####
  y <- colMeans(spls_res$Yloadings)
  U <- spls_res$projection
  UTP <- (dot(t(U),y)/dot(y,y))
  
  sPLS_bait_cands <- as.data.frame(UTP[order(abs(UTP), decreasing = TRUE)])
  RP <- prod(which(rownames(sPLS_bait_cands) %in% spikes))^(1/length(spikes))
  
  
  return(list(RP, sPLS_bait_cands,UTP))
}



Ranked_Coexp_HIGH <- function(X, baits, spikes, ncands){
  # baits_small <- baits[c(9:12)]
  
  ranked <- ranked_coexp(baits, X)
  
  # F1_scores <- F1_plot(ranked, spikes, ncands, TITLE = ": Ranked Correlation")
  
  RP <- prod(which(rownames(ranked) %in% spikes))^(1/length(spikes))
  
  
  return(list(RP, ranked))
  
}




ASCA_RES <- data.frame(matrix(nrow = ncands, ncol = nrow(pnum)))

COEXP_HIGH_RES <- data.frame(matrix(nrow = ncands, ncol = nrow(pnum)))

PLS_RES <- data.frame(matrix(nrow = ncands, ncol = nrow(pnum)))

MASCARA_RES  <- data.frame(matrix(nrow = ncands, ncol = nrow(pnum)))

WGCNA_RES <- data.frame(matrix(nrow = ncands, ncol = nrow(pnum)))




ASCA_CANDS <- data.frame(matrix(nrow = ncol(X), ncol = nrow(pnum)))

PLS_CANDS <- data.frame(matrix(nrow = ncol(X), ncol = nrow(pnum)))

COEXP_HIGH_CANDS <- data.frame(matrix(nrow = ncol(X), ncol = nrow(pnum)))

MASCARA_CANDS <- data.frame(matrix(nrow = ncol(X), ncol = nrow(pnum)))

WGCNA_CANDS <- data.frame(matrix(nrow = ncol(X), ncol = nrow(pnum)))





# i <- NULL
# for(i in 1:3){   #nrow(pnum)
# 
#   pid <- pnum$Var1[i]
#   transcripts_all <- merged$TRANSCRIPTID[which(merged$PATH == pid)]
#   size <- round(log2(length(transcripts_all)))
#   baits <- sample(transcripts_all,size)
#   spikes <- transcripts_all[-which(transcripts_all %in% baits)]
# 
# 
# 
# 
# 
# }


rownames(meta) <- meta[,1]
meta <- meta[,-c(1)]

bbaits <- list()
sspikes <- list()
med_cor <- c()
#Baits <- c(baits,spikes)

EI <- c()
i <- 1
while(i < nrow(pnum) + 1){   #

  tryCatch({


    pid <- pnum$Var1[i]
    transcripts_all <- merged$TRANSCRIPTID[which(merged$PATH == pid)]
    size <- round(log2(length(transcripts_all)))
    baits <- sample(transcripts_all,size)
    spikes <- transcripts_all[-which(transcripts_all %in% baits)]

    bbaits[[i]] <- baits
    sspikes[[i]] <- spikes

    ar <- ASCA(X, m = meta, Baits = baits, spikes = spikes, ncands = ncands, distance_calc = F, Return_Model = F)   ##need to save all metas?...
    ASCA_RES[,i] <- ar[[1]]
    ASCA_CANDS[1:length(rownames(ar[[2]])),i] <- rownames(ar[[2]])

    med_cor[i] <- median(cor(ar[[3]][,which(colnames(ar[[3]]) %in% baits)])[upper.tri(cor(ar[[3]][,which(colnames(ar[[3]]) %in% baits)]))])

    # cor(ar[[3]][,which(colnames(ar[[3]]) %in% baits)])

    MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
    MASCARA_RES[,i] <- MASCARAh[[1]]
    MASCARA_CANDS[1:length(rownames(MASCARAh[[2]])),i] <- rownames(MASCARAh[[2]])

    hcr <- Ranked_Coexp_HIGH(X, baits, spikes, ncands)
    COEXP_HIGH_RES[,i] <- hcr[[1]]
    COEXP_HIGH_CANDS[1:length(rownames(hcr[[2]])),i] <- rownames(hcr[[2]])

    hcr <- sPLSr(X, meta, baits, spikes, ncands)
    PLS_RES[,i] <- hcr[[1]]
    PLS_CANDS[1:length(rownames(hcr[[2]])),i] <- rownames(hcr[[2]])

    # MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
    # MASCARA_RES[,i] <- MASCARAh[[1]]
    # MASCARA_CANDS[1:length(rownames(MASCARAh[[2]])),i] <- rownames(MASCARAh[[2]])

    WGCNAh <- WGCNA4_test(sKME, baits, spikes = spikes)
    WGCNA_RES[,i] <- WGCNAh[[1]]
    WGCNA_CANDS[1:length(rownames(WGCNAh[[2]])),i] <- rownames(WGCNAh[[2]])



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

saveRDS(Noise_sim, "KEGG_Real_Path_Res_24_12_06_3.RDS")   #Noise_Sim_Coexp_24_05_24.RDS
####


s_B_res <- data.matrix(do.call(rbind,lapply(Noise_sim, function(x) x[1,])))

s_B_res <- cbind.data.frame("pid" = pnum$Var1, "Path_Len" = do.call(c,lapply(sspikes,length)), t(s_B_res))

colnames(s_B_res)[3:7] <- c("ASCA","PLS","Correlations","MASCARA","WGCNA")

s_B_res <- melt(s_B_res, id.vars = c("pid", "Path_Len"))
colnames(s_B_res)[3:4] <- c("Method","GMR")



i <- NULL
for(i in 1:nrow(s_B_res)){
  
  s_B_res$GMR[i] <- s_B_res$GMR[i]/ (prod(c(1:s_B_res$Path_Len[i]))^(1/s_B_res$Path_Len[i]))
  
}


s_B_res$GMR <- log2(s_B_res$GMR)

# 
# library(tidyverse)
# library(ggplot2)
# library(dplyr)
# 
# # Example list of data frames
# # df_list <- list(
# #   data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10)),
# #   data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10)),
# #   data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10))
# # )
# 
# # Extract the first row of each data frame and combine into a data frame
# first_rows <- lapply(Noise_sim, function(df) df[1, ])     #which(med_cor > 0.1)
# combined_df <- do.call(rbind, first_rows)
# 
# # Melt the data for ggplot2
# combined_long <- combined_df %>%
#   mutate(dataframe_id = factor(seq_along(Noise_sim))) %>%
#   pivot_longer(cols = -dataframe_id, names_to = "Variable", values_to = "Value")


s_B_res2 <- s_B_res[which(s_B_res$Path_Len < 50),]

# Plot using ggplot2
p <- ggplot(s_B_res2, aes(x = Method, y = GMR, group = pid, alpha = Path_Len)) +
  geom_line() +
  labs(x = "Method", y = "GMR", title = "KEGG_General") +
  theme_bw()

ggsave("KEGG_RICE7.png")






p <- ggplot(s_B_res2, aes(x = Method, y = GMR)) +  #, group = pid, alpha = Path_Len
  geom_boxplot() +
  labs(x = "Method", y = "GMR", title = "KEGG_General") +
  theme_bw()

ggsave("KEGG_RICE6.png")




# Initialize lists and counters
bbaits <- list()
sspikes <- list()
med_cor <- c()
EI <- c()
i <- 1

# Start the main loop over pnum
while (i <= nrow(pnum)) {

  pid <- pnum$Var1[i]
  transcripts_all <- merged$TRANSCRIPTID[which(merged$PATH == pid)]
  size <- round(log2(length(transcripts_all)))

  # Ensure there are enough transcripts to select baits and spikes
  if (length(transcripts_all) <= size) {
    message(sprintf("Not enough transcripts for pid %s to select baits.", pid))
    med_cor[i] <- NA
    EI <- c(EI, i)
    i <- i + 1
    next
  }

  tryCatch({

    # Get expression data for transcripts_all
    expr_data <- X[, colnames(X) %in% transcripts_all]

    # Compute correlation matrix
    cor_matrix <- cor(expr_data)

    # Compute average correlation for each transcript
    avg_cor <- rowMeans(cor_matrix)

    # Select top 'size' transcripts based on average correlation
    ranked_transcripts <- names(sort(avg_cor, decreasing = TRUE))
    baits <- ranked_transcripts[1:size]
    spikes <- setdiff(transcripts_all, baits)
    
    spikes <- ranked_transcripts[(size+1):(size + size)]  #sample(spikes, size)

    # Calculate median correlation among baits
    bait_indices <- which(colnames(expr_data) %in% baits)
    bait_cor_matrix <- cor_matrix[bait_indices, bait_indices]
    med_cor_value <- median(bait_cor_matrix[upper.tri(bait_cor_matrix)])
    med_cor[i] <- med_cor_value

    # Store baits and spikes
    bbaits[[i]] <- baits
    sspikes[[i]] <- spikes

    # Proceed with analyses using the selected baits and spikes
    ar <- ASCA(X, m = meta, Baits = baits, spikes = spikes, ncands = ncands,
               distance_calc = FALSE, Return_Model = FALSE)
    ASCA_RES[, i] <- ar[[1]]
    ASCA_CANDS[1:length(rownames(ar[[2]])), i] <- rownames(ar[[2]])

    MASCARAh <- MASCARA4_test(ar[[3]], baits, spikes = spikes)
    MASCARA_RES[, i] <- MASCARAh[[1]]
    MASCARA_CANDS[1:length(rownames(MASCARAh[[2]])), i] <- rownames(MASCARAh[[2]])

    hcr <- Ranked_Coexp_HIGH(X, baits, spikes, ncands)
    COEXP_HIGH_RES[, i] <- hcr[[1]]
    COEXP_HIGH_CANDS[1:length(rownames(hcr[[2]])), i] <- rownames(hcr[[2]])

    hcr <- sPLSr(X, meta, baits, spikes, ncands)
    PLS_RES[, i] <- hcr[[1]]
    PLS_CANDS[1:length(rownames(hcr[[2]])), i] <- rownames(hcr[[2]])

    WGCNAh <- WGCNA4_test(sKME, baits, spikes = spikes)
    WGCNA_RES[, i] <- WGCNAh[[1]]
    WGCNA_CANDS[1:length(rownames(WGCNAh[[2]])), i] <- rownames(WGCNAh[[2]])

  }, error = function(e) {
    # Handle errors
    message(sprintf("Error at i = %d: %s", i, e$message))
    EI <- c(EI, i)
    med_cor[i] <- NA
  })

  i <- i + 1  
  print(i)

}  
