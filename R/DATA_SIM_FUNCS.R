#######
Create_Core <- function(nreps, meta, irr_spikes = TRUE, struc_resid = FALSE, 
                        a_sigma = c(1,1), b_sigma = c(1.5,1.5), e_sigma = c(1,0.3),
                        noise_sd = 0.75, EffectSize = c(X_a_ab = 1, time = 1, E = 1, mu = 1), 
                        SCORE_SEED = 1000, plot = FALSE, score_plot = FALSE){
  
  ########### A_AB
  ##CREATE SCORES FOR PC1 top half
  set.seed(SCORE_SEED)
  a_ab_PC1 <- rnorm(length(unique(meta$time)))
  a_ab_PC1 <- abs(a_ab_PC1)
  a_ab_PC1 <- a_ab_PC1[order(a_ab_PC1, decreasing = FALSE)]
  a_ab_PC1_small <- a_ab_PC1
  
  ##CREATE SCORES FOR PC2 top half
  set.seed(1)
  a_ab_PC2 <- Create_Orthogonal(a_ab_PC1_small)
  a_ab_PC1 %*% a_ab_PC2 == 0
  a_ab_PC1 <- rep(a_ab_PC1, each = nreps)
  a_ab_PC2 <- rep(a_ab_PC2, each = nreps)
  
  ##CREATE FULL SCORES
  a_ab_PC1 <- matrix(c(a_ab_PC1, a_ab_PC1*(-1)))
  a_ab_PC2 <- matrix(c(a_ab_PC2, a_ab_PC2*(-1)))
  a_ab <- cbind.data.frame(meta[,1:2], a_ab_PC1, a_ab_PC2)
  colnames(a_ab)[c(3,4)] <- c("PC1","PC2")
  a_ab$PC1 %*% a_ab$PC2
  
  #CREATE LOADINGS 
  # a_ab_p1 <- c(rep(0,1988), rep(0.5,4), rep(1,4), rep(1.5,4))
  a_ab_p1 <- c(rep(0,1984), rep(0.5,6), rep(1,6), rep(1.5,4))
  
  
  if(irr_spikes == FALSE){
    a_ab_p2 <- c(rep(0,length(a_ab_p1)))
  }
  else{
    a_ab_p2 <- c(rep(0,500), rep(0.5,10), rep(0,1490))
  }
  P <- cbind(a_ab_p1, a_ab_p2)
  
  P_sigma <- P %*% diag(a_sigma)
  
  ############# B
  ##CREATE SCORES FOR PC1 TIME
  set.seed(2)
  b_PC1_small <- Create_Orthogonal(a_ab_PC1_small)
  b_PC1_small %*% a_ab_PC1_small
  
  ##CREATE SCORES FOR PC2 TIME
  b_PC2_small <- Create_Orthogonal(b_PC1_small)
  b_PC1 <- rep(b_PC1_small, each = nreps, times = 2)
  b_PC2 <- rep(b_PC2_small, each = nreps, times = 2)
  
  b <- cbind.data.frame(meta[,1:2],b_PC1,b_PC2)
  colnames(b)[3:4] <- c("PC1","PC2")
  
  #CREATE LOADINGS TIME
  b_p1 <- c(rep(0.1,1000), rep(0,1000))  #make dependent on length of a_ab_p1
  b_p2 <- c(rep(0,1000), rep(0.1,1000))  #make dependent on length of a_ab_p1
  Pb <- cbind(b_p1, b_p2)
  
  Pb_sigma <- Pb %*% diag(b_sigma)
  
  ################## MATRIX MULTIPLICATION OF SCORES AND LOADINGS PER EFFECT MATRIX
  X_a_ab <- EffectSize["X_a_ab"] * data.matrix(a_ab[c("PC1","PC2")]) %*% t(P_sigma)
  time <- EffectSize["time"] * data.matrix(b[c("PC1","PC2")]) %*% t(Pb_sigma)
  
  ##CREATE FEATURE MEAN VALUES and ensure they are above 0
  #set.seed here? 
  mu <- EffectSize["mu"] * matrix(rep(rnorm(dim(X_a_ab)[2], mean = mean(X_a_ab)), each = nrow(meta)), nrow = nrow(meta))
  mu <- mu + abs(min(mu))
  X_a_ab_mu <- X_a_ab + mu
  
  
  ####CREATE E
  
  if(struc_resid == TRUE){
    
    #make structured residuals instead of random noise:
    set.seed(125)
    e_PC1 <- Create_Orthogonal(a_ab_PC1) #also test orthogonality with other scores ?
    set.seed(126)
    e_PC2 <- Create_Orthogonal(e_PC1)
    
    #make loadings
    set.seed(127)
    # e_p1 <- c(rep(0,1988), rep(1,4), rep(0,8))  #make dependent on length of a_ab_p1
    e_p1 <- c(rep(0,1984), rep(1,12), rep(0,4))
    
    set.seed(128)
    # e_p2 <- c(rep(0,0), rnorm(2000))  #make dependent on length of a_ab_p1
    e_p2 <- c(rep(0.1,1984), rep(0,12), rep(0.1,4))
    
    Pe <- cbind(e_p1, e_p2)
    colnames(Pe) <- c("PC1","PC2")
    
    Pe_sigma <- Pe %*% diag(e_sigma)
    
    e <- cbind.data.frame(meta[,1:2], e_PC1, e_PC2)
    colnames(e)[c(3,4)] <- c("PC1","PC2")
    
    set.seed(1001)
    Struc_E <- data.matrix(e[c("PC1","PC2")]) %*% t(Pe_sigma)
    
    # E <- EffectSize["E"] * data.matrix(e[c("PC1","PC2")]) %*% t(Pe_sigma)
    E <- Struc_E + EffectSize["E"] * matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    
  }
  else{
    
    ## some random noise - i.e. differences between replicates
    set.seed(101)
    E <- EffectSize["E"] * matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    
    
  }
  
  
  ##ADD ALL TOGETHER TO CREATE "CORE"
  X_no.time <- X_a_ab_mu + E
  X <- X_no.time + time
  X <- t(((t(X) + abs(colMins(X)))))
  colnames(X) <- paste0("X_",c(1:dim(X)[2]))
  
  
  if(plot == TRUE){
    
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    t2 <- ggplot(a_ab, aes(x = time, y = PC2, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    P_anno <- cbind.data.frame(c(1:dim(P)[1]), P)
    colnames(P_anno) <- c("Feature", "PC1", "PC2")
    p1 <- ggplot(P_anno[which(P_anno$PC1 != 0),], aes(x = Feature, y = PC1)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    p2 <- ggplot(P_anno[which(P_anno$PC2 != 0),], aes(x = Feature, y = PC2)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    
    PLOT <- (t1 | t2) / (p1 | p2) + plot_annotation(title = "Simulated data, scores and loadings in PCs 1 and 2")
    
    return(list(X,P,Pb,a_ab,PLOT))
    
  }
  else if (score_plot == T){
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      ylim(-3,3) + 
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
      theme_classic() +
      theme(legend.position = "none")
    
    return(list(X,P,Pb,a_ab,t1))
  }
  else{
    return(list(X,P,Pb,a_ab))
  }
  
}


Create_Core_DEV <- function(nreps, meta, irr_spikes = TRUE, struc_resid = FALSE, 
                            a_sigma = c(1.5,0.75,0.6,0.6), b_sigma = c(1.5,1.5), e_sigma = c(1,0.8,0.5),
                            noise_sd = 0.75, EffectSize = c(X_a_ab = 1, time = 0.5, E = 0.5, mu = 1, Struc_E = 1), 
                            SCORE_SEED = 1000, plot = FALSE, score_plot = FALSE, Experiment_responders = 12, ts = 1234){
  
  ########### A_AB
  ##CREATE SCORES FOR PC1 top half
  # set.seed(SCORE_SEED)
  
  #    #     a_ab_PC1 <- c(0.001, 0.75, 0.85, 0.9)
  #    # 
  #    #   ET <- data.matrix(e[c("PC1","PC2","PC3")]) %*% t(Pe_sigma)
  #    # #svd(ET)$u as input alternatively
  #    # 
  #    # #matrix linear regression
  #    # A <- data.matrix(cbind.data.frame(a_ab[c("PC1","PC2")],b[c("PC1","PC2")]))
  #    # B <- solve(t(A)%*%A) %*% t(A) %*% ET
  #    # 
  #    # E <- ET - (A %*% B)  
  #    # 
  # a_ab_PC1 <- rnorm(length(unique(meta$time)))
  # a_ab_PC1 <- abs(a_ab_PC1)
  # a_ab_PC1 <- a_ab_PC1[order(a_ab_PC1, decreasing = FALSE)]
  # 
  # 
  
  
  
  # a_ab_PC1 <- c(0.001, 0.75, 0.85, 0.9)
  # 
  # a_ab_PC1_small <- a_ab_PC1
  # 
  # ##CREATE SCORES FOR PC2 top half
  # set.seed(1)
  # a_ab_PC2 <- Create_Orthogonal(a_ab_PC1_small)
  # a_ab_PC1 %*% a_ab_PC2 == 0
  # a_ab_PC1 <- rep(a_ab_PC1, each = nreps)
  # a_ab_PC2 <- rep(a_ab_PC2, each = nreps)
  # 
  # ##CREATE FULL SCORES
  # a_ab_PC1 <- matrix(c(a_ab_PC1, a_ab_PC1*(-1)))
  #   # a_ab_PC1 <- matrix(c(a_ab_PC1, rep(0,length(a_ab_PC1)))) #
  # 
  # 
  # a_ab_PC2 <- matrix(c(a_ab_PC2, a_ab_PC2*(-1)))
  # a_ab <- cbind.data.frame(meta[,1:2], a_ab_PC1, a_ab_PC2)
  # colnames(a_ab)[c(3,4)] <- c("PC1","PC2")
  # a_ab$PC1 %*% a_ab$PC2
  set.seed(ts)
  targets_a_ab <- abs(rnorm(12,mean = 0.25, sd = 0.25))
  
  
  # set.seed(1)
  set.seed(SCORE_SEED)
  
  a_ab <- as.data.frame(svd(matrix(rnorm(16), nrow = 4))$u)
  
  a_ab_PC1_small <- a_ab[,1]
  
  a_ab <- a_ab[rep(rownames(a_ab), each = nreps),]
  a_ab <- rbind.data.frame(a_ab, a_ab*-1)
  
  #CREATE LOADINGS 
  # a_ab_p1 <- c(rep(0,1988), rep(0.5,4), rep(1,4), rep(1.5,4))
  # a_ab_p1 <- c(rnorm(24,mean = -0.5, sd = 0.25),rep(0,1952), abs(rnorm(8,mean = 0.5, sd = 0.25)),abs(rnorm(16,mean = 1, sd = 0.25)) )#rep(1,8), rep(1,8), rep(1,8)
  
  
  # a_ab_p1 <- c(rnorm(24,mean = -0.5, sd = 0.25),                                #
  #              rep(0,2000 - (24+8+16+(Experiment_responders - 12))),            #
  #              abs(rnorm((Experiment_responders - 12),mean = 1, sd = 0.25)),    #
  #              abs(rnorm(12,mean = 0.5, sd = 0.25)),                            #
  #              abs(rnorm(12,mean = 1, sd = 0.25)))                              #                 rep(1,8), rep(1,8), rep(1,8)
  # 
  
  # set.seed(12)
  a_ab_p1 <- c(rep(0,24),                                                         #rnorm(24,mean = -0.5, sd = 0.25)
               rep(0,2000 - (24+4+12+(Experiment_responders))),
               abs(rnorm((Experiment_responders + 4),mean = 1, sd = 0.25)), 
               targets_a_ab)
  
  print(length(a_ab_p1))
  
  if(irr_spikes == FALSE){
    a_ab_p2 <- c(rep(0,length(a_ab_p1)))
    a_ab_p3 <- c(rep(0,length(a_ab_p1)))
    a_ab_p4 <- c(rep(0,length(a_ab_p1)))
  }
  else{
    a_ab_p2 <- c(rep(0,485),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1485))
    a_ab_p3 <- c(rep(0,185),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1785))
    a_ab_p4 <- c(rep(0,85),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1885))
  }
  P <- cbind(a_ab_p1, a_ab_p2, a_ab_p3, a_ab_p4)
  
  P_sigma <- P %*% diag(a_sigma)
  
  ############# B
  ##CREATE SCORES FOR PC1 TIME
  set.seed(22)
  b_PC1_small <- Create_Orthogonal(a_ab_PC1_small)/400
  b_PC1_small %*% a_ab_PC1_small
  
  ##CREATE SCORES FOR PC2 TIME
  b_PC2_small <- Create_Orthogonal(b_PC1_small)
  b_PC1 <- rep(b_PC1_small, each = nreps, times = 2)
  b_PC2 <- rep(b_PC2_small, each = nreps, times = 2)
  
  b <- cbind.data.frame(meta[,1:2],b_PC1,b_PC2)
  colnames(b)[3:4] <- c("PC1","PC2")
  
  #CREATE LOADINGS TIME
  b_p1 <- c(rep(0.1,1000), rep(0,1000)  )  #make dependent on length of a_ab_p1
  b_p2 <- c(rep(0,1000), rep(0.1,1000))  #make dependent on length of a_ab_p1
  Pb <- cbind(b_p1, b_p2)
  
  Pb_sigma <- Pb %*% diag(b_sigma)
  
  ################## MATRIX MULTIPLICATION OF SCORES AND LOADINGS PER EFFECT MATRIX
  X_a_ab <- EffectSize["X_a_ab"] * data.matrix(a_ab[]) %*% t(P_sigma)
  time <- EffectSize["time"] * data.matrix(b[c("PC1","PC2")]) %*% t(Pb_sigma)
  
  ##CREATE FEATURE MEAN VALUES and ensure they are above 0
  #set.seed here? 
  mu <- EffectSize["mu"] * matrix(rep(rnorm(dim(X_a_ab)[2], mean = mean(X_a_ab)), each = nrow(meta)), nrow = nrow(meta))
  mu <- mu + abs(min(mu))
  X_a_ab_mu <- X_a_ab + mu
  
  
  ####CREATE E
  
  if(struc_resid == TRUE){
    
    #make structured residuals instead of random noise:
    set.seed(127)
    
    # t_l <- svd(matrix(rnorm(2000*24), nrow = 24))$u
    
    t_l <- matrix(rnorm(24*3), nrow = 24)
    
    e_PC1 <- t_l[,1]
    e_PC2 <- t_l[,2]
    e_PC3 <- t_l[,3]
    
    
    
    #create random correlated effect in t_l instead of the random matrix from above
    l_t <- c(1.5,1)
    l_p <- abs(rnorm(16,sd = 1))
    l_x <- l_t %*% t(l_p)
    
    set.seed(1277)
    l_t2 <- Create_Orthogonal(l_t)
    l_p2 <- abs(rnorm(16,sd =0.6))
    l_x2 <- l_t2 %*% t(l_p2)
    
    l_X <- l_x + l_x2
    l_X <- l_X/norm(l_X, type = "F")
    # plot(l_X[1,],l_X[2,])
    
    
    ###make one strong relationship between baits and then some random weaker ones 
    #
    
    
    # e_p1 <- c(rep(0,1976), abs(abs(rnorm(16,sd = 0.2, mean = 0.5))) , rep(0,8)) # rnorm(8, sd = sd(l_X[1,]))l_X[1,] rep(0.5,16)
    # # rtruncnorm(12, a = 0.1, b = 0.8, mean = 0.3)
    # set.seed(128)
    # # e_p2 <- c(rep(0,0), rnorm(2000))  #make dependent on length of a_ab_p1
    # e_p2 <- c(rep(0,1976), abs(rnorm(16,sd = 0.2, mean = 0.5)),rep(0,8) ) #rnorm(8, sd = sd(l_X[2,])) l_X[2,] rep(0.5,16)
    # # rep(0,12)
    # # e_p3 <- c(rep(0,1968),rnorm(32,mean = -1, sd = 0.2)) # rnorm(8,sd =0.1) #baits + spikes
    # # e_p3 <- c(rep(0,1976), abs(rnorm(16,sd = 0.2, mean = 0.5)),rep(0,8) ) # extra source for baits
    # # e_p3 <- c(rep(0,1000), abs(rnorm(16,sd = 0.2, mean = 0.5)), rep(1976,8),
    # # abs(rnorm(16,sd = 0.2, mean = 0.5)),rep(0,8) ) # baits + other
    # e_p3 <- c(rep(0,1976), abs(rnorm(16,sd = 0.2, mean = 0.5)),rep(0,8) ) # extra some baits plus spikes
    # 
    
    set.seed(1288)
    e_p1 <- c(rep(0,1984), abs(abs(rnorm(16,sd = 0.2, mean = 0.5)))) # rnorm(8, sd = sd(l_X[1,]))l_X[1,] rep(0.5,16)
    # rtruncnorm(12, a = 0.1, b = 0.8, mean = 0.3)
    set.seed(128)
    # e_p2 <- c(rep(0,0), rnorm(2000))  #make dependent on length of a_ab_p1
    e_p2 <- c(rep(0,1984), abs(rnorm(16,sd = 0.2, mean = 0.5))) #rnorm(8, sd = sd(l_X[2,])) l_X[2,] rep(0.5,16)
    e_p3 <- c(rep(0,1984), abs(rnorm(16,sd = 0.2, mean = 0.5))) # extra some baits plus spikes
    
    
    Pe <- cbind(e_p1, e_p2, e_p3)
    colnames(Pe) <- c("PC1","PC2", "PC3")
    
    Pe_sigma <- Pe %*% diag(e_sigma)
    e <- cbind.data.frame(meta[,1:2], e_PC1, e_PC2, e_PC3)
    colnames(e)[c(3,4,5)] <- c("PC1","PC2","PC3")
    
    
    #################
    
    ET <- data.matrix(e[c("PC1","PC2","PC3")]) %*% t(Pe_sigma)
    #svd(ET)$u as input alternatively
    
    #matrix linear regression
    A <- data.matrix(cbind.data.frame(a_ab[],b[c("PC1","PC2")]))
    B <- solve(t(A)%*%A) %*% t(A) %*% ET
    
    E <- ET - (A %*% B)  
    
    set.seed(10021)
    Struc_E <- EffectSize["Struc_E"] * E
    # Struc_E <- EffectSize["Struc_E"] * data.matrix(e[c("PC1","PC2","PC3")]) %*% t(Pe_sigma)
    
    dim(Struc_E)
    print(E[,1] %*% time[,1])
    print(time[,1] %*% X_a_ab[,1])
    
    # E <- EffectSize["E"] * data.matrix(e[c("PC1","PC2")]) %*% t(Pe_sigma)
    
    Non_Struc_E <- matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    
    #matrix linear regression
    # A <- data.matrix(cbind.data.frame(X_a_ab, time, Struc_E))
    # B <- solve(t(A)%*%A) %*% t(A) %*% Non_Struc_E 
    # 
    # Non_Struc_E <- ET - (A %*% B)  
    # 
    
    
    E <- Struc_E + (EffectSize["E"] * Non_Struc_E)
    
    
    ############
    
  }
  else{
    
    ## some random noise - i.e. differences between replicates
    set.seed(101)
    E <- EffectSize["E"] * matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    
    
  }
  
  
  ##ADD ALL TOGETHER TO CREATE "CORE"
  X_no.time <- X_a_ab_mu + E
  X <- X_no.time + time
  X <- t(((t(X) + abs(colMins(X)))))
  colnames(X) <- paste0("X_",c(1:dim(X)[2]))
  
  
  Effects <- list("main" = X_a_ab, "cofactor" = time, "residual" = E, "ambient" = Struc_E, "noise" = Non_Struc_E)
  
  
  if(plot == TRUE){
    
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    t2 <- ggplot(a_ab, aes(x = time, y = PC2, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    P_anno <- cbind.data.frame(c(1:dim(P)[1]), P)
    colnames(P_anno) <- c("Feature", "PC1", "PC2")
    p1 <- ggplot(P_anno[which(P_anno$PC1 != 0),], aes(x = Feature, y = PC1)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    p2 <- ggplot(P_anno[which(P_anno$PC2 != 0),], aes(x = Feature, y = PC2)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    
    PLOT <- (t1 | t2) / (p1 | p2) + plot_annotation(title = "Simulated data, scores and loadings in PCs 1 and 2")
    
    return(list(X,P,Pb,a_ab,PLOT, X_a_ab,time, Struc_E ))
    
  }
  else if (score_plot == T){
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      ylim(-3,3) + 
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
      theme_classic() +
      theme(legend.position = "none")
    
    return(list(X,P,Pb,a_ab,t1))
  }
  else{
    return(list(X,P,Pb,a_ab,Effects))  #Effect inputs here?..
  }
  
}

Core2NB <- function(X){
  a <- X
  
  set.seed(41)
  mean <- 0.04282106      #mean(dds2@rowRanges@elementMetadata$dispersion)
  sd <- 0.08265171        #sd(dds2@rowRanges@elementMetadata$dispersion)
  shape <- mean^2/sd^2
  scale <- sd^2/mean
  
  set.seed(42)
  disp <- rgamma(dim(a) [2], shape = shape, scale = scale)
  
  set.seed(6)
  a <- t(a)
  d <- matrix(rnbinom(dim(a)[1]*dim(a)[2], mu = 2^a, size = 1/disp), nrow = length(disp), byrow = FALSE)  #NB counts
  X_nb <- t(d)
  
  colnames(X_nb) <- colnames(X)
  
  return(X_nb)
}


# set.seed(1000)
# X <- Create_Core(3)
# X_nb <- Core2NB(X)
#  

set.seed(1000) 
Simulate <- function(nreps){
  
  Data <- Core2NB(Create_Core(nreps))
  return(Data)
  
}

# Y <- Simulate(5) #not working - need to change meta first

# X2 <- Simulate(3)
# sum(X_nb != X2)  #should be 0 as same seed used


VS_transform <- function(X2, meta){
  
  dds <- DESeqDataSetFromMatrix(countData = t(X2), 
                                colData = meta, 
                                design = ~ ID)
  vsd <- varianceStabilizingTransformation(object = dds, 
                                           blind = TRUE,           
                                           fitType = "parametric")
  
  X_vst <- t(assay(vsd))
  return(X_vst)
  
}





PCA <- function(X, baits, spikes = NULL, ncands, distance_calc = FALSE){
  PCD <- prcomp(X)
  ASCA_cands <- get_ASCA_cands(PCD, baits = baits, distance_calc = distance_calc)
  
  if(distance_calc == T){
    F1_scores <- F1_plot(ASCA_cands, baits = spikes, ncands)
    
  }else{
    F1_scores <- F1_plot(ASCA_cands, baits = baits, ncands)
    
  }
  
  
  return(list(F1_scores, ASCA_cands))
}

# PCA(X_vst, baits, ncands)



ASCA <- function(X, m, Baits, spikes, ncands, distance_calc = FALSE, Return_Model = FALSE){
  
  # X <- X[,-which(colnames(X) %in% Baits)]
  
  m <- m
  m[,1:2] <- lapply(m[,1:2],factor)
  res_ASCAplus <- ASCA_decompose(d = m[,1:2], x = X,
                                 f = "growth_condition + time + growth_condition:time")
  minT <-  res_ASCAplus$decomposition$growth_condition + res_ASCAplus$decomposition$`growth_condition:time`
  PCD <- prcomp(minT)
  print(Baits)
  ASCA_cands <- get_ASCA_cands(PCD, meta, distance_calc= distance_calc, baits = NULL, 
                               spikes = NULL, ret_candN = nrow(PCD$rotation))
  
  #PCD = PCD,  meta = m, baits = baits, distance_calc = distance_calc
  
  # PCD, meta, distance_calc= FALSE, baits = NULL, spikes = NULL, ret_candN = nrow(PCD$rotation)
  
  resids <- res_ASCAplus$residuals
  
  ASCA_cands <- ASCA_cands[-which(rownames(ASCA_cands) %in% baits),]
  
  if(distance_calc == TRUE){
    # F1_scores <- F1_plot(ASCA_cands, baits = spikes, ncands)
    RP <- prod(which(rownames(ASCA_cands) %in% spikes))^(1/length(spikes))
    
    
  }else{
    # F1_scores <- F1_plot(ASCA_cands, baits = Baits, ncands)
    RP <- prod(which(rownames(ASCA_cands) %in% spikes))^(1/length(spikes))
    
    
  }
  
  if(Return_Model == TRUE){
    
    return(list(RP, PCD))
    
  }else{return(list(RP, ASCA_cands, resids, res_ASCAplus))}
  
  
  
}


ASCA_notime <- function(X, m, baits, spikes, ncands, distance_calc = FALSE, Return_Model = FALSE){
  
  m <- m
  m <- apply(m,2,factor)
  res_ASCAplus <- ASCA_decompose(d = m, x = X,
                                 f = "growth_condition")
  minT <-  res_ASCAplus$decomposition$growth_condition 
  PCD <- prcomp(minT)
  print(baits)
  ASCA_cands <- get_ASCA_cands(PCD,  baits = baits, distance_calc = distance_calc)
  
  
  if(distance_calc == TRUE){
    F1_scores <- F1_plot(ASCA_cands, baits = spikes, ncands)
    
  }else{
    F1_scores <- F1_plot(ASCA_cands, baits = baits, ncands)
    
  }
  
  if(Return_Model == TRUE){
    
    return(list(F1_scores, PCD))
    
  }else{return(list(F1_scores, ASCA_cands))}
  
  
  
}


# ASCA(X_vst, meta, baits, ncands)


#to calculate VE per factor
ASCA_calc <- function(X,Y){
  
  # levels <- apply(Y, 2, function(x) y <- as.numeric(as.factor(x)))
  
  M_ASCA_rg_scaled <- ASCA.Calculate(X, levels, equation.elements = "1,2,12")
  
  minT <- M_ASCA_rg_scaled$'1'$means.matrix + M_ASCA_rg_scaled$'12'$means.matrix
  
  rownames(minT) <- rownames(X)
  colnames(minT) <- colnames(X)
  PCD <- prcomp(minT)
  
  ##project '1' + '12' + E onto PCA space
  newdata <- M_ASCA_rg_scaled$'1'$means.matrix + M_ASCA_rg_scaled$'12'$means.matrix + M_ASCA_rg_scaled$remainder
  colnames(newdata) <- colnames(minT)
  
  #this predict method does matrix multiplication of the input (new data) and the rotation (loading) matrix 
  PCDE <- predict(PCD, newdata)
  PCDE <- cbind.data.frame(Y$growth_condition, Y$time, PCDE)
  colnames(PCDE)[1:2] <- c("day","growth_condition")
  
  return(list(PCDE,M_ASCA_rg_scaled,PCD))
  
  
  
}


PLS <- function(X, meta, baits, ncands){
  
  y <- meta$growth_condition
  y <- gsub("-","n",y)
  
  nCore=12   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='PLS'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  
  PLSModel <- invisible(MUVR(X=X, y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method))
  
  PLS_cands <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
  
  if(nrow(PLS_cands) < ncol(X)){
    extra <- matrix(rep(ncol(X),3 * (ncol(X) - nrow(PLS_cands))), ncol = 3)
    PLS_cands <- rbind(PLS_cands, extra)
    
  }
  
  
  stopCluster(cl)
  F1_scores <- F1_plot(PLS_cands, baits, ncands, TITLE = ": PLS with MUVR")
  
  return(list(F1_scores, PLS_cands))
  
  
}

RF <- function(X, meta, baits, ncands){
  
  y <- meta$growth_condition
  y <- gsub("-","n",y)
  
  nCore=12  # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='RF'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  methParam <- customParams(NZV = TRUE)
  
  RFModel = invisible(MUVR(X=X, y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method, methParam = methParam))
  
  RF_cands <- RFModel$VIP[order(RFModel$VIP[,2]),]
  
  if(nrow(RF_cands) < ncol(X)){
    extra <- matrix(rep(ncol(X),3 * (ncol(X) - nrow(RF_cands))), ncol = 3)
    RF_cands <- rbind(RF_cands, extra)
    
  }
  
  stopCluster(cl)
  F1_scores <- F1_plot(RF_cands, baits, ncands, TITLE = ": RF with MUVR")
  
  return(list(F1_scores, RF_cands))
  
  
}

# RF(X_vst, meta, baits, ncands)

Ranked_Coexp <- function(X, baits, spikes, ncands){
  
  ranked <- ranked_coexp(baits, X)
  
  F1_scores <- F1_plot(ranked, spikes, ncands, TITLE = ": high effect baits")
  
  return(list(F1_scores, ranked))
  
}


Ranked_Coexp_HIGH <- function(X, baits, ncands){
  baits_small <- baits[c(9:12)]
  
  ranked <- ranked_coexp(baits_small, X)
  
  F1_scores <- F1_plot(ranked, baits[-which(baits %in% baits_small)], ncands, TITLE = ": high effect baits")
  
  return(list(F1_scores, ranked))
  
}

# Ranked_Coexp_HIGH(X_vst, baits, ncands)


Ranked_Coexp_LOW <- function(X, baits, ncands){
  baits_small <- baits[c(1:4)]
  
  ranked <- ranked_coexp(baits_small, X)
  
  F1_scores <- F1_plot(ranked, baits[-which(baits %in% baits_small)], ncands, TITLE = ": low effect baits")
  
  return(list(F1_scores, ranked))
  
}

# Ranked_Coexp_LOW(X_vst, baits, ncands)




ASCA_calc <- function(x,y){
  
  levels <- apply(x, 2, function(x) y <- as.numeric(as.factor(x)))
  
  M_ASCA_rg_scaled <- ASCA.Calculate(y, levels, equation.elements = "1,2,12")
  ##rgs
  
  minT <- M_ASCA_rg_scaled$'1'$means.matrix + M_ASCA_rg_scaled$'12'$means.matrix
  
  rownames(minT) <- rownames(y)
  colnames(minT) <- colnames(y)
  PCD <- prcomp(minT)
  ##project '1' + '12' + E onto PCA space
  newdata <- M_ASCA_rg_scaled$'1'$means.matrix + M_ASCA_rg_scaled$'12'$means.matrix + M_ASCA_rg_scaled$remainder
  colnames(newdata) <- colnames(minT)
  
  #this predict method does matrix multiplication of the input (new data) and the rotation (loading) matrix 
  PCDE <- predict(PCD, newdata)
  x <- as.data.frame(x)
  PCDE <- cbind.data.frame(x$time, x$growth_condition, PCDE)
  colnames(PCDE)[1:2] <- c("time","growth_condition")
  
  return(list(PCDE,M_ASCA_rg_scaled,PCD))
  
}

##distribution comparisons materials 
##without normal noise test... - better but not ideal
##noise b slope per feature

ASCA_pred <- function(d,x){
  
  res_ASCAplus <- ASCA_decompose(d = d, x = x,
                                 f = "growth_condition + time + growth_condition:time",
                                 glm_par = list(family = quasipoisson(link = "log")))
  
  #find the terms of interest here
  terms <- names(res_ASCAplus$decomposition)
  #> terms
  #[1] "mu"    "cond"  "error"
  minT <-  res_ASCAplus$decomposition[,2,] + res_ASCAplus$decomposition[,4,]
  PCD <- prcomp(minT)
  
  newdata <- res_ASCAplus$decomposition[,2,] + res_ASCAplus$decomposition[,4,] + res_ASCAplus$decomposition[,5,]
  
  #ASCA+ decomposition "error" is not as expected
  #needs to be calculate by X-estimateX (i.e. sum of all effect matrices)
  
  #this predict method does matrix multiplication of the input (new data) and the rotation (loading) matrix
  #projects the data into the PCA space
  PCDE <- predict(PCD, newdata)
  PCDE <- cbind(levels, PCDE)
  
  return(list(PCDE,res_ASCAplus,PCD))
  
}



get_ASCA_cands <- function(PCD, meta, distance_calc= FALSE, baits = NULL, spikes = NULL, ret_candN = nrow(PCD$rotation)){
  #############
  
  if(distance_calc== TRUE){
    
    cands <- ranked_dist(baits,PCD)
    
  }else{
    absload <- abs(data.matrix(PCD$rotation[,1:2])) %*% diag(summary(PCD)$importance[2,1:2])
    combscore <- rowSums(absload[,1:2])
    
    
    orderedload <- cbind(combscore, absload)
    cands <- as.data.frame(orderedload[order(orderedload[,1], decreasing = TRUE),])
    colnames(cands) <- c("VIP", "PC1","PC2")
    
    cands <- round(cands[1:ret_candN,], 4)
    
    
  }
  
  
  
  
  
  ##############
  #calculate loadings * variance explained
  #ASCA_cands <- PCD$rotation[order(PCD$rotation[,1]),]
  #
  #load <- PCD$rotation[,1:3] %*% diag(summary(PCD)$importance[2,1:3])
  #
  #this makes no sense as the different directions can cancel out the overall measure
  #combscore_dir <- rowSums(load[,1:3])
  #
  #
  #ol <- cbind(combscore_dir, load)
  #
  #makes no sense to order like this as the direction of importance needs to be factored in
  #ol <- as.data.frame(ol[order(ol[,1], decreasing = TRUE),])
  
  #get "eigen position" of each condition from PCA on ASCA effect matrix
  #in this case just the scores calculated from "cond" effect matrix
  
  # head(PCD$x[,1:5])
  # head(ASCA_res_concat[[1]][,1:5])
  
  # range(PCD$rotation[,1])
  
  
  
  #find features that are most indicative of certain conditions (and those that are more indicative of the difference)
  
  
  
  #also try the g_b ranking - see loadings_plot_notes
  
  
  
  return(cands)
  
  
}

Create_Orthogonal_bin_fact <- function(x){
  #not finished
  a_ab_PC2_incomplete <- rnorm(length(x) - 1)
  a_ab_PC2_incomplete_last <- - (x[1:(length(x)-1)] %*% a_ab_PC2_incomplete)/x[length(x)]
  a_ab_PC2 <- c(a_ab_PC2_incomplete, a_ab_PC2_incomplete_last)
  
  return(a_ab_PC2)
  
}

Create_Orthogonal <- function(x){
  
  a_ab_PC2_incomplete <- rnorm(length(x) - 1)
  a_ab_PC2_incomplete_last <- - (x[1:(length(x)-1)] %*% a_ab_PC2_incomplete)/x[length(x)]
  a_ab_PC2 <- c(a_ab_PC2_incomplete, a_ab_PC2_incomplete_last)
  
  return(a_ab_PC2)
  
}



M_SD <- function(data, meta, ...){
  samples2 <- data
  rownames(samples2) <- paste0(meta$growth_condition, "_", meta$time)
  samples_melt <- reshape2::melt(samples2)
  samples_melt$Var2 <- as.character(samples_melt$Var2)
  
  samples_melt$FULL_ID <- paste(samples_melt$Var1, samples_melt$Var2, sep = "_")
  colnames(samples_melt)[2:3] <- c("Feature","Expression")
  
  #FULL_ID contains info on metabolite and sample condition
  #data summarised by this column
  data_summary <- samples_melt %>% group_by(FULL_ID) %>%
    summarise(mean=mean(Expression), sd=sd(Expression))
  
  data_summary <- cbind.data.frame(data_summary, do.call(rbind, strsplit(data_summary$FULL_ID, split = "_")))
  colnames(data_summary)[4:6] <- c("growth_condition", "time","Feature")
  
  
  #might need to change aes depending on experiment
  plot_normal <- ggplot(data_summary, aes(x = log2(mean), y = log2(sd), shape = time, colour = growth_condition)) +
    geom_point() +
    theme_bw() +
    theme(...)
  
  return(list(data_summary, samples_melt, plot_normal))
  
}


Power_Transform <- function(data,meta){
  
  m_sd <- M_SD(data, meta, legend.position = "none")
  
  data_summary <- m_sd[[1]]
  #calculate slope of mean sd relationship
  data_summary$log2_sd <- log2(data_summary$sd)
  data_summary$log2_mean <- log2(data_summary$mean)
  
  data_summary[is.na(data_summary) | data_summary =="Inf" | data_summary =="-Inf"] <- NA
  
  model <- lm(log2_sd ~ log2_mean, data = data_summary)
  
  ##transform the data based on regression slope
  samples_melt <- m_sd[[2]]
  samples_melt$transform <- samples_melt$Expression^(1-model$coefficients[2])
  
  #samples_melt$transform[is.na(samples_melt$transform) | samples_melt$transform =="Inf" | samples_melt$transform =="-Inf"] <- NA
  #now replot the heteroscedasticity
  data_summary <- samples_melt %>% group_by(FULL_ID) %>% # make a dataframe to save the values
    summarise(mean=mean(transform, na.rm = TRUE), sd=sd(transform, na.rm = TRUE)) 
  
  data_summary <- cbind.data.frame(data_summary, do.call(rbind, strsplit(data_summary$FULL_ID, split = "_")))
  colnames(data_summary)[4:6] <- c("growth_condition", "time","Feature")
  
  #data_summary[is.na(data_summary) | data_summary =="Inf" | data_summary =="-Inf"] <- NA
  
  plot1 <- m_sd[[3]]
  plot2 <- ggplot(data_summary, aes(x = log2(mean), y = log2(sd), shape = time, colour = growth_condition)) +
    geom_point() + 
    theme_bw() 
  
  return(list(samples_melt, plot1, plot2, model))
}


Power_Transform_FW <- function(data,meta){
  
  m_sd <- M_SD(data, meta, legend.position = "none")
  
  data_summary <- m_sd[[1]]
  #calculate slope of mean sd relationship
  data_summary$log2_sd <- log2(data_summary$sd)
  data_summary$log2_mean <- log2(data_summary$mean)
  
  data_summary[is.na(data_summary) | data_summary =="Inf" | data_summary =="-Inf"] <- NA
  
  model <- lm(log2_sd ~ log2_mean, data = data_summary)
  
  ##transform the data based on regression slope
  samples_melt <- m_sd[[2]]
  samples_melt$transform <- samples_melt$Expression^(1-model$coefficients[2])
  
  #samples_melt$transform[is.na(samples_melt$transform) | samples_melt$transform =="Inf" | samples_melt$transform =="-Inf"] <- NA
  #now replot the heteroscedasticity
  data_summary <- samples_melt %>% group_by(FULL_ID) %>% # make a dataframe to save the values
    summarise(mean=mean(transform, na.rm = TRUE), sd=sd(transform, na.rm = TRUE)) 
  
  data_summary <- cbind.data.frame(data_summary, do.call(rbind, strsplit(data_summary$FULL_ID, split = "_")))
  colnames(data_summary)[4:6] <- c("growth_condition", "time","Feature")
  
  #data_summary[is.na(data_summary) | data_summary =="Inf" | data_summary =="-Inf"] <- NA
  
  plot1 <- m_sd[[3]]
  plot2 <- ggplot(data_summary, aes(x = log2(mean), y = log2(sd), shape = time, colour = growth_condition)) +
    geom_point() + 
    theme_bw() 
  
  return(list(samples_melt, plot1, plot2, model))
}



SL_plot <- function(PCD,meta){
  
  a_c.PCD <- PCD
  meta <- as.data.frame(meta)
  a_c.PCD$x <- cbind.data.frame(time = meta$time, growth_condition = meta$growth_condition, a_c.PCD$x)
  
  a_c.PCD$rotation <- cbind.data.frame(paste0("X_",c(1:dim(X)[2])), a_c.PCD$rotation)
  colnames(a_c.PCD$rotation)[1] <- "Feature"
  a_c.PCD$rotation$Feature <- factor(a_c.PCD$rotation$Feature, levels = a_c.PCD$rotation$Feature)
  
  ## *Recreated* simulated data, scores and loadings in PCs 1 and 2
  
  t1_rec <- ggplot(a_c.PCD$x, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
    geom_line() +
    theme_classic() +
    theme(legend.position = "none")
  
  t2_rec <- ggplot(a_c.PCD$x, aes(x = time, y = PC2, group=growth_condition, colour=growth_condition)) +
    geom_line() +
    theme_classic() +
    theme(legend.position = "none")
  
  ###
  
  p1_rec <- ggplot(a_c.PCD$rotation[which(abs(a_c.PCD$rotation$PC1) > 0.05),], aes(x = Feature, y = PC1)) +     #CHANGE to include candidate selection here
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_blank())
  
  p2_rec <- ggplot(a_c.PCD$rotation[which(abs(a_c.PCD$rotation$PC2) > 0.05),], aes(x = Feature, y = PC2)) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_blank())
  
  
  plot <- (t1_rec | t2_rec) / (p1_rec | p2_rec)
  return(plot)
  
}





scoreplot <- function(df, meta, ve, ...){
  
  #df <- as.data.frame(df)
  df <- cbind.data.frame(meta,df)
  df$time <- factor(df$time, levels = c(1:4))
  df$growth_condition <- factor(df$growth_condition)
  
  
  ggplot(df, ...)+
    geom_point() +
    xlab(paste0("PC1 ",ve[1])) +
    ylab(paste0("PC2 ",ve[2])) + 
    theme_bw()
  
}



loadingplot <- function(df, meta, baits = ref[which(ref$Baits == "Pathway")], ve, ref, PCs = c(1,2), ...){
  #df is PCD$rotation
  
  #rematch REF, ref, df
  df <- as.data.frame(df)
  df <- cbind.data.frame(ref,df)
  
  # df$Feature <- rownames(df)
  df$Feature[-which(df$Feature %in% baits)] <- "Other"
  df$Feature[which(df$Feature %in% baits)] <- "Pathway"
  
  df$Feature <- factor(df$Feature, levels = c("Other","Pathway"))
  
  #pdf(paste0(Sys.Date(),"_feature_loadings_project_full2.pdf"))
  ggplot(df, ...) + 
    geom_point(alpha = 0.1) +
    xlab(paste0("PC", PCs[1]," ",ve[1])) +
    ylab(paste0("PC", PCs[2]," ",ve[2])) + 
    #geom_label_repel(data = df[which(rownames(df) %in% baits),], aes(label = rownames(df[which(rownames(df) %in% baits),])), nudge_y = 0.003, nudge_x = 0.001) +
    geom_point(data = df[which(df$feature %in% baits),]) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme_bw() 
  
}


####only works with baits for PC1 of a+ab - generalise 
heatmap <- function(data, meta){
  
  
  meta <- cbind.data.frame(meta, rep(c(1:3), times = 8))
  colnames(meta)[5] <- "reps"
  
  data <- t(data)
  
  colnames(data) <- paste(meta$growth_condition, meta$time, meta$reps)
  
  rownames(meta) <- colnames(data)
  
  #create annotation for clarity
  anno <- cbind(meta$growth_condition, as.character(meta$time))
  anno <- as.data.frame(anno)
  rownames(anno) <- rownames(meta)
  colnames(anno) <- c("growth_condition","time")
  #This needs to be a factor for the heatmap call
  anno$time <- factor(anno$time, level = c("1","2","3","4"))
  
  annor <- as.data.frame(c(rep(0.5,4), rep(1,4), rep(1.5,4)))
  colnames(annor) <- "Effect_Level"
  annor$Effect_Level <- factor(annor$Effect_Level)
  rownames(annor) <- baits
  
  heat <- pheatmap::pheatmap(data.matrix(data), 
                             cluster_rows = TRUE,
                             cluster_cols = FALSE, 
                             show_rownames = TRUE, 
                             show_colnames = FALSE,
                             annotation_col = anno,
                             annotation_row = annor,
                             main = "Coexpressed Features",silent = TRUE)
  
  heat
  
}





qqplot <- function(data, meta, ...){
  
  df <- cbind.data.frame(meta,data)
  
  df$time <- factor(df$time, levels = c(1:4))
  df$growth_condition <- factor(df$growth_condition)
  
  
  ggplot(df, ...) + 
    geom_point() +
    theme_bw()
  
  
}

narm <- function(df, keep_col1 = FALSE){
  df[which(is.na(df), arr.ind = TRUE)] <- 0
  if(keep_col1 == FALSE){
    df <- df[,-1]
  }else{df <- df}
}


corclean <- function(cor_mat){
  
  ####get rid of upper triangle of correlation matrix 
  #this line is problematic for downstream analysis
  #cor_mat[upper.tri(cor_mat)] <- 999
  
  melt_adj <- reshape2::melt(cor_mat)
  
  #melt_adj <- melt_adj[-which(melt_adj$value == 999),]
  melt_adj <- melt_adj[-which(melt_adj$Var1 == melt_adj$Var2),]
  
  return(melt_adj)
  
  
}

com_filt <- function(res,common, dist_metric = "cor"){
  res <- res[which(res$Feature %in% common),]
  
  if(dist_metric == "dist"){
    res <- cbind(res, c(dim(res)[1]:1))
    
    
  }else{
    res <- cbind(res, c(1:dim(res)[1]))
    
    
  }
  res <- res[order(res$Feature),]
  colnames(res)[4] <- "Rank"
  
  return(res)
}


#take baits and return top n (user defined) coexpression candidates
#make function that takes each bait, and it's top coexpression candidates
baiter <- function(x, melt_adj){
  
  if (x %in% melt_adj$Var1 == FALSE){
    # print(paste0(x," doesn't pass filters"))
    CoExpCands <- data.frame(matrix(c(x,"doesn't pass filters", 0), ncol = 3))
  } else{
    
    CoExpCands <- melt_adj[which(melt_adj$Var1 == x),]
    CoExpCands <- CoExpCands[order(CoExpCands$value, decreasing= TRUE),]
    
    
    #write.table(CoExpCands, paste0(x,"_Top_CoExpCands.txt"), row.names = FALSE, quote = FALSE)
    
    
    
  }
  rownames(CoExpCands) <- NULL
  colnames(CoExpCands) <- c("Bait","Feature","PCC")
  
  return(CoExpCands)
  
}


ranked_coexp <- function(baits, data){
  
  adj_matrix <- cor(data)
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



ranked_dist <- function(baits, PCD, ncomp = 5){
  
  sum <- as.data.frame(summary(PCD)$importance)
  data <- PCD$rotation %*% diag(sum[2,])
  
  data <- data[,c(1:ncomp)]
  
  adj_matrix <- as.matrix(dist(data))
  melt_adj <- corclean(adj_matrix)
  rm(adj_matrix)
  
  i <- NULL
  RES <- list()
  for(i in 1:length(baits)){
    RES[[i]] <- baiter(baits[i], melt_adj)
  }
  
  rm(melt_adj)
  names(RES) <- baits
  RES <- RES[c(which(do.call(rbind, lapply(RES,dim))[,1] > 1))]
  l <- lapply(RES, "[[",2)
  l <- lapply(l, as.character)
  
  common <- Reduce(intersect, l)
  rres <- lapply(RES, com_filt, dist_metric = "dist", common = common)
  
  
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




F1 <- function(ASCA_cands, baits, ncands){
  
  actual <- c(rownames(ASCA_cands) %in% baits)
  actual[which(actual == FALSE)] <- 0
  actual[which(actual == TRUE)] <- 1
  actual <- factor(actual)
  
  
  pred <- c(rep(1,ncands),rep(0,nrow(ASCA_cands) - ncands))
  pred <- factor(pred)
  
  CM <- caret::confusionMatrix(pred,actual, mode = "everything", positive = "1")
  
  return(CM$byClass["F1"])
  
}



F1_ALL <- function(ASCA_cands, baits){
  
  baits <- baits[which(baits %in% rownames(ASCA_cands))]
  
  actual <- c(rownames(ASCA_cands) %in% baits)
  actual[which(actual == FALSE)] <- 0
  actual[which(actual == TRUE)] <- 1
  actual <- factor(actual)
  
  ncands <- which(rownames(ASCA_cands) %in% baits)[length(baits)]
  
  pred <- c(rep(1,ncands),rep(0,nrow(ASCA_cands) - ncands))
  pred <- factor(pred)
  
  CM <- confusionMatrix(pred,actual, mode = "everything", positive = "1")
  
  return(CM$byClass["F1"])
  
}




F1_plot <- function(ASCA_cands, baits, n, TITLE = NULL){
  
  RES <- c()
  i <- NULL
  
  for(i in 1:n){
    
    f1 <- F1(ASCA_cands, baits, i)
    RES <- c(RES,f1)
    
  }
  
  
  
  df <- cbind(c(1:n), RES)
  colnames(df) <- c("N_Selected_Features", "F1_score")
  
  grob <- grobTree(textGrob(paste0("Max F1 score = ",round(max(RES, na.rm = TRUE), digits = 4),"\n N_selected = ",which(RES == max(RES, na.rm = TRUE)))))
  
  plot <- ggplot(df) + 
    geom_point(aes(x = N_Selected_Features, y = F1_score)) +
    geom_vline(xintercept = which(RES == max(RES, na.rm = TRUE))) +
    ggtitle(paste0("F1 scores",TITLE)) +
    annotation_custom(grob) + 
    theme_bw()
  
  
  
  return(list(RES,plot))
  
}


F1_plot_no_plot <- function(ASCA_cands, baits, n, TITLE = NULL){
  
  RES <- c()
  i <- NULL
  
  for(i in 1:n){
    
    f1 <- F1(ASCA_cands, baits, i)
    RES <- c(RES,f1)
    
  }
  
  
  return(RES)
  
}





PLSr_HIGH <- function(X, meta, baits, ncands){
  
  baits_small <- baits[c(9:12)]  #change this to make one single function, also for coexp and randomforest
  baits <- baits[-which(baits %in% baits_small)]
  
  nCore=1   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='PLS'            # Selected core modelling algorithm
  cl=makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits_small)){
    y <- X[,which(colnames(X) %in% baits_small[i])]
    
    
    PLSModel = invisible(MUVR::MUVR(X=X[,-which(colnames(X) %in% baits_small)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method))#, methParam = methParam)
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    
    
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    VIPS[[i]] <- VIPS[[i]][,c("Feature","Rank")]
    
    
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, dplyr::left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  
  F1_scores <- F1_plot(RES_all, baits, ncands, TITLE = ": PLS-R High")
  
  return(list(F1_scores,RES_all))
  
  
}

PLSr_LOW <- function(X, meta, baits, ncands){
  
  baits_small <- baits[c(1:4)]  #change this to make one single function, also for coexp and randomforest
  baits <- baits[-which(baits %in% baits_small)]
  
  nCore=12   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='PLS'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits_small)){
    y <- X[,which(colnames(X) %in% baits_small[i])]
    
    
    PLSModel = invisible(MUVR(X=X[,-which(colnames(X) %in% baits_small)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method))#, methParam = methParam)
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    
    
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    VIPS[[i]] <- VIPS[[i]][,c("Feature","Rank")]
    
    
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  F1_scores <- F1_plot(RES_all, baits, ncands, TITLE = ": PLS-R Low")
  
  
  return(list(F1_scores,RES_all))
  
  
}


PLSr <- function(X, meta, baits, spikes, ncands){
  
  # baits_small <- baits[c(9:12)]  #change this to make one single function, also for coexp and randomforest
  # baits <- baits[-which(baits %in% baits_small)]
  # 
  nCore=1   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='PLS'            # Selected core modelling algorithm
  cl=makeCluster(nCore)
  doParallel::registerDoParallel(cl)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits)){
    y <- X[,which(colnames(X) %in% baits[i])]
    
    
    PLSModel = invisible(MUVR::MUVR(X=X[,-which(colnames(X) %in% baits)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method))#, methParam = methParam)
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    
    
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    VIPS[[i]] <- VIPS[[i]][,c("Feature","Rank")]
    
    
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, dplyr::left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  
  # F1_scores <- F1_plot(RES_all, spikes, ncands, TITLE = ": PLS-R")
  RP <- prod(which(rownames(RES_all) %in% spikes))^(1/length(spikes))
  
  
  return(list(RP,RES_all))
  
  
}



RFr <- function(X, meta, baits, spikes, ncands){
  
  # baits_small <- baits[c(9:12)]  #change this to make one single function, also for coexp and randomforest
  # baits <- baits[-which(baits %in% baits_small)]
  # 
  nCore=12   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='RF'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  methParam <- customParams(NZV = TRUE)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits)){
    y <- X[,which(colnames(X) %in% baits[i])]
    
    
    PLSModel = invisible(MUVR(X=X[,-which(colnames(X) %in% baits)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method, methParam = methParam))
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    
    
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    VIPS[[i]] <- VIPS[[i]][,c("Feature","Rank")]
    
    
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  F1_scores <- F1_plot(RES_all, spikes, ncands, TITLE = ": RF-R")
  
  return(list(F1_scores,RES_all))
  
  
}


RFr_HIGH <- function(X, meta, baits, ncands){
  
  baits_small <- baits[c(9:12)]  #change this to make one single function, also for coexp and randomforest
  baits <- baits[-which(baits %in% baits_small)]
  
  nCore=12   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='RF'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  methParam <- customParams(NZV = TRUE)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits_small)){
    y <- X[,which(colnames(X) %in% baits_small[i])]
    
    
    PLSModel = invisible(MUVR(X=X[,-which(colnames(X) %in% baits_small)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method, methParam = methParam))
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    
    
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    VIPS[[i]] <- VIPS[[i]][,c("Feature","Rank")]
    
    
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  F1_scores <- F1_plot(RES_all, baits, ncands, TITLE = ": RF-R High")
  
  return(list(F1_scores,RES_all))
  
  
}

RFr_LOW <- function(X, meta, baits, ncands){
  
  baits_small <- baits[c(1:4)]  #change this to make one single function, also for coexp and randomforest
  baits <- baits[-which(baits %in% baits_small)]
  
  nCore=12   # Number of processor threads to use
  nRep=nCore              # Number of MUVR repetitions
  nOuter=3                # Number of outer cross-validation segments
  varRatio=0.8            # Proportion of variables kept per iteration 
  method='RF'            # Selected core modelling algorithm
  cl=makeCluster(nCore)   
  registerDoParallel(cl)
  methParam <- customParams(NZV = TRUE)
  
  i <- NULL
  RES <- list()
  VIPS <- list()
  for(i in 1:length(baits_small)){
    y <- X[,which(colnames(X) %in% baits_small[i])]
    
    
    PLSModel = invisible(MUVR(X=X[,-which(colnames(X) %in% baits_small)], y , nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method, methParam = methParam))
    RES[[i]] <- PLSModel
    
    VIPS[[i]] <- PLSModel$VIP[order(PLSModel$VIP[,2]),]
    VIPS[[i]] <- cbind.data.frame(Feature = rownames(VIPS[[i]]),VIPS[[i]],Rank = c(1:nrow(VIPS[[i]])))
    
  }
  stopCluster(cl)
  
  RES_all <- purrr::reduce(VIPS, left_join, by = "Feature")
  RES_all$comb_rank <- rowMeans(RES_all[,grep("Rank",colnames(RES_all))])
  RES_all <- RES_all[order(RES_all$comb_rank),]
  rownames(RES_all) <- RES_all$Feature
  
  if(nrow(RES_all) < ncol(X)){
    extra <- matrix(rep(ncol(X),ncol(RES_all) * (ncol(X) - nrow(RES_all))), ncol = ncol(RES_all))
    colnames(extra) <- colnames(RES_all)
    RES_all <- rbind(RES_all, extra)
    
  }
  
  F1_scores <- F1_plot(RES_all, baits, ncands, TITLE = ": RF-R Low")
  
  return(list(F1_scores,RES_all))
  
  
}


sPLSr <- function(X, meta, baits, spikes, ncands, TITLE = NULL){
  
  #mean centre X
  # X <- scale(X)
  
  splsR <- spls(X[,-which(colnames(X) %in% baits)],X[,which(colnames(X) %in% baits)], ncomp = 2)
  sPLS_cands <- cbind.data.frame(abs(splsR$loadings$X) %*% splsR$prop_expl_var$Y, splsR$loadings$X)
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
  
  # F1_scores <- F1_plot(sPLS_cands, spikes, ncands, TITLE = paste0(": ",TITLE))
  # F1_scores <- F1_plot(RES_all, spikes, ncands, TITLE = ": PLS-R")
  RP <- prod(which(rownames(sPLS_cands) %in% spikes))^(1/length(spikes))
  
  return(list(RP,sPLS_cands,splsR))
  
}


MASCARA4_test <- function(resids,ref, baits, spikes, ncands){
  
  Fref <- ref
  
  spls_res <- spls(resids[,-which(colnames(resids) %in% baits)], 
                   resids[,which(colnames(resids) %in% baits)], all.outputs = TRUE)
  
  sPLS_cands <- cbind.data.frame(abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X, spls_res$loadings$X) #  * spls_res$prop_expl_var$X
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
  
  importance <- paste0(round(spls_res$prop_expl_var$X * 100, digits = 2), "%")
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_cands),Fref2$Feature),]
  
  #####
  
  
  #calculating direction of baits to determine direction of candidate selection spls loading space
  
  
  angles <- rbind.data.frame(spls_res$loadings$Y, colMeans(spls_res$loadings$Y))
  
  # x <- angles[4,]
  y <- angles[nrow(angles),]
  
  #adapted from https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r
  angle <- function(x,y){
    
    x <- c(as.matrix(x))
    y <- c(as.matrix(y))
    
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    return(as.numeric(theta))
  }
  
  angles_loadings <- apply(sPLS_cands[,-1],1,angle, y = y)
  
  
  #define this better
  sPLS_cands <- cbind.data.frame(sPLS_cands, "angle" = angles_loadings, "score" = 3 * sPLS_cands$`abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X` * (1-angles_loadings))
  
  sPLS_cands <- sPLS_cands[order(sPLS_cands$score, decreasing = TRUE),]
  
  #########
  
  avg_point <- angles[nrow(angles),]/nrow(angles)
  
  # rot.mat <- matrix(c(0,-1,1,0),nrow = 2, byrow = TRUE)  #90 degree clockwise rotation
  # 
  # rotated_point <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat)
  # 
  # ang <- angle(avg_point,c(1,0))
  # 
  # rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  # rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  #######
  

  
  #rotated loadings
  ang <- angle(avg_point,c(1,0))
  
  rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  
  rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  rotate <- function(x,y){
    x <- c(as.matrix(x))
    x %*% y
  }
  
  loadings_r_baits <- t(apply(sPLS_cands[,c("comp1","comp2")],1,rotate, y = rot.mat2))
  colnames(loadings_r_baits) <- c("comp1","comp2")
  
  
  ######
  
  sPLS_bait_cands <- loadings_r_baits[order(loadings_r_baits[,1], decreasing = TRUE),]
  
  
  # #####
  # 
  # #loop over rotation matrices from -pi/4 -> pi/4
  # #resolution of 100 (?)
  # 
  # range.ang <- seq(from = -pi/6, to = pi/6, length.out = 60)
  # 
  # ###
  # 
  # i <- NULL
  # L_R <- data.frame()
  # LOADINGS_R <- list()
  # 
  # for(i in 1:length(range.ang)){
  #   
  #   
  #   ang <- range.ang[i]
  #   
  #   rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  #   loadings_r <- t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat))
  #   colnames(loadings_r) <- c("comp1","comp2")
  #   L_R <- rbind(L_R,loadings_r)
  #   LOADINGS_R[[i]] <- loadings_r
  #   
  #   
  # }
  # 
  # # colnames(LOADINGS_R) <- c("comp1","comp2")
  # L_R$rad <- rep(range.ang, each = nrow(sPLS_cands))
  # L_R$frame <- rep(c(1:length(range.ang)), each = nrow(sPLS_cands))
  # 
  # 
  # 
  # #get sd of original comp1 distribution
  # sigma <- sd(sPLS_cands$comp1)
  # 
  # big.T.sig <- function(x,y){
  #   big <- length(which(x > 2.5*y))
  #   return(big)
  # }
  # 
  # i <- NULL
  # res <- c()
  # for(i in 1:length(LOADINGS_R)){
  #   
  #   r <- big.T.sig(LOADINGS_R[[i]][,1],sigma)
  #   res <- c(res,r)
  # }
  # 
  # 
  # RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  # 
  # density <- ggplot(RES, aes(x = Rad, y = Density))+
  #   geom_smooth() + theme_bw() + ggtitle("Number of points above 2*sigma of original rotation")
  # 
  # density
  # 
  # 
  # RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  # 
  # #taken from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
  # inflect <- function(x, threshold = 1){
  #   up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  #   down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  #   a    <- cbind(x,up,down)
  #   list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  # }
  # 
  # 
  # ggres <- ggplot_build(density)$data[[1]]
  # infs <- inflect(ggres$y)
  # rads <- ggres$x[infs$maxima]
  # 
  # #######
  # 
  # maxima.ang <- rads
  # 
  # ###
  # 
  # i <- NULL
  # # L_R <- data.frame()
  # sPLS_maxima_cands <- list()
  # 
  # for(i in 1:length(maxima.ang)){
  #   
  #   
  #   ang <- maxima.ang[i]
  #   
  #   rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  #   loadings_r <- as.data.frame(t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat)))
  #   colnames(loadings_r) <- c("comp1","comp2")
  #   # L_R <- rbind(L_R,loadings_r)
  #   sPLS_maxima_cands[[i]] <- loadings_r[order(loadings_r$comp1, decreasing = TRUE),]
  #   
  # }
  
  # sPLS_maxima_cands <- list(sPLS_bait_cands)
  sPLS_maxima_cands <- as.data.frame(sPLS_bait_cands)
  
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_maxima_cands[[1]]),Fref2$Feature),]
  
  # loadings <- loadingplot(sPLS_maxima_cands[[1]], meta, baits = Fref2$Feature[which(Fref2$Baits %in% c("Pathway","PSI"))],
  #                      importance, PCs = c(1,2), ref = Fref2, aes(x = comp1, y = comp2, colour = Baits, label = Description))
  # 
  
  # F1_scores <- F1_plot(sPLS_maxima_cands, n = 500, baits = spikes)
  
  RP <- prod(which(rownames(sPLS_maxima_cands) %in% spikes))^(1/length(spikes))
  
  
  # Results <- list("Candidates" = sPLS_maxima_cands, "Plot" = loadings, "Density" = density)
  
  return(list(RP, sPLS_maxima_cands))
}

MASCARA <- function(resids,ref, baits){
  
  Fref <- ref
  
  spls_res <- spls(resids[,-which(colnames(resids) %in% baits)], 
                   resids[,which(colnames(resids) %in% baits)], all.outputs = TRUE)
  
  sPLS_cands <- cbind.data.frame(abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X, spls_res$loadings$X) #  * spls_res$prop_expl_var$X
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
  
  importance <- paste0(round(spls_res$prop_expl_var$X * 100, digits = 2), "%")
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_cands),Fref2$Feature),]
  
  #####
  
  
  #calculating direction of baits to determine direction of candidate selection spls loading space
  
  
  angles <- rbind.data.frame(spls_res$loadings$Y, colMeans(spls_res$loadings$Y))
  
  # x <- angles[4,]
  y <- angles[nrow(angles),]
  
  #adapted from https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r
  angle <- function(x,y){
    
    x <- c(as.matrix(x))
    y <- c(as.matrix(y))
    
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    return(as.numeric(theta))
  }
  
  angles_loadings <- apply(sPLS_cands[,-1],1,angle, y = y)
  
  
  #define this better
  sPLS_cands <- cbind.data.frame(sPLS_cands, "angle" = angles_loadings, "score" = 3 * sPLS_cands$`abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X` * (1-angles_loadings))
  
  sPLS_cands <- sPLS_cands[order(sPLS_cands$score, decreasing = TRUE),]
  
  #########
  
  avg_point <- angles[nrow(angles),]/nrow(angles)
  
  # rot.mat <- matrix(c(0,-1,1,0),nrow = 2, byrow = TRUE)  #90 degree clockwise rotation
  # 
  # rotated_point <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat)
  # 
  # ang <- angle(avg_point,c(1,0))
  # 
  # rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  # rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  #######
  
  
  #rotated loadings
  ang <- angle(avg_point,c(1,0))
  
  rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  
  rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  rotate <- function(x,y){
    x <- c(as.matrix(x))
    x %*% y
  }
  
  loadings_r_baits <- t(apply(sPLS_cands[,c("comp1","comp2")],1,rotate, y = rot.mat2))
  colnames(loadings_r_baits) <- c("comp1","comp2")
  
  
  ######
  
  sPLS_bait_cands <- loadings_r_baits[order(loadings_r_baits[,1], decreasing = TRUE),]
  
  
  # #####
  # 
  # #loop over rotation matrices from -pi/4 -> pi/4
  # #resolution of 100 (?)
  # 
  # range.ang <- seq(from = -pi/6, to = pi/6, length.out = 60)
  # 
  # ###
  # 
  # i <- NULL
  # L_R <- data.frame()
  # LOADINGS_R <- list()
  # 
  # for(i in 1:length(range.ang)){
  #   
  #   
  #   ang <- range.ang[i]
  #   
  #   rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  #   loadings_r <- t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat))
  #   colnames(loadings_r) <- c("comp1","comp2")
  #   L_R <- rbind(L_R,loadings_r)
  #   LOADINGS_R[[i]] <- loadings_r
  #   
  #   
  # }
  # 
  # # colnames(LOADINGS_R) <- c("comp1","comp2")
  # L_R$rad <- rep(range.ang, each = nrow(sPLS_cands))
  # L_R$frame <- rep(c(1:length(range.ang)), each = nrow(sPLS_cands))
  # 
  # 
  # 
  # #get sd of original comp1 distribution
  # sigma <- sd(sPLS_cands$comp1)
  # 
  # big.T.sig <- function(x,y){
  #   big <- length(which(x > 2.5*y))
  #   return(big)
  # }
  # 
  # i <- NULL
  # res <- c()
  # for(i in 1:length(LOADINGS_R)){
  #   
  #   r <- big.T.sig(LOADINGS_R[[i]][,1],sigma)
  #   res <- c(res,r)
  # }
  # 
  # 
  # RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  # 
  # density <- ggplot(RES, aes(x = Rad, y = Density))+
  #   geom_smooth() + theme_bw() + ggtitle("Number of points above 2*sigma of original rotation")
  # 
  # density
  # 
  # 
  # RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  # 
  # #taken from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
  # inflect <- function(x, threshold = 1){
  #   up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  #   down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  #   a    <- cbind(x,up,down)
  #   list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  # }
  # 
  # 
  # ggres <- ggplot_build(density)$data[[1]]
  # infs <- inflect(ggres$y)
  # rads <- ggres$x[infs$maxima]
  # 
  # #######
  # 
  # maxima.ang <- rads
  # 
  # ###
  # 
  # i <- NULL
  # # L_R <- data.frame()
  # sPLS_maxima_cands <- list()
  # 
  # for(i in 1:length(maxima.ang)){
  #   
  #   
  #   ang <- maxima.ang[i]
  #   
  #   rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  #   loadings_r <- as.data.frame(t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat)))
  #   colnames(loadings_r) <- c("comp1","comp2")
  #   # L_R <- rbind(L_R,loadings_r)
  #   sPLS_maxima_cands[[i]] <- loadings_r[order(loadings_r$comp1, decreasing = TRUE),]
  #   
  # }
  
  sPLS_maxima_cands <- list(sPLS_bait_cands)
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_maxima_cands[[1]]),Fref2$Feature),]
  
  loadings <- loadingplot(sPLS_maxima_cands[[1]], meta, baits = Fref2$Feature[which(Fref2$Baits %in% c("Pathway","PSI"))],
                          importance, PCs = c(1,2), ref = Fref2, aes(x = comp1, y = comp2, colour = Baits, label = Description))
  
  Results <- list("Candidates" = sPLS_maxima_cands, "Plot" = loadings, "Density" = density)
  
  return(Results)
}



MASCARA_DEV <- function(resids,ref, baits){
  
  Fref <- ref
  
  spls_res <- spls(resids[,-which(colnames(resids) %in% baits)], 
                   resids[,which(colnames(resids) %in% baits)], all.outputs = TRUE)
  
  sPLS_cands <- cbind.data.frame(abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X, spls_res$loadings$X) #  * spls_res$prop_expl_var$X
  sPLS_cands <- sPLS_cands[order(sPLS_cands[,1], decreasing = TRUE),]
  
  importance <- paste0(round(spls_res$prop_expl_var$X * 100, digits = 2), "%")
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_cands),Fref2$Feature),]
  
  #####
  
  
  #calculating direction of baits to determine direction of candidate selection spls loading space
  
  
  angles <- rbind.data.frame(spls_res$loadings$Y, colMeans(spls_res$loadings$Y))
  
  # x <- angles[4,]
  y <- angles[nrow(angles),]
  
  #adapted from https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r
  angle <- function(x,y){
    
    x <- c(as.matrix(x))
    y <- c(as.matrix(y))
    
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    return(as.numeric(theta))
  }
  
  angles_loadings <- apply(sPLS_cands[,-1],1,angle, y = y)
  
  
  #define this better
  sPLS_cands <- cbind.data.frame(sPLS_cands, "angle" = angles_loadings, "score" = 3 * sPLS_cands$`abs(spls_res$loadings$X) %*% spls_res$prop_expl_var$X` * (1-angles_loadings))
  
  sPLS_cands <- sPLS_cands[order(sPLS_cands$score, decreasing = TRUE),]
  
  #########
  
  avg_point <- angles[nrow(angles),]/nrow(angles)
  
  # rot.mat <- matrix(c(0,-1,1,0),nrow = 2, byrow = TRUE)  #90 degree clockwise rotation
  # 
  # rotated_point <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat)
  # 
  # ang <- angle(avg_point,c(1,0))
  # 
  # rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  # rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  #######
  
  
  #rotated loadings
  ang <- angle(avg_point,c(1,0))
  
  rot.mat2 <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
  
  rotated_point2 <- as.data.frame(c(as.matrix(avg_point)) %*% rot.mat2)
  
  rotate <- function(x,y){
    x <- c(as.matrix(x))
    x %*% y
  }
  
  loadings_r_baits <- t(apply(sPLS_cands[,c("comp1","comp2")],1,rotate, y = rot.mat2))
  colnames(loadings_r_baits) <- c("comp1","comp2")
  
  
  ######
  
  sPLS_bait_cands <- loadings_r_baits[order(loadings_r_baits[,1], decreasing = TRUE),]
  
  
  #####
  
  #loop over rotation matrices from -pi/4 -> pi/4
  #resolution of 100 (?)
  
  range.ang <- seq(from = -pi/6, to = pi/6, length.out = 60)
  
  ###
  
  i <- NULL
  L_R <- data.frame()
  LOADINGS_R <- list()
  
  for(i in 1:length(range.ang)){
    
    
    ang <- range.ang[i]
    
    rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
    loadings_r <- t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat))
    colnames(loadings_r) <- c("comp1","comp2")
    L_R <- rbind(L_R,loadings_r)
    LOADINGS_R[[i]] <- loadings_r
    
    
  }
  
  # colnames(LOADINGS_R) <- c("comp1","comp2")
  L_R$rad <- rep(range.ang, each = nrow(sPLS_cands))
  L_R$frame <- rep(c(1:length(range.ang)), each = nrow(sPLS_cands))
  
  
  
  #get sd of original comp1 distribution
  sigma <- sd(sPLS_cands$comp1)
  
  big.T.sig <- function(x,y){
    big <- length(which(x > 2.5*y))
    return(big)
  }
  
  i <- NULL
  res <- c()
  for(i in 1:length(LOADINGS_R)){
    
    r <- big.T.sig(LOADINGS_R[[i]][,1],sigma)
    res <- c(res,r)
  }
  
  
  RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  
  density <- ggplot(RES, aes(x = Rad, y = Density))+
    geom_smooth() + theme_bw() + ggtitle("Number of points above 2*sigma of original rotation")
  
  density
  
  
  RES <- cbind.data.frame("Rad" = range.ang, "Density" = res)
  
  #taken from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
  inflect <- function(x, threshold = 1){
    up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
    down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
    a    <- cbind(x,up,down)
    list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  }
  
  
  ggres <- ggplot_build(density)$data[[1]]
  infs <- inflect(ggres$y)
  rads <- ggres$x[infs$maxima]
  
  #######
  
  maxima.ang <- rads
  
  ###
  
  i <- NULL
  # L_R <- data.frame()
  sPLS_maxima_cands <- list()
  
  for(i in 1:length(maxima.ang)){
    
    
    ang <- maxima.ang[i]
    
    rot.mat <- matrix(c(cos(ang),-sin(ang),sin(ang),cos(ang)),nrow = 2, byrow = TRUE) #clockwise rotation
    loadings_r <- as.data.frame(t(apply(loadings_r_baits[,c("comp1","comp2")],1,rotate, y = rot.mat)))
    colnames(loadings_r) <- c("comp1","comp2")
    # L_R <- rbind(L_R,loadings_r)
    sPLS_maxima_cands[[i]] <- loadings_r[order(loadings_r$comp1, decreasing = TRUE),]
    
  }
  
  
  Fref2 <- Fref[-which(Fref$Feature %in% baits),]
  Fref2 <- Fref2[match(rownames(sPLS_maxima_cands[[1]]),Fref2$Feature),]
  
  loadings <- loadingplot(sPLS_maxima_cands[[1]], meta, baits = Fref2$Feature[which(Fref2$Baits %in% c("Pathway","PSI"))],
                          importance, PCs = c(1,2), ref = Fref2, aes(x = comp1, y = comp2, colour = Baits, label = Description))
  
  Results <- list("Candidates" = sPLS_maxima_cands, "Plot" = loadings, "Density" = density)
  
  return(Results)
}




# # ATTEMPT at average bait cosine similarity calculations
# bait_avg <- colMeans(PCD$rotation[1997:2000,1:4])
# barplot(t(rbind(PCD$rotation[1989:2000,1:4], bait_avg)))
# 
# 
# loadings_4_cosine <- t(rbind(PCD$rotation[,1:4], bait_avg))
# 
# cos_sim <- cosine(loadings_4_cosine)
# 
# #calculate cosine similarity of baits and other features
# m_c <- melt(cos_sim)
# colnames(m_c)[3] <- "cosine_sim"
# 
# 
# m_c <- m_c[which(m_c$Var1 == "bait_avg"),]
# m_c <- m_c[order(m_c$cosine_sim, decreasing = TRUE),]


#modified from MetStat
# ASCA.GetSummary <- {
#   function (asca, quietly = FALSE) 
#   {
#     summaryPCA <- matrix(NA, nrow = length(asca$ee.names) + 1, 
#                          ncol = 10)
#     rownames(summaryPCA) <- c("data", asca$ee.names)
#     colnames(summaryPCA) <- paste("PC", 1:10, sep = "")
#     relevant.pcs <- which(asca$svd$var.explained[1:10] > 0.01)
#     max.pc <- 0
#     for (i in relevant.pcs) {
#       max.pc <- max(max.pc, i)
#       summaryPCA[1, i] <- asca$svd$var.explained[i]
#     }
#     ee.index <- 2
#     for (ee in asca$ee.names) {
#       relevant.pcs <- which(asca[[ee]]$svd$var.explained[1:10] > 
#                               0.01)
#       for (i in relevant.pcs) {
#         max.pc <- max(max.pc, i)
#         summaryPCA[ee.index, i] <- asca[[ee]]$svd$var.explained[i]
#       }
#       ee.index <- ee.index + 1
#     }
#     summaryPCA <- summaryPCA[, 1:max.pc, drop = FALSE]
#     summarySSQ <- matrix(NA, nrow = 1, ncol = length(asca$ee.names) + 
#                            2)
#     rownames(summarySSQ) <- "Contribution to ssq"
#     colnames(summarySSQ) <- c("Overall means", asca$ee.names, 
#                               "Residuals")
#     summarySSQ[1, 1] <- asca$ssq.mean/asca$ssq
#     ee.index <- 2
#     for (ee in asca$ee.names) {
#       summarySSQ[1, ee.index] <- asca[[ee]]$ssq/asca$ssq
#       ee.index <- ee.index + 1
#     }
#     summarySSQ[1, dim(summarySSQ)[2]] <- asca$ssq.remainder/asca$ssq
#     summarySSQCentered <- matrix(NA, nrow = 1, ncol = length(asca$ee.names) + 
#                                    1)
#     rownames(summarySSQCentered) <- "Contribution to ssq"
#     colnames(summarySSQCentered) <- c(asca$ee.names, "Residuals")
#     ee.index <- 1
#     for (ee in asca$ee.names) {
#       summarySSQCentered[1, ee.index] <- asca[[ee]]$ssq/(asca$ssq - 
#                                                            asca$ssq.mean)
#       ee.index <- ee.index + 1
#     }
#     summarySSQCentered[1, dim(summarySSQCentered)[2]] <- asca$ssq.remainder/(asca$ssq - 
#                                                                                asca$ssq.mean)
#     if (!quietly) {
#       cat("Variance explained per principal component (if >1%):\n")
#       for (r in 1:dim(summaryPCA)[1]) {
#         if (r == 1) {
#           cat(paste("Whole data set \t", sep = ""))
#         }
#         else if (nchar(rownames(summaryPCA)[r]) == 1) {
#           cat(paste("Factor ", rownames(summaryPCA)[r], 
#                     "     \t", sep = ""))
#         }
#         else {
#           cat(paste("Interaction ", rownames(summaryPCA)[r], 
#                     "\t", sep = ""))
#         }
#         for (c in 1:dim(summaryPCA)[2]) {
#           varExplained <- paste(colnames(summaryPCA)[c], 
#                                 ": ", formatC(summaryPCA[r, c] * 100, digits = 2, 
#                                               format = "f"), "% ", sep = "")
#           varExplained <- paste(varExplained, paste(rep(" ", 
#                                                         (13 - nchar(varExplained))), collapse = ""))
#           cat(varExplained)
#         }
#         cat("\n")
#       }
#       cat("\n")
#       cat("Percentage each effect contributes to the total sum of squares:\n")
#       for (c in 1:dim(summarySSQ)[2]) {
#         if (c == 1) {
#           cat(paste("Overall means  \t", sep = ""))
#         }
#         else if (c == dim(summarySSQ)[2]) {
#           cat(paste("Residuals      \t", sep = ""))
#         }
#         else if (nchar(colnames(summarySSQ)[c]) == 1) {
#           cat(paste("Factor ", colnames(summarySSQ)[c], 
#                     "     \t", sep = ""))
#         }
#         else {
#           cat(paste("Interaction ", colnames(summarySSQ)[c], 
#                     "\t", sep = ""))
#         }
#         cat(paste(formatC(summarySSQ[1, c] * 100, digits = 2, 
#                           format = "f"), "%\n", sep = ""))
#       }
#       cat("\n")
#       cat("Percentage each effect contributes to the sum of squares of the centered data:\n")
#       for (c in 1:dim(summarySSQCentered)[2]) {
#         if (c == dim(summarySSQCentered)[2]) {
#           cat(paste("Residuals      \t", sep = ""))
#         }
#         else if (nchar(colnames(summarySSQCentered)[c]) == 
#                  1) {
#           cat(paste("Factor ", colnames(summarySSQCentered)[c], 
#                     "     \t", sep = ""))
#         }
#         else {
#           cat(paste("Interaction ", colnames(summarySSQCentered)[c], 
#                     "\t", sep = ""))
#         }
#         cat(paste(formatC(summarySSQCentered[1, c] * 100, 
#                           digits = 2, format = "f"), "%\n", sep = ""))
#       }
#       cat("\n")
#     }
#     list(summary.pca = summaryPCA, summary.ssq = summarySSQ, summary.centred.ssq = summarySSQCentered)
#   }
# }




Create_Core_DEV_2 <- function(nreps, meta, irr_spikes = TRUE, struc_resid = FALSE, 
                              a_sigma = c(1.5,0.75,0.6,0.6), b_sigma = c(1,0.8), e_sigma = c(1,0.8,0.5),
                              noise_sd = 0.75, EffectSize = c(X_a_ab = 1, time = 0.5, E = 0.5, mu = 1, Struc_E = 1), 
                              SCORE_SEED = 1000, plot = FALSE, score_plot = FALSE, Experiment_responders = 12, ts = 1234, struc_seed = 127){
  
  if(struc_resid == TRUE){
    set.seed(ts)
    targets_a_ab <- abs(rnorm(12,mean = 0.25, sd = 0.25))
  }
  else{
    targets_a_ab <- rep(0,12)
  }
  
  # set.seed(1)
  set.seed(SCORE_SEED)
  
  a_ab <- as.data.frame(svd(matrix(rnorm(16), nrow = 4))$u)
  
  a_ab_PC1_small <- a_ab[,1]
  
  a_ab <- a_ab[rep(rownames(a_ab), each = nreps),]
  a_ab <- rbind.data.frame(a_ab, a_ab*-1)
  
  a_ab_p1 <- c(rep(0,24),                                                         
               rep(0,2000 - (24+4+12+(Experiment_responders))),
               abs(rnorm((Experiment_responders + 4),mean = 1, sd = 0.25)), 
               targets_a_ab)
  
  print(length(a_ab_p1))
  
  if(irr_spikes == FALSE){
    a_ab_p2 <- c(rep(0,length(a_ab_p1)))
    a_ab_p3 <- c(rep(0,length(a_ab_p1)))
    a_ab_p4 <- c(rep(0,length(a_ab_p1)))
  }
  else{
    a_ab_p2 <- c(rep(0,485),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1485))
    a_ab_p3 <- c(rep(0,185),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1785))
    a_ab_p4 <- c(rep(0,85),rnorm(15,mean = -0.5, sd = 0.25), rnorm(15,mean = 0.5, sd = 0.25), rep(0,1885))
  }
  P <- cbind(a_ab_p1, a_ab_p2, a_ab_p3, a_ab_p4)
  
  P_sigma <- P %*% diag(a_sigma)
  
  ############# B
  ##CREATE SCORES FOR PC1 TIME
  set.seed(22)
  b_PC1_small <- Create_Orthogonal(a_ab_PC1_small)#/400
  b_PC1_small %*% a_ab_PC1_small
  
  ##CREATE SCORES FOR PC2 TIME
  b_PC2_small <- Create_Orthogonal(b_PC1_small)
  b_PC1 <- rep(b_PC1_small, each = nreps, times = 2)
  b_PC2 <- rep(b_PC2_small, each = nreps, times = 2)
  
  b <- cbind.data.frame(meta[,1:2],b_PC1,b_PC2)
  colnames(b)[3:4] <- c("PC1","PC2")
  
  #CREATE LOADINGS TIME
  b_p1 <- c(rep(0.1,1000), rep(0,1000)  )  #make dependent on length of a_ab_p1
  b_p2 <- c(rep(0,1000), rep(0.1,1000))  #make dependent on length of a_ab_p1
  Pb <- cbind(b_p1, b_p2)
  
  Pb_sigma <- Pb %*% diag(b_sigma)
  
  ################## MATRIX MULTIPLICATION OF SCORES AND LOADINGS PER EFFECT MATRIX
  X_a_ab <- data.matrix(a_ab[]) %*% t(P_sigma)
  X_a_ab <- EffectSize["X_a_ab"] * (X_a_ab/norm(X_a_ab, type = "F"))
  norm_X_a_ab <- norm(X_a_ab, type = "F")
  
  time <- data.matrix(b[c("PC1","PC2")]) %*% t(Pb_sigma)
  time <- EffectSize["time"] * (time/norm(time, type = "F"))
  norm_time <- norm(time, type = "F")
  
  
  
  ##CREATE FEATURE MEAN VALUES and ensure they are above 0
  #set.seed here? 
  mu <- EffectSize["mu"] * matrix(rep(rnorm(dim(X_a_ab)[2], mean = mean(X_a_ab)), each = nrow(meta)), nrow = nrow(meta))
  mu <- mu + abs(min(mu))
  X_a_ab_mu <- X_a_ab + mu
  
  
  ####CREATE E
  
  if(struc_resid == TRUE){
    
    #make structured residuals instead of random noise:
    set.seed(struc_seed)
    
    # t_l <- svd(matrix(rnorm(2000*24), nrow = 24))$u
    
    t_l <- matrix(rnorm(24*3), nrow = 24)
    
    e_PC1 <- t_l[,1]
    e_PC2 <- t_l[,2]
    e_PC3 <- t_l[,3]
    
    
    
    #create random correlated effect in t_l instead of the random matrix from above
    l_t <- c(1.5,1)
    l_p <- abs(rnorm(16,sd = 1))
    l_x <- l_t %*% t(l_p)
    
    set.seed(1277)
    l_t2 <- Create_Orthogonal(l_t)
    l_p2 <- abs(rnorm(16,sd =0.6))
    l_x2 <- l_t2 %*% t(l_p2)
    
    l_X <- l_x + l_x2
    l_X <- l_X/norm(l_X, type = "F")
    
    
    set.seed(1288)
    e_p1 <- c(rep(0,1984), abs(abs(rnorm(16,sd = 0.2, mean = 0.5)))) # rnorm(8, sd = sd(l_X[1,]))l_X[1,] rep(0.5,16)
    # rtruncnorm(12, a = 0.1, b = 0.8, mean = 0.3)
    set.seed(128)
    # e_p2 <- c(rep(0,0), rnorm(2000))  #make dependent on length of a_ab_p1
    e_p2 <- c(rep(0,1984), abs(rnorm(16,sd = 0.2, mean = 0.5))) #rnorm(8, sd = sd(l_X[2,])) l_X[2,] rep(0.5,16)
    e_p3 <- c(rep(0,1984), abs(rnorm(16,sd = 0.2, mean = 0.5))) # extra some baits plus spikes
    
    
    Pe <- cbind(e_p1, e_p2, e_p3)
    colnames(Pe) <- c("PC1","PC2", "PC3")
    
    Pe_sigma <- Pe %*% diag(e_sigma)
    e <- cbind.data.frame(meta[,1:2], e_PC1, e_PC2, e_PC3)
    colnames(e)[c(3,4,5)] <- c("PC1","PC2","PC3")
    
    
    #################
    
    ET <- data.matrix(e[c("PC1","PC2","PC3")]) %*% t(Pe_sigma)
    #svd(ET)$u as input alternatively
    
    #matrix linear regression
    A <- data.matrix(cbind.data.frame(a_ab[],b[c("PC1","PC2")]))
    B <- solve(t(A)%*%A) %*% t(A) %*% ET
    
    E <- ET - (A %*% B)  
    
    set.seed(10021)
    Struc_E <- E
    Struc_E <- EffectSize["Struc_E"] * (Struc_E/norm(Struc_E, type = "F"))
    norm_Struc_E <- norm(Struc_E, type = "F")
    
    
    dim(Struc_E)
    print(E[,1] %*% time[,1])
    print(time[,1] %*% X_a_ab[,1])
    
    
    Non_Struc_E <- matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    Non_Struc_E <- EffectSize["E"] * (Non_Struc_E/norm(Non_Struc_E, type = "F"))
    norm_Non_Struc_E <- norm(Non_Struc_E, type = "F")
    
    
    E <- Struc_E + Non_Struc_E
    
    norms <- c("main" = norm_X_a_ab, "cofactor" = norm_time, "struc_E" = norm_Struc_E, "Non_struc_E" = norm_Non_Struc_E)
    
    ############
    
  }
  else{
    
    ## some random noise - i.e. differences between replicates
    set.seed(1001)  #previously 101
    E <- matrix(rnorm(dim(X_a_ab)[1] * dim(X_a_ab)[2], mean = 0, sd = noise_sd), nrow = nrow(meta))
    
    
    E <- EffectSize["E"] * (E/norm(E, type = "F"))
    norm_Non_Struc_E <- norm(E, type = "F")
    
    
    #change this 
    norms <- c("main" = norm_X_a_ab, "cofactor" = norm_time, "Non_struc_E" = norm_Non_Struc_E)
    
  }
  
  ##ADD ALL TOGETHER TO CREATE "CORE"
  X_no.time <- X_a_ab_mu + E
  X <- X_no.time + time
  X <- t(((t(X) + abs(colMins(X)))))
  colnames(X) <- paste0("X_",c(1:dim(X)[2]))
  
  
  Effects <- list("main" = X_a_ab, "cofactor" = time, "residual" = E, "ambient" = Struc_E, "noise" = Non_Struc_E)
  # norms <- c("main" = norm_X_a_ab, "cofactor" = norm_time, "struc_E" = norm_Struc_E, "Non_struc_E" = norm_Non_Struc_E)
  
  if(plot == TRUE){
    
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    t2 <- ggplot(a_ab, aes(x = time, y = PC2, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      theme_classic() +
      theme(legend.position = "none")
    P_anno <- cbind.data.frame(c(1:dim(P)[1]), P)
    colnames(P_anno) <- c("Feature", "PC1", "PC2")
    p1 <- ggplot(P_anno[which(P_anno$PC1 != 0),], aes(x = Feature, y = PC1)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    p2 <- ggplot(P_anno[which(P_anno$PC2 != 0),], aes(x = Feature, y = PC2)) +
      geom_bar(stat="identity") +
      theme_classic() +
      # scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) 
      theme(axis.text.x = element_blank())
    
    PLOT <- (t1 | t2) / (p1 | p2) + plot_annotation(title = "Simulated data, scores and loadings in PCs 1 and 2")
    
    return(list(X,P,Pb,a_ab,PLOT, X_a_ab,time, Struc_E ))
    
  }
  else if (score_plot == T){
    t1 <- ggplot(a_ab, aes(x = time, y = PC1, group=growth_condition, colour=growth_condition)) +
      geom_line() +
      ylim(-3,3) + 
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
      theme_classic() +
      theme(legend.position = "none")
    
    return(list(X,P,Pb,a_ab,t1))
  }
  else{
    return(list(X,P,Pb,a_ab,"Effects" = Effects, norms))  #Effect inputs here?..
  }
  
}


