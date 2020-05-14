# File: functions.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Defines functions for different weighted swapping algorithms


########################
####### Sampling #######
########################

swap <- function(A, W){
  # Weighted Checkerboard Swapping
  # A: binary matrix
  # W: weight matrix (positive elements, same size as A)
  
  # sample two 1's at random
  ones <- which(A == 1, arr.ind=T)
  ones_smp <- sample(1:nrow(ones),2)  
  o1 <- ones[ones_smp[[1]],]
  o2 <- ones[ones_smp[[2]],]
 
  # construct coords of opposing corners
  c1 <- c(o1[1], o2[2])
  c2 <- c(o2[1], o1[2])
  
  # Determine if a swap is allowed
  is_checkerboard <- (A[c1[1],c1[2]]==0) & (A[c2[1],c2[2]]==0)
  is_not_struct_zeros <- (W[c1[1],c1[2]]!=0) & (W[c2[1],c2[2]]!=0)
  can_swap <- is_checkerboard & is_not_struct_zeros
 
  did_swap <- p <- NA
  if(can_swap){
    # compute probability of swap and flip coin
    p_swap <- W[c1[1],c1[2]]*W[c2[1],c2[2]]/(W[c1[1],c1[2]]*W[c2[1],c2[2]] + W[o1[1],o1[2]]*W[o2[1],o2[2]])
    if(runif(1) < p_swap){
      # perform the swap
      A[o1[1],o1[2]] = A[o2[1],o2[2]] = 0
      A[c1[1],c1[2]] = A[c2[1],c2[2]] = 1
    }
  }
  return(A)
}


curveball_trade <- function(A, W){
  # Weighted curveball swapping
  # A: binary matrix
  # W: weight matrix (positive elements, same size as A)
  
  # sample two rows at random
  rows <- sample(1:nrow(A))
  r1 <- A[rows[1],]
  r2 <- A[rows[2],]
  w1 <- W[rows[1],]
  w2 <- W[rows[2],]
  # get columns shared and not shared
  A12 <- which(r1+r2==2)  # Shared between r1 and r2
  A1m2 <- which(r1-r2>0 & w2>0)  # in r1 but not r2 (and not a struct zero in 2)
  A2m1 <- which(r2-r1>0 & w1>0)  # in r2 but not r1 (and not a struct zero in 1)
  
  A1m2sz <- which(r1-r2>0 & w2==0)  # in r1 but not r2 (and a struct zero in 2)
  A2m1sz <- which(r2-r1>0 & w1==0)  # in r2 but not r1 (and a struct zero in 1)
  
  # Check if can trade
  can_trade <- !(length(A1m2)==0 | length(A2m1)==0)
  if(!can_trade){
    return(A)
  }
  # shuffle the tradable indices
  B <- sample(c(A1m2, A2m1))
  B1m2 <- B[1:length(A1m2)]
  B2m1 <- B[(length(A1m2)+1):length(B)]
  # compute trade probability
  r1_01 <- B1m2[is.na(match(B1m2, A1m2))]  # cols that went 0 -> 1 in row 1
  r2_01 <- B2m1[is.na(match(B2m1, A2m1))]  # cols that went 0 -> 1 in row 2
  p_trade <- prod(W[rows[1],r1_01])*prod(W[rows[2],r2_01])/(prod(W[rows[1],r1_01])*prod(W[rows[2],r2_01]) + prod(W[rows[2],r1_01])*prod(W[rows[1],r2_01]))
  if(runif(1) < p_trade){
    # construct traded rows
    r1B <- rep(0, length(r1))
    r1B[c(A12, A1m2sz, B1m2)] <- 1  # shared 1's, sz 1's, and traded 1's
    r2B <- rep(0, length(r2))
    r2B[c(A12, A2m1sz, B2m1)] <- 1  # shared 1's, sz 1's, and traded 1's
    
    A[rows[1],] <- r1B
    A[rows[2],] <- r2B
  }
  return(A)
}

##########################
####### Statistics #######
##########################

diag_stat <- function(X){
  # Diagonal Divergence
  ones <- which(X==1, arr.ind=T)
  return(sum(abs(ones[,1] - ones[,2])))
}

cscore <- function(X){
  # From Gotelli 2000, Null Model Analysis of Species Cooccurrence Patterns
  S <- rowSums(X)  # Islands per species
  R <- length(S)  # number of rows
  CO <- X %*% t(X)  # species cooccurrence
  D1 <- matrix(S, nrow(CO), ncol(CO)) - CO           # differences
  D2 <- matrix(S, nrow(CO), ncol(CO), byrow=T) - CO  # differences
  CB <- as.data.frame(D1*D2) %>% rownames_to_column() %>% 
    pivot_longer(row.names(X)) %>%
    filter(rowname < name)
  C <- sum((D1*D2)[upper.tri(CO)])/(R*(R-1)/2)  # c score
  return(list(C=C, CB=CB))
}