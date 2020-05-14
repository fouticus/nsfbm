# File: 3_1_small_example.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Verify for a small example that we've sampled correctly.


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")  # Load swapping functions
library(itsmr)         # For computing standard errors

# sampling parameters
N <- 1000         # number of samples
burnin <- 10000   # burn in before sampling
thin <- 10000     # thinning

# cross reference for state names
letters <- list("1,0,0,0,1,0,0,0,1" = "A",
                "0,1,0,1,0,0,0,0,1" = "B",
                "1,0,0,0,0,1,0,1,0" = "C", 
                "0,0,1,0,1,0,1,0,0" = "D",
                "0,0,1,1,0,0,0,1,0" = "E",
                "0,1,0,0,0,1,1,0,0" = "F")

########################
####### Sampling #######
########################


# Build A, W
A <- diag(3)
W <- matrix(c(1,2,1,
              2,1,2,
              1,2,1), 3,3)

# Function to compute standard errors:
TH_se <- function(x){
  # Tukey Hanning Standard Error Estimates
  wn <- function(k, bn){
    # Blackman-Tukey window
    a <- 1/4  # Makes this the Tukey-Hanning window
    return((1 - 2*a + 2*a*cos(pi*abs(k)/bn)) * (1*(abs(k)<bn)))
  }
  m <- length(x)
  nu <- 0.5  # convenient choice to satisfy theorem assumptions
  gammas <- itsmr::acvf(x, h=m-1)
  ns <- 1:m
  sig2_hat <- numeric(m)
  for(i in 1:length(ns)){
    n <- ns[i]
    bn <- floor(n^nu)
    ws <- wn(0:(m-1), bn)
    sig2_hat[[i]] <- sum(ws*gammas)
  }
  return(sqrt(sig2_hat))
}

# set the seed
set.seed(utf8ToInt("toyland, toyland, little girl and boy land"))
  
# Do the sampling
A2 <- A
As <- array(0, dim=c(3, 3, N))   # remember each state
for(k in 1:burnin){
  A2 <- swap(A2, W)
}
for(j in 1:N){
  for(k in 1:thin){
    A2 <- swap(A2, W)
  }
  As[,,j] <- A2
}

# get all the states that we've encountered
states <- new.env()
for(i in 1:dim(As)[3]){
  state <- paste0(as.character(As[,,i]), collapse=",")   # represent each state as a character vector
  assign(state, 1, envir=states)
}
state_names <- ls(states)
ns <- length(state_names)
df <- data.frame(counts=numeric(ns), unnormed_likelihood=numeric(ns), probs=numeric(ns), emp_probs=numeric(ns))
rownames(df) <- state_names

# count states at each step
S <- matrix(0, N, ns)
colnames(S) <- state_names
for(i in 1:dim(As)[3]){
  state <- paste0(as.character(As[,,i]), collapse=",")
  S[i,state] <- 1
}

# cumulative count for each state
SC <- apply(S, 2, function(x){cumsum(x)})

# relative probability for each state 
SP <- t(apply(SC, 1, function(x){x/sum(x)}))

# compute empirical probability
df$counts <- SC[N,]
df$emp_probs <- df$counts/N

# compute theoretical probability (Assume we've visited all states)
for(name in state_names){
  A <- matrix(as.numeric(strsplit(name, ",")[[1]]), 3, 3)
  prod(W^A)
  df[name, "unnormed_likelihood"] <- prod(W^A)
}
df$probs <- df$unnormed_likelihood/sum(df$unnormed_likelihood)

# compute KL-Divergence at each step
KLD <- apply(SP, 1, function(x){sum(df$probs*log(df$probs/x))})

# Compute Tukey Hanning Standard errors:
for(name in state_names){
  df[name, "SE"] <- round(TH_se(S[,name])[N]/sqrt(N), 4)
}

# add state labels
df$letter <- unlist(lapply(rownames(df), function(x) letters[[x]]))

print(df)
print(paste("KLD:", tail(KLD, 1)))