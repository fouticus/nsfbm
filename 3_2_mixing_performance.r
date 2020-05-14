# File: 3_2_mixing_performance.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Examing mixing for different weights and structural zeros


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")  # Load swapping functions
  
library(dplyr)    # for manipulating data
library(ggplot2)  # for plots
theme_set(theme_minimal() + theme(axis.ticks=element_line()))

# output directory
output_dir <- file.path("output")
dir.create(output_dir, showWarnings=F)

# size of matrix
n <- m <- 20

# fill densities of A
ps <- c(0.1, 0.25, 0.5) 

# weights for W's
wss <- c("uniform", "runif", "exp", "runif2")

# struct zero schemes
zss <- c("runif_10", "runif_25", "runif_50", "tri_10", "tri_25", "tri_50")


# MCMC stuff
N <- 10000      # iterations per case
seeds <- c(7272, 7128, 9972, 819233, 1112, 5342, 6104, 9712, 2287, 71293)
lm <- 5000  # Max lag in autocorrelation plots

########################
####### Sampling #######
########################

# Create A's, W's, and Z's
set.seed(utf8ToInt("Root, root for the home team!"))

# Make W's
W_list <- list()
for(ws in wss){
  if(ws == "uniform"){
    W <- matrix(1, n, m)
  } else if(ws == "runif"){
    W <- matrix(runif(n*m), n, m)
  } else if(ws == "exp"){
    W <- matrix(rexp(n*m), n, m)
  } else if(ws == "runif2"){
    W <- matrix(runif(n*m, 0.5), n, m)
  }
  W_list[[ws]] <- W
}

# Make Z's
Z_list <- list()
for(zs in zss){
  if(grepl("tri", zs)){
    nz <- sqrt(as.numeric(strsplit(zs, "_")[[1]][[2]])/100 *n*m)
    Z <- matrix(1*(outer(1:n, 1:m, function(x,y){x+y})>nz), n, m)
  } else if(grepl("runif", zs)){
    pct <- as.numeric(strsplit(zs, "_")[[1]][[2]])/100
    Z <- matrix(1*(runif(n*m)>pct), n, m)
  }
  Z_list[[zs]] <- Z
}

# Make A's
A_list <- list()
for(p in ps){
  A <- matrix(runif(n*m)<p, n, m) * 1;
  A_list[[paste0("rand_", p)]] <- A
}
# take only 0.25 case:
A_list[["rand_0.1"]] <- NULL
A_list[["rand_0.5"]] <- NULL
ass <- names(A_list)



# Function for running the simulation on each case
sim <- function(as, ws, zs){
  dfs <- list()
  case <- paste(as, ws, zs, sep="_")
  for(seed in seeds){
    set.seed(seed)
    
    print(paste(case, seed))
    # get relevant matrices
    A <- A_list[[as]]
    W <- W_list[[ws]]
    Z <- Z_list[[zs]]
    A <- A*Z   # force A to be zero if there are structural zeros
    W <- W*Z   # force W to be zero if there are structural zeros
  
    # Do checkerboard swaps
    A1 <- A
    A1s <- array(0, dim=c(n, m, N))
    A1s[,,1] <- A
    for(i in 2:N){
      A1 <- swap(A1, W)
      A1s[,,i] <- A1
    }
  
    # Do curveball trades
    A2 <- A
    A2s <- array(0, dim=c(n, m, N))
    A2s[,,1] <- A
    for(i in 2:N){
      A2 <- curveball_trade(A2, W)
      A2s[,,i] <- A2
    }
  
    # compute diag statistic for each matrix
    ds1 <- apply(A1s, 3, diag_stat)
    ds2 <- apply(A2s, 3, diag_stat)
    
    # make dataframe
    zsv <- strsplit(zs, "_")[[1]]
    asv <- strsplit(as, "_")[[1]]
    df <- data.frame(method=c(rep("Swap",N), rep("Curveball",N)), iter=rep(1:N,2), stat=c(ds1, ds2), 
                     seed = seed, as=as.factor(as), ws=as.factor(ws), zs=as.factor(zs), 
                     zs1=zsv[1], zs2=zsv[2], as1=asv[1], as2=asv[2])
   
    # compute effective sample size 
    df$ess <- 0
    for(method in c("Swap", "Curveball")){
      # Effective sample size at each step
      stat <- df[df$method==method,c("stat")]
      N <- length(stat)
      ess <- rep(0, N) 
      for(j in 2:N){ # step j
        ess[j] <- coda::effectiveSize(stat[1:j])
      }
      df[df$method==method,c("ess")] <- ess
    }
    
    dfs[[seed]] <- df
  }
  df <- do.call(rbind, dfs)
  return(df)
}

# Run each simulation. 
params <- expand.grid(ass, wss, zss, stringsAsFactors=F)
colnames(params) <- c("as", "ws", "zs")
dfs <- lapply(1:nrow(params), function(i) {do.call(sim, as.list(params[i,]))})
df <- do.call(rbind, dfs); rm(dfs)
df$zs2 <- factor(df$zs2, levels=c("0", "10", "25", "50"))
df$zs1 <- factor(df$zs1, levels=c("runif", "tri"))

########################
####### Plotting #######
########################

# Labels for the plot
weight_labels <- as_labeller(c(exp="Weights:\nExp(1)", runif="Weights:\nUniform(0,1)", runif2="Weights:\nUniform(0.5,1)", uniform="Weights:\nAll 1's"))
make_percent <- function(ps){
  ps <- as.numeric(ps)
  ps <- vapply(ps, function(p){if(p<=1){p <- p*100}else{p}}, 1)
  paste0(round(ps), "%")
}
fill_perc <- function(ps){paste("Fill:", make_percent(ps), sep="\n")}
zero_perc <- function(ps){paste("Structural Zeros:", make_percent(ps), sep="\n")}

# Plot 1
df_ave <- df %>% 
  group_by(method, iter, as, ws, zs, zs1, zs2, as1, as2) %>%
  summarize(ess=mean(ess))

df_ave %>% 
  filter(iter %% 50 == 0) %>%
  ggplot(map=aes(x=iter, y=ess, color=zs1, linetype=method)) + 
  facet_grid(ws~zs2, scales="fixed", 
             labeller=labeller(ws=weight_labels, zs2=zero_perc)) + 
  geom_line(size=0.8) + 
  scale_linetype_manual(name="Method:", values=c("solid", "dashed"), labels=c("Curveball", "Swap")) + 
  scale_color_manual("Structural Zeros Type:", labels=c("Random", "Monotonic"), values=c("#4275f5", "#ff2424")) +
  labs(x="Iteration", y="Effective Sample Size") + 
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000)) + 
  scale_y_continuous(breaks=c(0, 25, 50), limits=c(0,50)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # rotate x labels
        strip.text.y=element_text(size=12, angle=0),  # rotate facets
        strip.text.x = element_text(size=12),
        legend.justification = c(0,1),
        legend.position= c(0.01,0.99),
        legend.box = "horizontal",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12),
        legend.background = element_rect(fill="white", size=0.4),
        panel.spacing = unit(0.1, "cm"), 
        panel.border=element_rect(fill=NA),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        ) + 
  guides(linetype = guide_legend(order = 1), 
         color = guide_legend(order = 2)); p
ggsave(file.path(output_dir, "3_2_ess_aves.png"), height=8, width=8, dpi="retina")


# Plot 2
df_ex <- df %>% 
  filter(as2=="0.25", zs2=="50", ws=="uniform", iter %% 10 == 0)

df_ex_ave <- df_ex %>% 
  group_by(as1, as2, zs1, zs2, ws, method, iter) %>%
  summarize(ess=mean(ess))


ggplot(map=aes(x=iter, y=ess, color=zs1, linetype=method)) + 
  geom_line(data=df_ex, aes(group=paste(zs1, method, seed)), size=0.50, alpha=0.3) + 
  geom_line(data=df_ex_ave, size=1.1) + 
  scale_color_manual("Structural Zeros Type:", labels=c("Random", "Monotonic"), values=c("#4275f5", "#ff2424")) +
  scale_linetype_manual(name="Method:", values=c("solid", "dashed")) + 
  labs(x="Iteration", y="Effective Sample Size") + 
  scale_x_continuous(breaks=c(0, 2500, 5000, 7500, 10000)) + 
  scale_y_continuous(breaks=c(0, 25, 50), limits=c(0,50)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # rotate x labels
        strip.text.y=element_text(size=12, angle=0),  # rotate facets
        strip.text.x = element_text(size=12),
        legend.justification = c(0,1),
        legend.position= c(0.01,0.99),
        legend.box = "horizontal",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12),
        legend.background = element_rect(fill="white", size=0.4),
        panel.spacing = unit(0.1, "cm"), 
        panel.border=element_rect(fill=NA),
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        ) + 
  guides(linetype = guide_legend(order = 1), 
         color = guide_legend(order = 2))
ggsave(file.path(output_dir, "3_2_ess_example.png"), height=5, width=8, dpi="retina")



