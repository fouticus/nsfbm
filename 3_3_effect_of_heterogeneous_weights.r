# File: 3_3_effect_of_heterogeneous_weights.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: test a statistic under different weighting schemes


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")
library(ggplot2)        # for plotting
library(reshape2)       # for melt
library(grid)           # for newpage
theme_set(theme_minimal() + theme(axis.ticks=element_line()))


output_dir <- file.path("output")
dir.create(output_dir, showWarnings=F)


# MCMC stuff
N <- 5000   # number of samples per weight
burnin <- 5000  # burn in before sampling
thin <- 1000   # thinning
ps <- c(-10, -4, -2, 0, 2, 4, 10)  # which powers to evaluate


########################
####### Sampling #######
########################


# Build A, W
n <- m <- 50

zeros <- matrix(0, n/2, m/2)
ones <- matrix(1, n/2, m/2)
A <- rbind(cbind(zeros, ones), cbind(ones, zeros))
W <- outer(1:n, 1:m, FUN=function(x, y){n+m-abs(x-y)})/(n+m)

# sample:
sim <- function(p){
  A2 <- A
  W2 <- W^p/max(W^p)
  As <- array(0, dim=c(n, m, N))
  for(k in 1:burnin){
    A2 <- swap(A2, W2)
  }
  for(j in 1:N){
    for(k in 1:thin){
      A2 <- swap(A2, W2)
    }
    As[,,j] <- A2
  }
  
  # compute summary stats
  stats <- apply(As, 3, diag_stat)
  
  df <- data.frame(p=as.factor(p), iter=1:N, stat=stats)
  return(df)
}

dfs <- lapply(ps, sim)
df <- do.call(rbind, dfs)



########################
####### Plotting #######
########################

# plot sampling distributions
denom <- n*sum(A)

xs = c(15300, 17300, 19000, 20800, 22600, 24300, 26700)/denom
ys = c(850, 750, 700, 700, 750, 800, 850)
texts = c("p==10", "p==4", "p==2", "p==0", "p==-2", "p==-4", "p==-10")

df %>% 
ggplot(aes(x=stat/denom, y=..count.., fill=p)) + 
  geom_histogram(alpha=0.8, color="black", bins=200, size=0.2) +
  scale_fill_viridis_d(name="p:") + 
  scale_x_continuous(breaks=c(0.25, 0.3, 0.35, 0.4, 0.45)) + 
  labs(x="Diagonal Divergence (T)", y="Frequency") + 
  annotate("text", x=xs, y=ys, label=texts, parse=TRUE, size=5) + 
  theme(legend.position="none",
        axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        )
ggsave("3_3_diag_sampling_dist.png", path=output_dir, width=10, height=4, dpi="retina")

# Plot A and W
eb <- element_blank()
plot_bare <- function(p){
  gt <- ggplotGrob(p)
  gt <-  gtable::gtable_filter(gt, "panel")
  grid.newpage()
  grid.draw(gt)
  class(gt) <-  c("Panel", class(gt))
  print.Panel <- function(x){
    grid.newpage()
    grid.draw(x)
  }
}
plot_01 <- function(X){
  p <- ggplot(melt(X), aes(x=Var1, y=Var2)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient(low="white", high="black", limits=c(0,max(X))) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_reverse(expand=c(0,0)) + 
    coord_fixed(ratio=1)+
    theme(axis.title = eb,
          axis.text = eb,
          axis.ticks = eb,
          legend.position = "none",
          panel.grid=eb,
          aspect.ratio=1.0,
          panel.border=element_rect(fill=NA))
  plot_bare(p)
}

plot_01(A)
ggsave(file.path(output_dir, "3_3_A.png"), height=4, width=4, dpi="retina")

plot_01(W)
ggsave(file.path(output_dir, "3_3_W.png"), height=4, width=4, dpi="retina")
