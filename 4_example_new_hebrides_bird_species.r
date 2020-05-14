# File: 4_example_new_hebrides_bird_species.r
# Author: Alex Fout (fout@colostate.edu)
# Purpose: Perform weighted sampling on the new hebrides data


#####################
####### Setup #######
#####################

rm(list=ls())

source("functions.r")
library(reshape2)     # For melt
library(ggplot2)      # for plotting
theme_set(theme_minimal() + theme(axis.ticks=element_line()))
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")


output_dir <- file.path("output")
dir.create(output_dir, showWarnings=F)


# MCMC stuff
N <- 5000   # number of samples per weight
burnin <- 1000  # burn in before sampling
thin <- 500   # thinning

########################
####### Sampling #######
########################

# Prep A, W
A <- read.csv("cooccurrence.csv", row.names=1)
A <- 1*(A==1)

W1 <- W2 <- matrix(1, nrow(A), ncol(A))
k <- 2
w <- 15
r <- ceiling(nrow(A)/k)
c <- ceiling(ncol(A)/k)
for(i in 1:(k-1)){
  W2[((i-1)*r):((i*r)-1), ((i-1)*c):((i*c)-1)] <- w
}
W2[((k-1)*r):nrow(A), ((k-1)*c):ncol(A)] <- w


scores0 <- cscore(A)
C0 <- scores0$C
CB0 <- scores0$CB

# Sample:
{
  set.seed(utf8ToInt("birds aren't real"))
  # Uniform sampling
  A1 <- A
  Cs1 <- array(N)
  CBs1 <- list()
  for(k in 1:burnin){
    A1 <- curveball_trade(A1, W1)
  }
  for(j in 1:N){
    for(k in 1:thin){
      A1 <- curveball_trade(A1, W1)
    }
    stats <- cscore(A1)
    Cs1[j] <- stats$C
    CBs1[[j]] <- cbind(stats$CB, iter=j)
  }
  CB1 <- do.call(rbind, CBs1)
}
{
  # Weighted Sampling
  A2 <- A
  Cs2 <- array(N)
  CBs2 <- list()
  for(k in 1:burnin){
    A2 <- curveball_trade(A2, W2)
  }
  for(j in 1:N){
    for(k in 1:thin){
      A2 <- curveball_trade(A2, W2)
    }
    stats <- cscore(A2)
    Cs2[j] <- stats$C
    CBs2[[j]] <- cbind(stats$CB, iter=j)
  }
  CB2 <- do.call(rbind, CBs1)
}

########################
####### Plotting #######
########################

# Plot A, W
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
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(expand=c(0,0)) + 
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

legend1 <-  theme(legend.position="bottom")
rotatexlabs <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
scalex1 <- scale_x_continuous(breaks=seq(0, 10000, 2500))
scaley2 <- scale_y_continuous(breaks=seq(0, 1250, 500))

spec <- read.csv("species.csv")
isla <- read.csv("islands.csv")

dfA <- melt(A)
dfA$Species <- factor(dfA$Var1, levels=rev(spec$abbr), labels=rev(spec$name))
dfA$Island <- factor(dfA$Var2, levels=isla$abbr, labels=isla$name)
ggplot(dfA, aes(x=Island, y=Species)) + 
  geom_tile(aes(fill=as.factor(value)), color="grey60") + 
  scale_fill_manual(name="Fill:", values=c("white", "grey80")) + 
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_fixed(ratio=1)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "top",
        aspect.ratio=1.0,
        panel.border=element_rect(fill=NA))
ggsave(file.path(output_dir, "4_nh_A_2.png"), height=6, width=6, dpi="retina")


rownames(W2) <- rev(spec$name)
colnames(W2) <- isla$name
dfW2 <- melt(W2)
dfW2$Species <- dfW2$Var1
dfW2$Island <- dfW2$Var2
ggplot(dfW2, aes(x=Island, y=Species)) + 
  geom_tile(aes(fill=as.factor(value)), color="grey60") + 
  scale_fill_manual(name="Fill:", values=c("white", "grey80")) + 
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_fixed(ratio=1)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "top",
        aspect.ratio=1.0,
        panel.border=element_rect(fill=NA))
ggsave(file.path(output_dir, "4_nh_W2_2.png"), height=6, width=6, dpi="retina")


# Plot sampling distributions with p-values
cs <- data.frame(cscore=c(Cs1, Cs2), weight=rep(c("uniform", "weighted"), each=N))
ps <- cs %>% group_by(weight) %>%
  mutate(p=cscore>C0) %>% 
  summarize(p=mean(p))
p1 <- (ps %>% filter(weight=="uniform") %>% select(p))[[1]]
p2 <- (ps %>% filter(weight=="weighted") %>% select(p))[[1]]

xs = c(8.6, 8.6)
ys = c(750, 500)
texts = c(paste0("p[Uniform]==",p1), paste0("p[Weighted]==",p2))
cs %>%
  ggplot(aes(cscore, fill=weight)) + 
  geom_histogram(position="identity", bins=40, alpha=0.5, color="black") + 
  geom_vline(xintercept=C0, linetype="dashed") + 
  scale_fill_manual(name="", labels=c("Uniform", "Weighted"), values=c("uniform"="red", "weighted"="blue")) + 
  annotate("text", x=xs, y=ys, label=texts, parse=TRUE, size=4.5, hjust=0) + 
  labs(x="C score", y="Frequency") + 
  theme(legend.position="bottom")
ggsave(file.path(output_dir, "4_nh_cscore_hist.png"), height=4, width=7, dpi="retina")
