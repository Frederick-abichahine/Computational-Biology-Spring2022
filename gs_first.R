library(isma)
library(mND)
library(limma)
library(igraph)
library(pcaPP)
library(SciViews)
library(tidyverse)
#library(rsample)
setwd("Desktop/LAU/Spring 2022/BIF 498HG (CSC 613)/Assignments/Final Project/Implementation")
source("GeneSurrounder.R")

data(X0)
data(A)

###########################
### mND only      #########
###########################
W <- normalize_adj_mat(A)
X0_perm <- perm_X0(X0, r = 50, W, seed_n = 2)

Xs <- ND(X0_perm, W, cores = 2)
#data(Xs)

ind_adj <- neighbour_index(W)

mND_score <- mND(Xs, ind_adj, k=3, cores = 2)

mND_score <- signif_assess(mND_score)

###########################
### GeneSurrounder ########
###########################

## Preprocessing
ge <- X0[,2]
plot(sort(ge))
int_network <- graph_from_adjacency_matrix(A, mode = "undirected")
ge_dist <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
ge_cor <- calcCorMatrix(ge, corMethod = "pearson", exprName = "ge", useMethod = "everything")
diam <- diameter(int_network)
genes.assayedETnetwork <- intersect(
  rownames(ge_dist),
  rownames(ge_cor))

## Resampling
#ge_boot <- bootstraps(as.vector(ge), times = 1000)

ge_resampled <- data.frame(matrix(ncol = nrow(A), nrow = 1000))
colnames(ge_resampled) <- rownames(A)

for(i in 1:1000){
  ge_resampled[i,] <- sample(ge, length(ge), replace = TRUE)
}

## Gene Surrounder
gs_results <- data.frame()
for(i in 1:nrow(A)){
gs <- geneNIDG(distance.matrix = ge_dist, cor.matrix = ge_cor, geneStats.observed = ge,
               perm.geneStats.matrix = as.matrix(ge_resampled), genes.assayedETnetwork = genes.assayedETnetwork, diameter = diam,
               num.Sphere.resamples = 1000, gene.id = rownames(A)[i]
               )
gs <- gs %>% mutate(p.Fisher = -2*(ln(p.Sphere) + ln(p.Decay)))
gs_results <- rbind(gs_results, gs[which.min(gs$p.Fisher),])
}


######### NOTES ##########
# unique(idx) takes so much time i couldn't wait for it. How to subset the matrix by conditions to give a smaller matrix (not a list!)
# Check GeneSurrounder.R line 325
# Error thrown: Error in cor.fk(abs(ge[igenes.names]), igenes.distances) : 
#                       x and y must have same length.
#### UPDATE
# GS only runs on 1 gene at a time, only outputs p.Sphere and p.Decay not p.Fisher
# Created for loop to run gs on each gene, calculate p.Fisher, and keep results for min(p.Fisher) only
# -2(ln(p.Sphere) + ln(p.Decay)) does not give a p-value. It gives X2 distribution with 4 degrees of freedom
# How to convert to p-value?



