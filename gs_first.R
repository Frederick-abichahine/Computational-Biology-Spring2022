library(TCGAbiolinks)
library(isma)
library(mND)
library(limma)
library(igraph)
library(pcaPP)
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
gs <- geneNIDG(distance.matrix = ge_dist, cor.matrix = ge_cor, geneStats.observed = ge,
               perm.geneStats.matrix = as.matrix(ge_resampled), genes.assayedETnetwork = genes.assayedETnetwork, diameter = diam,
               num.Sphere.resamples = 1000, gene.id = rownames(A)
               )

### Troubleshooting
distances <- ge_dist[rownames(A),genes.assayedETnetwork]

for(i in 1:diam){
  igenes.distances <- distances[distances <= i & distances > 0]
  igenes.names <- names(igenes.distances)
  print(length(abs(ge[igenes.names])))
  print(length(igenes.distances))
  return(
    
    cor.fk(abs(ge[igenes.names]), igenes.distances)
    
  )
}

idx <- which(distances <= 3 & distances > 0, arr.ind = TRUE)
distances[idx]

######### NOTES ##########
# unique(idx) takes so much time i couldn't wait for it. How to subset the matrix by conditions to give a smaller matrix (not a list!)
# Check GeneSurrounder.R line 325
# Error thrown: Error in cor.fk(abs(ge[igenes.names]), igenes.distances) : 
#                       x and y must have same length.
        