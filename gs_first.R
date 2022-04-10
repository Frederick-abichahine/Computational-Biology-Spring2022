## ==========================================================================
## Installing Libraries
## ==========================================================================
## 
## Install Several Required packages:
install.packages(c("devtools", "gplots",
                   "ggplot2", "igraph",
                   "lattice", "knitr","rsample",
                   "RColorBrewer", "rmarkdown",
                   "stringr", "UpSetR", "vcfR",
                   "pcaPP", "SciViews", "tidyverse"))
## ---------------------------------------------------------------------------
## TCGAbiolinks package (from github):
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
## ---------------------------------------------------------------------------
## Bioconductor packages: (R >= 3.5)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("IRanges", "GenomicRanges", "GenomicFeatures",
                       "org.Hs.eg.db", "Rsamtools",
                       "SummarizedExperiment",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "VariantAnnotation"))
## ---------------------------------------------------------------------------
## Install ISMA - For Integrative retrieval
## and analysis of Mutations
install.packages("isma/isma_0.1.4.tar.gz", repos=NULL)
## ---------------------------------------------------------------------------
## Install mND - For Multi Network Diffusion
install.packages("mND/mND_0.1.7.tar.gz", repos = NULL)
## ---------------------------------------------------------------------------
## Get GeneSurrounder.R and calc_p.R from github and Source it
source("genesurrounder/GeneSurrounder.R")
source("genesurrounder/run_geneSurrounder.R")
source("mND/calc_p.R")

## ==========================================================================
## Loading Libraries
## ==========================================================================

library(isma)
library(mND)
library(limma)
library(igraph)
library(pcaPP)
library(SciViews)
library(tidyverse)
#library(rsample)

## ==========================================================================
## Get Data Inputs
## ==========================================================================

## Now we are ready to start our analysis
data(X0)
data(A)

## ==========================================================================
## First we Try mND Only - Without GeneSurrounder
## ==========================================================================

## Normalize the Adjacency Matrix
W <- normalize_adj_mat(A)
## Permute the Layers Matrix
X0_perm <- perm_X0(X0, r = 50, W, seed_n = 2)

## Perform Network Diffusion - Non-Windows
Xs <- ND(X0_perm, W, cores = 2)

## Perform Network Diffusion - Windows
Xs <- ND(X0_perm, W)

## Get Indices of Adjacent Neighbours
ind_adj <- neighbour_index(W)

## Perform mND considering sets of 3-neighbors - Non-Windows
mND_score <- mND(Xs, ind_adj, k = 3, cores = 2)

## Perform mND considering sets of 3-neighbors - Windows
mND_score <- mND(Xs, ind_adj, k = 3)

mND_score <- signif_assess(mND_score)

## ==========================================================================
## GeneSurrounder (Trial)
## ==========================================================================

## Preprocessing
ge <- X0[,2]
plot(sort(ge))
ggplot(data.frame(ge = ge), aes(ge)) + geom_density()
ge_filtered <- ge[ge >= 10] # gives 438 genes, roughly 2 hours to run GeneSurrounder

int_network <- graph_from_adjacency_matrix(A, mode = "undirected")
ge_dist <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")

ge_cor <- calcCorMatrix(data.frame(ge_filtered), corMethod = "pearson", exprName = "ge_filtered", useMethod = "everything")

diam <- diameter(int_network)

genes.assayedETnetwork <- intersect(rownames(ge_dist), rownames(ge_cor))

## Resampling
#ge_boot <- bootstraps(as.vector(ge), times = 1000)

ge_resampled <- data.frame(matrix(ncol = length(ge_filtered), nrow = 1000))
colnames(ge_resampled) <- rownames(data.frame(ge_filtered))
set.seed(123)
for(i in 1:1000){
  ge_resampled[i,] <- sample(ge_filtered, length(ge_filtered), replace = TRUE)
}

## GeneSurrounder
gs_results <- data.frame()
time <- vector()

for(i in 1:length(ge_filtered)){
  start_time <- Sys.time()
  print(paste("Run", i))
  
  gs <- geneNIDG(distance.matrix = ge_dist, 
               cor.matrix = ge_cor, 
               geneStats.observed = ge_filtered,
               perm.geneStats.matrix = as.matrix(ge_resampled), 
               genes.assayedETnetwork = genes.assayedETnetwork, 
               diameter = diam, # diameter >= 8 # diameter < 8 gives an error due to geneid.d line 376 of GeneSurrounder.R
               num.Sphere.resamples = 1000, 
               gene.id = rownames(data.frame(ge_filtered))[i] 
               )
  
  print(gs[which.min(gs$p.Fisher),])
  gs_results <- rbind(gs_results, gs[which.min(gs$p.Fisher),])
  end_time <- Sys.time()
  time[i] <- end_time - start_time
  print(paste("Time:", time[i], ""))
}
gs_results <- gs_results %>% mutate(time = time)
write.csv(gs_results, "gs_result_1.csv")

ggplot(gs_results, aes(p.Sphere)) +
  geom_density()
ggplot(gs_results, aes(p.Fisher)) +
  geom_density()
ggplot(gs_results, aes(time)) +
  geom_density()

## ==========================================================================
## GeneSurrounder on mND (Trial)
## ==========================================================================

mND_sig_genes <- rownames(mND_score$mND)[mND_score$mND$mNDp <= 0.05]

ge_mND_filtered <- ge[mND_sig_genes]
ge_mND_filtered <- ge_mND_filtered[ge_mND_filtered != 0]

ge_cor_mND <- calcCorMatrix(ge_mND_filtered, corMethod = "pearson", exprName = "ge_mND_filtered", useMethod = "everything")

ge_resampled_mND <- data.frame(matrix(ncol = length(ge_mND_filtered), nrow = 1000))
colnames(ge_resampled_mND) <- names(ge_mND_filtered)
set.seed(123)
for(i in 1:1000){
  ge_resampled_mND[i,] <- sample(ge_mND_filtered, length(ge_mND_filtered), replace = TRUE)
}

gs_new_results <- run_geneSurrounder(distance.matrix = ge_dist, 
                                     cor.matrix = ge_cor_mND, 
                                     geneStats.observed = ge_mND_filtered,
                                     perm.geneStats.matrix = as.matrix(ge_resampled_mND),
                                     diameter = diam, 
                                     num.Sphere.resamples = 1000, 
                                     gene.id = names(ge_mND_filtered),
                                     decay_only = TRUE
)

ggplot(gs_new_results, aes(p.Decay)) +
  geom_density()
ggplot(gs_new_results, aes(radius)) +
  geom_bar()
ggplot(gs_new_results, aes(time)) +
  geom_density()
sum(gs_new_results$time)/3600

######### P values need to be based on ND score not mND (excluding top k neighbors)
NDp <- data.frame(calc_p(Xs))

NDp_filtered_genes <- rownames(NDp %>% filter(L1 <= 0.05 | L2 <= 0.05))
ge_NDp_filtered <- X0[NDp_filtered_genes,2]
ge_NDp_filtered <- ge_NDp_filtered[ge_NDp_filtered > 0]
ge_cor_NDp <- calcCorMatrix(ge_NDp_filtered, corMethod = "pearson", exprName = "ge_NDp_filtered", useMethod = "everything")

ge_resampled_NDp <- data.frame(matrix(ncol = length(ge_NDp_filtered), nrow = 1000))
colnames(ge_resampled_NDp) <- names(ge_NDp_filtered)
set.seed(123)
for(i in 1:1000){
  ge_resampled_NDp[i,] <- sample(ge_NDp_filtered, length(ge_NDp_filtered), replace = TRUE)
}

gs_NDp_filtered_results <- run_geneSurrounder(distance.matrix = ge_dist, 
                                     cor.matrix = ge_cor_NDp, 
                                     geneStats.observed = ge_NDp_filtered,
                                     perm.geneStats.matrix = as.matrix(ge_resampled_NDp),
                                     diameter = diam, 
                                     num.Sphere.resamples = 1000, 
                                     gene.id = names(ge_NDp_filtered),
                                     decay_only = TRUE
)

ggplot(gs_NDp_filtered_results, aes(p.Decay)) +
  geom_density()
ggplot(gs_NDp_filtered_results, aes(radius)) +
  geom_bar()
ggplot(gs_NDp_filtered_results, aes(time)) +
  geom_density()
sum(gs_NDp_filtered_results$time)/60
### These results can then be used to adjust top k neighbors in mND

## ==========================================================================
## mND on GeneSurrounder
## ==========================================================================

# Step 1: GeneSurrounder

X0_new <- data.frame(X0) %>%
  mutate(sum = L1 + L2) %>%
  filter(sum > 0) %>%
  select(-sum)

ge_no_0 <- data.frame(ge[ge > 0])
ge_cor <- calcCorMatrix(ge_no_0, corMethod = "pearson", exprName = "ge_no_0", useMethod = "everything")

ge_resampled <- data.frame(matrix(ncol = length(ge[ge > 0]), nrow = 1000))
colnames(ge_resampled) <- names(ge[ge > 0])
set.seed(123)
for(i in 1:1000){
  ge_resampled[i,] <- sample(ge[ge > 0], length(ge[ge > 0]), replace = TRUE)
}

gs_results_4 <- run_geneSurrounder(distance.matrix = ge_dist, 
                                    cor.matrix = ge_cor, 
                                    geneStats.observed = ge[ge > 0],
                                    perm.geneStats.matrix = as.matrix(ge_resampled),
                                    diameter = diam, 
                                    num.Sphere.resamples = 1000, 
                                    gene.id = names(ge[ge > 0][3001:4000]),
                                    decay_only = TRUE,
                                    file_name = "gs_all_results_4.csv"
                                   ) # Took too much time, running chunks, didn't run underneath

X0_adjusted <- X0_new %>% mutate(L2 = -log10(gs_results_all$p.Decay)*L2)

W_new <- normalize_adj_mat(A[rownames(A)%in%rownames(ge_no_0),colnames(A)%in%rownames(ge_no_0)])
X0_perm_new <- perm_X0(X0_adjusted, r = 50, W_new, seed_n = 2)

Xs_new <- ND(X0_perm_new, W_new, cores = 2)

ind_adj_new <- neighbour_index(W_new)

mND_score_new <- mND(Xs_new, ind_adj_new, k=3, cores = 2)

mND_score_new <- signif_assess(mND_score_new)


######### NOTES ##########
# unique(idx) takes so much time i couldn't wait for it. How to subset the matrix by conditions to give a smaller matrix (not a list!)
# Check GeneSurrounder.R line 325
# Error thrown: Error in cor.fk(abs(ge[igenes.names]), igenes.distances) : 
#                       x and y must have same length.
#### UPDATE 1
# GS only runs on 1 gene at a time, only outputs p.Sphere and p.Decay not p.Fisher
# Created for loop to run gs on each gene, calculate p.Fisher, and keep results for min(p.Fisher) only
# -2(ln(p.Sphere) + ln(p.Decay)) does not give a p-value. It gives X2 distribution with 4 degrees of freedom
# How to convert to p-value?
#### UPDATE 2
# Converted to p-value with pchisq in GeneSurrounder.R
# geneNIDG() now outputs p.Fisher
# GeneSurrounder requires scores per sample, modified it to calculate p.Decay only with L2 values
# Still time intensive, couldn't run genesurrounder on all genes to have -log10(p.Decay) as input to mND. The code is written in any case
# Tried filtering genes by mNDp <= 0.05 to shortlist genesurrounder input. Also filter for gene expression > 0 --> 3753 genes estimated to take 5.2 hours
# Currently running
#### UPDATE 3
# Sample output was generated (gs_result_NDp_filter.csv with genes filtered based on mND scores)
# Fixed so that genes can be filtered based on ND p-values (using calc_p() from calc_p.R taken from mND package)
# Sample output for genes filtered by GE >= 10 are in gs_result_GE_filter.csv (done without parallelization)
# Added parallelization to run_genesurrounder() to run on all genes (last chunk)
# Couldn't execute because system storage got full
#### UPDATE 4
# run_genesurrounder() fixed to output df
# Function was moved to genesurrounder/run_geneSurrounder.R, in case you are using Windows set cores = 1
# currently running it on 1000 genes at a time to get complete output
# Saved ge_resampled.csv so that we use same resampled data for all genes
# Couldn't push to repo due to size > 100 mb, tell me if you need it and I'll send or set.seed(123)
# Seeds on R were problematic before, but output of ge_resampled seems consitent here
# gs_results_first_3000.csv saved (6397 genes remaining at a rate of 16.25 sec/gene on average --> estimated 4hrs 30mins/1000 genes using 2 cores)
# --> a bit less than 29 hours of runtime remaining
#### UPDATE 5
# gs_results_first_4000.csv saved (5397 genes remaining at an updated rate of 16.4 sec/gene on average --> estimated 4hrs 30mins/1000 genes using 2 cores)
# --> a bit more than 24 hours 30 minutes of runtime remaining

## next
## Hopefully to run the remaining genes
## Countdown: 4000/9397
