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
## Get GeneSurrounder.R from github and Source it
source("genesurrounder/GeneSurrounder.R")
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
library(rsample)
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

## Perform mND considering sets of 3-neighbors
mND_score <- mND(Xs, ind_adj, k=3)

mND_score <- signif_assess(mND_score)

## ==========================================================================
## GeneSurrounder
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

## mND on GeneSurrounder
X0_new <- data.frame(X0) %>%
  mutate(sum = L1 + L2) %>%
  filter(sum > 0) %>%
  select(-sum)

ge_new <- X0_new %>% filter(L2 > 0) %>% select(L2)
ge_cor_new <- calcCorMatrix(setNames(ge_new$L2, rownames(ge_new)), corMethod = "pearson", exprName = "ge_new", useMethod = "everything")

run_geneSurrounder <- function(distance.matrix, 
                               cor.matrix, 
                               geneStats.observed, 
                               perm.geneStats.matrix,
                               diameter,
                               num.Sphere.resamples = 1,
                               gene.id,
                               decay_only = TRUE,
                               file_name = "gs_results.csv",
                               cores = 2){
  genes.assayedETnetwork <- intersect(rownames(ge_dist), rownames(cor.matrix))
  gs_results <- data.frame()
  time <- vector()
  if(!decay_only){
    parallel::mclapply(1:length(gene.id), function(i){
    start_time <- Sys.time()
    print(paste("Run", i))
    
    gs <- geneNIDG(distance.matrix = distance.matrix, 
                   cor.matrix = cor.matrix, 
                   geneStats.observed = geneStats.observed,
                   perm.geneStats.matrix = perm.geneStats.matrix, 
                   genes.assayedETnetwork = genes.assayedETnetwork, 
                   diameter = diam, # diameter >= 8 # diameter < 8 gives an error due to geneid.d line 376 of GeneSurrounder.R
                   num.Sphere.resamples = num.Sphere.resamples, 
                   gene.id = gene.id[i] 
    )
    
    print(gs[which.min(gs$p.Fisher),])
    gs_results <- rbind(gs_results, gs[which.min(gs$p.Fisher),])
    end_time <- Sys.time()
    time[i] <- end_time - start_time
    print(paste("Time:", time[i], ""))
  }, mc.cores = cores)
  gs_results <- gs_results %>% mutate(time = time)
  write.csv(gs_results, file_name)
  gs_results
  }else{
    gs_results <- parallel::mclapply(1:length(gene.id), function(i){
      start_time <- Sys.time()
      print(paste("Run", i))
      distances <- distance.matrix[gene.id[i],
                                   genes.assayedETnetwork]
      
      
      sizes <- vapply(1:diameter,function(RADIUS){
        
        igenes.distances <- distances[distances <= RADIUS
                                      & distances > 0]
        
        length(igenes.distances)
        
        
      },
      numeric(1))
      observed.tau_b <- Observed.DecayDE(distance.matrix,
                                         gene.id[i],
                                         genes.assayedETnetwork,
                                         diameter,
                                         geneStats.observed)
      
      
      null.tau_b <- Resample.DecayDE(distance.matrix,
                                     gene.id[i],
                                     genes.assayedETnetwork,
                                     diameter,
                                     perm.geneStats.matrix,
                                     sizes)
      
      
      num.Decay.resamples <- nrow(perm.geneStats.matrix)
      
      # proportion of null taub \LEQ observed i.e more discordant 
      p.Decay <- vapply(1:diameter,function(distance){
        
        observed.p <- observed.tau_b[distance]
        
        null.p <- null.tau_b[distance,1:num.Decay.resamples]
        
        return(length(null.p[null.p <= observed.p])/length(null.p))
        
      },
      numeric(1))
      
      p.Decay[p.Decay == 0] <- 1/(num.Decay.resamples+1)
      
      gs <- data.frame(gene.id = rep(gene.id[i],diameter),
                 radius = 1:diameter,
                 size = sizes,
                 observed.tau_b = observed.tau_b,
                 p.Decay = p.Decay)
      
      print(gs[which.min(gs$p.Decay),])
      gs_results <- gs[which.min(gs$p.Decay),]
      #gs_results <- rbind(gs_results, gs[which.min(gs$p.Decay),])
      end_time <- Sys.time()
      #time[i] <- end_time - start_time
      gs_results$time <- end_time - start_time
      print(paste("Time:", time[i], ""))
      return(gs_results)
    }, mc.cores = cores)
    #gs_results <- gs_results %>% mutate(time = time)
    write.csv(gs_results, file_name)
    return(gs_results)
  }
}

ge_resampled_new <- data.frame(matrix(ncol = nrow(ge_new), nrow = 1000))
colnames(ge_resampled_new) <- rownames(ge_new)
set.seed(123)
for(i in 1:1000){
  ge_resampled_new[i,] <- sample(ge_new[,1], nrow(ge_new), replace = TRUE)
}

gs_new_results <- run_geneSurrounder(distance.matrix = ge_dist, 
                                     cor.matrix = ge_cor_new, 
                                     geneStats.observed = setNames(ge_new$L2, rownames(ge_new)),
                                     perm.geneStats.matrix = as.matrix(ge_resampled_new),
                                     diameter = 8, 
                                     num.Sphere.resamples = 1000, 
                                     gene.id = rownames(ge_new),
                                     decay_only = TRUE
                                     ) # Took too much time didn't run underneath

  X0_adjusted <- X0_new %>% mutate(L2 = -log10(gs_new_results$p.Decay))
  
  W_new <- normalize_adj_mat(A[rownames(A)%in%rownames(ge_new),colnames(A)%in%rownames(ge_new)])
  X0_perm_new <- perm_X0(X0_adjusted, r = 50, W_new, seed_n = 2)
  
  Xs_new <- ND(X0_perm_new, W_new, cores = 2)
  
  ind_adj_new <- neighbour_index(W_new)
  
  mND_score_new <- mND(Xs_new, ind_adj_new, k=3, cores = 2)
  
  mND_score_new <- signif_assess(mND_score_new)

## GeneSurrounder on mND
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

### GS on all genes (trial)
ge_no_0 <- data.frame(ge[ge > 0])
ge_cor <- calcCorMatrix(ge_no_0, corMethod = "pearson", exprName = "ge_no_0", useMethod = "everything")

ge_resampled <- data.frame(matrix(ncol = length(ge[ge > 0]), nrow = 1000))
colnames(ge_resampled) <- names(ge[ge > 0])
set.seed(123)
for(i in 1:1000){
  ge_resampled[i,] <- sample(ge[ge > 0], length(ge[ge > 0]), replace = TRUE)
}

gs_results <- run_geneSurrounder(distance.matrix = ge_dist, 
                                              cor.matrix = ge_cor, 
                                              geneStats.observed = ge[ge > 0],
                                              perm.geneStats.matrix = as.matrix(ge_resampled),
                                              diameter = diam, 
                                              num.Sphere.resamples = 1000, 
                                              gene.id = names(ge[ge > 0]),
                                              decay_only = TRUE,
                                              file_name = "gs_all_results.csv"
)
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

## next
## Someone should try to run the last chunk
## CSV output for run_genesurrounder() should be fixed because mclapply() for parallelization returns list (has to be df)

