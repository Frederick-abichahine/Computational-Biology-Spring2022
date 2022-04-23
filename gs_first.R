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
#Xs <- ND(X0_perm, W, cores = 2)

## Perform Network Diffusion - Windows
#Xs <- ND(X0_perm, W)
#saveRDS(Xs, "Data/Xs.rds")
data(Xs)

## Get Indices of Adjacent Neighbours
ind_adj <- neighbour_index(W)

## Perform mND considering sets of 3-neighbors - Non-Windows
mND_score <- mND(Xs, ind_adj, k = 3, cores = 2)

## Perform mND considering sets of 3-neighbors - Windows
mND_score <- mND(Xs, ind_adj, k = 3)

mND_score <- signif_assess(mND_score)

#saveRDS(mND_score, "Data/mND_scores.rds")

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

#### Step 1: GeneSurrounder ####

#X0_new <- data.frame(X0) %>%
#  mutate(sum = L1 + L2) %>%
#  filter(sum > 0) %>%
#  select(-sum)

data(X0)
ge <- X0[,2]

ge_no_0 <- data.frame(ge[ge > 0])
ge_cor <- calcCorMatrix(ge_no_0, corMethod = "pearson", exprName = "ge_no_0", useMethod = "everything")

ge_resampled <- data.frame(matrix(ncol = length(ge[ge > 0]), nrow = 1000))
colnames(ge_resampled) <- names(ge[ge > 0])
set.seed(123)
for(i in 1:1000){
  ge_resampled[i,] <- sample(ge[ge > 0], length(ge[ge > 0]), replace = TRUE)
}

gs_results_all <- run_geneSurrounder(distance.matrix = ge_dist, 
                                    cor.matrix = ge_cor, 
                                    geneStats.observed = ge[ge > 0],
                                    perm.geneStats.matrix = as.matrix(ge_resampled),
                                    diameter = diam, 
                                    num.Sphere.resamples = 1000, 
                                    gene.id = names(ge[ge > 0]), #subset for chunks
                                    decay_only = TRUE,
                                    file_name = "gs_results_all.csv"
                                   ) # Took too much time, ran in chunks

#write.csv(gs_results_all, "Data/gs_results_all.csv")

#### Step 2: mND on GeneSurrounder - Adjusted scores ####
library(lubridate)
gs_results_all <- read.csv("Data/gs_results_all.csv") %>% mutate(time = duration(time)) %>% select(-X)
rownames(gs_results_all) <- gs_results_all$gene.id

X0_gs_adjusted <- data.frame(X0) %>%
  rownames_to_column('rn')  %>%
  mutate(L2 = if_else(rn %in% gs_results_all$gene.id, -log10(gs_results_all[rn,]$p.Decay)*L2, L2)) %>%
  column_to_rownames('rn')

X0_gs_adjusted <- as.matrix(X0_gs_adjusted)
#saveRDS(X0_gs_adjusted, "Data/X0_gs_adjusted.rds")

ggplot(data.frame(X0), aes(L2)) + geom_density()
ggplot(data.frame(X0_gs_adjusted), aes(L2)) + geom_density()

data(A)
W <- normalize_adj_mat(A)
ind_adj <- neighbour_index(W)

X0_perm_new <- perm_X0(X0_gs_adjusted, r = 50, W, seed_n = 2)

Xs_new <- ND(X0_perm_new, W, cores = 2)
#saveRDS(Xs_new, "Data/Xs_new.rds")

mND_score_new_k2 <- mND(Xs_new, ind_adj, k=2, cores = 2)
#mND_score_new <- mND(Xs_new, ind_adj, k=3, cores = 2)

mND_score_new_k2 <- signif_assess(mND_score_new_k2)
#mND_score_new <- signif_assess(mND_score_new)

#saveRDS(mND_score_new_k2, "Data/mND_gs_adjusted_scores_k2.rds")
#saveRDS(mND_score_new, "Data/mND_gs_adjusted_scores.rds")


## ==========================================================================
## Results
## ==========================================================================

#Load Data
data(A)
W <- normalize_adj_mat(A)
ind_adj <- neighbour_index(W)

data(X0)
data(Xs)
X0_gs_adjusted <- readRDS("Data/X0_gs_adjusted.rds")
Xs_new <- readRDS("Data/Xs_new.rds")

#data(mND_score) #not same format as mND_score.rds (i think old format)
mND_score <- readRDS("Data/mND_scores.rds")
#mND_score_new_k2 <- readRDS("Data/mND_gs_adjusted_scores_k2.rds")
mND_score_new <- readRDS("Data/mND_gs_adjusted_scores.rds")

##### mND Alone results #####

#H1: genes with a mutation frequency greater than zero;
#H2: top 1200 differentially expressed genes (FDR < 10^-7).
#Further, we set the cardinalities of gene sets N1 and N2, containing the genes with the highest scoring neighborhoods, as |H1|=|N1| and |H2|=|N2|
Hl <- list(l1 = rownames(X0[X0[,1]>0,]), 
           l2 = names(X0[order(X0[,2], decreasing = T),2][1:1200])
)
top_Nl <- unlist(lapply(Hl, function(x) length(x)))
top_Nl

class_res <- classification(mND_score, X0, Hl, top = top_Nl)

#Classification of genes in every layer
head(class_res$gene_class)

#Occurrence of (M; L; I; NS) for each gene across layers
head(class_res$occ_labels)

plot_results(mND_score, class_res, W, n = 150, directory = "Results/mND_results/")

#Optimizing k (Mac only)
k_val <- seq(1,6,1)
k_results <- optimize_k(Xs, X0, k_val, ind_adj, W, top = 200, cores = 2)
k_results <- data.frame(k_results)
colnames(k_results) <- k_val
#write.csv(k_results, "Data/k_results.csv")
# can also be loaded through data(k_results) but is a list not data.frame

##### GS-Adjusted mND results #####


#H1: genes with a mutation frequency greater than zero;
#H2: top 1200 differentially expressed genes (FDR < 10^-7).
#Further, we set the cardinalities of gene sets N1 and N2,
#containing the genes with the highest scoring neighborhoods, as |H1|=|N1| and |H2|=|N2|

Hl_new <- list(l1 = rownames(X0[X0[,1]>0,]), 
           l2 = names(X0_gs_adjusted[order(X0_gs_adjusted[,2], decreasing = T),2][1:1200])
)
top_Nl_new <- unlist(lapply(Hl_new, function(x) length(x)))
top_Nl_new
class_res_new <- classification(mND_score_new, X0_gs_adjusted, Hl_new, top = top_Nl_new)
#class_res_new <- classification(mND_score_new, X0_gs_adjusted, Hl_new, top = top_Nl_new)

#Classification of genes in every layer
head(class_res_new$gene_class)

#Occurrence of (M; L; I; NS) for each gene across layers
head(class_res_new$occ_labels)

plot_results(mND_score_new, class_res_new, W, n = 150, directory = "Results/mND_results_new_k2/")
#plot_results(mND_score_new, class_res_new, W, n = 150, directory = "Results/mND_results_new/")

#Optimizing k (Mac only)
k_results_new <- optimize_k(Xs_new, X0_gs_adjusted, k_val, ind_adj, W, top = 200, cores = 2)
k_results_new <- data.frame(k_results_new)
colnames(k_results_new) <- k_val
#write.csv(k_results_new, "Data/k_results_new.csv")
#k = 2 seems like a better choice in this case

## ==========================================================================
## Analysis - Ghadi
## ==========================================================================

##### % Label change #####

#Sanity check:
class_res_new$gene_class <- class_res_new$gene_class[match(rownames(class_res$gene_class), rownames(class_res_new$gene_class)),]
sum(class_res_new$gene_class[,1] != class_res$gene_class[,1])

shift_sm <- table(mND = class_res$gene_class[,1], GS_adjusted_mND = class_res_new$gene_class[,1])

shift_ge <- table(mND = class_res$gene_class[,2], GS_adjusted_mND = class_res_new$gene_class[,2])

# Percentages over mND genes
results_ge_mND <- shift_ge
for(i in 1:nrow(shift_ge)){
  results_ge_mND[i,] <- (shift_ge[i,]/sum(shift_ge[i,]))*100
}
results_ge_mND
# Total selected mND genes (I, L, M)
11796 - sum(shift_ge[4,]) #2297 k = 3

# Percentages over GS_adjusted genes
results_ge_GS_adjusted <- shift_ge
for(i in 1:ncol(shift_ge)){
  results_ge_GS_adjusted[,i] <- (shift_ge[,i]/sum(shift_ge[,i]))*100
}
results_ge_GS_adjusted
# Total selected GS-adjusted mND genes (I, L, M)
11796 - sum(shift_ge[,4]) #2109 k = 3
11796 - sum(shift_ge[,4]) - (11796 - sum(shift_ge[4,])) # --> -188 genes k = 3, -205 k = 2
# Difference in M
sum(shift_ge[,3]) #GS-adjusted (314 with k = 3) --> +37 modules
sum(shift_ge[3,]) #mND (277 with k = 3) 

# Percent over all genes
(shift_ge/11796)*100

##### Network Visualization #####
library(visNetwork)
library(RCy3)
ls("package:igraph", pattern = "^layout_.")

## mND graph
genes_subset_mND <- class_res$gene_class %>% filter(L2 != "NS")
A_subset_mND <- A[rownames(A) %in% rownames(genes_subset_mND), 
                  colnames(A) %in% rownames(genes_subset_mND)]
graph_mND <- graph_from_adjacency_matrix(A_subset_mND) %>% 
             set_vertex_attr("class", value = genes_subset_mND[,2])
vis_mND <- toVisNetworkData(graph_mND)
visNetwork(nodes = vis_mND$nodes, edges = vis_mND$edges) %>%
  visIgraphLayout(layout = "layout_") %>%
  visOptions(nodesIdSelection = TRUE, selectedBy = "class") %>%
  visNodes(color = vis_mND$nodes$class)
#cytoscapePing()
#createNetworkFromIgraph(graph_mND,new.title='mND')
#cytoscape was taking time idk

#1 = Isolated, 2 = Linker, 3 = Module
## GS-adjusted mND graph
genes_subset_adjusted_mND <- class_res_new$gene_class %>% filter(L2 != "NS")
A_subset_adjusted_mND <- A[rownames(A) %in% rownames(genes_subset_adjusted_mND), 
                  colnames(A) %in% rownames(genes_subset_adjusted_mND)]
graph_adjusted_mND <- graph_from_adjacency_matrix(A_subset_adjusted_mND) %>%
                      set_vertex_attr("class", value = genes_subset_adjusted_mND[,2])
vis_adjusted_mND <- toVisNetworkData(graph_adjusted_mND)
visNetwork(nodes = vis_adjusted_mND$nodes, edges = vis_adjusted_mND$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") %>%
  visOptions(nodesIdSelection = TRUE, selectedBy = "class")

## ================================
## Running Analysis - Bilal
## ================================
#saveRDS(mND_score_new, "Data/mND_gs_adjusted_scores.rds")
data(A)
W <- normalize_adj_mat(A)
ind_adj <- neighbour_index(W)

data(X0)
data(Xs)
X0_gs_adjusted <- readRDS("Data/X0_gs_adjusted.rds")
Xs_new <- readRDS("Data/Xs_new.rds")

#data(mND_score) #not same format as mND_score.rds (i think old format)
mND_score <- readRDS("Data/mND_scores.rds")
mND_score_new_k2 <- readRDS("Data/mND_gs_adjusted_scores_k2.rds")
mND_score_new <- readRDS("Data/mND_gs_adjusted_scores.rds")

#get new and old genes
new_genes = rownames(mND_score_new$mND)[which(mND_score_new$mND$mNDp < 0.05)]
old_genes = rownames(mND_score$mND)[which(mND_score$mND$mNDp < 0.05)]

## Ratio - we have a 7% change in the genes 
#   Ghadi: I think the way you calculated it means:
# " 7.7% of the total genes were additionally significant with a 0.05 mNDp threshold "
# Ratio would be length(new_genes)/length(old_genes) (diff lengths means hard to get a single percentage to reflect change)
sum(!(new_genes %in% old_genes))/nrow(mND_score$mND)

length(new_genes)
length(old_genes)

## ================================
## Old
data(X0)
Hl_old <- list(l1 = rownames(X0[X0[,1]>0,]), 
           l2 = names(X0[order(X0[,2], decreasing = T),2][1:1200])
)
top_Nl_old <- unlist(lapply(Hl_old, function(x) length(x)))
top_Nl_old
class_res_old <- classification(mND_score, X0, Hl_old, top = top_Nl_old)
head(class_res_old$gene_class)


## New 
Hl_new <- list(l1 = rownames(X0_gs_adjusted[X0_gs_adjusted[,1]>0,]), 
               l2 = names(X0_gs_adjusted[order(X0_gs_adjusted[,2], decreasing = T),2][1:1200])
)
top_Nl_new <- unlist(lapply(Hl_new, function(x) length(x)))
top_Nl_new
class_res_new <- classification(mND_score_new, X0_gs_adjusted, Hl_new, top = top_Nl_new)
#class_res_new <- classification(mND_score_new_k2, X0_gs_adjusted, Hl_new, top = top_Nl_new)

#Classification of genes in every layer
head(class_res_new$gene_class)
##-------------------------
## After obtaining the classification of each gene
## we wish to check the cumulative_decay_score on target genes
## Target Genes are new genes that are obtained
## Cumulative Decay Score = GE_neighbor/d_neighbor_target
## Stategy:
## For ever node in the graph:
##  for all modules in the list of old modules:
##    distance = get_distance (old_module, node)
##    score[node] += Measure(old_module)/distance
##------------------------------------------------
## Now, this Measure could be taken as GE
## or the adjusted -log_10(p-value)
##-------------------------------------------------
## First lets convert the Adjacency Matrix to Distance Matrix
## install.packages("netmeta")
library(netmeta)

int_network <- graph_from_adjacency_matrix(A, mode = "undirected")
ge_dist <- calcAllPairsDistances(int_network, directionPaths = "all", networkName = "int_network")
## Second, lets get the decay effect of each gene from gs_results_all
gs_results_all <- read.csv("Data/gs_results_all.csv")
## Get List of all old nodes classification
old_class = class_res_old$gene_class[,2]
## Get list of all new nodes classification
new_class = class_res_new$gene_class[,2]
## Remember - adjusted scores of genes
X0_gs_adjusted

## initialize empty vector
results_c_decay_score = rep(0,nrow(ge_dist))
new_module = rep(0,nrow(ge_dist))
## Now start the strategy:
## For every gene

for(i in seq(1,nrow(gs_results_all))){
  ## Progress
  ## install.packages("svMisc")
  ## library(svMisc)
  ##progress(i)
  ## ----
  gene = gs_results_all[i,]
  g_id = gene$gene.id
  g_rd = gene$radius
  g_old_class = old_class[g_id]
  g_new_class = new_class[g_id]
  ## Get the measure of the module node
  ## we choose measyre as -log_10_pvalue
  ## Other choice is original GE
  #measure = X0_gs_adjusted[i,2]
  measure <- -log10(gene$p.Decay)
  ##print(measure)
  ## If the gene is an old module
  if(g_old_class == "M"){
    ## Traverse other genes in the surrounding radius
    ## and update their cumulative_decay_score
    for(j in seq(1, nrow(ge_dist))){
      if (j == i){
        next
      }
      ## Get distance between the two nodes
      g_dist = ge_dist[i,j]
      if (g_dist <= g_rd){
        ## Update cumulative score of gene_j
        results_c_decay_score[j] = as.numeric(results_c_decay_score[j]) + measure/g_dist
      }
    }
  }
  ## Not old module - check if new, if yes put 1 in the vector
  else{
    if (g_new_class == "M"){
      new_module[i] = 1
    }
  }
}

## We have 69 new modules - validate this
# Ghadi: 72 with k = 2, 69 with k = 3
sum(new_module==1)
## How many modules are in new but are not in old in old but they are in new
sum(class_res_old$gene_class[,2]=="M")
sum(class_res_new$gene_class[,2]=="M")

## Get md genes in new
new_mdg = rownames(class_res_new$gene_class)[class_res_new$gene_class[,2] == "M"]
old_mdg = rownames(class_res_old$gene_class)[class_res_old$gene_class[,2] == "M"]  
sum(new_mdg %in% old_mdg)
sum(!(old_mdg %in% new_mdg)) #32 (k=2) modules discarded by GS
## True - we have 69 genes as new modules previously not as modules
# Ghadi: 72 with k = 2, 69 k = 3

sum(!(new_mdg %in% old_mdg))
distinct_mdg = new_mdg[which(!(new_mdg %in% old_mdg))]
discarded_mdg <- new_mdg[which(!(old_mdg %in% new_mdg))]
## length(distinct_mdg)

## Indices of distinct mdg
dmdg_indices = which(rownames(X0_gs_adjusted)%in%distinct_mdg)
discarded_mdg_indices <- which(rownames(X0_gs_adjusted)%in%discarded_mdg)

## Mean Results
mean(results_c_decay_score[dmdg_indices])
mean(results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%distinct_mdg))])

mean(results_c_decay_score[discarded_mdg_indices])
mean(results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%discarded_mdg))])

## Distinct Density Plots
ggplot(as.data.frame(as.matrix(results_c_decay_score[dmdg_indices])), aes(V1)) + geom_density()
ggplot(as.data.frame(as.matrix(results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%distinct_mdg))])), aes(V1)) + geom_density()

## Overlaying the Density Plots
# Distinct
a = data.frame(cumulative_decay_score = results_c_decay_score[dmdg_indices])
b = data.frame(y = results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%distinct_mdg))])
ggplot() + 
  geom_density(data = a, aes(x = cumulative_decay_score), 
               fill = "#DAF7A6", color = "black", alpha = 0.7) + 
  geom_density(data = b, aes(x = y),
               fill = "#C70039", color = "black", alpha = 0.7)
#Discarded
x = data.frame(cumulative_decay_score = results_c_decay_score[discarded_mdg_indices])
y = data.frame(y = results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%discarded_mdg))])
ggplot() + 
  geom_density(data = x, aes(x = cumulative_decay_score), 
               fill = "#DAF7A6", color = "black", alpha = 0.7) + 
  geom_density(data = y, aes(x = y),
               fill = "#C70039", color = "black", alpha = 0.7)

## Next: 
## Significance Test between two vectors which are (distinct):
## 1. results_c_decay_score[dmdg_indices]
## 2. results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%distinct_mdg))]
## ----------------------------------------------------------------------------

vector1 <- results_c_decay_score[dmdg_indices]
vector2 <- results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%distinct_mdg))]

## Comparing Variances:
#Using var.test() in order to compare the variances
#Running a two-sided variance test on both vectors with alpha = 5% (0.05) and null hypothesis (H0) -> variances are equal
#ratio (default is 1) is the expected ratio for H0, and the alternative hypothesis (H1) displays the other side (default is two.sided)

var.test(x = vector1, y = vector2, ratio = 1, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)

#we obtain a p-value = 0.002163 with a cut-off of 0.05 which is statistically significant
#we obtain a ratio = 0.554 for x/y; This means that the variance of the distinctive genes that were selected as new modules was significantly lower than the old modules for the calculated cumulative decay score
#the 95% confidence interval is [0.406, 0.801]

## Comparing Means:
#Using t.test() in order to compare the means
#mu is either the population's expected mean or the expected difference in the means of the two samples
#mu = 0 as the null hypothesis.
#paired indicates if the data is paired or not
#var.equal indicates if both the samples are treated as equal or not equal w.r.t their variances
#the above var.equal would be based on the var.test done before (in our case = FALSE)
#since our data is unpaired => we perform the following comparison:

t.test(x = vector1, y = vector2, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#we obtain a p-value = 5.896e-09 with a cut-off of 0.05 which is statistically significant
#the 95% confidence interval is [29.82, 55.46]
#the means of vector1 and vector2 are 168.97 and 126.33 respectively.

## Next: 
## Significance Test between the other two vectors which are (discarded):
## 3. results_c_decay_score[discarded_mdg_indices]
## 4. results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%discarded_mdg))]
## ----------------------------------------------------------------------------

vector3 <- results_c_decay_score[discarded_mdg_indices]
vector4 <- results_c_decay_score[which(!(rownames(X0_gs_adjusted)%in%discarded_mdg))]

#Now we perform the exact same steps as before and analyze the output for the discarded genes...

## Comparing Variances:
var.test(x = vector3, y = vector4, ratio = 1, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)

#we obtain a p-value = 0.08254 with a cut-off of 0.05 which implies that it is not statistically significant
#we obtain a ratio = 0.605 for x/y; This means that the variance of the discarded genes that were selected as new modules was not significantly lower than the old modules for the calculated cumulative decay score
#the 95% confidence interval is [0.388, 1.070]

## Comparing Means:
t.test(x = vector3, y = vector4, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#we obtain a p-value = 5.055e-06 with a cut-off of 0.05 which is statistically significant
#the 95% confidence interval is [33.95, 74.00]
#the means of vector3 and vector4 are 180.41 and 126.44 respectively.
                           
## ----------------------------------------------------------------------------
#Sanity check:
sum(class_res_new$gene_class[,1] != class_res_old$gene_class[,1])

shift_sm <- table(mND = class_res_old$gene_class[,1], GS_adjusted_mND = class_res_new$gene_class[,1])

shift_ge <- table(mND = class_res_old$gene_class[,2], GS_adjusted_mND = class_res_new$gene_class[,2])



# Percentages over mND genes
results_ge_mND <- shift_ge
for(i in 1:nrow(shift_ge)){
  results_ge_mND[i,] <- (shift_ge[i,]/sum(shift_ge[i,]))*100
}
results_ge_mND
# Total selected mND genes (I, L, M)
11796 - sum(shift_ge[4,])

# Percentages over GS_adjusted genes
results_ge_GS_adjusted <- shift_ge
for(i in 1:ncol(shift_ge)){
  results_ge_GS_adjusted[,i] <- (shift_ge[,i]/sum(shift_ge[,i]))*100
}
results_ge_GS_adjusted
# Total selected GS-adjusted mND genes (I, L, M) --> -205 genes
11796 - sum(shift_ge[,4])
11796 - sum(shift_ge[,4]) - (11796 - sum(shift_ge[4,]))

# Percent over all genes
(shift_ge/11796)*100

## ================================
# Analyzing coverage of Breast Cancer Related Genes
## ================================

# Getting Disease-related genes with AutoSeed (includes eDGAR, DrugBank, and MalaCards)
library(Autoseed)
bc_genes_search <- AutoSeed("breast cancer")
bc_genes_search2 <- AutoSeed("breast carcinoma")
bc_genes <- Reduce(union, c(bc_genes_search$edgar,bc_genes_search$malacards, bc_genes_search$drugbank,
                            bc_genes_search2$edgar,bc_genes_search2$malacards, bc_genes_search2$drugbank))
bc_genes_clean <- unique(vapply(bc_genes, function(x){return(unlist(strsplit(x, split='::', fixed=TRUE))[1])}, c("x")))
bc_genes_clean <- bc_genes_clean[bc_genes_clean%in%rownames(X0)]
# Matched distinct genes
sum(distinct_mdg %in% bc_genes_clean)
# Matched discarded genes
sum(discarded_mdg %in% bc_genes_clean)
# Matched old M genes
sum(old_mdg %in% bc_genes_clean)
# Matched new M genes
sum(new_mdg %in% bc_genes_clean)

genes_selected_old <- rownames(class_res_old$gene_class %>% filter(L2 != "NS"))
genes_selected_new <- rownames(class_res_new$gene_class %>% filter(L2 != "NS"))

# Matched selected old
sum(genes_selected_old %in% bc_genes_clean)/length(bc_genes_clean)*100
# Matched selected new
sum(genes_selected_new %in% bc_genes_clean)/length(bc_genes_clean)*100

# Analyzing coverage by mNDp cutoffs
cutoffs <- seq(0.0001, 0.06, 0.0001)
genes_selected_cutoffs <- data.frame(mNDp_cutoffs = 1:length(cutoffs), old = 1:length(cutoffs), new = 1:length(cutoffs), new_k2 = 1:length(cutoffs))
for(i in 1:length(cutoffs)){
  genes_cutoff_old <- rownames(mND_score$mND %>% filter(mNDp <= cutoffs[i]))# %>% filter(L2 != "NS"))
  genes_cutoff_new <- rownames(mND_score_new$mND %>% filter(mNDp <= cutoffs[i]))# %>% filter(L2 != "NS"))
  genes_cutoff_new_k2 <- rownames(mND_score_new_k2$mND %>% filter(mNDp <= cutoffs[i]))
  print(paste("Cutoff:", i))
  # Matched selected by mNDp cutoff old
  old <- sum(genes_cutoff_old %in% bc_genes_clean)/length(bc_genes_clean)*100
  print("Old")
  print(old)
  # Matched selected by mNDp cutoff new
  new <- sum(genes_cutoff_new %in% bc_genes_clean)/length(bc_genes_clean)*100
  print("New (k = 3)")
  print(new)
  # Matched selected by mNDp cutoff new k = 2
  new_k2 <- sum(genes_cutoff_new_k2 %in% bc_genes_clean)/length(bc_genes_clean)*100
  print("New (k = 2)")
  print(new_k2)
  
  genes_selected_cutoffs[i,] <- c(cutoffs[i], old, new, new_k2)
}
# Interesting how we get more coverage no matter the cutoff although converges 
# slightly at xtrmly low cutoffs. k = 2 shows slightly less coverage but at lower 
# cutoffs coverage is almost the same (seems that lower k prioritizes important genes)
ggplot() +
  geom_line(data = genes_selected_cutoffs, aes(cutoffs, old), color = "black") +
  geom_line(data = genes_selected_cutoffs, aes(cutoffs, new), color = "red") +
  geom_line(data = genes_selected_cutoffs, aes(cutoffs, new_k2), color = "blue", alpha = 0.5) +
  ylab("Cumulative Coverage (%)") +
  xlab("mNDp Cutoff") +
  ggtitle("Cumulative Coverage of disease-related genes by mNDp cutoffs", 
          subtitle = "mND (k = 3) in black, GS-adjusted mND (k = 2) in blue, GS-adjusted mND (k = 3) in red")

## ================================
# Enrichment Analysis
## ================================

require(KEGGREST)

pathways.list <- keggList("pathway", "hsa") #Needs Internet connection
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
#saveRDS(genes.by.pathway, "Data/genes.by.pathway.rds")
genes.by.pathway <- readRDS("Data/genes.by.pathway.rds")

enrich <- function(mND_score, pathways.list, pathway.codes, genes.by.pathway)
{
  geneList <- mND_score$mND$mNDp
  names(geneList) <- rownames(mND_score$mND)

  # Wilcoxon test for each pathway
  pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                             }
  ))
  # Assemble output table
  outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
  outdat$pathway.name <- pathways.list[paste("path:",outdat$pathway.code, sep = "")]
  outdat$p.value <- pVals.by.pathway[,"p.value"]
  outdat$Annotated <- pVals.by.pathway[,"Annotated"]
  return(outdat[order(outdat$p.value),])
}

enrichment_old <- enrich(mND_score, pathways.list, pathway.codes, genes.by.pathway)
enrichment_new <- enrich(mND_score_new, pathways.list, pathway.codes, genes.by.pathway)
#write.csv(enrichment_old, "Results/enrichment_old.csv")
#enrichment_old <- read.csv("Results/genes.by.pathway.csv")
#write.csv(enrichment_new, "Results/enrichment_new.csv")
#enrichment_new <- read.csv("Results/genes.by.pathway.csv")

comparative_enrich <- function(X0, mND_scores = list(), pathways.list, pathway.codes, genes.by.pathway){
  result <- list()
  for(i in 1:length(genes.by.pathway)){
    pathway_result <- data.frame(matrix(nrow = length(seq (50, 500, 1)), 
                                        ncol = length(mND_scores)+2))
    colnames(pathway_result) <- append(names(mND_scores), c("DE_scores","cutoff"))
    count <- 0
    for(mND_score in mND_scores){
      count <- count + 1
      count2 <- 0
      genes <- mND_score$mND %>% arrange(desc(mND))
      control_genes <- sort(X0[,2], decreasing = TRUE) 
      for(cutoff in seq (50, 500, 1)){
        count2 <- count2 + 1
        pathway_result[count2,count] <- sum(rownames(genes[1:cutoff,]) %in% genes.by.pathway[[i]])/length(genes.by.pathway[[i]])*100
        pathway_result[count2,ncol(pathway_result)] <- cutoff
        
        if(count == 1){
          pathway_result[count2,ncol(pathway_result)-1] <- sum(names(control_genes[1:cutoff]) %in% genes.by.pathway[[i]])/length(genes.by.pathway[[i]])*100
        }
      }
    }
    result[[pathways.list[[i]]]] <- pathway_result
  }
  return(result)
}

comparative_enrichment <- comparative_enrich(X0, list(mND_old_k3 = mND_score, mND_new_k3 = mND_score_new, mND_new_k2 = mND_score_new_k2), pathways.list, pathway.codes, genes.by.pathway)

bc_pathways <- c("Breast cancer - Homo sapiens (human)", 
                 "Epstein-Barr virus infection - Homo sapiens (human)", 
                 "Pathways in cancer - Homo sapiens (human)", 
                 "Viral carcinogenesis - Homo sapiens (human)")
other_pathways <- c("Cell cycle - Homo sapiens (human)",
                    "FoxO signaling pathway - Homo sapiens (human)",
                    "Hippo signaling pathway - Homo sapiens (human)",
                    "Oocyte meiosis - Homo sapiens (human)",
                    "p53 signaling pathway - Homo sapiens (human)",
                    "PI3K-Akt signaling pathway - Homo sapiens (human)",
                    "Progesterone-mediated oocyte maturation - Homo sapiens (human)",
                    "Proteasome - Homo sapiens (human)",
                    "Rap1 signaling pathway - Homo sapiens (human)",
                    "Ras signaling pathway - Homo sapiens (human)",
                    "Sphingolipid signaling pathway - Homo sapiens (human)",
                    "TGF-beta signaling pathway - Homo sapiens (human)")

for(pathway in append(bc_pathways, other_pathways)){

  plot <- ggplot(data = comparative_enrichment[[pathway]]) +
            geom_line(aes(cutoff, mND_old_k3), color = "black") +
            geom_line(aes(cutoff, mND_new_k3), color = "red") +
            geom_line(aes(cutoff, mND_new_k2), color = "blue") +
            geom_line(aes(cutoff, DE_scores), color = "magenta") +
            ylab("Cumulative Coverage (%)") +
            xlab("Cutoff of sorted gene list") +
            ggtitle(pathway) 
              #subtitle = "Black: mND (k = 3), Blue: GS-adjusted mND (k = 2), Red: GS-adjusted mND (k = 3), Magenta: original DE scores")
  print(plot)
}

## ================================
# Visualization of M & L (Looks like neighbors are not included in class M)
## ================================

## mND Graph
library(visNetwork)
genes_subset_mND <- class_res_old$gene_class %>% filter(!(L2 %in% c("NS", "I")))
A_subset_mND <- A[rownames(A) %in% rownames(genes_subset_mND), 
                  colnames(A) %in% rownames(genes_subset_mND)]
graph_mND <- graph_from_adjacency_matrix(A_subset_mND) %>% 
  set_vertex_attr("class", value = genes_subset_mND[,2])
vis_mND <- toVisNetworkData(graph_mND)
visNetwork(nodes = vis_mND$nodes, edges = vis_mND$edges) %>%
  visIgraphLayout(layout = "layout_nicely") %>%
  visOptions(nodesIdSelection = TRUE, selectedBy = "class") %>%
  visNodes(color = vis_mND$nodes$class)

## GS-adjusted mND graph
genes_subset_adjusted_mND <- class_res_new$gene_class %>% filter(!(L2 %in% c("NS", "I")))
A_subset_adjusted_mND <- A[rownames(A) %in% rownames(genes_subset_adjusted_mND), 
                           colnames(A) %in% rownames(genes_subset_adjusted_mND)]
graph_adjusted_mND <- graph_from_adjacency_matrix(A_subset_adjusted_mND) %>%
  set_vertex_attr("class", value = genes_subset_adjusted_mND[,2])
vis_adjusted_mND <- toVisNetworkData(graph_adjusted_mND)
visNetwork(nodes = vis_adjusted_mND$nodes, edges = vis_adjusted_mND$edges) %>%
  visIgraphLayout(layout = "layout_nicely") %>%
  visOptions(nodesIdSelection = TRUE, selectedBy = "class")

## ================================
## ================================
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
#### UPDATE 6
# gs_results_first_5000.csv saved (4397 genes remaining at an updated rate of 16.65 sec/gene on average --> estimated 4hrs 36mins/1000 genes using 2 cores)
# --> a bit more than 20 hours 20 minutes of runtime remaining
#### UPDATE 7
# gs_results_first_6000.csv saved (3397 genes remaining at an updated rate of 16.62 sec/gene on average --> estimated 4hrs 36mins/1000 genes using 2 cores)
# --> a bit more than 15 hours 30 minutes of runtime remaining
#### UPDATE 8
# gs_results_first_7000.csv saved (2397 genes remaining at an updated rate of 16.54 sec/gene on average --> estimated 4hrs 36mins/1000 genes using 2 cores)
# --> a bit more than 11 hours 00 minutes of runtime remaining
#### UPDATE 9
# gs_results_first_8000.csv saved (1397 genes remaining at an updated rate of 16.47 sec/gene on average --> estimated 4hrs 30mins/1000 genes using 2 cores)
# --> a bit more than 6 hours 36 minutes of runtime remaining
#### UPDATE 10
# gs_results_all.csv saved --> rate: 16.39 secs/gene using 2 cores
## Countdown: 9397/9397
# Ran mND with initial and gs adjusted scores (saved as mND_scores.rds and mND_gs_adjusted_scores.rds in Data folder)
# We have to see how to validate now: visualizations, metrics, maybe checking if genes that were enhanced by gs were new modules in mND
#### UPDATE 11 
# Added results for mND alone (k = 3) and GS adjusted mND (k = 2, k = 3) as indicated by k optimization
# created .rds objects or .csv files for needed info to be loaded to run results, saved plots as well
#### UPDATE 12
# Added shift from class to class percentages
# Added network visualization for selected genes

#### next ####
# Evaluation
# Comparison with mND alone
# 1. We could probably check % transition from class to another (Isolated, Linker, Module, NotSignificant)
# 2. Dr. Hanna, Frederick, and Ghadi discussed quantifying correlation between mND score of genes not selected 
#    before GS (but selected after; lets say decayed genes) and connection of a gene to nodes labeled as M by mND alone (w/out GS)
#    We discussed doing this at many distance thresholds to see if the correlation drops when M genes are further
#    We would do mND score / distance to normalize (or maybe -log10(mNDp)/distance) and add up for all M genes that
#    contain the gene of interest in their radius and are at a distance == threshold to gene of interest
#    This score is calculated for decayed genes and we check correlation to mND adjusted (or -log10(mNDp adjusted))
#    --> From here we can either take cut-off for this score (mND score / distance) to consider a gene decayed
#        maybe trying many distance thresholds where genes selected are distance =< threshold and we check 
#        where correlation stabilizes?? we need to try it and see the trend
#       OR
#        We could just leave the score as is (mND score / distance) to describe that the gene is 
#        influenced by decay this much
# 3. We could visualize the 1200 genes for mND alone and GS-adjusted mND colored but I, L, M (3 colors)
#    We could check for connectivity and overall topology
