run_geneSurrounder <- function(distance.matrix, 
                               cor.matrix, 
                               geneStats.observed, 
                               perm.geneStats.matrix,
                               diameter,
                               num.Sphere.resamples = 1,
                               gene.id,
                               decay_only = TRUE,
                               file_name = "gs_results.csv",
                               cores = 2 # Set to 1 for Windows
                               ){
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
    gs_results <- bind_rows(gs_results)
    write.csv(gs_results, file_name)
    return(gs_results)
  }
}