# Simulations for the number of focal and fixed units 

source("setup.R")

source('functions/data_generators.R')
source('functions/covariate_functions.R')
source('functions/response_functions.R')
source('functions/existing_estimators.R')
source('functions/proposed_estimators.R')
source('functions/variance_estimators.R')
source('functions/precompute_matrices.R')
source('functions/permutation_tests.R')


generate_network_sm <- function(n) sample_smallworld(dim = 1, 
                                                    size = n, # number of nodes will be size^dim
                                                    nei = 10,
                                                    p = 0.1) # probability of relinking to a different node)
generate_network_pl <- function(n){
  net = sample_fitness_pl(no.of.nodes = n, no.of.edges = (10*n), exponent.out = 3)
  components = components(net)
  if(components$no > 1){
    vert_ids <- V(net)[components$membership == which.max(components$csize)]
    induced_subgraph(net, vert_ids)
  } # Avoid having unconnected nodes
  else net
} 

runsim = function(n, gen, prob){
  g = gen(n)
  
  v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
  Nf <- which(v == max(v)); rm(v)
  
  
  F_epsnet2 <- random_focals_epsnet(g, eps = 2)
  F_random2 <- random_focals(g, Nf, eps = 2)
  F_epsnet3 <- random_focals_epsnet(g, eps = 3)
  F_random3 <- random_focals(g, Nf/mean(degree(g)), eps = 3)
  
  c(Nfoc_epsnet2 = length(F_epsnet2$focals),
    Nfoc_random2 = length(F_random2$focals),
    Naux_epsnet2 = n - length(F_epsnet2$fixed),
    Naux_random2 = n - length(F_random2$fixed),
    
    Nfoc_epsnet3 = length(F_epsnet3$focals),
    Nfoc_random3 = length(F_random3$focals),
    Naux_epsnet3 = n - length(F_epsnet3$fixed),
    Naux_random3 = n - length(F_random3$fixed))
}

params = purrr::cross(list(
  n = c(100, 1000),
  gen = c(generate_network_sm, generate_network_pl),
  prob = c(0.5, 0.2)
))

estimates = sapply(params, function(x){
  foreach(rep = 1:100, .combine=rbind,
          .export = lsf.str(envir = .GlobalEnv, all.names = TRUE),
          .packages = c('dplyr', 'broom', 'spdep', 'ade4', 'igraph')) %do% {
            runsim(x$n, x$gen, x$prob)
          } %>% data.frame %>% apply(2, mean)
})


load("../data/netCaltech.RData")
prob = 0.5
v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
Nf <- which(v == max(v)); rm(v)
n = 762

estimates_caltech = foreach(rep = 1:100, .combine=rbind,
                            .export = lsf.str(envir = .GlobalEnv, all.names = TRUE),
                            .packages = c('dplyr', 'broom', 'spdep', 'ade4', 'igraph')) %do% {
                              F_epsnet2 <- random_focals_epsnet(g, eps = 2)
                              F_random2 <- random_focals(g, Nf, eps = 2)
                              F_epsnet3 <- random_focals_epsnet(g, eps = 3)
                              F_random3 <- random_focals(g, Nf/mean(degree(g)), eps = 3)
                              
                              c(Nfoc_epsnet2 = length(F_epsnet2$focals),
                                Nfoc_random2 = length(F_random2$focals),
                                Naux_epsnet2 = n - length(F_epsnet2$fixed),
                                Naux_random2 = n - length(F_random2$fixed),
                                
                                Nfoc_epsnet3 = length(F_epsnet3$focals),
                                Nfoc_random3 = length(F_random3$focals),
                                Naux_epsnet3 = n - length(F_epsnet3$fixed),
                                Naux_random3 = n - length(F_random3$fixed))
                              
                            } %>% data.frame %>% apply(2, mean)
