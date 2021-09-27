# Simulation script for the second model variable, linear response, prob = 0.2 with frac2 vs num2
# Parameters in cofig: network and network_size

source("setup.R")

source('functions/data_generators.R')
source('functions/covariate_functions.R')
source('functions/response_functions.R')
source('functions/existing_estimators.R')
source('functions/proposed_estimators.R')
source('functions/variance_estimators.R')
source('functions/precompute_matrices.R')
source('functions/permutation_tests.R')

# Define static parameters
n_cores = 48
nsim = 100 # number of simulations of the data
nsim_save = 20
nperm = 500 # number of permutations for permutation tests
prob = 0.2

# Define changing parameters: network, prob, 
source(as.character(commandArgs(trailingOnly=TRUE)))

output_file <- paste0("sim2linear02new", network, network_size, ".csv")

covariate_fns_for_response=list(
  frac_nbh = fraction_trt_nbrs,
  frac_nbh2 = fraction_trt_nbrs2_old,
  num_nbh2 = number_trt_nbrs2_old
)

vars = list(
  vars1 = c('frac_nbh', 'frac_nbh2'),
  vars2 = c('frac_nbh', 'num_nbh2')
)

params = purrr::cross(list(
  beta0 = list(c(0, 0.5, 0, 0), c(0, 0.5, 0.1, 0), c(0, 0.5, 0.3, 0)),
  beta1 = list(c(1, 0.5, 0, 0), c(1, 0.5, 0.1, 0), c(1, 0.5, 0.3, 0)),
  noise_sd = list(0.1, 0.25)
))

if(network == "sm") generate_network <- function(n) sample_smallworld(dim = 1, 
                                                                      size = n, # number of nodes will be size^dim
                                                                      nei = 10,
                                                                      p = 0.1) # probability of relinking to a different node)
if(network == "pl") generate_network <- function(n){
  net = sample_fitness_pl(no.of.nodes = n, no.of.edges = (10*n), exponent.out = 3)
  components = components(net)
  if(components$no > 1){
    vert_ids <- V(net)[components$membership == which.max(components$csize)]
    induced_subgraph(net, vert_ids)
  } # Avoid having unconnected nodes
  else net
} 

#### Start where we left
if(output_file %in% list.files("results_final")){
  result <- read.csv(file=paste0("results_final/", output_file), header = F)
  init <- nrow(result)/nsim_save + 1
  rm(result)
}else{ init <- 1 }

#### Simulation definition ####

run_sim = function(param, covariate_fns_for_response, vars, prob, nperm, generate_network, network_size, pid = NULL) {
  # param: row of parameters
  
  g = generate_network(network_size)
  glistw <- nb2listw(neig2nb(neig(edges = get.edgelist(g))), zero.policy = TRUE)  # if there are no neighbors we asign zero weights
  # Precompute adjacency matrix and products
  adj = as_adj(g)
  adj2 = adj %*% adj
  # diag(adj2) <- 0 # remove ego as two-step neighbor
  rm(adj)
  
  # create list with covariates
  data = generate_covariate_data(g, covariate_fns_for_response, prob) # also generates treatment assignment
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param)
  
  # Baseline model residuals
  res0 = linear_residuals(data, vars = c("frac_nbh"))
  
  
  # Geary's test
  geary0 = geary.mc(res0, listw = glistw, nsim = nperm)$p.value < 0.05
  
  # Athey's test with random or epsnet based focal units and Geary's
  # Optimal number of focal units
  v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
  Nf <- which(v == max(v))/mean(degree(g)); rm(v) # Modified this
  
  F_epsnet <- random_focals_epsnet(g, eps = 3)
  F_random <- random_focals(g, Nf, eps = 3) # Aronow's estimate Nf makes sense for SUTVA testing but not so much for higher order interference
  
  vars2 = sapply(vars, function(y) y[2])
  tests = sapply(vars2,
                 function(x){
                   res = linear_residuals(data, c("frac_nbh", x))
                   other = vars2 [- which(vars2 == x)]
                   c(geary = geary.mc(res, listw = nb2listw(neig2nb(neig(edges = get.edgelist(g)))), nsim =nperm)$p.value < 0.05,
                     athey_epsnet = cond_perm_test(res0, w = data$w, g = g, focal = F_epsnet, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm,
                                                   adj2 = adj2, conf.level = 0.05 ),
                     athey_random = cond_perm_test(res0, w = data$w, g = g, focal = F_random, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm,
                                                   adj2 = adj2, conf.level = 0.05 ),
                     athey_compare_epsnet1 = cond_perm_test(res, w = data$w, g = g, focal = F_epsnet, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]]]],
                                                            method = "spearman", alternative = "greater", nperm = nperm,
                                                            adj2 = adj2, conf.level = 0.05 ),
                     athey_compare_random1 = cond_perm_test(res, w = data$w, g = g, focal = F_random, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]]]],
                                                            method = "spearman", alternative = "greater", nperm = nperm,
                                                            adj2 = adj2, conf.level = 0.05 ))
                 })
  
  # Precompute variance
  vf = lapply(vars, function(x){
    tryCatch({
      precompute_variance(
        g, 
        list(frac_nbh = fraction_trt_nbrs,
             frac_nbh2 = covariate_fns_for_response[[x[2]]]), 
        n_boot_reps = 200, 
        n_cores = n_cores
      )
    }, error = function(e) NULL)
  })
  
  vf$vars0 = precompute_variance(
    g,
    list(frac_nbh = fraction_trt_nbrs),
    n_boot_reps = 200,
    n_cores = n_cores
  )
  
  return(data.frame(
    pid=pid,
    
    adj0=data %>% linear_adjustment(vars=c("frac_nbh")),
    adj1=data %>% linear_adjustment(vars=vars$vars1),
    adj2=data %>% linear_adjustment(vars=vars$vars2),

    var0 = if(is.null(vf$vars0)) NA else data %>% linear_variance_estimate(vf$vars0, vars=c("frac_nbh")),
    var1 = if(is.null(vf$vars1)) NA else data %>% linear_variance_estimate(vf$vars1, vars=vars$vars1),
    var2 = if(is.null(vf$vars2)) NA else data %>% linear_variance_estimate(vf$vars2, vars=vars$vars2),

    geary0=geary0,
    geary1=tests[1,1],
    geary2=tests[1,2],

    athey1_epsnet = tests[2,1],
    athey2_epsnet = tests[2,2],

    athey1_random = tests[3,1],
    athey2_random = tests[3,2],

    athey12_epsnet = tests[4,1],
    athey21_epsnet = tests[4,2],

    athey12_random = tests[5,1],
    athey21_random = tests[5,2]
    
  ))
}

run_sim(params[[1]], covariate_fns_for_response, vars, prob, nperm, generate_network, network_size, pid = 1)


#### Run simulation ####
registerDoParallel(cores=n_cores)
print('Running simulation...')

start = proc.time()
foreach(k = init:(length(params)*(nsim/nsim_save)), .combine=rbind,
        .export = lsf.str(envir = .GlobalEnv, all.names = TRUE),
        .packages = c('dplyr', 'broom', 'spdep', 'ade4', 'igraph')) %do% {
          # Set seed for reproducibility, even if Euler kicks you out
          set.seed(k)
          i = ceiling(k/(nsim/nsim_save))
          param = params[[i]]
          print(unlist(param))
          estimates = foreach(rep = 1:nsim_save, .combine=rbind,
                              .export = lsf.str(envir = .GlobalEnv, all.names = TRUE),
                              .packages = c('dplyr', 'broom', 'spdep', 'ade4', 'igraph')) %dopar% {
                                run_sim(param, covariate_fns_for_response, vars, prob, nperm, generate_network, network_size, pid = i)
                              } %>% data.frame
          write.table(estimates, file=paste0("results_final/", output_file), append=TRUE, col.names=FALSE, row.names=FALSE, sep=',')
          print(proc.time())
        }