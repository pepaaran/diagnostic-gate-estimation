# Simulation script for the first model variable and prob = 0.5 (without exposure mappings) and dynamic response
# Parameters: network type and size

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
prob = 0.5

# Define changing parameters: network size, prob, response
source(as.character(commandArgs(trailingOnly=TRUE)))

output_file <- paste0("sim1dynamic05", network, network_size,".csv")

covariate_fns_for_response=list(
  frac_nbh = fraction_trt_nbrs,
  num_nbh = number_trt_nbrs
)

vars = list(
  vars1 = c('frac_nbh'),
  vars2 = c('num_nbh')
)

params = purrr::cross(list(
  b_intercept = 0,
  b_direct = 1,
  b_spill = c(0, 0.5, 1),
  max_t = c(2, 4),
  is_probit = FALSE,
  noise_sd = c(0.5, 1)
))

# Real GATE empirically computed for each simulated network and response parameters
compute_real_gate <- function(g, param, n_reps = 100){
  foreach(rep = 1:n_reps, .combine=rbind, .inorder=FALSE,
          .export = lsf.str(envir = .GlobalEnv, all.names = TRUE)) %dopar% {
            mean( dynamic_time_response(w=rep(1, vcount(g)), g, param)) - mean( dynamic_time_response(w=rep(0, vcount(g)), g, param))
          } %>% mean()
}

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
  glistw <- nb2listw(neig2nb(neig(edges = get.edgelist(g))), zero.policy = TRUE)  
  # if there are no neighbors we asign zero weights but other zero policies can be used for the considered interference features
  
  # create list with covariates
  data = generate_covariate_data(g, covariate_fns_for_response, prob) # also generates treatment assignment
  
  # generate response
  data$y = dynamic_time_response(data$w, g, param)
  
  # Baseline model residuals
  res0 <- data$y - mean(data$y[data$w==0])
  res0[data$w==1] <- data$y[data$w==1] - mean(data$y[data$w==1])
  res0 <- as.vector(res0)
  
  # Geary's test
  geary0 = geary.mc(res0, listw = glistw, nsim = nperm)$p.value < 0.05
  
  # Athey's test with random or epsnet based focal units and Geary's
  # Optimal number of focal units
  v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
  Nf <- which(v == max(v)); rm(v)
  
  F_epsnet <- random_focals_epsnet(g, eps = 2)
  F_random <- random_focals(g, Nf, eps = 2) # Aronow's estimate Nf makes sense for SUTVA testing but not so much for higher order interference
  
  tests = sapply(vars,
                 function(x){
                   res = linear_residuals(data, x)
                   other = vars [- which(vars == x)]
                   c(geary = geary.mc(res, listw = glistw, nsim =nperm)$p.value < 0.05,
                     athey_epsnet = cond_perm_test(res0, w = data$w, g = g, focal = F_epsnet, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm, conf.level = 0.05),
                     athey_random = cond_perm_test(res0, w = data$w, g = g, focal = F_random, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm, conf.level = 0.05),
                     athey_compare_epsnet1 = cond_perm_test(res, w = data$w, g = g, focal = F_epsnet, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]] ]],
                                                            method = "spearman", alternative = "greater", nperm = nperm, conf.level = 0.05),
                     athey_compare_random1 = cond_perm_test(res, w = data$w, g = g, focal = F_random, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]] ]],
                                                            method = "spearman", alternative = "greater", nperm = nperm, conf.level = 0.05))
                 })
  
  # Precompute variance
  vf = lapply(vars, function(x){
    tryCatch({
      precompute_variance(
        g, 
        list(frac_nbh = covariate_fns_for_response[[x]]), 
        n_boot_reps = 200, 
        n_cores = n_cores,
        prob = prob
      )
    }, error = function(e) NULL)
  })
  
  return(data.frame(
    pid=pid,
    gate = compute_real_gate(g, param, n_reps = 100), 
    
    dm=data %>% difference_in_means,
    adj1=data %>% linear_adjustment(vars=vars$vars1),
    adj2=data %>% linear_adjustment(vars=vars$vars2),

    var0 = dm_variance_estimate(data),
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


# Trial simulation run
run_sim(params[[1]], covariate_fns_for_response, vars, prob, nperm, generate_network, network_size, pid = 1)


#### Run simulation ####
registerDoParallel(cores=n_cores)
print('Running simulation...')

start = proc.time()
# results <- data.frame()
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
          
          # results <- rbind(results, estimates)
          print(proc.time())
        }