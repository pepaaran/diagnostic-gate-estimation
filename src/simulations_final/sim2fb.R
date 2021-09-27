# Simulation script for the second model variable and facebook network, passing parameters as a config file
# This scripts works for linear or dynamic response and pi=0.2 or 0.5 
# NOTE: for dynamic response, the real GATE is estimated locally when doing the analysis

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

# Define changing parameters: network, prob, 
source(as.character(commandArgs(trailingOnly=TRUE)))

output_file <- paste0("sim2",response, "_", network, "_pi", prob, ".csv")

load(paste0(network, ".RData"))

covariate_fns_for_response=list(
  frac_nbh = fraction_trt_nbrs,
  frac_nbh2 = fraction_trt_nbrs2,
  frac_nbh2w = fraction_trt_nbrs2_weighted
)

vars = list(
  vars1 = c('frac_nbh', 'frac_nbh2'),
  vars2 = c('frac_nbh', 'frac_nbh2w')
)

if(response == "linear") params = purrr::cross(list(
  beta0 = list(c(0, 0.5, 0, 0), c(0, 0.5, 0.1, 0), c(0, 0.5, 0.3, 0)),
  beta1 = list(c(1, 0.5, 0, 0), c(1, 0.5, 0.1, 0), c(1, 0.5, 0.3, 0)),
  noise_sd = list(0.1)
))
if(response == "dynamic") params = purrr::cross(list(
  b_intercept = 0,
  b_direct = 1,
  b_spill = c(0, 0.5, 1),
  max_t = c(2, 4),
  is_probit = FALSE,
  noise_sd = c(0.5)
))


#### Start where we left
if(output_file %in% list.files("results_final")){
  result <- read.csv(file=paste0("results_final/", output_file), header = F)
  init <- nrow(result)/nsim_save + 1
  rm(result)
}else{ init <- 1 }

#### Simulation definition ####

run_sim = function(param, g, covariate_fns_for_response, vars, prob, nperm, response, As, Nf, vf, glistw, pid = NULL) {
  # param: row of parameters
  
  # create list with covariates
  data = generate_covariate_data(g, covariate_fns_for_response, prob) # also generates treatment assignment
  
  # generate response
  if(response == "linear") data$y = linear_response(data$w, data$x_obs, param)
  if(response == "dynamic") data$y = dynamic_time_response(data$w, g, param)
  
  # Baseline model residuals
  res0 = linear_residuals(data, vars = c("frac_nbh"))
  
  # Geary's test
  geary0 = geary.mc(res0, listw = glistw, nsim = nperm)$p.value < 0.05
  
  # Athey's test with random or epsnet based focal units and Geary's
  F_epsnet <- random_focals_epsnet(g, eps = 3)
  F_random <- random_focals(g, Nf, eps = 3) # Aronow's estimate Nf makes sense for SUTVA testing but not so much for higher order interference
  
  vars2 = sapply(vars, function(y) y[2])
  tests = sapply(vars2,
                 function(x){
                   res = linear_residuals(data, c("frac_nbh", x))
                   other = vars2 [- which(vars2 == x)]
                   c(geary = geary.mc(res, listw = glistw, nsim =nperm)$p.value < 0.05,
                     athey_epsnet = cond_perm_test(res0, w = data$w, g = g, focal = F_epsnet, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm,
                                                   A = As[[x]], conf.level = 0.05 ),
                     athey_random = cond_perm_test(res0, w = data$w, g = g, focal = F_random, covariate_fn = covariate_fns_for_response[[x]],
                                                   method = "spearman", alternative = "greater", nperm = nperm,
                                                   A = As[[x]], conf.level = 0.05 ),
                     athey_compare_epsnet1 = cond_perm_test(res, w = data$w, g = g, focal = F_epsnet, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]]]],
                                                            method = "spearman", alternative = "greater", nperm = nperm,
                                                            A = As[[other[[1]]]], conf.level = 0.05 ),
                     athey_compare_random1 = cond_perm_test(res, w = data$w, g = g, focal = F_random, 
                                                            covariate_fn = covariate_fns_for_response[[other[[1]]]],
                                                            method = "spearman", alternative = "greater", nperm = nperm,
                                                            A = As[[other[[1]]]], conf.level = 0.05 ))
                 })
  
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

# Precompute adjacency matrix and products
adj = as_adj(g)
adj2 = adj %*% adj
diag(adj2) <- 0 # remove ego as two-step neighbor

As <- list(
  frac_nbh2 = adj2 > 0,
  frac_nbh2w = adj2,
  frac_nbh2o = (adj2 > 0) & (!adj) 
)
rm(adj); rm(adj2)

glistw <- nb2listw(neig2nb(neig(edges = get.edgelist(g))), zero.policy = TRUE)  # if there are no neighbors we asign zero weights

# Optimal number of focal units
v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
Nf <- which(v == max(v))/mean(degree(g)); rm(v)

# Precompute variance
vf = lapply(vars, function(x){
  tryCatch({
    precompute_variance(
      g, 
      list(frac_nbh = fraction_trt_nbrs,
           frac_nbh2 = covariate_fns_for_response[[x[2]]]), 
      n_boot_reps = 200, 
      n_cores = n_cores, prob = prob
    )
  }, error = function(e) NULL)
})

vf$vars0 = precompute_variance(
  g,
  list(frac_nbh = fraction_trt_nbrs),
  n_boot_reps = 200,
  n_cores = n_cores, prob = prob
)

# Trial simulation
run_sim(params[[1]], g, covariate_fns_for_response, vars, prob, nperm, response, As, Nf, vf, glistw, pid = 1)


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
                                run_sim(param, g, covariate_fns_for_response, vars, prob, nperm, response, As, Nf, vf, glistw, pid = i)
                              } %>% data.frame
          write.table(estimates, file=paste0("results_final/", output_file), append=TRUE, col.names=FALSE, row.names=FALSE, sep=',')
          print(proc.time())
        }