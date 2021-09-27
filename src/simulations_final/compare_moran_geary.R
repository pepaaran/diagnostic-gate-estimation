# Simulation exploring Geary's C vs Moran's I in the simplest model

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
nsim = 1000 # number of simulations of the data
nsim_save =100 # (with 100 it runs out of memory with mem=20000 per core)
nperm = 500 # number of permutations for permutation tests
prob = 0.5

output_file <- "explore_geary_vs_moran_1000.csv"

load("netCaltech.RData")

covariate_fns_for_response=list(
  frac_nbh = fraction_trt_nbrs#,  num_nbh = number_trt_nbrs
)

params = purrr::cross(list(
  beta0 = list(c(0, 0), c(0, 0.2), c(0, 0.5)),
  beta1 = list(c(1, 0), c(1, 0.2), c(1, 0.5)),
  noise_sd = list(0.1, 0.25)
))



#### Simulation definition ####

run_sim = function(param, g, covariate_fns_for_response, prob, nperm, pid = NULL) {
  # param: row of parameters
  
  # create list with covariates
  data = generate_covariate_data(g, covariate_fns_for_response, prob) # also generates treatment assignment
  glistw = nb2listw(neig2nb(neig(edges = get.edgelist(g))))
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param)

  # Baseline model residuals
  res0 <- data$y - mean(data$y[data$w==0])
  res0[data$w==1] <- data$y[data$w==1] - mean(data$y[data$w==1])
  res0 <- as.vector(res0)
  
  # Geary's test
  geary0 = geary.mc(res0, listw = glistw, nsim = nperm)$p.value
  moran0 = moran.mc(res0, listw = glistw, nsim = nperm)$p.value
  
  res1 = linear_residuals(data, vars = "frac_nbh")
  
  # Geary's test
  geary1 = geary.mc(res1, listw = glistw, nsim = nperm)$p.value
  moran1 = moran.mc(res1, listw = glistw, nsim = nperm)$p.value
  
  
  return(data.frame(
    pid=pid,
    
    dm=data %>% difference_in_means,
    adj1=data %>% linear_adjustment(vars="frac_nbh"),
    
    geary0=geary0,
    geary1=geary1,
    
    moran0 = moran0,
    moran1 = moran1
  ))
}

# Trial simulation run
run_sim(params[[1]], g, covariate_fns_for_response, prob, nperm, pid = 1)


#### Run simulation ####
registerDoParallel(cores=n_cores)
print('Running simulation...')

start = proc.time()
# results <- data.frame()
foreach(k = 1:(length(params)*(nsim/nsim_save)), .combine=rbind,
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
                                run_sim(param, g, covariate_fns_for_response, prob, nperm, pid = i)
                              } %>% data.frame
          write.table(estimates, file=paste0("results_final/", output_file), append=TRUE, col.names=FALSE, row.names=FALSE, sep=',')
          
          # results <- rbind(results, estimates)
          print(proc.time())
        }