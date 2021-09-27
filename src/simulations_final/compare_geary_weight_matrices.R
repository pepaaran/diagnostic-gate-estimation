# Simulation comparing three specifications of the weight matrix for Geary's test of 2-step interference

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
nsim_save =50 # (with 100 it runs out of memory with mem=20000 per core)
nperm = 500 # number of permutations for permutation tests
prob = 0.5

output_file <- "explore_geary_matrices_1000.csv"

load("netCaltech.RData")

covariate_fns_for_response=list(
  frac_nbh = fraction_trt_nbrs,
  frac_nbh2 = fraction_trt_nbrs2#,  num_nbh = number_trt_nbrs
)

params = purrr::cross(list(
  beta0 = list(c(0, 0.5, 0), c(0, 0.5, 0.1), c(0, 0.5, 0.3)),
  beta1 = list(c(1, 0.5, 0), c(1, 0.5, 0.1), c(1, 0.5, 0.3)),
  noise_sd = list(0.1)
))



#### Simulation definition ####

run_sim = function(param, glistw, g2listw, covariate_fns_for_response, prob, nperm, pid = NULL) {
  # param: row of parameters
  
  # create list with covariates
  data = generate_covariate_data(g, covariate_fns_for_response, prob) # also generates treatment assignment
  
  # generate response
  data$y = linear_response(data$w, data$x_obs, param)
  
  # Baseline model residuals
  # res0 <- data$y - mean(data$y[data$w==0])
  # res0[data$w==1] <- data$y[data$w==1] - mean(data$y[data$w==1])
  # res0 <- as.vector(res0)
  
  # Geary's test
  # geary0 = geary.mc(res0, listw = glistw, nsim = nperm)$p.value
  # geary20 = geary.mc(res0, listw = g2listw, nsim = nperm)$p.value
  
  res1 = linear_residuals(data, vars = "frac_nbh")
  res2 = linear_residuals(data, vars = c("frac_nbh", "frac_nbh2"))
  
  # Geary's test
  geary1 = geary.mc(res1, listw = glistw, nsim = nperm)$p.value
  geary21 = geary.mc(res1, listw = g2listw, nsim = nperm)$p.value
  geary2 = geary.mc(res2, listw = glistw, nsim = nperm)$p.value
  geary22 = geary.mc(res2, listw = g2listw, nsim = nperm)$p.value
  
  
  return(data.frame(
    pid=pid,
    
    adj1=data %>% linear_adjustment(vars="frac_nbh"),
    adj2=data %>% linear_adjustment(vars=c("frac_nbh", "frac_nbh2")),
    
    geary1=geary1,
    geary2=geary2,
    
    geary21 = geary21,
    geary22 = geary22
  ))
}

# Weight matrices in spdep format
glistw = nb2listw(neig2nb(neig(edges = get.edgelist(g))))
g2 <- graph_from_adjacency_matrix(as_adj(g) %*% as_adj(g), mode = 'undirected')
g2listw <- nb2listw(neig2nb(neig(edges=get.edgelist(g2))))

# Trial simulation run
run_sim(params[[1]], glistw, g2listw, covariate_fns_for_response, prob, nperm, pid = 1)

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
                                run_sim(param, glistw, g2listw, covariate_fns_for_response, prob, nperm, pid = i)
                              } %>% data.frame
          write.table(estimates, file=paste0("results_final/", output_file), append=TRUE, col.names=FALSE, row.names=FALSE, sep=',')
          
          # results <- rbind(results, estimates)
          print(proc.time())
        }