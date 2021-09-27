# Computes Gamma and Delta for variance estimation

precompute_variance = function(g, covariate_fns, n_boot_reps, n_cores, prob = 0.5) {
  registerDoParallel(cores=n_cores) 
  moments = foreach(i = 1:n_boot_reps,
                    .packages = 'igraph',
                    .export = lsf.str(envir = .GlobalEnv, all.names = TRUE)) %dopar% {
    data = generate_covariate_data(g, covariate_fns, prob_treat = prob)
    w = data$w
    
    x0 = cbind(1, data$x_obs[w==0,]) %>% as.matrix
    x1 = cbind(1, data$x_obs[w==1,]) %>% as.matrix
    
    # return sample covariances
    list(
      S1inv = solve(t(x1) %*% x1),
      S0inv = solve(t(x0) %*% x0)
    )
  }
  
  B = length(moments)

  S0inv = foreach(m = moments) %do% {m$S0inv} %>% (function(M) {Reduce(`+`, M) / B})
  S1inv = foreach(m = moments) %do% {m$S1inv} %>% (function(M) {Reduce(`+`, M) / B})

  
  data = generate_covariate_data(g, covariate_fns, prob_treat = prob)
  omega0 = c(1, colMeans(data$x_ctrl))
  omega1 = c(1, colMeans(data$x_trt))
  
  (t(omega0) %*% S0inv %*% omega0 + t(omega1) %*% S1inv %*% omega1) %>% as.vector
}
