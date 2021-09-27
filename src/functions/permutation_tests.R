# Residual permutation test functions

linear_residuals <- function(data, vars=NULL){
  if (is.null(vars)) vars = names(data$x_obs)
  w = data$w
  y0 = data$y[w==0]
  y1 = data$y[w==1]
  x0 = data$x_obs %>% select(one_of(vars)) %>% filter(w==0)
  x1 = data$x_obs %>% select(one_of(vars)) %>% filter(w==1)
  
  res0 = lm(y0 ~ ., data=x0) %>% residuals
  res1 = lm(y1 ~ ., data=x1) %>% residuals
  res <- c()
  res[w==0] <- res0
  res[w==1] <- res1
  res
}

exposure_residuals <- function(data){
  w = data$w
  y0 = data$y[w==0]
  y1 = data$y[w==1]
  x0 = data$x_obs$exp_nbh[w == 0]
  x1 = data$x_obs$exp_nbh[w == 1]
  
  if(length(unique(x0)) <2){
    res0 = y0 - mean(y0)
  }else{
    res0 = lm(y0 ~ ., data=x0) %>% residuals
  }
  
  if(length(unique(x1)) <2){
    res1 = y1 - mean(y1)
  }else{
    res1 = lm(y1 ~ ., data=x1) %>% residuals
  }
  
  res <- c()
  res[w==0] <- res0
  res[w==1] <- res1
  res
}

geary_C_test <- function(resid, g, nsim = 1000){
  g_nb <- neig2nb(neig(edges=get.edgelist(g)))
  
  t = geary.test(resid,
                listw = nb2listw(g_nb))
  c(C = t$estimate[1],
    pvalue = t$p.value,
    mc_pvalue = geary.mc(resid, listw = nb2listw(g_nb), nsim = nsim)$p.value)
}

geary_C_test_nbh2<- function(resid, g, nsim = 1000){
  g2 <- graph_from_adjacency_matrix(as_adj(g) %*% as_adj(g) > 0, mode = 'undirected')
  g_nb <- neig2nb(neig(edges=get.edgelist(g2)))

  t = geary.test(resid,
                 listw = nb2listw(g_nb))
  c(C = t$estimate[1],
    pvalue = t$p.value,
    mc_pvalue = geary.mc(resid, listw = nb2listw(g_nb), nsim = nsim)$p.value)
}

geary_C_test_nbh2only<- function(resid, g, nsim = 1000){
  g2 <- graph_from_adjacency_matrix((as_adj(g) %*% as_adj(g) - as_adj(g)) > 0, mode = 'undirected')
  g_nb <- neig2nb(neig(edges=get.edgelist(g2)))

  t = geary.test(resid,
                 listw = nb2listw(g_nb))
  c(C = t$estimate[1],
    pvalue = t$p.value,
    mc_pvalue = geary.mc(resid, listw = nb2listw(g_nb), nsim = nsim)$p.value)
}

geary_C_test_A2<- function(resid, g, nsim = 1000){
  g2 <- graph_from_adjacency_matrix(as_adj(g) %*% as_adj(g), mode = 'undirected')
  g_nb <- neig2nb(neig(edges=get.edgelist(g2)))
  
  t = geary.test(resid,
                 listw = nb2listw(g_nb))
  c(C = t$estimate[1],
    pvalue = t$p.value,
    mc_pvalue = geary.mc(resid, listw = nb2listw(g_nb), nsim = nsim)$p.value)
} # This is the same test as before but without normalizing the matrix to 0 and 1 

.cor_permuted <- function(W, g, res_focal, covariate_fn, focal, fixed, method, prob_treat = 0.5, ...){
  # focal is a vector indicating which nodes are focal
  # fixed is a logical vector indicating which nodes are fixed for the randomization (includes the focals)
  w <- W
  w[!fixed] <- rbinom(sum(!fixed), size=1, prob = prob_treat) # draw from the conditional distribution
  
  x_obs = covariate_fn(g, w, ...)
  cor(x = x_obs[focal], y = res_focal, method = method, use = "complete.obs")
}

cond_perm_test <- function(res, w, focal, g, covariate_fn, method, alternative = "greater", 
                           prob = 0.5, nperm = 1000, seed = NULL, conf.level = NULL, ...){
  # Conditional permutation test from Athey et al. 2018
  # focal is a list with elements focal (indicating the focal nodes) and fixed (indicating the 
  # neighboring nodes over which we don't randomize)
  # method: indicates which correlation method to use("spearman" or "pearson")
  # conf.level: indicates the confidence level for which we accept or reject (if NULL the whole test is returned, 
  #             otherwise only TRUE for rejecting and FALSE for accepting)
  # ... is used for passing the matrix as an argument and speeding up simulations

  # NOTE: maybe I should pass focal and fixed as vectors, or maybe I should compute fixed internally
  if(! is.null(seed)) set.seed(seed)
  
  # v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob) * (gorder(g)-x)*prob)
  # n_fixed <- which(v == max(v)) # This is approximate based on the treatment randomization
  # 
  # v_fixed = sample(1:gorder(g), n_fixed)
  
  res_focal <- res[focal$focals]
  
  observed <- cor(x = covariate_fn(g, w, ...)[focal$focals],
                  y = res_focal, method = method, use = "complete.obs")
  names(observed) <- "statistic"
  
  fixed <- 1:gorder(g) %in% focal$fixed
  
  foreach(i = 1:nperm, .combine=rbind,
          .export = lsf.str(envir = .GlobalEnv, all.names = TRUE),
          .packages = c('dplyr', 'broom', 'ade4', 'igraph')) %dopar% {
            
            .cor_permuted(W = w, g = g,
                          res_focal = res_focal, 
                          covariate_fn = covariate_fn, focal = focal$focals, fixed = fixed,
                          method = method, ...)
            
          } %>% as.vector() -> permutations
  
  # return(list(
  #   greater = mean(permutations > observed),
  #   less = mean(permutations < observed),
  #   cor_obs = observed,
  #   cor_perms = permutations
  # ))
  
  resrank <- sum(permutations < observed)
  names(resrank) <- "observed rank"
  
  if(alternative == "greater"){
    pval <- punif((resrank+1)/(nperm+1), lower.tail = FALSE)
  }else if(alternative == "less"){
    pval <- punif((resrank+1)/(nperm+1))
  }
  names(pval) = "p-value"
  
  if(! is.null(conf.level)) return(pval < conf.level)
  
  tres = list(
    statistic = observed,
    parameter = resrank,
    p.value = pval,
    alternative = alternative,
    method = "Conditional permutation test of interference correlation",
    data.name = data.name <- paste(deparse(substitute(res)), "\nnumber of focal units:", 
                                   length(focal$focals), "\nnumber of simulations:", 
                                   nperm , "\ninterference feature:",
                                   deparse(substitute(covariate_fn)), "\n"),
    res = permutations
  )
  
  
  class(tres) <- c("htest")
  tres
}

random_focals_epsnet <- function(g, eps = 2, seed = NULL){
  # Note: 2-net looks at neighbors (test SUTVA), 3-net and neighbors of neighbors (test CTR assumption)
  
  if(! is.null(seed)) set.seed(seed)
  
  if(eps == 2){
    D <- as.matrix(get.adjacency(g))
  } else if(eps == 3) {
    d <- as.matrix(get.adjacency(g))
    D <- d %*% d
  }else {
    d <- as.matrix(get.adjacency(g)) %^% (eps - 2)
    D <- d %*% as.matrix(get.adjacency(g))
  }
  
  unmarked <- 1:gorder(g)
  focals <- c()
  while(length(unmarked) > 0){
    j <- sample(unmarked, 1)
    
    unmarked <- unmarked[unmarked != j]
    focals <- c(focals, j)
    
    unmarked <- unmarked[! (unmarked %in% which(D[j, ] > 0))] # this indicates the vertices at distance <= eps
    
  }
  
  fixed <- focals
  if(eps > 2){
    for(f in focals){
      fixed <- c(fixed, which(d[f, ] > 0))
    }
    fixed <- unique(fixed)
  }
  
  list(focals = focals, fixed = fixed)
}

random_focals <- function(g, Nf = NULL, eps = 2, prob_treat = 0.5, seed = NULL){
  # Note: eps=2 looks at neighbors (test SUTVA),and eps=3 at neighbors of neighbors (test CTR assumption)
  
  if(! is.null(seed)) set.seed(seed)
  
  if(is.null(Nf)){
    v = sapply(1:gorder(g), function(x) x * (gorder(g)-x)*(1-prob_treat) * (gorder(g)-x)*prob_treat)
    Nf <- which(v == max(v))
  }

  focals <- sample(1:gorder(g), Nf)
  
  
  if(eps == 3) {
    d <- as.matrix(get.adjacency(g))
  }else if(eps > 3){
    d <- as.matrix(get.adjacency(g)) %^% (eps - 2)
  }
  
  fixed <- focals
  if(eps > 2){
    for(f in focals){
      fixed <- c(fixed, which(d[f, ] > 0))
    }
    fixed <- unique(fixed)
  }
  
  list(focals = focals, fixed = fixed)
}
