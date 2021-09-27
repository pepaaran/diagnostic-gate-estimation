# Defines commonly used covariate generators

# Fraction of treated neighbors
fraction_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w) / degree(g)}

# Number of treated neighbors
number_trt_nbrs = function(g, w) {as.vector(as_adj(g) %*% w)}

# Fraction of treated individuals in the 2-hop neighborhood.
fraction_trt_nbrs2 = function(g, w, A = NULL) {
  if(is.null(A)){
    adj = as_adj(g)
    adj2 = adj %*% adj
    diag(adj2) <- 0 # remove ego as two-step neighbor
    A <- adj2 > 0 # transform into indicator
  }
  # This allows to pass adj2 as argument and avoid repeated computations
  as.vector(A %*% w) / apply(A, 2, sum) 
}

fraction_trt_nbrs2_weighted = function(g, w, A = NULL) {
  if(is.null(A)){
    adj = as_adj(g)
    A = adj %*% adj
    diag(A) <- 0 # remove ego as two-step neighbor
  }
  # This allows to pass adj2 as argument and avoid repeated computations
  as.vector(A %*% w) / apply(A, 2, sum) 
}

fraction_trt_nbrs2_only <- function(g, w, A = NULL){
  if(is.null(A)){
    adj = as_adj(g)
    adj2 = adj %*% adj
    diag(adj2) <- 0
    A = ((adj2 > 0) - adj) > 0
  }
  as.vector(A %*% w) / apply(A, 2, sum)
}

# Fraction of treated neighbors in the 3-hop nieghborhood
fraction_trt_nbrs3 = function(g, w, adj3 = NULL) {
  if(is.null(adj3)){
    adj = as_adj(g)
    adj3 = adj %*% adj %*% adj
  }
  
  as.vector(adj3 %*% w) / apply(adj3, 2, sum)
}

number_trt_nbrs2 = function(g, w, adj2 = NULL) {
  if(is.null(adj2)){
    adj = as_adj(g)
    adj2 = adj %*% adj
  }
  # This allows to pass adj2 as argument and avoid repeated computations
  as.vector(adj2 %*% w)
}

# exposure_trt_nbrs = function(g, w){
#   # Four different exposure models: 00, 01, 10, 00
#   num_nebrs = as.vector(as_adj(g) %*% w)
#   exposure = 1 * (num_nebrs > 0)
#   paste0(w, exposure)
# }
exposure_trt_nbrs = function(g, w){
  1 * ( as.vector(as_adj(g) %*% w) > 0)
}

qfrac = function(g, w, q = 0.75){
  frac <- as.vector(as_adj(g) %*% w) / degree(g)
  as.vector(frac > q)
} # Turning this into 0 and 1 I can easily fit the regression model afterwards, but it doesn't allow for more than 2 categories

# The old functions still count the ego
fraction_trt_nbrs2_old = function(g, w, adj2 = NULL) {
  if(is.null(adj2)){
    adj = as_adj(g)
    adj2 = adj %*% adj
  }
  as.vector(adj2 %*% w) / apply(adj2, 2, sum) 
}

number_trt_nbrs2_old = function(g, w, adj2 = NULL) {
  if(is.null(adj2)){
    adj = as_adj(g)
    adj2 = adj %*% adj
  }
  as.vector(adj2 %*% w)
}