# This script allows to use the facebook100 dataset available at https://archive.org/details/oxford-2005-facebook-matrix

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages 
source('setup.R')

# Read directly from MatLab
data <- readMat("../data/facebook100/facebook100/Caltech36.mat")

# str(data) # check structure of dataset
net <- graph_from_adjacency_matrix(data$A)

# Check for connected components
components <- components(net, mode = "strong")

# Select only largest connected component
vert_ids <- V(net)[components$membership == which.max(components$csize)]
g <- as.undirected(induced_subgraph(net, vert_ids))


save(g, file = "../data/netCaltech.RData")


#### Check size of networks ####
size_networks = rbindlist(lapply(list.files("../data/facebook100/facebook100")[-c(22, 23, 58)], function(f){
  data.frame(name = f, size = nrow(readMat(paste0("../data/facebook100/facebook100/", f))$A))
}))
