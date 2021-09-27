# setup R packages
dependencies <- c(
  'igraph', 'data.table', 'expm', 'R.matlab', 'coin', 'pbmcapply', 'corrplot', 'mc2d', 'spdep',
  'ade4', 'gam', 'foreign', 'foreach', 'doParallel', 'dplyr', 'tidyr', 'broom', 'devtools', 'data.table',
  'plot.matrix', 'RColorBrewer', 'optparse'
)
new_packages <- dependencies[!(dependencies %in% installed.packages()[, "Package"])]

if (length(new_packages)) install.packages(new_packages)

# load packages
lapply(dependencies, require, character.only = TRUE)

rm(dependencies, new_packages)

# Install Puelz's library to test causal effects under general interference
# install_github("dpuelz/CliqueRT")
# library(CliqueRT)

# Set global option to avoid error "zero length neighbour vector"
set.ZeroPolicyOption(TRUE)
