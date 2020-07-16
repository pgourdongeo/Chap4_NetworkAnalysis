## Test insert python chunk for community detection

library(reticulate)
use_python("/opt/python/3.6/bin/python3.6")

virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "cmake")
virtualenv_install("r-reticulate", "boost")


import("biSBM")
