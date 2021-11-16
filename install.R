## General dependencies:
install.packages("tidyverse")

## BioC dependencies
install.packages("BiocManager")
BiocManager::install(c("MOFA2"))
system("pip install mofapy2")
library(reticulate)
use_python(system("which python",intern=T), required=TRUE)
