library(tidyverse)
library(MOFA2)
# setup reticulate with correct python install
library(reticulate)
use_python(system("which python3.7",intern=T), required=TRUE)

# Q1
#a) read data
Data.Clin = read_tsv("medical_genomics_course/Practical2/Data/Data.Clin.txt")
RNA = read.table("medical_genomics_course/Practical2/Data/RNA.tsv",h=T,row.names=1)
DNAMeth_promoter = read.table("medical_genomics_course/Practical2/Data/DNAMeth_promoter.tsv",h=T,row.names=1)
DNAMeth_genebody = read.table("medical_genomics_course/Practical2/Data/DNAMeth_genebody.tsv",h=T,row.names=1)
DNAMeth_enhancer = read.table("medical_genomics_course/Practical2/Data/DNAMeth_enhancer.tsv",h=T,row.names=1)

#b) 
head(Data.Clin)
head(RNA)
dim(RNA)
head(DNAMeth_promoter)
dim(DNAMeth_promoter)
head(DNAMeth_genebody)
dim(DNAMeth_genebody)
head(DNAMeth_enhancer)
dim(DNAMeth_enhancer)

#Q2
##a) densities
ggplot( data=tibble(value=unlist(RNA)), mapping = aes(x=value) ) + geom_histogram()
ggplot( data=tibble(value=unlist(DNAMeth_enhancer)), mapping = aes(x=value) ) + geom_histogram()
ggplot( data=tibble(value=unlist(DNAMeth_genebody)), mapping = aes(x=value) ) + geom_histogram()
ggplot( data=tibble(value=unlist(DNAMeth_promoter)), mapping = aes(x=value) ) + geom_histogram()
### not really Gaussian but not too widespread and globally symmetric
##b) create mofa object
mofa_untrained = create_mofa( list(RNA=as.matrix(RNA),Meth_pro=as.matrix(DNAMeth_promoter),Meth_bod=as.matrix(DNAMeth_genebody),Meth_enh=as.matrix(DNAMeth_enhancer)) )

#Q3
##a)
opt_train = get_default_training_options(mofa_untrained)
opt_train$convergence_mode ="slow"
opt_model = get_default_model_options(mofa_untrained)
opt_model$num_factors = 5

#b)
mofa_untrained = prepare_mofa(mofa_untrained,model_options = opt_model,training_options = opt_train)

#c) 
mofa_trained = run_mofa(mofa_untrained)

#Q4
##a) 
plot_variance_explained(mofa_trained)

##b)
plot_factors(mofa_trained)

##c) 
plot_factor_cor(mofa_trained)

## clean session
rm(mofa_untrained)
rm(list = c("DNAMeth_enhancer",  "DNAMeth_genebody",  "DNAMeth_promoter","RNA","file","opt_model","opt_train"))
gc()
