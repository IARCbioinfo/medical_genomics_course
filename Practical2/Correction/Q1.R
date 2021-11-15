library(tidyverse)

# read data
Data.Clin = read_tsv("medical_genomics_course/Practical2/Data/Data.Clin.txt")
RNA = read_tsv("medical_genomics_course/Practical2/Data/RNA.tsv")
DNAMeth_promoter = read_tsv("medical_genomics_course/Practical2/Data/DNAMeth_promoter.tsv")
DNAMeth_genebody = read_tsv("medical_genomics_course/Practical2/Data/DNAMeth_genebody.tsv")
DNAMeth_enhancer = read_tsv("medical_genomics_course/Practical2/Data/DNAMeth_enhancer.tsv")

