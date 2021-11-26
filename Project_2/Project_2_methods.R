
### Project 2 methods ###

### 1. Download IDAT files, sample sheet and script


### 2. Process the data into usable m and beta tables using the scripts provided here: https://github.com/IARCbioinfo/Methylation_analysis_scripts/blob/master/Methylation_pre-processing.R
# Did any samples fail quality control? 
# Did any samples have an incorrect sex label?
# How many probes were dropped at each step? 
# What is the size of the resulting dataset?
# Do you notice any effects on the data from technical or clinical variables?


### 3. Create datasets for analysis
## a) separate methylation probes into gene body, promoter, and enhancer regions
# Read in beta and m-tables, use data.table::fread() with option data.table=FALSE for faster loading

# Filter the beta and m tables to retain only probes with an absolute beta value difference of > 0.1
beta.diff <- apply(beta, 1, max, na.rm=TRUE) - apply(beta, 1, min, na.rm=TRUE)
names(beta.diff) <- rownames(beta)
length(beta.diff) 
beta.diff <- beta.diff[which(beta.diff > 0.1)]
length(beta.diff) 
# How many probes remain? What proportion of the original probes passed this filter?
m_filt <- m[which(rownames(m) %in% names(beta.diff)), ]
beta_filt <- beta[which(rownames(beta) %in% names(beta.diff)), ]

# Load EPIC array annotation file and label probes by promoter, enhancer and gene body status
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann$class <- sapply(1:nrow(ann), function(x)
  if(ann$Regulatory_Feature_Group[x] == "Promoter_Associated"){"Promoter"}
  else if(ann$Phantom5_Enhancers[x] != ""){"Enhancer"}
  else if(grepl("Body|Exon", ann$UCSC_RefGene_Group[x])){"Body"}
  else if(ann$Regulatory_Feature_Group[x] != "Promoter_Associated" & ann$Phantom5_Enhancers[x] == "" & !grepl("Body|Exon", ann$UCSC_RefGene_Group[x])){NA})

# Subset the m_filt tables into three: m_pro, m_enh, and m_bod using the class labels for each probe you generated above

## b) select the 5,000 most variable sites in each and combine into a single dataset
# Shrink the m_xxx tables to retain only the 5,000 most variable probes, create matching beta_xxx tables, follow the method shown here for m_pro
vv <- apply(m_pro,1,var)
cv <- cumsum(sort(vv,decreasing = T))/sum(vv)
m_pro.red <- m_pro[order(match(rownames(m_pro),names(cv))),]
m_pro.red <- m_pro.red[1:5000, ] 
rm(vv, cv)
beta_pro.red <- beta_filt[rownames(beta_filt) %in% rownames(m_pro.red), ]

# combine into single datasets
df_list <- list(m_pro.red, m_enh.red, m_bod.red)
m_red <- Reduce(function(x,y) rbind(x,y), df_list)
rm(df_list)
df_list <- list(beta_pro.red, beta_enh.red, beta_bod.red)
beta_red <- Reduce(function(x,y) rbind(x,y), df_list)


### 4. Perform consensus clustering
# Generate consensus clusters
library(ConsensusClusterPlus)
d_temp <- as.matrix(m_red)
d <- sweep(d_temp,1,apply(d_temp,1,median,na.rm=T))
dir <- # set an ouptut directory
results <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1, title=dir,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

# Create a data frame of samples and their clusters
Meth_clusters <- as.data.frame(colnames(m_red))
colnames(Meth_clusters)[1] <- "sample_id"
for(i in 2:6){
  Meth_clusters.temp <- as.data.frame(results[[i]][["consensusClass"]])
  Meth_clusters.temp$sample_id <- rownames(Meth_clusters.temp)
  colnames(Meth_clusters.temp)[1] <- paste0("K",i,"_group")
  Meth_clusters <- merge(Meth_clusters, Meth_clusters.temp, by.x="sample_id", by.y="sample_id")
}


### 5. Create tumour maps with UMAP (R package umap) using the 15,000 probe m table, colour tumour samples by different features of interest including consensus cluster groups, sex, age, and tumour type
# To generate umap coordinates
library(umap)
model_n100 <- umap(t(m_red), n_neighbors=100)
coordinates_df <- as.data.frame(model_n100$layout)
# What happens when you change the number of nearest neighbors?

# Use ggplot2 to make tumour maps, an example is shown below
library(ggplot2)
ggplot(coordinates_df, aes(x=V1, y=V2, color=K3_group)) +
  geom_point(size=3, alpha = 0.9) +
  scale_color_manual(values=c("#ffd27f","#de425b", "#FF8B61")) +
  theme(legend.position="bottom") +
  labs(title="UMAP dimension 1 vs 2: 15,000 probes", x ="Dim1", y = "Dim2") 











