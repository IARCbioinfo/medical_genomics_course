library(MOFA2)

# a)
plot_factors(mofa_trained,factors = 1:5)

# b)
plot_factors(mofa_trained,factors = 1:5,color_by = Data.Clin$Histopathology)
plot_factors(mofa_trained,factors = 1:5,color_by = Data.Clin$Molecular_clusters)

## alternative 
LFs = get_factors(mofa_trained)
plot(LFs$group1[,1],LFs$group1[,2],col=factor(Data.Clin$Histopathology))

## c) 
mofa_trained_umap = run_umap(mofa_trained)

ggplot( mofa_trained_umap@dim_red$UMAP , aes(x=UMAP1,y=UMAP2,shape=Data.Clin$Histopathology, col=Data.Clin$Molecular_clusters) ) + 
  geom_point() + theme_classic()

# Q6
## a) 
### plot
plot_weights(mofa_trained)
plot_top_weights(mofa_trained,view = "RNA",factors = 1)
### get names
WRNA      = get_weights(mofa_trained,views = "RNA")$RNA
### get top 10
RNA_LF1_top10 = names(sort(abs(WRNA[,1]),decreasing = T))[1:10]
### read annotations
RNA_annot = read_tsv("medical_genomics_course/Practical2/Data/RNA_gene_annotation.tsv")
lapply( RNA_LF1_top10 , function (x) RNA_annot[RNA_annot$gene_id==x,])

## b)
plot_top_weights(mofa_trained,view = "Meth_enh",factors = 1)
WMeth_enh = get_weights(mofa_trained,views = "Meth_enh")$Meth_enh
Methenh_LF1_top10 = names(sort(abs(WMeth_enh[,1]),decreasing = T))[1:10]
Methenh_annot = read_tsv("medical_genomics_course/Practical2/Data/DNAMeth_enhancer_annotation.tsv")
lapply( Methenh_LF1_top10 , function (x) Methenh_annot[Methenh_annot$Ilmid==x,])

## c) 
library(MOFAdata)
data("MSigDB_v6.0_C2_human") 
head(rownames(MSigDB_v6.0_C2_human), n=3)
enrichment.parametric <- run_enrichment(mofa_trained,view = "RNA", factors = 1:3,
                                        feature.sets = MSigDB_v6.0_C2_human,
                                        sign = "negative",
                                        statistical.test = "parametric")
plot_enrichment(enrichment.parametric, factor = 1, max.pathways = 15)
plot_enrichment(enrichment.parametric, factor = 3, max.pathways = 15)

# clean session
rm(MSigDB_v6.0_C2_human)
gc()

# Q7
LFs = get_factors(mofa_trained)
summary( lm(LFs$group1[,1]~Sex + Age + Histopathology, data=Data.Clin) )
summary( lm(LFs$group1[,2]~Sex + Age + Histopathology, data=Data.Clin) )
summary( lm(LFs$group1[,3]~Sex + Age + Histopathology, data=Data.Clin) )