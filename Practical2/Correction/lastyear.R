Data.Clin=read.table("Data.Clin.txt", header = TRUE, sep = "\t")
Manifest.hg38 = read.delim("Manifest_hg38.txt") 3
Manifest.hg38$UCSC_RefGene_Name = as.character(Manifest.hg38$UCSC_RefGene_Name)```
D_met=read.table("NormalisedFilteredMTable_noInf.csv", header = TRUE, sep = ",")
rownames(D_met)=D_met$X
D_met=D_met[,-1]
colnames(D_met)=sapply(colnames(D_met), function(x) as.character(Data.Clin$Sample_ID[which(as.character(Data.Clin$Meth_ID) == x)]))

D_beta=read.table("NormalisedFilteredBetaTable.csv",header = TRUE, sep = ",")
rownames(D_beta)=D_beta$X
D_beta=D_beta[,-1]
colnames(D_beta)=sapply(colnames(D_beta), function(x) as.character(Data.Clin$Sample_ID[which(as.character(Data.Clin$Meth_ID) == x)]))

Gene_exp = read.csv("gene_FPKM_matrix.csv",row.names = 1)
colnames(Gene_exp) = sapply(strsplit(colnames(Gene_exp),"FPKM."),"[[",2)
expr_all = read.csv("gene_count_matrix_1pass.csv",row.names = 1)

# 1 
`FPKM.diff = apply(Gene_exp,1,max,na.rm=TRUE) - apply(Gene_exp,1,min,na.rm=TRUE)`
- `names(FPKM.diff) = rownames(Gene_exp)`
- `FPKM.diff = FPKM.diff[which(FPKM.diff >= 1)]`
- `expr_all = expr_all[which(rownames(expr_all) %in% names(FPKM.diff)),]`

- `library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)`
- `sexcpgs = rownames(Locations)[Locations$chr %in% c("chrX","chrY","chrM")]`
- `expr_annot = read.table("SRR7646229_pass1_gene_abund.tab",header = T,sep = "\t")[,1:6]`
- `expr_annot$Reference` gives for each gene the related chromosome

e.g.
- `deseqexpr = DESeqDataSetFromMatrix(expr_all, colData = data.frame(colnames(expr_all)), design = ~1, tidy = F)`
- `annot_ordered = expr_annot[sapply( rownames(deseqexpr), function(x) which(expr_annot$Gene.ID==x)[1] ),]`
- `deseqexpr_nosex = deseqexpr[!annot_ordered$Reference %in% c("chrM","chrX","chrY"),]`
- `vstexpr_nosex = varianceStabilizingTransformation(deseqexpr_nosex,blind = T)`
- `vstexpr_nosex = assay(vstexpr_nosex)`

vv = apply(vstexpr_nosex,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
vstexpr_nosex = vstexpr_nosex[order(match(rownames(vstexpr_nosex),names(cv))),]
D_expr_redB = vstexpr_nosex[1:5000,]

`Beta.diff = apply(D_beta,1,max,na.rm=TRUE) - apply(D_beta,1,min,na.rm=TRUE)`
- `names(Beta.diff) = rownames(D_beta)`
- `Beta.diff = Beta.diff[which(Beta.diff >= 0.1)]`
- `D_met_red = D_met[which(rownames(D_met) %in% names(Beta.diff)),]`


 `D_met.pro = D_met_red[which(rownames(D_met_red) %in% as.character(Manifest.hg38$IlmnID[which(Manifest.hg38$class == "Promoter")])),]`
- `D_met.bod = D_met_red[which(rownames(D_met_red) %in% as.character(Manifest.hg38$IlmnID[which(Manifest.hg38$class == "Body")])),]`
- `D_met.enh = D_met_red[which(rownames(D_met_red) %in% as.character(Manifest.hg38$IlmnID[which(Manifest.hg38$class == "Enhancer")])),]`

`vv = apply(D_met.pro,1,var)`
- `cv = cumsum(sort(vv,decreasing = T))/sum(vv)`
- `D_met.pro = D_met.pro[order(match(rownames(D_met.pro),names(cv))),]`
- `D_met.proB = D_met.pro[1:ifelse(nrow(D_met.pro) > 5000, 5000, nrow(D_met.pro)),]`

RNA.Meth.samples = sort(unique(c(colnames(D_expr_redB),colnames(D_met.proB)))) #make a vector with all the cohort samples
D_expr_redB = D_expr_redB[,order(match(colnames(D_expr_redB),RNA.Meth.samples))]
D_exprB_MOFA = as.matrix(D_expr_redB)
dim(D_exprB_MOFA)

D_met.proB_MOFA = matrix(NA,nrow(D_met.proB),length(which(!RNA.Meth.samples%in%colnames(D_met.proB))),dimnames = list(rownames(D_met.proB),RNA.Meth.samples[which(!RNA.Meth.samples%in%colnames(D_met.proB))]))
D_met.proB_MOFA = as.data.frame(D_met.proB_MOFA)
D_met.proB_MOFA = cbind(D_met.proB,D_met.proB_MOFA)
D_met.proB_MOFA = D_met.proB_MOFA[,order(match(colnames(D_met.proB_MOFA),RNA.Meth.samples))]
D_met.proB_MOFA = as.matrix(D_met.proB_MOFA)

D_met.bodB_MOFA = matrix(NA,nrow(D_met.bodB),length(which(!RNA.Meth.samples%in%colnames(D_met.bodB))),dimnames = list(rownames(D_met.bodB),RNA.Meth.samples[which(!RNA.Meth.samples%in%colnames(D_met.bodB))]))
D_met.bodB_MOFA = as.data.frame(D_met.bodB_MOFA)
D_met.bodB_MOFA = cbind(D_met.bodB,D_met.bodB_MOFA)
D_met.bodB_MOFA = D_met.bodB_MOFA[,order(match(colnames(D_met.bodB_MOFA),RNA.Meth.samples))]
D_met.bodB_MOFA = as.matrix(D_met.bodB_MOFA)

D_met.enhB_MOFA = matrix(NA,nrow(D_met.enhB),length(which(!RNA.Meth.samples%in%colnames(D_met.enhB))),dimnames = list(rownames(D_met.enhB),RNA.Meth.samples[which(!RNA.Meth.samples%in%colnames(D_met.enhB))]))
D_met.enhB_MOFA = as.data.frame(D_met.enhB_MOFA)
D_met.enhB_MOFA = cbind(D_met.enhB,D_met.enhB_MOFA)
D_met.enhB_MOFA = D_met.enhB_MOFA[,order(match(colnames(D_met.enhB_MOFA),RNA.Meth.samples))]
D_met.enhB_MOFA = as.matrix(D_met.enhB_MOFA)

par(mfrow=c(2,2),mar = rep(2, 4))
plot(density(as.matrix(D_expr_redB)), xlab="Gene expression (nrc)", ylab="Frequency", main=paste0("Density on RNA-seq data (n=", ncol(D_expr_redB)-1, " ; N=", nrow(D_expr_redB), ")"))
plot(density(as.matrix(D_met.proB), na.rm=T, from=min(D_met.proB, na.rm=T)), xlab="CpGs methylation (M-values)", ylab="Frequency", main=paste0("Density on MetPro data (n=", ncol(D_met.proB), " ; N=", nrow(D_met.proB), ")"))
plot(density(as.matrix(D_met.bodB), na.rm=T, from=min(D_met.bodB, na.rm=T)), xlab="CpGs methylation (M-values)", ylab="Frequency", main=paste0("Density on MetBod data (n=", ncol(D_met.bodB), " ; N=", nrow(D_met.bodB), ")"))
plot(density(as.matrix(D_met.enhB), na.rm=T, from=min(D_met.enhB, na.rm=T)), xlab="CpGs methylation (M-values)", ylab="Frequency", main=paste0("Density on MetEnh data (n=", ncol(D_met.enhB), " ; N=", nrow(D_met.enhB), ")"))
dev.off()

MOFAobject = create_mofa(list("RNA" = D_exprB_MOFA,"MethPro" = D_met.proB_MOFA,"MethBod" = D_met.bodB_MOFA,"MethEnh" = D_met.enhB_MOFA))
plot_data_overview(MOFAobject)

data_optsB = get_default_data_options(MOFAobject)
model_optsB = get_default_model_options(MOFAobject)
model_optsB$num_factors = 5
train_opts = get_default_training_options(MOFAobject)
train_opts$convergence_mode = "slow"
stochastic_opts = get_default_stochastic_options(MOFAobject)

MOFAobject = prepare_mofa(object=MOFAobject, data_options=data_optsB, model_options=model_optsB, training_options=train_opts )

MOFAobject.trained = run_mofa(MOFAobject, save_data = T, outfile = "MOFAobject.hdf5")

plot_variance_explained(MOFAobject.trained)

plot_factors(MOFAobject.trained, 1:10)

plot_factor_cor(MOFAobject.trained)

library(ade4)
require(cowplot)
pca_run = dudi.pca(t(D_expr_redB), scannf = F, nf = 10, center = T, scale = F)
PCs.RNA = pca_run$li
pca_run = dudi.pca(t(rbind(D_met.proB,D_met.enhB,D_met.bodB)), scannf = F, nf = 10, center = T, scale = F)
PCs.Meth = pca_run$li

LFs = get_factors(MOFAobject.trained)$group1

LFs.Meth = LFs[which(rownames(LFs) %in% rownames(PCs.Meth)),]
PCs.Meth = PCs.Meth[order(match(rownames(PCs.Meth),rownames(LFs))),]

M = cor(LFs.Meth,PCs.Meth)
corrplot(M, type="full", order="original",
         col=brewer.pal(n=8, name="RdYlBu"))

LFs.RNA = LFs[which(rownames(LFs) %in% rownames(PCs.RNA)),]
PCs.RNA = PCs.RNA[order(match(rownames(PCs.RNA),rownames(LFs))),]

R = cor(LFs.RNA,PCs.RNA)
corrplot(R, type="full", order="original",
         col=brewer.pal(n=8, name="RdYlBu"))

sample = rownames(LFs)
classe = unlist(lapply(sample, function(x){
  if(x %in% Data.Clin$Sample_ID[which(Data.Clin$Molecular_clusters == "LC1")]) {"#009975ff"}
  else if(x %in% Data.Clin$Sample_ID[which(Data.Clin$Molecular_clusters == "LC2")]) {"#ff0000ff"}
  else if(x %in% Data.Clin$Sample_ID[which(Data.Clin$Molecular_clusters == "LC3")]) {"#e88200ff"}
}))

plot(LFs[,1:2],type="p",pch=16, col = classe)
grid()

require(umap)
umapMOFA = umap(LFs)
plot(umapMOFA$layout,col = classe, pch=16,
     xlab= "UMAP1",
     ylab= "UMAP2",
     main=paste0("UMAP from MOFA based on RNA-seq and 850K data (n=", nrow(LFs), ")"))
grid()

M = cor(LFs,umapMOFA$layout)
corrplot(M, type="full", order="original",
         col=brewer.pal(n=8, name="RdYlBu"))

plot_top_weights(MOFAobject.trained,view = "RNA",factor = 1, nfeatures = 10, scale = T)
plot_weights(MOFAobject.trained,view = "RNA",factor = 1, nfeatures = 10, scale = T)

w1.rna = plot_weights(MOFAobject.trained,view = "RNA",factor = 1, nfeatures = 10, scale = T, return_data = T)
w1.rna$feature = as.character(w1.rna$feature)
expr_annot$Gene.ID = as.character(expr_annot$Gene.ID)
expr_annot$Gene.Name = as.character(expr_annot$Gene.Name)
all(w1.rna$feature %in% expr_annot$Gene.ID)
any(duplicated(expr_annot$Gene.ID))
w1.rna$gene = sapply(w1.rna$feature, function(x) paste0(unique(expr_annot$Gene.Name[which(expr_annot$Gene.ID == x)]), collapse = ";"))
head(w1.rna)
w1.rna$abs_w = abs(w1.rna$value)
w1.rna$gene[order(w1.rna$abs_w, decreasing = T)][1:10]

w2.rna = plot_weights(MOFAobject.trained,view = "RNA",factor = 2, nfeatures = 10, scale = T, return_data = T)
w2.rna$feature = as.character(w2.rna$feature)
w2.rna$gene = sapply(w2.rna$feature, function(x) paste0(unique(expr_annot$Gene.Name[which(expr_annot$Gene.ID == x)]), collapse = ";"))
w2.rna$abs_w = abs(w2.rna$value)
w2.rna$gene[order(w2.rna$abs_w, decreasing = T)][1:10]

Manifest_hg38$UCSC_RefGene_Name = as.character(Manifest_hg38$UCSC_RefGene_Name)
w1.enh = plot_weights(MOFAobject.trained,view = "MethEnh",factor = 1, nfeatures = 10, scale = T, return_data = T)
w1.enh$feature = as.character(w1.enh$feature)
w1.enh$gene = sapply(w1.enh$feature, function(x) paste0(unique(Manifest_hg38$UCSC_RefGene_Name[which(Manifest_hg38$IlmnID == x)]), collapse = ";"))
w1.enh$abs_w = abs(w1.enh$value)
w1.enh$gene[order(w1.enh$abs_w, decreasing = T)][1:10]

w2.enh = plot_weights(MOFAobject.trained,view = "MethEnh",factor = 2, nfeatures = 10, scale = T, return_data = T)
w2.enh$feature = as.character(w2.enh$feature)
w2.enh$gene = sapply(w2.enh$feature, function(x) paste0(unique(Manifest_hg38$UCSC_RefGene_Name[which(Manifest_hg38$IlmnID == x)]), collapse = ";"))
w2.enh$abs_w = abs(w2.enh$value)
w2.enh$gene[order(w2.enh$abs_w, decreasing = T)][1:10]


write.table(w1.rna$gene,"w1.rna.txt",row.names = T,col.names = T,quote = F, sep="\t")
write.table(w2.rna$gene,"w2.rna.txt",row.names = T,col.names = T,quote = F, sep="\t")
write.table(unlist(strsplit(w1.enh$gene[which(w1.enh$gene != "")],";")),"w1.enh.txt",row.names = T,col.names = T,quote = F, sep="\t")
write.table(unlist(strsplit(w2.enh$gene[which(w2.enh$gene != "")],";")),"w2.enh.txt",row.names = T,col.names = T,quote = F, sep="\t")

LFs = LFs[order(match(rownames(LFs),Data.Clin$Sample_ID)),]
all(Data.Clin$Sample_ID == rownames(LFs))

cli = colnames(Data.Clin)[-c(1,ncol(Data.Clin))] #Clinical data

var_cli = cbind(Data.Clin[,-c(1,ncol(Data.Clin))],LFs)
rownames(var_cli) = rownames(LFs)

p_cli = data.frame()
q_cli = data.frame()
for(j in 1:length(cli)){
  for(i in 1:ncol(LFs)){
    fact = colnames(LFs)[i]
    fit = lm(var_cli[,fact] ~ as.character(var_cli[,j]), data=var_cli)
    p_cli[j,i] = lmp(fit)
    colnames(p_cli)[i] = fact
  }
  rownames(p_cli)[j] = colnames(var_cli)[j]
  for(i in 1:ncol(LFs)){
    fact = colnames(LFs)[i]
    q_cli[j,i] = p.adjust(p_cli[j,],method = "BH")[i]
    colnames(q_cli)[i] = fact
  }
  rownames(q_cli)[j] = colnames(var_cli)[j]
}
q_cli = as.data.frame(matrix(p.adjust(as.vector(as.matrix(p_cli)), method="BH"),ncol=10))
colnames(q_cli) = colnames(p_cli)
rownames(q_cli) = rownames(p_cli)

# correlation matrix
qcor_all = matrix(t(q_cli),ncol(LFs),length(cli))
colnames(qcor_all) = rownames(q_cli)
rownames(qcor_all) = colnames(q_cli)


install.packages("factoextra")
library(factoextra)

fviz_nbclust(LFs[,1:2], kmeans, method = "wss")
km.res = kmeans(LFs[,1:2], 3, nstart = 25)
plot(LFs[,1:2], col = km.res$cluster, pch = 19)


cosmic = read.table("Census_allThu Jan 31 16_52_52 2019.csv",quote="\"",stringsAsFactors=F, na.strings = c("NA","."),sep=",",header=T)
w1.rna$gene[which(w1.rna$gene %in% cosmic$Gene.Symbol)]
w2.rna$gene[which(w2.rna$gene %in% cosmic$Gene.Symbol)]

Data.Clin$cluster = factor(km.res$cluster)
rownames(Data.Clin) = Data.Clin$Sample_ID

dds = DESeqDataSetFromMatrix(countData = expr_all, colData = Data.Clin, design = ~ cluster)
design(dds)
dds = DESeq(dds)
res = results(dds)
head(res)


