library(ade4)
library(corrplot)
library(ggpubr)
library(patchwork)
library(MOFA2)

# run PCAs
pcaRNA = dudi.pca(t(RNA),nf = 5,scannf = F)
pcaMethprom = dudi.pca(t(DNAMeth_promoter[,!is.na(DNAMeth_promoter[1,])]),nf = 5,scannf = F)
pcaMethgene = dudi.pca(t(DNAMeth_genebody[,!is.na(DNAMeth_genebody[1,])]),nf = 5,scannf = F)
pcaMethenh  = dudi.pca(t(DNAMeth_enhancer[,!is.na(DNAMeth_enhancer[1,])]),nf = 5,scannf = F)

# plot variance explained by each PCA
## create variance table
R2tab = tibble(axis=c(1:29, rep(1:17,3)), 
               Omic=c(rep("RNA",29),rep("Meth_pro",17),rep("Meth_bod",17),rep("Meth_enh",17)), 
               var=c(pcaRNA$eig/sum(pcaRNA$eig),
                     pcaMethprom$eig/sum(pcaMethprom$eig),
                     pcaMethgene$eig/sum(pcaMethgene$eig),
                     pcaMethenh$eig/sum(pcaMethenh$eig)) )
## order omic levels
R2tab$Omic = factor(R2tab$Omic,levels = c("RNA","Meth_pro","Meth_bod","Meth_enh"))

## plot
ggplot( data= R2tab %>% filter(axis<=5), aes(x=axis,y=var,fill=var)) + geom_bar(stat="identity") + 
  facet_grid(~Omic) + theme_classic() + scale_fill_gradient2() + ylab("Proportion of variance")

# plot comparisons between PCA axes and MOFA factors
## create table with all axes and factors
pca_join = full_join(left_join(bind_cols(Sample=rownames(pcaMethgene$li), pcaMethgene$li), #%>% pivot_longer(cols=Axis1:Axis5) %>% mutate(Omic="Meth_bod")
bind_cols(Sample=rownames(pcaMethprom$li), pcaMethprom$li), by="Sample", suffix = c(".Meth_bod",".Meth_pro") ), #%>% pivot_longer(cols=Axis1:Axis5) %>% mutate(Omic="Meth_bod")
left_join(bind_cols(Sample=rownames(pcaRNA$li), pcaRNA$li) , 
    bind_cols(Sample=rownames(pcaMethenh$li), pcaMethenh$li), by="Sample", suffix = c(".RNA",".Meth_enh")) , by="Sample" )  #%>% pivot_longer(cols=Axis1:Axis5) %>% mutate(Omic="Meth_bod")

pca_join = left_join(pca_join , 
                     bind_cols(Sample=rownames(mofa_trained@expectations$Z$group1), mofa_trained@expectations$Z$group1), 
                     by ="Sample" )

## plot some comparisons
gg1bp <- ggplot(data = pca_join, aes(y=Axis1.Meth_bod, x= Axis1.Meth_pro)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)
gg1be <- ggplot(data = pca_join, aes(y=Axis1.Meth_bod, x= Axis1.Meth_enh)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)
gg1br <- ggplot(data = pca_join, aes(y=Axis1.Meth_bod, x= Axis1.RNA)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)
gg1bM <- ggplot(data = pca_join, aes(y=Axis1.Meth_bod, x= Factor1)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)
gg1rM <- ggplot(data = pca_join, aes(y=Axis1.RNA, x= Factor1)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)
gg13rM <- ggplot(data = pca_join, aes(y=Axis1.RNA, x= Factor3)) + geom_point() + theme_classic() + geom_smooth(method = "lm", se = FALSE)

## combine plots
gg1bp+gg1be+gg1br +gg1bM

gg13rM

# compute correlations between MOFA factors and PCA axes
rRNA <- abs(cor(x = do.call(rbind, get_factors(mofa_trained)), 
             y = pcaRNA$li, use = "complete.obs"))
rMp <- abs(cor(x = do.call(rbind, get_factors(mofa_trained))[rownames(pcaMethprom$li),], 
                y = pcaMethprom$li, use = "complete.obs"))
rMg <- abs(cor(x = do.call(rbind, get_factors(mofa_trained))[rownames(pcaMethprom$li),], 
               y = pcaMethgene$li, use = "complete.obs"))
rMe <- abs(cor(x = do.call(rbind, get_factors(mofa_trained))[rownames(pcaMethprom$li),], 
               y = pcaMethenh$li, use = "complete.obs"))

## plot correlation matrices
corrplot(rRNA, tl.col = "black")
corrplot(rMp, tl.col = "black")
corrplot(rMg, tl.col = "black")
corrplot(rMe, tl.col = "black")
