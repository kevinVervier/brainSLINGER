---
title: "pharmacogenomics_bipolar"
author: "Kevin Vervier"
date: "April 16, 2018"
output: html_document
---

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. 

In this document, we summarized work done on brain-specific transcriptome imputation and its applications to identifying potential new treatments in psychiatry.

# **Brain transcriptome inference models**

In this section, we descrided how we adapted the SLINGER framework, which initially applies to whole blood transcriptome, and created gene-level predictive models. To inject brain-related information, we relied on TiSAn scores for each genotyped SNP, found in the DGN study.

![](../pic/brain_slinger_scheme.png)

## **TiSAn score extraction**

We provide the scores for the DGN SNPs (list found here: *../data/DGN_hg19_locations.bed*) used in this project: *../data/tisan_brain_SNP.Rdata*. To extract TiSAn brain scores for the SNPs found in a different project, please refer to http://github.com/kevinVervier/TiSAn. 

We used the source in **train_model.R**, combined with **run_model.sh** and **main.sh** to make it parallel on a cluster. Each of the jobs returned a set of gene-level elastic-net models for brain tissues. It is also possible to run all the gene models at the same time without parallelization (can take days).

```{r,eval=F,echo=T}
# merge all the brainSlinger models in one file
brain.models = NULL
for(i in 1:100){
  x = read.table(paste('../data/01.DGN-data.Rdata/brain_elasticNet_alpha0.5_weights_chunk',i,'.txt',sep=''),header=TRUE,sep='\t')
  brain.models = rbind(brain.models, x)
}

write.table(brain.models,file='../data/singler_models_Brain.txt',row.names = FALSE,quote = FALSE,sep='\t')
```

## **Brain models performance versus SLINGER**

In this section, we present several metrics made to compare brain-SLINGER models with regular whole-blood SLINGER models found here (http://github.com/kevinVervier/SLINGER/tree/master/model/AllSNPs-elasticNet.txt).
```{r}

# load control: tissue-agnostic SLINGER models
control = read.table('../data/slinger_models_WB.txt',header=TRUE)
control.genes = as.character(unique(control$gene)) # 11,959

# load models trained using TiSAn-brain scores
x = read.table('../data/singler_models_Brain.txt',header=TRUE,sep='\t')
brain.genes = unique(as.character(x$gene))
# 9,747 genes with at least one 'brain-SNP'

newGenes = brain.genes[which(!(brain.genes %in% control.genes))]
```

```{r,eval=F,echo=T}
# find genes in brain-slinger but not in regular slinger
write.table(brain.genes[which(!(brain.genes %in% control.genes))],file='../data/new_brain_genes.txt',quote = FALSE,row.names = FALSE,col.names = FALSE) # 1,337 new genes
write.table(control.genes[which(!(control.genes %in% brain.genes))],file='../data/missing_brain_genes.txt',quote = FALSE,row.names = FALSE,col.names = FALSE) # 3,549 lost genes
write.table(control.genes,file='../data/control_genes.txt',quote = FALSE,row.names = FALSE,col.names = FALSE)

```

We compared genes with a model in SLINGER and in brain-SLINGER and found `r length(newGenes)` genes only accessible through the brain enriched procedure, and not in SLINGER. Moreover, we used this list of new genes in Panther website (versus all the remaining genes) for GO terms enrichment and found the following significant enrichments:

* neurological system process (BP)
* signal transducer activity (MF)
* receptor activity (MF)

Interestingly, in the process of training brain SLINGER models, we also lost the ability to predict transcriptome for `r  length(which(!(control.genes %in% brain.genes)))` genes. When compared with the newly acquired genes, we did not found any GO term enrichment for the missing genes, and found *stress response* as the only enriched term in the new genes.

We also look at the accuracy of the lost models in terms of their cross-validated r-squared values in SLINGER, compared to the rest of the SLINGER models.

```{r}
perfs <- read.csv('slinger_accuracy.csv',header=TRUE)
tmp = perfs[which(perfs$gene %in% control.genes[which(!(control.genes %in% brain.genes))]),]

t.test(tmp$R2_all_SNPs,perfs$R2_all_SNPs)
```
This test tells us that, on average, the genes we lost tend to have a significantly lower accuracy, suggesting that using TiSAn also lead to less spurious association between genetic variation and gene expression.

```{r, fig.height = 8, fig.width = 12, fig.align = "center"}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(VennDiagram,quietly=T)
grid.newpage()

venn.plot <- draw.pairwise.venn(9747,11959, 9747-1337, c("brain-SLINGER", "WB-SLINGER"),lwd=5,cex=3,fill=cbPalette[3:2],scaled=TRUE,
                                alpha = 0.6,fontface=4,cat.col=cbPalette[3:2],cat.cex=3,cat.dist=c(0.04,0.05),cat.pos=c(-150,45))
grid.draw(venn.plot)

```

Following findings from *Functional Architectures of Local and Distal Regulation of Gene Expression in Multiple Human Tissues* (Liu,2017), we investigated if the tissue-specific models tend to be enriched in distal regulatory elements.

First, we compute local and distal regulatory elements found in SLINGER models.

```{r, echo =T, eval =F}

# load models
models = read.table('/wdata/kvervier/PrediXcan/output/01.DGN-data.Rdata/AllSNPs-elasticNet.txt',header=TRUE)
# load gene info (genecode)
genes = read.table('/wdata/kvervier/PrediXcan/input/gencode',header=FALSE)

# get snps info
load('/wdata/kvervier/PrediXcan/input/01.DGN-data.Rdata')
# need to convert to GrChr37 (liftover)
# load hg19 location
loc = read.table('/wdata/kvervier/pharmacogenomics/input/DGN_hg19_locations.bed',sep = ':')
loc[,1] = gsub(loc[,1],pattern = 'chr',replacement = '')
loc[,3] = gsub(loc[,2],pattern='.*-',replacement = '')
loc[,2] = loc[,3]
# read deleted positions
del_pos = read.table('/wdata/kvervier/pharmacogenomics/input/DGN_deleted_pos.txt',header=FALSE) # 349 positions
# remove them from snp_locations
idx = match(paste(paste(paste('chr',snp_locations[,2],sep=''),snp_locations[,3],sep=':'),snp_locations[,3],sep='-'), del_pos[,1])
snp_locations = snp_locations[-which(!is.na(idx)),]

# update snp_locations
snp_locations[,2] = as.numeric(loc[,1])
snp_locations[,3] = as.numeric(loc[,2])
snp_locations[,4] = as.numeric(loc[,2])


##########################################
# analysis on interactions
snps = models[which(models$SNP %in% snp_locations$refsnp_id),]
snps$chr.snp = snp_locations$chr_name[match(snps$SNP,snp_locations$refsnp_id)]
snps$pos.snp = snp_locations$chrom_start[match(snps$SNP,snp_locations$refsnp_id)]

snps$chr.gene = gsub(pattern = 'chr',replacement = '',genes$V1[match(snps$gene,genes$V6)])
snps$start.gene = genes$V3[match(snps$gene,genes$V6)]
snps$end.gene =   genes$V4[match(snps$gene,genes$V6)]
snps$ens.id =  gsub(pattern = '\\..*',replacement='',x=genes$V5[match(snps$gene,genes$V6)])

########################################
# separate local and distal interactions

idx = which(snps$chr.snp != snps$chr.gene | abs(snps$pos.snp - snps$start.gene) > 1000000 |  abs(snps$pos.snp - snps$end.gene) > 1000000 )

distal = snps[idx,]
local = snps[-idx,]

save(distal,local,file='/wdata/kvervier/PrediXcan/output/01.DGN-data.Rdata/18.data_split.Rdata')

```

```{r}
load('/wdata/kvervier/PrediXcan/output/01.DGN-data.Rdata/18.data_split.Rdata')
# get overall proportion of local/distal
nrow(local)/(nrow(local)+nrow(distal)) 
```
```{r, echo=T,eval=F}
# get proportion of local/distal per gene
tmp = rbind(local,distal)
count = sapply(unique(tmp$gene),function(g) nrow(local[which(local$gene==g),])/(nrow(local[which(local$gene==g),])+nrow(distal[which(distal$gene==g),])))
save(count,distal,local,file='/wdata/kvervier/PrediXcan/output/01.DGN-data.Rdata/18.data_split.Rdata')
```
```{r}
load('/wdata/kvervier/PrediXcan/output/01.DGN-data.Rdata/18.data_split.Rdata')
mean(count) 
```

Then, we applied a similar methodology to brain-SLINGER models:
```{r, echo =T, eval =F}
# load models
models = read.table('/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/brain_models.txt',header=TRUE,sep='\t')
# load gene info (genecode)
genes = read.table('/wdata/kvervier/PrediXcan/input/gencode',header=FALSE)

# get snps info
load('/wdata/kvervier/PrediXcan/input/01.DGN-data.Rdata')
# need to convert to GrChr37 (liftover)
# load hg19 location
loc = read.table('/wdata/kvervier/pharmacogenomics/input/DGN_hg19_locations.bed',sep = ':')
loc[,1] = gsub(loc[,1],pattern = 'chr',replacement = '')
loc[,3] = gsub(loc[,2],pattern='.*-',replacement = '')
loc[,2] = loc[,3]
# read deleted positions
del_pos = read.table('/wdata/kvervier/pharmacogenomics/input/DGN_deleted_pos.txt',header=FALSE) # 349 positions
# remove them from snp_locations
idx = match(paste(paste(paste('chr',snp_locations[,2],sep=''),snp_locations[,3],sep=':'),snp_locations[,3],sep='-'), del_pos[,1])
snp_locations = snp_locations[-which(!is.na(idx)),]

# update snp_locations
snp_locations[,2] = as.numeric(loc[,1])
snp_locations[,3] = as.numeric(loc[,2])
snp_locations[,4] = as.numeric(loc[,2])

##########################################
# analysis on interactions
snps = models[which(models$SNP %in% snp_locations$refsnp_id),]
snps$chr.snp = snp_locations$chr_name[match(snps$SNP,snp_locations$refsnp_id)]
snps$pos.snp = snp_locations$chrom_start[match(snps$SNP,snp_locations$refsnp_id)]

snps$chr.gene = gsub(pattern = 'chr',replacement = '',genes$V1[match(snps$gene,genes$V6)])
snps$start.gene = genes$V3[match(snps$gene,genes$V6)]
snps$end.gene =   genes$V4[match(snps$gene,genes$V6)]
snps$ens.id =  gsub(pattern = '\\..*',replacement='',x=genes$V5[match(snps$gene,genes$V6)])

########################################
# separate local and distal interactions

idx = which(snps$chr.snp != snps$chr.gene | abs(snps$pos.snp - snps$start.gene) > 1000000 |  abs(snps$pos.snp - snps$end.gene) > 1000000 )

distal.brain = snps[idx,]
local.brain = snps[-idx,]

save(distal.brain,local.brain,file='/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/18.data_split.Rdata')
```
```{r}
# get overall proportion of local/distal
load('/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/18.data_split.Rdata')
nrow(local.brain)/(nrow(local.brain)+nrow(distal.brain))
```
```{r, echo =T, eval =F}
# get proportion of local/distal per gene
tmp = rbind(local.brain,distal.brain)
count.brain = sapply(unique(tmp$gene),function(g) nrow(local.brain[which(local.brain$gene==g),])/(nrow(local.brain[which(local$gene==g),])+nrow(distal.brain[which(distal.brain$gene==g),])))
save(count.brain,local.brain,distal.brain,file='/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/18.data_split.Rdata')
```
```{r}
load('/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/18.data_split.Rdata')
mean(count.brain)
```

We tested if the two distribtions were different using Kolmogorow-Smirnof test:
```{r, fig.align = "center"}
ks.test(count,count.brain,alternative = 'less')
ks.test(count,count.brain,alternative = 'less')$p.val

hist(count,col=rgb(1,0,0,0.5),main='Proportion of Local Enrichment in gene models')
hist(count.brain,add=TRUE,col=rgb(0,0,1,0.5))
legend(x='topright',legend = c('SLINGER','brain SLINGER'),fill =c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)))
```

## **Brain models performance versus PrediXcan-GTEx**

In this section, we proposed to compare the brain-SLINGer models trained using a large collection of whole blood expressions and computationally enriched in brain-related elements by TiSAn, with models trained using PrediXcan on a limited number of brain expression samples from GTEx database.

To compare these approaches, we measured the correlation between their predicted gene expression and the actual measured expression in GTEx 10 brain regions.

PrediXcan models performances were obtained using */wdata/kvervier/pharmacogenomics/src/PRX/04.get_gtex_v6_preds.R*. 

```{r, echo=T,eval=F}
tissues = c('Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex',
            'Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus',
            'Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia')

PVAL.sum = NULL
for(tissue in tissues){
  #load PrediXcan predictions
  load(paste('/wdata/kvervier/PrediXcan/output/GTEx/PRX_gtex',tissue,'_expr.Rdata',sep=''))
  info = data.frame(gene=names(R2),pred.perf.R2 = R2, pred.perf.pval = PVAL,cor=COR, pred.perf.qval = p.adjust(PVAL, "fdr"))
  # convert the Ensembl gene name to symbol
  library(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  map_genes <- getBM(attributes=c('ensembl_gene_id',
                                  'hgnc_symbol'),filters = 'ensembl_gene_id', values = info$gene ,mart = ensembl)
  info$symbol = map_genes$hgnc_symbol[match(info$gene,map_genes$ensembl_gene_id)]
  min.prx = min(info$pred.perf.R2)

  # brain SLINGER model predictions
  load(file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_expr.Rdata',sep=''))
  
  # use adjusted R2 value
  idx = sort(intersect(names(R2),info$symbol),decreasing=FALSE)
  sub.info = info[which(info$symbol %in% idx),]
  sub.info = sub.info[order(sub.info$symbol,decreasing=FALSE),]
  sub.cor = COR[idx]  
  # remove pairs of NAs
  idx.filt = which(is.na(sub.cor) & is.na(sub.info$cor))  
  sub.cor = sub.cor[-idx.filt]
  sub.info = sub.info[-idx.filt,]
   
  cat(tissue,':',t.test(sub.cor,sub.info$cor,alternative = 'greater',paired=TRUE)$p.val,'\n')
  PVAL.sum = c(PVAL.sum,t.test(sub.cor,sub.info$cor,alternative = 'greater',paired=TRUE)$p.val)
}

names(PVAL.sum) = tissues
save(PVAL.sum,file='/wdata/kvervier/pharmacogenomics/output/GTEx/pval_comparison_slinger_prx.Rdata')
```

```{r, echo=T, eval=F}


# produce plot (boxplot or barplot) for comparison between Brain-slinger and GTEX PRX

DF = NULL

tissues = c('Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex',
            'Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus',
            'Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia')[order(PVAL.sum,decreasing=FALSE)]

for(tissue in tissues){

  load(paste('/wdata/kvervier/PrediXcan/output/GTEx/PRX_gtex',tissue,'_expr.Rdata',sep=''))
  info = data.frame(gene=names(R2),pred.perf.R2 = R2, pred.perf.pval = PVAL,cor=COR, pred.perf.qval = p.adjust(PVAL, "fdr"))
  # convert the Ensembl gene name to symbol
  library(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  map_genes <- getBM(attributes=c('ensembl_gene_id',
                                  'hgnc_symbol'),filters = 'ensembl_gene_id', values = info$gene ,mart = ensembl)
  info$symbol = map_genes$hgnc_symbol[match(info$gene,map_genes$ensembl_gene_id)]
  min.prx = min(info$pred.perf.R2)
  #################################
  # brain SLINGER model predictions
  load(file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_expr.Rdata',sep=''))
  
  # use adjusted R2 value
  idx = sort(intersect(names(R2),info$symbol),decreasing=FALSE)

  sub.info = info[which(info$symbol %in% idx),]

  sub.info = sub.info[order(sub.info$symbol,decreasing=FALSE),]

  sub.cor = COR[idx]  
  # remove pairs of NAs
  idx.filt = which(is.na(sub.cor) & is.na(sub.info$cor))  
  sub.cor = sub.cor[-idx.filt]
  names(sub.cor) = NULL
  sub.info = sub.info[-idx.filt,]

 DF = rbind(DF,cbind(rep(tissue,2*length(sub.cor)),c(rep('PrediXcan',length(sub.cor)),rep('brain-SLINGER',length(sub.cor))),c(sub.info$cor,sub.cor)))
}

df = as.data.frame(DF)
names(df) = c('tissue','method','cor')
df$cor = as.numeric(as.character(df$cor))
df$tissue = factor(df$tissue,levels = unique(df$tissue))  
save(df,file='/wdata/kvervier/pharmacogenomics/output/GTEx/pval_comparison_slinger_prx_cor.Rdata')
```

```{r, fig.height = 10, fig.width = 12, fig.align = "center"}
load('/wdata/kvervier/pharmacogenomics/output/GTEx/pval_comparison_slinger_prx.Rdata')
sort(PVAL.sum,decreasing=FALSE)
load('/wdata/kvervier/pharmacogenomics/output/GTEx/pval_comparison_slinger_prx_cor.Rdata')
############
# Comparison plot (only for tissues with a significant difference)

idx = names(which(PVAL.sum < 0.05))
sub.df = df[which(df$tissue %in% idx),]
sub.df$tissue = factor(sub.df$tissue,levels = unique(sub.df$tissue))  

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <-ggplot(sub.df, aes(x=tissue,y= cor, fill=method))

p <- p + geom_boxplot()+ ylab('Spearman correlation between predicted and \nobserved gene expressions') + 
  theme_grey(base_size = 30) + theme(axis.text.x=element_text(colour="black")) +
  theme(legend.title = element_text(size=22, face="bold")) + 
  theme(legend.text = element_text(size = 22)) + xlab('GTEx tissue')+ #+ ggtitle('Comparison between brain-SLINGER and GTEx-PrediXcan') 
  scale_fill_manual(values=cbPalette[3:2],name="method") +
  theme(legend.key.size = unit(1.5, "cm")) #+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
p 

```

This study suggests that brain-SLINGER performances are comparable with models trained on actual brain expression data (PrediXcan-GTEx).

Even better, for two brain regions, brain-SLINGER shows a stronger correlation between estimated gene expression and measured expression, across thousands of genes.

Moreover, we split all the considered genes in two groups: the ones with a stronger correlation found in brain SLINGER, and the ones for which regular SLINGER was better at estimating GTEX brain expressions.

**TODO: need to do the same analysis with PRX models**

```{r,eval=F,echo=T}
tissues = c('Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex',
            'Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus',
            'Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia')

for(tissue in tissues){
  cat(tissue,'\n')
  # slinger
  load(file=paste('/wdata/kvervier/PrediXcan/output/GTEx/SLINGER_',tissue,'_expr.Rdata',sep=''))
  slinger.cor = COR
  #brain-slinger
  load(file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_expr.Rdata',sep=''))
  brain.cor = COR
  idx = intersect(names(brain.cor),names(slinger.cor))
  brain.cor = brain.cor[idx]
  slinger.cor = slinger.cor[idx]
  #length(brain.cor[which(brain.cor > slinger.cor)])
  tmp = brain.cor[which(brain.cor > slinger.cor)]
  write.table(names(tmp),file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_higher_cor.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
  tmp2 = tmp[tmp>0]
  write.table(names(tmp2),file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_higher_cor_positive_only.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
  tmp3 = slinger.cor[which(slinger.cor > brain.cor)]
  write.table(names(tmp3),file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_null_model.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
#get null model when only brain_cor > 0 are kept:
  tmp4 = c(tmp3,tmp[tmp<0])
  write.table(names(tmp4),file=paste('/wdata/kvervier/pharmacogenomics/output/GTEx/brain_model_',tissue,'_null_model_positive_only.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
  
}

```

We ran GO term enrichment analysis (Panther.DB) for each brain region using the set of genes found to be better predicted with regular SLINGER as the background set of genes. This table summarizes the top GO term findings (one example per region).

GTEx tissue  | GO term  | category | corrected pvalue 
------------ | ------- | -------- | ---------------  
Anterior cingulate cortex (BA24)    | regulation of axogenesis | BP |   5.13E-03      
Caudate basal ganglia       | regulation of histone acetylation | BP |  1.77E-02   
Cerebellar hemisphere   | muscle contraction | BP |   4.56E-03   
Cerebellum       | locomotion | BP |  1.15E-03
Cortex    | transmembrane receptor | BP |   3.42E-02    
Frontal cortex     | cell differentiation | BP |  1.59E-02 
Hippocampus    | neuron differentiation | BP |  1.64E-02  
Hypothalamus     | positive regulation of neurogenesis | BP |  8.84E-06
Anterior cingulate cortex (BA24)    | regulation of axogenesis | BP |   5.13E-03      
Putamen basal ganglia      | calcium channel activity | MF |  4.10E-04

# **Candidate genes for bipolar disorder**

In this section, we applied brain transcriptome inference on the WTCCC bipolar disorder cohort.

## **Brain transcriptome imputation**

Here, we applied brain-SLINGER models to predict gene expression for WTCCC cohorts: 2 controls (numbers 1 and 2) and 7 diseases (including bipolar disorder, number 3).

```{r, echo =T, eval = F}

# get brain-slinger predictions (with exclusion lists)

for(i in 1:9){
  cat(i,'\n')
  #load WTCCC data
  load(paste('/wdata/kvervier/PrediXcan/input/WTCCC/EGAD0000000000',i,'_filtered_indiv_excluded_SNP.Rdata',sep=''))
  
  ###############################################################################
  # init gene expression matrix (output)
  EXPR2 = NULL #will be a indiv x gene_expr matrix (obtained with non-excluded SNPs)
  NAMES = NULL
  # loop over all models
  for(chunk in 1:100){
    cat('Chunk',chunk,'\n')
    #read models in chunk i
    models = read.delim(paste('/wdata/kvervier/pharmacogenomics/output/01.DGN-data.Rdata/brain_elasticNet_alpha0.5_weights_chunk',chunk,'.txt',sep=''),header=TRUE)
    for(gene in unique(models$gene)){
      #store names
      NAMES = c(NAMES,gene)
      #subset beta weights
      idx = which(models$gene == gene)
      tmp2 = sub.GENO[,which(colnames(sub.GENO) %in% models$SNP[idx])]
      beta = models$beta[idx]
      names(beta) = models$SNP[idx]
      
      if(!is.null(dim(tmp2))){
        EXPR2 = cbind(EXPR2, tmp2%*%beta[colnames(tmp2)])
      }else{
        EXPR2 = cbind(EXPR2, tmp2*beta[colnames(sub.GENO)[which(colnames(sub.GENO) %in% models$SNP[idx])]])
      }  
    }
  }
  colnames(EXPR2) = NAMES
  save(EXPR2,file=paste('/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD0000000000',i,'_expr_filtered_indiv_excluded_SNP.Rdata',sep=''))
}

```

## **Differentially Expressed Genes in Bipolar Disorder**
```{r, echo=T, eval=F}
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000001_expr_filtered_indiv_excluded_SNP.Rdata')
controls = EXPR2
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000002_expr_filtered_indiv_excluded_SNP.Rdata')
controls = rbind(controls,EXPR2)
# load bipolar disorder
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000003_expr_filtered_indiv_excluded_SNP.Rdata')
case = EXPR2
# get all genes with diff expression:
#init labels
Y = c(rep(1,nrow(case)),rep(-1,nrow(controls)))
Y = factor(Y)

#significancy threshold
thresh= 0.05/ncol(case)

#init output vectors
RES.all = NULL
RES.brain = NULL
RES.up = NULL
RES.down = NULL

# get associations for every gene
for(i in 1:ncol(case)){
  if(i%%100 == 0) cat('Gene ', i, '\n')
  #both controls
  expr = c(case[,i],controls[,i])
  model <- glm(Y ~ expr, family=binomial(link='logit'))
  pvalue = coef(summary(model))[,4][2]
  if(!is.na(pvalue)){ if(pvalue < thresh){
    RES.brain = rbind(RES.brain, c(colnames(case)[i],pvalue,model$coefficients['expr']))
    RES.all = rbind(RES.all, c(colnames(case)[i],pvalue,model$coefficients['expr']))
    if(model$coefficients['expr'] > 0){
      RES.up = rbind(RES.up,c(colnames(case)[i],pvalue,model$coefficients['expr']))
    }else{
      RES.down = rbind(RES.down,c(colnames(case)[i],pvalue,model$coefficients['expr']))
    }
  }else{
    RES.all = rbind(RES.all, c(colnames(case)[i],pvalue,model$coefficients['expr']))
  }
  }
}

save(RES.brain,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/06.brain_bp_signif.Rdata')

```
```{r, echo=F}
load('/wdata/kvervier/pharmacogenomics/output/WTCCC/06.brain_bp_signif.Rdata')
```

From this analysis, we found `r nrow(RES.brain)` genes with a significant different gene expression between cases and controls.
We also plot two examples of such genes (one up- and one down-regulated in cases).
```{r, fig.height = 10, fig.width = 12, fig.align = "center"}
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000001_expr_filtered_indiv_excluded_SNP.Rdata')
controls = EXPR2
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000002_expr_filtered_indiv_excluded_SNP.Rdata')
controls = rbind(controls,EXPR2)
# load bipolar disorder
load(file='/wdata/kvervier/pharmacogenomics/output/WTCCC/brain_EGAD00000000003_expr_filtered_indiv_excluded_SNP.Rdata')
case = EXPR2
# get all genes with diff expression:
#init labels
Y = c(rep(1,nrow(case)),rep(-1,nrow(controls)))
Y = factor(Y)
Y = factor(ifelse(Y==1,'bipolar','control'))

##################################################
# get one example for up and down regulated genes --> smaller pvalues
down.gene = "DLG3" #RES.down[which.min(RES.down[,2]),1]
up.gene = "BCL11B"#RES.up[which.min(RES.up[,2]),1]

# init df
expr = c(case[,down.gene],controls[,down.gene])
df = data.frame(group=Y,expr = scale(expr),gene = rep(down.gene,length(expr)))
df = rbind(df,data.frame(group=Y,expr = scale(c(case[,up.gene],controls[,up.gene])),gene = rep(up.gene,length(expr))))

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <-ggplot(df, aes(x=gene,y= expr, fill=group))

p <- p + geom_boxplot()+ ylab('predicted gene expression') + 
  theme_grey(base_size = 30) + theme(axis.text.x=element_text(colour="black")) +
  theme(legend.title = element_text(size=20, face="bold")) + 
  scale_fill_manual(values=cbPalette[5:6],guide = guide_legend(title = "group")) +
  theme(legend.text = element_text(size = 25)) +
  theme(legend.key.size = unit(1.5, "cm")) + xlab('') + theme(axis.text.x = element_text(angle = 60, hjust = 1))
p 

```

Here is the full distribution of significant genes by their direction of effect between cases and controls.

```{r, fig.height = 6, fig.width = 20, fig.align = "center"}

# load list of genes and estimated effect
load('/wdata/kvervier/pharmacogenomics/output/WTCCC/06.brain_bp_signif.Rdata')

RES.brain = RES.brain[order(as.numeric(RES.brain[,2]),decreasing=TRUE),]

# look at known 'brain genes'
genes = read.table('/home/kvervier/varann/aim1/output/brain-gene-analysis/brain_gene.db',header=TRUE)
#
idx.lit = which(RES.brain[,1] %in% as.character(genes[,3]))

# pubmed for bipolar
source('/home/kvervier/STUDIES/programing_with_literature/functions.R')
cit.l = NULL
for(i in 1: nrow(RES.brain)){
  tmp = getIDs(paste(RES.brain[i,1],"AND bipolar [abstract]",sep=' '))
  cit.l = c(cit.l,length(tmp))
}

idx.lit = sort(unique(c(idx.lit,which(cit.l > 0))))


tmp = -log10(as.numeric(RES.brain[,2])) * sign(as.numeric(RES.brain[,3]))
names(tmp) = RES.brain[,1]

tmp = sort(tmp)

# from 06.wtccc_tisan.R
thresh=5.129784e-06

library(RColorBrewer)
colors <- colorRampPalette(c("#67a9cf", "#f7f7f7","#ef8a62"))(length(tmp))

#png('~/Documents/talks/comp_psych_2018/pic/deg.png',width =2050,height = 400)
x<-barplot(tmp, xaxt="n",col=colors)
title(ylab='Association strength (signed pvalue)', cex.lab=1.5)
labs <- names(tmp)
labs[-idx.lit] = ''

text(cex=1.3, x=x-1.7, y=-37.5, labs, xpd=TRUE, srt=45)
abline(h=-log10(thresh),lwd=2,lty=2)
abline(h=log10(thresh),lwd=2,lty=2)

```

We also include an interactive table to help browsing the list of genes in temrs of co-citations in Pubmed for brain, neuron and psychiatric terms:

```{r}
library(DT)

load('/wdata/kvervier/pharmacogenomics/output/WTCCC/06.brain_bp_signif_anno.Rdata')
rownames(DB) = NULL
colnames(DB)[1] = 'gene.name'
datatable(DB)

```

# **Pharmacogenomics for bipolar disorder**

In this section, we propose to combine the list of DDE genes found in the previous section with the CMAP (clue.io) database.

Given a list of genes, one can found drugs with comparable effects on expression.

Here, the idea is to identify potnetial new treatments for bipolar disorder. Therefore, we mostly focus on drugs that show the strongest dissimilarity in temrs of gene signatures.

```{r, fig.align = "center"}
knitr::include_graphics('/home/kvervier/Documents/talks/comp_psych_2018/pic/deg_drug.png')
```

Results from CMAP can be found in */wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/all_up_down_bp.gct*, where the last column represents the correlation between dde genes in bipolar disorder and effect of drugs.

In order to validate the strongest candidate drugs, we proposed to use literature mining on drugs already associated with bipolar disorder.

```{r}

source('/home/kvervier/STUDIES/programing_with_literature/functions.R')

# load list of drugs
l1000_db = readLines('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/all_up_down_bp.gct')
tmp = strsplit(l1000_db,split = '\t')
merge.tmp=do.call('rbind',tmp)

l1000_db = merge.tmp[-1,]
colnames(l1000_db) = merge.tmp[1,]
l1000_db = as.data.frame(l1000_db)

l1000_db$`UPTAG-summary` = as.numeric(as.character(l1000_db$`UPTAG-summary` ))
```

We actually get co-citations count for each drug name and three key words: 'bipolar disorder', 'psychiatr', and 'neuro'.
```{r, eval =F, echo =T}

#get pubmed papers with co-citations
BP = list()
PSYCH = list()
NEURO = list()
CALC = list()
for(i in 1:nrow(l1000_db)){
  # get cocitations with bipolar disorder
  tmp = getIDs(paste(l1000_db$name[i],"AND bipolar disorder [abstract]",sep=' '))
  if(!is.null(tmp)){
    BP[[as.character(l1000_db$name[i])]] = tmp
  }else{
    BP[[as.character(l1000_db$name[i])]] = NA # length will be one, but given we require at least 2 citations, it will automatically filter them
  }
  # get cocitations with psychiatry
  tmp = getIDs(paste(l1000_db$name[i],"AND psychiatr [abstract]",sep=' '))
  if(!is.null(tmp)){
    PSYCH[[as.character(l1000_db$name[i])]] = tmp
  }else{
    PSYCH[[as.character(l1000_db$name[i])]] = NA
  }
  # get cocitations with neurology
  tmp = getIDs(paste(l1000_db$name[i],"AND neuro [abstract]",sep=' '))
  if(!is.null(tmp)){
    NEURO[[as.character(l1000_db$name[i])]] = tmp
  }else{
    NEURO[[as.character(l1000_db$name[i])]] = NA
  }
  # get cocitations with calcium
  tmp = getIDs(paste(l1000_db$name[i],"AND calcium [abstract]",sep=' '))
  if(!is.null(tmp)){
    CALC[[as.character(l1000_db$name[i])]] = tmp
  }else{
    CALC[[as.character(l1000_db$name[i])]] = NA
  }
}
# need to keep drug names for matching
save(BP,PSYCH,NEURO,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query.Rdata')
save(CALC,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_calcium.Rdata')

```

We also need to compute the number of co-citations, and also check for drugs that could be missing (failed Pubmed queries).

```{r, echo=T , eval=F}
load('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query.Rdata')
BP.l = sapply(BP,length)
PSYCH.l = sapply(PSYCH,length)
NEURO.l = sapply(NEURO,length)
CALC.l = sapply(CALC,length)
# if some PubMed queries failed, filter them
if(length(BP.l) != length(PSYCH.l) | length(BP.l) != length(NEURO.l) ){
  idx = intersect(names(NEURO),intersect(names(BP),names(PSYCH)))
  cat(length(idx),'out of',nrow(l1000_db),'drugs were successfully queried for the 3 keywords.\n')
  idx.bp = which(names(BP)%in%idx)
  idx.neuro = which(names(NEURO)%in%idx)
  idx.psych = which(names(PSYCH)%in%idx)
  BP.l = BP.l[idx.bp]
  PSYCH.l = PSYCH.l[idx.psych]
  NEURO.l = NEURO.l[idx.neuro]
}

DB = rbind(BP.l,PSYCH.l,NEURO.l)
save(DB,BP.l,PSYCH.l,NEURO.l,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_length.Rdata')
save(CALC.l,CALC,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_calcium.Rdata')

```

Using a binning strategy (20 groups), we derive a bin enrichment score based on the number of drugs with at least 25 co-citations with bipolar disorder.

```{r}
load('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_length.Rdata')
# is there more bp-drugs in the bin, than in the rest of the drugs
nbins = 20
idx.pval = cut(l1000_db$`UPTAG-summary`,right = TRUE,breaks = quantile(l1000_db$`UPTAG-summary`,probs = seq(0,1,by = 1/nbins)))
l1000_db$bins = idx.pval
#co-citations threshold below which a drug is not considered as related to bipolar disorder
thresh = 25
# init working bin
work.bin = NULL
rem.dat = l1000_db
rem.l = BP.l
#
PVAL = NULL
PVAL.non.cumul = NULL
REL.ENR = NULL
ENR = NULL
for(bin in levels(l1000_db$bins)){
  work.bin = c(work.bin,BP.l[l1000_db$name[which(l1000_db$bins == bin)]])
  idx.del = which(rem.dat$bins == bin)
  rem.dat = rem.dat[-idx.del,]
  rem.l = rem.l[-idx.del]
  mat = matrix(c(sum(work.bin>thresh) , sum(work.bin < thresh),sum(rem.l>thresh),sum(rem.l<thresh)),ncol=2)
  mat.non.cumul = matrix(c(sum(BP.l[l1000_db$name[which(l1000_db$bins == bin)]]>thresh) , sum(BP.l[l1000_db$name[which(l1000_db$bins == bin)]] < thresh),
                           sum(BP.l[l1000_db$name[which(l1000_db$bins != bin)]]>thresh),sum(BP.l[l1000_db$name[which(l1000_db$bins != bin)]]<thresh)),ncol=2)
  PVAL = c(PVAL,fisher.test(mat)$p.val)
  PVAL.non.cumul = c(PVAL.non.cumul,fisher.test(mat.non.cumul)$p.val)
  REL.ENR = c(REL.ENR,(mat[1,1]/mat[2,1])/(mat[1,2]/mat[2,2]))
  ENR = c(ENR,mat[1,1]/mat[2,1])
}

bp.drug.name = names(which(BP.l[l1000_db$name[which(l1000_db$bins == levels(l1000_db$bins)[1])]] > thresh))

save(PVAL,ENR,REL.ENR,file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/02.rel_enrichment.Rdata')
#get top drugs (not only BP related)
write.table(names(BP.l[l1000_db$name[which(l1000_db$bins == levels(l1000_db$bins)[1])]]),file='/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/top_drugs.txt',quote = FALSE,row.names = FALSE,col.names = FALSE)

```

```{r}
library(DT)

load('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_length.Rdata')
load('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/pubmed_query_calcium.Rdata')

top_drugs = read.table('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/top_drugs.txt',header=F)

DB = data.frame('drug.name'=top_drugs)
rownames(DB) = DB[,1]

# get drug families

DB$description = rep(0,nrow(DB))
idx = match(rownames(DB),l1000_db$name)
if(sum(is.na(idx)) == 0){
  DB$description = l1000_db$description[idx] 
  }else{
    DB$description[which(!is.na(idx))] = l1000_db$description[idx[-which(is.na(idx))]] 
  }


DB$bp.citations = rep(0,nrow(DB))
idx = match(rownames(DB),names(BP.l))
if(sum(is.na(idx)) == 0){
  DB$bp.citations = BP.l[idx] 
  }else{
    DB$bp.citations[which(!is.na(idx))] = BP.l[idx[-which(is.na(idx))]] 
  }

  DB$psych.citations = rep(0,nrow(DB))
idx = match(rownames(DB),names(PSYCH.l))
if(sum(is.na(idx)) == 0){
  DB$psych.citations = PSYCH.l[idx] 
  }else{
    DB$psych.citations[which(!is.na(idx))] = PSYCH.l[idx[-which(is.na(idx))]] 
    }

 DB$calcium.citations = rep(0,nrow(DB))
idx = match(rownames(DB),names(CALC.l))
if(sum(is.na(idx)) == 0){
  DB$calcium.citations = CALC.l[idx] 
  }else{
    DB$calcium.citations[which(!is.na(idx))] = CALC.l[idx[-which(is.na(idx))]] 
  }

DB$neuro.citations = rep(0,nrow(DB))
idx = match(rownames(DB),names(NEURO.l))
if(sum(is.na(idx)) == 0){
  DB$neuro.citations = NEURO.l[idx] 
  }else{
    DB$neuro.citations[which(!is.na(idx))] = NEURO.l[idx[-which(is.na(idx))]] 
    }

rownames(DB) = NULL
colnames(DB)[1] = 'drug.name'
datatable(DB)

```


Plot the enrichment for bipolar drugs with know bipolar drugs in the first top 5% bin.

```{r,  fig.align = "center",fig.height = 12, fig.width = 20}
#fig.height = 12, fig.width = 20,
#load('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/02.rel_enrichment.Rdata')

library(ggplot2)
# relative enrichment
df4 = data.frame('asso'=REL.ENR,'bin'=levels(l1000_db$bins))
df4$bin = factor(df4$bin,levels = unique(df4$bin))

library(RColorBrewer)
library(gplots)
color.palette  <- rev(colorRampPalette(c("#67a9cf", "#f7f7f7","#ef8a62"))(nlevels(l1000_db$bins)-1))

#png('/wdata/kvervier/pharmacogenomics/output/WTCCC/CMAP/literature_enrichment_bp_rel_enrichment.png',width = 800,height = 650)
par(mar=c(4, 7, 4, 5) + 0.1)
dodge <- position_dodge(width = 0.9)
d <- ggplot() +
  geom_bar(data=df4, aes(x=bin,y=asso,fill=bin),stat='identity') + ylab('literature enrichment \nfor bipolar disorder drugs') +
  theme_grey(base_size = 30) + theme(axis.text.x=element_text(colour="black")) +
  theme(legend.position="none") +
  scale_fill_manual(name = "bins",values=color.palette) + 
  scale_x_discrete(name='drug dissimilarity with BP genes (quantile)',breaks=df4$bin[c(1,6,11,16,nrow(df4))],labels=rev(c('100%','75%','50%','25%','top 5%')))+ 
  geom_hline(yintercept = 1,linetype = 2)  + annotate("text", size = 7,x = 4, y = seq(1.2,1.9,length.out = length(bp.drug.name)), label = rev(bp.drug.name),hjust = 0)+
  geom_segment(aes(xend=rep(2,length(bp.drug.name)), x=rep(3.6,length(bp.drug.name)), yend=seq(1.3,1.6,length.out = length(bp.drug.name)), y=seq(1.2,1.90,length.out = length(bp.drug.name))), arrow = arrow(length = unit(0.2, "cm")))
d

```