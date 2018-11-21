#-----------------------------------------------------------#
# Source for predictive models based on PrediXcan framework #
#                   on Neon, with a subset of genes         #
#-----------------------------------------------------------#
args = commandArgs(trailingOnly = T)
alpha=as.numeric(args[1])
chunk=as.numeric(args[2])

require(glmnet)
library(methods)

# load DGN data for training
data.path= '../data/DGN-data.Rdata'
system(paste('mkdir -p ../output/',tail(unlist(strsplit(data.path,split = '/')),1),sep=''))
en.dir = paste('../output/',tail(unlist(strsplit(data.path,split = '/')),1),'/',sep='')

load(data.path)

#gene info
gencodefile <- '../data/gencode'
gencode <- read.table(gencodefile)

#paste function
"%&%" = function(a,b) paste(a,b,sep="")

resultsarray <- array(0,c(ncol(DGN.trans),8))
dimnames(resultsarray)[[1]] <- colnames(DGN.trans)
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- en.dir %&% "working_" %&% tis %&% "_exp_10-foldCV_elasticNet_alpha" %&% alpha %&% "_chunk" %&% chunk %&% ".txt"
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene","SNP","chr","start","end","referenceAllele","effectAllele","beta")
workingweight <- en.dir %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_weights" %&% "_chunk" %&% chunk %&% ".txt"
write(weightcol,file=workingweight,ncol=8,sep="\t")

########################################
# load tissue-specific data (from TiSAn)
load('../data/tisan_brain_SNP.Rdata')

# check if order is the same
tmp1= tisan_anno[,1]
tmp2 = colnames(allele.count)
tmp3 = match(tmp1,tmp2)

allele.count = allele.count[,tmp3]
# find which SNPs need to be filtered (TiSAN =0)
idx.filt = which(tisan_anno[,5] == 0)
allele.count = allele.count[,-idx.filt] # TO TEST ?


  
set.seed(1001)

start = (chunk-1)*round(ncol(DGN.trans)/100) + 1
end = (chunk)*round(ncol(DGN.trans)/100)
if(chunk == 100) end =  ncol(DGN.trans)

for(i in start:end){
  if(i %% 10 == 0) cat(i,'/',ncol(DGN.trans),'\n')
  #expression data
  exppheno = DGN.trans[,i]
  exppheno <- scale(exppheno, center=T, scale=T)
  exppheno[is.na(exppheno)] <- 0
  #rownames(exppheno) <- rownames(exp.w.geno)
  ##run Cross-Validation over alphalist
  fit <- cv.glmnet(allele.count,exppheno,nfolds=10,alpha=alpha,keep=T) 
  fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
  best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, ...)
  cvm.best = best.lam[,1]
  lambda.best = best.lam[,2]
  nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
  
  ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
  ret[ret == 0.0] <- NA
  bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
  names(bestbetas) = rownames(ret)[which(!is.na(ret))] 
  
  pred.mat <- fit$fit.preval[,nrow.best]
  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    genename <- as.character(colnames(DGN.trans)[i])
    rsq <- res$r.squared
    pval <- res$coef[2,4]
    
    resultsarray[i,] <- c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)
    
    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- snp_locations[which(snp_locations$refsnp_id %in% bestbetalist),]
    betatable<-as.matrix(cbind(genename,bestbetainfo,bestbetas))
    #betafile<-cbind(genename,betatable[,1],betatable[,5],betatable[,6],betatable[,7]) ##output "gene","SNP","refAllele","effectAllele","beta"
    write(t(betatable),file=workingweight,ncolumns=8,append=T,sep="\t") # t() necessary for correct output from write() function
    
  }else{
    genename <- as.character(colnames(DGN.trans)[i])
    resultsarray[i,1] <- genename
    resultsarray[i,2:8] <- c('no-feature',NA,NA,NA,0,NA,NA)
  }
  write(resultsarray[i,],file=workingbest,ncolumns=8,append=T,sep="\t")
}
