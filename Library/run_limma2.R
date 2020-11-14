########################################
# limma with Library Size as Covariate #
########################################

# Install or Load Required Packages
if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'car', 'reshape2')
if(! require("limma")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
}
library(limma)

##########################
# Fit limma To A Dataset #
##########################

fit.limma2 = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {
  
  if (!transformation %in% c('AST', 'LOG', 'LOGIT', 'NONE')) {
    stop ('Transformation should be one of AST/LOG/LOGIT/NONE for limma.')
  }
  
  design <- model.matrix(~., data=cbind.data.frame(metadata, libSize))
  
  # Transformation
  if (transformation =='LOG')   {
    x<-apply(features, 2, LOG)
  }

  if (transformation =='LOGIT')   {
    x<-apply(features, 2, LOGIT)
  }
  
  if (transformation =='AST')   {
    x<-apply(features, 2, AST)
  }
  
  if (transformation =='NONE')   {
    x<-features
  }
  
  # Convert to matrix, and transpose
  y<-t(as.matrix(x)) 
  
  # Fit limma 
  # Random Effect Adjustment
  if (!length(ID)==length(unique(ID))){
    dupcor <-  limma::duplicateCorrelation(y, design, block = ID)
    fit <- limma::lmFit(y,design, block = ID, correlation = dupcor$cor)
  } else{ 
    fit <- limma::lmFit(y,design)}
  
  # Empirical Bayes Adjustment
  fit <- limma::eBayes(fit)
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-c(1, ncol(fit$coefficients))], 'coef')
    pvalue.vector<-rename.features(fit$p.value[,-c(1, ncol(fit$coefficients))], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  } 
  else{
    coef<-fit$coefficients[,-c(1, ncol(fit$coefficients))]
    pval<-fit$p.value[,-c(1, ncol(fit$coefficients))]
    DD<-cbind.data.frame(coef,pval)
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  return(DD)
}

###################################
# Fit limma To A List of Datasets #
###################################

list.limma2<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.limma2", "rename.features", "LOG", "LOGIT", "AST"),
          .packages = c("limma", "dplyr", "car", "reshape2"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.limma2(features, metadata, libSize, ID, transformation, MultTestCorrection)
            DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
            wh.TP = intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
            newname = paste0(DD$pairwiseAssociation[wh.TP], "_TP")
            DD$pairwiseAssociation[wh.TP] <- newname;
            DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
            stop.time <- Sys.time()
            time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
            DD$time<-time
            return(DD)
          }
}


