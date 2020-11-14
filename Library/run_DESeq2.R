###########################
# DESeq2 Default Pipeline #
###########################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  }
library(DESeq2)

###########################
# Fit DESeq2 To A Dataset #
###########################

fit.DESeq2<-function(features, metadata, libSize, ID, transformation, MultTestCorrection){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DESeq2 model. Use NONE.')
  
  
  # Random Effect Adjustment
  if(!length(ID)==length(unique(ID))){
    stop('edgeR random effect model is currently not implemented.')
  }
  
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  
  
  # Fit Model
  x <- DESeqDataSetFromMatrix(countData = t(features), colData = metadata, design = formula)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(x), 1, gm_mean)
  x = estimateSizeFactors(x, geoMeans = geoMeans)
  fit <- DESeq(x)
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(coef(fit)[,-1], 'coef')
    pvalMatrix<-get_pval_DESeq2(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  else{
    coef<-coef(fit)[,-1]
    pval<-results(fit,name=resultsNames(fit)[2])$pvalue
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

# Get P-values from DESeq2 Fit
get_pval_DESeq2<-function(fit){
  List <- list()
  for(i in 1:length(resultsNames(fit))){
    List[[i]] <- results(fit,name=resultsNames(fit)[i])$pvalue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-names(fit)
  colnames(Matrix)<-resultsNames(fit)
  return(Matrix)
}

####################################
# Fit DESeq2 To A List of Datasets #
####################################

list.DESeq2<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.DESeq2", "rename.features", "get_pval_DESeq2"),
          .packages = c("DESeq2", "dplyr", "reshape2"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.DESeq2(features, metadata, libSize, ID, transformation, MultTestCorrection)
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



