##########################
# edgeR Default Pipeline #
##########################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("edgeR")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
}
library(edgeR)

##########################
# Fit edgeR To A Dataset #
##########################

fit.edgeR = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default edgeR model. Use NONE.')
  
  d <- DGEList(counts=t(features))
  d <- edgeR::calcNormFactors(d, method='TMM')
  
  # Random Effect Adjustment
  if(!length(ID)==length(unique(ID))){
    stop('edgeR random effect model is currently not implemented.')
  }
  design <- model.matrix(~., data=metadata)
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  fit <- glmFit(d,design)
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalMatrix<-get_pval_edgeR(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  else{
    fit<-glmLRT(fit, 2)
    coef<-fit$coefficients[,-1]
    pval<-fit$table$PValue
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

# Get P-values from An edgeR Fit

get_pval_edgeR<-function(fit){
  mat<-as.matrix(fit$coefficients)
  rownames<-rownames(mat)
  colnames<-colnames(mat)
  n<-dim(fit$coefficients)[2]
  List <- list()
  for(i in 1:n){
    List[[i]] <- glmLRT(fit, i)$table$PValue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-rownames
  colnames(Matrix)<-colnames
  return(Matrix)
}

###################################
# Fit edgeR To A List of Datasets #
###################################

list.edgeR<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.edgeR", "rename.features", "get_pval_edgeR"),
          .packages = c("edgeR", "dplyr", "reshape2"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.edgeR(features, metadata, libSize, ID, transformation, MultTestCorrection)
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
