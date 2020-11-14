##############################
# limmaVOOM Default Pipeline #
##############################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("limma")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
}
library(limma)

##############################
# Fit limmaVOOM To A Dataset #
##############################

fit.limmaVOOM = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default limmaVOOM model. Use NONE.')
  
  # Fit limmaVOOM
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  design <- model.matrix(~., data=metadata)
  y <- voom(x,design,plot=FALSE)
  
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
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalue.vector<-rename.features(fit$p.value[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  } 
  else{
    coef<-fit$coefficients[,-1]
    pval<-fit$p.value[,-1]
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

#######################################
# Fit limmaVOOM To A List of Datasets #
#######################################

list.limmaVOOM<-function(physeq, transformation='NONE', MultTestCorrection= 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.limmaVOOM", "rename.features"),
          .packages = c("limma", "dplyr", "reshape2"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.limmaVOOM(features, metadata, libSize, ID, transformation, MultTestCorrection)
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


