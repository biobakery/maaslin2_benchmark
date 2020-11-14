##################################
# metagenomeSeq Default Pipeline #
##################################

# Install or Load Required Packages
if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'devtools', 'reshape2')
if(! require("metagenomeSeq")) {
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("metagenomeSeq")
  devtools::install_github('jnpaulson/metagenomeSeq', force=TRUE)
}
library(metagenomeSeq)

##################################
# Fit metagenomeSeq To A Dataset #
##################################

fit.metagenomeSeq <- function(features, metadata, libSize, ID, transformation, MultTestCorrection){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default metagenomeSeq model. Use NONE.')
  
  design <- model.matrix(~., data=metadata)

  # Collect data and normalize
  count_table <- t(as.matrix(features)) 
  mgsdata <- newMRexperiment(counts = count_table)
  mgsp <- cumNormStat(mgsdata)
  mgsdata <- cumNorm(mgsdata, mgsp)
  
  # Fit model
  # Random Effects Adjustment
  if(!length(ID)==length(unique(ID))){
    fit <- fitZig(obj=mgsdata,mod=design, useMixedModel=TRUE,block=ID)
  } else{
    fit <- fitZig(obj=mgsdata,mod=design)
  }
  
  # Collect Output
  if(dim(metadata)[2]>1){
    coef<-fit$fit$coefficients[,!colnames(fit$fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    pval<-fit$eb$p.value[,!colnames(fit$fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    coef.vector<-rename.features(coef, 'coef')
    pvalue.vector<-rename.features(pval, 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    rownames(DD)<-NULL;
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
  } 
  else{
    coef<-fit$fit$coefficients[,!colnames(fit$fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    pval<-fit$eb$p.value[,!colnames(fit$fit$coefficients) %in% c("(Intercept)", "scalingFactor")]
    DD<-cbind.data.frame(coef,pval)
    colnames(DD)<-c('coef', 'pval')
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    rownames(DD)<-NULL;
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
  }
  return(DD)
}

###########################################
# Fit metagenomeSeq To A List of Datasets #
###########################################

list.metagenomeSeq<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, .export=c("fit.metagenomeSeq", "rename.features", "LOG", "LOGIT", "AST"),
          .packages = c("metagenomeSeq", "car", "dplyr", "reshape2"), .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.metagenomeSeq(features, metadata, libSize, ID, transformation, MultTestCorrection)
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


