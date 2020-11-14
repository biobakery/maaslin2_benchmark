#############################################
# Vanilla LM with Library Size as Covariate #
#############################################

##################################### 
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
library('dplyr')
pacman::p_load('pbapply', 'car', 'nlme')

#################################
# Fit Linear Model To A Dataset #
#################################

fit.LM2 <- function(features, metadata, libSize, ID, transformation, MultTestCorrection){
  
  # Transformation
  if (!transformation %in% c('AST', 'LOG', 'LOGIT', 'NONE')) {
    stop ('Transformation should be one of AST/LOG/LOGIT/NONE for LM')
  }
  
  if (transformation =='LOG')   {
    features<-apply(features, 2, LOG)
  }

  if (transformation =='LOGIT')   {
    features<-apply(features, 2, LOGIT)
  }
  
  if (transformation =='AST')   {
    features<-apply(features, 2, AST)
  }
  
  if (transformation =='NONE')   {
    features<-features
  }
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    featuresVector <- features[, x]
    
    # Scrap All-Zero Features
    
    if(sum(featuresVector!=0)<1){
      print(paste("Cannot fit model to all zeroes for feature", x, "returning NA"))
      para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
      colnames(para)<-c('coef', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
      rownames(para)<-NULL
    }
    else {
      
      # Fit Model
      dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize, ID)
      dat_sub2<- dat_sub[, !colnames(dat_sub) %in% c('expr', 'ID')]
      formula<-as.formula(paste("expr ~ ", paste(colnames(dat_sub2), collapse= "+")))
    
      # Random effect adjustment
      if(!length(ID)==length(unique(ID))){
        fit <- tryCatch({
          fit1 <- lme(formula, data = dat_sub, random= ~ 1 | ID)
        }, error=function(err){
          fit1 <- try({lme(formula, data = dat_sub, random= ~ 1 | ID)}) 
          return(fit1)
        })
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$tTable)[-1,-c(2:4)]
          para<-para[-nrow(para),] 
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
        else{
          print(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$features<-colnames(features)[x]
          rownames(para)<-NULL
        } 
      }
      
      else{
        fit <- tryCatch({
          fit1 <- glm(formula, data = dat_sub, family='gaussian')
        }, error=function(err){
          fit1 <- try({glm(formula, data = dat_sub, family='gaussian')}) 
          return(fit1)
        })
        
        if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
          para<-para[-nrow(para),] 
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        } 
        else{
          print(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$features<-colnames(features)[x]
          rownames(para)<-NULL
        }
      }
    }
    return(para)
  })    
   
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL;
  return(paras)  
}

##########################################
# Fit Linear Model To A List of Datasets #
##########################################

list.LM2<-function(physeq, transformation='NONE', MultTestCorrection = "BH"){
  foreach(physeq=physeq, .export=c("fit.LM2", "AST", "LOG", "LOGIT"),
          .packages=c("dplyr", "pbapply", "car", "nlme"), .errorhandling = 'remove') %dopar% 
  {
    start.time <- Sys.time()
    features<-physeq$features; 
    metadata<-physeq$metadata;
    libSize <- physeq$libSize;
    ID<-physeq$ID;
    DD<-fit.LM2(features, metadata, libSize, ID, transformation, MultTestCorrection)
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

