########################################
# Compound Poisson Linear Model (CPLM) #
########################################
 
#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'pbapply', 'cplm', 'glmmTMB')

#########################
# Fit CPLM To A Dataset #
#########################

fit.CPLM <- function(features, 
                     metadata, 
                     libSize, 
                     ID, 
                     transformation, 
                     MultTestCorrection){
 
 #########################
 # Transformation if any #
 #########################
 
 if (transformation!='NONE') stop ('Transformation currently not supported for a default CPLM model. Use NONE.')
 
 #####################
 # Per-feature model #
 #####################
 
 paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
   
   ###############################
   # Extract features one by one #
   ###############################
   
   featuresVector <- features[, x]
   
   #################################
   # Create per-feature input data #
   #################################
   
   dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize, ID)
   formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
   
   ##############################################
   # Automatic library size adjustment for GLMs #
   ##############################################
   
   if(length(unique(libSize)) > 1){ # To prevent offsetting with TSS-normalized data 
     formula<-update(formula, . ~ . - offset(log(libSize)))
     }
   
   #######################
   # Random effect model #
   #######################
   
   if(!length(ID) == length(unique(ID))){ 
     formula<-update(formula, . ~ . + (1|ID))
     fit <- tryCatch({
       fit1 <- glmmTMB::glmmTMB(formula = formula,  
                                data = dat_sub, 
                                family = glmmTMB::tweedie(link = "log"), 
                                ziformula = ~0)
       }, error=function(err){
         fit1 <- try({glmmTMB::glmmTMB(formula = formula,  
                                       data = dat_sub, 
                                       family = glmmTMB::tweedie(link = "log"), 
                                       ziformula = ~0)}) 
         return(fit1)
         })
     
     ###################################
     # Summarize Coefficient Estimates #
     ###################################
     
     if (class(fit) != "try-error"){
       para<-as.data.frame(coef(summary(fit))$cond)[-1,-3]
       } else{
         print(paste("Fitting problem for feature", x, "returning NA"))
         para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=3))
         }
     colnames(para)<-c('coef', 'stderr', 'pval')
     para$metadata<-colnames(metadata)
     para$feature<-colnames(features)[x]
     }
   
   #######################
   # Fixed effects model #
   #######################
   
   else{ 
     fit <- tryCatch({
       fit1 <- cplm::cpglm(formula, 
                           data = dat_sub)
       }, error=function(err){
         fit1 <- try({cplm::cpglm(formula, 
                                  data = dat_sub)}) 
         return(fit1)
         })
     
     ###################################
     # Summarize Coefficient Estimates #
     ###################################
     
     if (class(fit) != "try-error"){
       para<-as.data.frame(summary(fit)$coefficients)[-1,-3]
       } else{
         print(paste("Fitting problem for feature", x, "returning NA"))
         para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=3))  
         }
     colnames(para)<-c('coef', 'stderr', 'pval')
     para$metadata<-colnames(metadata)
     para$feature<-colnames(features)[x]
     }
   return(para)
   })
 
 #################
 # Return output #
 #################
 
 paras<-do.call(rbind, paras)
 paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
 paras<-paras[order(paras$qval, decreasing=FALSE),]
 paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
 rownames(paras)<-NULL
 return(paras)   
}

##################################
# Fit CPLM To A List of Datasets #
##################################

list.CPLM<-function(physeq, transformation = 'NONE', MultTestCorrection = "BH"){
  foreach(physeq = physeq, 
          .export = c("fit.CPLM"), 
          .packages = c("dplyr", "pbapply", "cplm", "glmmTMB"),
          .errorhandling = "remove") %dopar% 
  {
    start.time<-Sys.time()
    features<-physeq$features
    metadata<-physeq$metadata
    libSize<-physeq$libSize
    ID<-physeq$ID
    DD<-fit.CPLM(features, metadata, libSize, ID, transformation, MultTestCorrection)
    DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
    wh.TP<-intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
    newname<-paste0(DD$pairwiseAssociation[wh.TP], "_TP")
    DD$pairwiseAssociation[wh.TP]<-newname
    DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
    stop.time<-Sys.time()
    time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
    DD$time<-time
    return(DD)
  }
}


#######################################
# ALTERNATIVE APPROACH NOW DEPRECATED #
#######################################

##############################################################################
# Fitting the GLMM (Profile Likelihood for Selecting Power Parameter + GLMM) #
##############################################################################

# MASS::glmmPQL(fixed = formula, random = list(~1|ID), family=statmod::tweedie(var.power = tweedie::tweedie.profile(formula, p.vec=seq(1.0,2.0,by=0.1), data=dat_sub)$p.max, link.power=0), data= dat_sub) (VERY VERY SLOW!)
# MASS::glmmPQL(fixed = formula, random = list(~1|ID), family=statmod::tweedie(var.power = 1.5, link.power=0), data= dat_sub) # Only runs for fixed p
