############
# Spearman #
############

##################################### 
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr')

#############################
# Fit Spearman To A Dataset #
#############################

fit.Spearman  = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default Spearman model. Use NONE.')
  if (length(unique(metadata$Metadata1_TP))!=length((metadata$Metadata1_TP))) stop ('Spearman only valid for continuous metadata.')
  
  
  spe <- function(x){
    tryCatch({
    fit1 <- cor.test(x, as.vector(t(metadata)), method='spearman')
  }, error=function(err){
    fit1 <- try({cor.test(x, as.vector(t(metadata)), method='spearman')}) 
    return(fit1)
  })}
  
    # Collect results
   spes <- apply(features, 2, spe)
   output_df <- data.frame(pval = sapply(spes, function(x) x$p.value))
   output_df$coef <- sapply(spes, function(x) x$estimate)
   output_df$metadata<-colnames(metadata)
   output_df$feature<-colnames(features)
   output_df$qval<-as.numeric(p.adjust(output_df$pval, method = MultTestCorrection))
   output_df<-output_df[order(output_df$qval, decreasing=FALSE),]
   output_df<-dplyr::select(output_df, c('feature', 'metadata'), everything())
   rownames(output_df)<-NULL
   return(output_df)
}

######################################
# Fit Spearman To A List of Datasets #
######################################

list.Spearman <-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.Spearman"),
          .packages = "dplyr", 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.Spearman(features, metadata, libSize, ID, transformation, MultTestCorrection)
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
