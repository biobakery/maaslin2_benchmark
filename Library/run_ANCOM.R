#########
# ANCOM #
#########

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('devtools', 'tibble', 'exactRankTests', 'openxlsx', 'DT', 'dplyr', 'coin', 'compositions')
devtools::source_url("https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R?raw=TRUE")

###########################
# Fit ANCOM To A Dataset #
###########################

fit.ANCOM = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {

  if (transformation!='NONE') stop ('Transformation currently not supported for a default ANCOM model. Use NONE.')

  #############################################################################################
  # ANCOM standard pipeline for DA (adopted from https://github.com/FrederickHuangLin/ANCOM) #
  #############################################################################################
  
  ###############################
  # Step 1: ANCOM preprocessing #
  ###############################

  features_ancom<-as.data.frame(t(features))
  metadata_ancom<-rownames_to_column(metadata, 'ID')
  preprocess.ancom = feature_table_pre_process(feature_table = features_ancom, 
                                     meta_data = metadata_ancom, 
                                     sample_var = "ID", 
                                     group_var = NULL,
                                     out_cut = 0.05,
                                     zero_cut = 0.90,
                                     lib_cut = 1000,
                                     neg_lb = FALSE)
  feature_table = preprocess.ancom$feature_table 
  meta_data = preprocess.ancom$meta_data 
  struc_zero = preprocess.ancom$structure_zeros 
  
  
  #################
  # Step 2: ANCOM #
  #################

  if(!length(ID)==length(unique(ID))){
    res <- ANCOM(feature_table = feature_table, 
                 meta_data = meta_data, 
                 struc_zero = struc_zero,
                 main = colnames(metadata),
                 alpha = 0.05,
                 adj_formula = NULL,
                 rand_formula = "~ 1 | ID",
                 p_adj_method = MultTestCorrection) 
  } else{
    res <- ANCOM(feature_table = feature_table, 
                 meta_data = meta_data, 
                 struc_zero = struc_zero,
                 main = colnames(metadata),
                 alpha = 0.05,
                 adj_formula = NULL,
                 rand_formula = NULL,
                 p_adj_method = MultTestCorrection)
    }

  #########################################################
  # Step 3: Extract results and enforce meaningful format #
  #########################################################
  
  df <- data.frame(coef = res$out)
  df$pval<-1 # Fake p-values
  df$feature = colnames(features)
  df$metadata = names(metadata)
  df$pval[df$coef.detected_0.7] <- 0 # Fake p-values
  df$qval<-df$pval # Fake q-values
  df<-df[order(df$qval, decreasing=FALSE),]
  df<-dplyr::select(df, c('feature', 'metadata', 'pval', 'qval'))
  rownames(df)<-NULL;
  return(df)
}

###################################
# Fit ANCOM To A List of Datasets #
###################################

list.ANCOM<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq,
          .export="fit.ANCOM",
          .packages=c('tibble', 'exactRankTests', 'openxlsx', 'DT', 'dplyr', 'coin', 'compositions'),
          .errorhandling = 'remove') %dopar%
          {
            start.time <- Sys.time()
            features<-physeq$features;
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.ANCOM(features, metadata, libSize, ID, transformation, MultTestCorrection)
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



