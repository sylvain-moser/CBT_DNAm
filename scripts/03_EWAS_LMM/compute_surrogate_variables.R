#' Calculate surrogate variables for biological variable
#' 
#' @param dataset A dataframe where to find the variable of interest and the possible surrogates
#' @param feature The name of the column defining the variable of interest
#' @param surrogates The name of the columns defining the possible surrogate
#' @param scale indicate if features and surrogates have to be scale. Should be TRUE if scaling not done before
#' @return A dataframe with the surrogate variables for the complete cases of the input dataset
library(dplyr)
library(sva)
library(tidyr)

compute_surrogates=function(dataset,feature,surrogates_matrix,surrogate_df_path,scale=F){
#creeate the phenodata
  pheno=dataset %>% dplyr::select(all_of(feature))
# remove feature from the surrogate_matrix:
  surrogate_matrix=surrogates_matrix[rownames(surrogates_matrix) != feature,]
  #define models
  mod_g = model.matrix(as.formula(paste0("~",feature)),data=pheno)
  mod0_g = model.matrix(~1,data=pheno)
  svobj=sva(surrogates_matrix,mod_g,mod0_g)
  sv_df=as.data.frame(svobj$sv)
  save(sv_df,file=paste0(surrogate_df_path,feature,"_surrogates.RData"))
  return (sv_df)
}


get_surrogates_matrix=function(dataset,scale=F,n_var){
  library(dplyr)
  library(tidyr)
  variances= dataset %>% summarise_all(var)
  high_variance_surrogates=names(sort(variances,decreasing = T)[1:n_var])
  surrogates=dataset %>% dplyr::select(all_of(high_variance_surrogates)) %>% drop_na()
  if (scale==T){
    surrogates=dataset %>% mutate_all(.funs = ~scale(.,center = TRUE,scale = TRUE))
  }
  surrogates=as.matrix(surrogates)
  surrogates=WGCNA::transposeBigData(surrogates)
}

