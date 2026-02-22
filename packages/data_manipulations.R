###############################################################################
#' The following script contains various functions for common data manipulations
#'    of microbiome-related and metadata feature tables.
#'  
#' Last updated: 01/01/2023
#'  
#' Implemented functions:  
#' - remove_rare_features
#' - clr_transform
###############################################################################

#' Remove rare features (columns) based on either a prevalence cutoff (i.e. % of
#'    non-zero values) or a mean-abundance cutoff, or both (in which case they 
#'    will be applied consecutively). 
#'    
#' @param feat_table Any feature table representing abundances (taxonomic/
#'    functional...)
#' @param filter_method One of 'prevalence', 'avg_abundance', 'both'
#' @param prevalence_cutoff The minimum portion of samples that should have a 
#'    non-zero value, to qualify the feature as non-rare (between 0 and 1, 
#'    default: 0.1). 
#' @param avg_abundance_cutoff Features with an average abundance (over samples)
#'    less than this threshold are considered rare and discarded.
#' 
#' @return An updated feature table 
remove_rare_features <- function(feat_table, 
                                 filter_method = 'both', 
                                 prevalence_cutoff = 0.1, 
                                 avg_abundance_cutoff = 0.005) {
  # Required libraries loading
  require(dplyr)
  
  # Argument verifications 
  if (prevalence_cutoff < 0 | prevalence_cutoff > 1) 
    error('Provided an invalid prevalence_cutoff value')
  if (avg_abundance_cutoff < 0) 
    error('Provided an invalid avg_abundance_cutoff value')
  
  # Initialize table with all current features
  new_feat_table <- feat_table
  n_samples <- nrow(feat_table)
  n_taxa_before_filter <- ncol(feat_table)-1
  
  # Prevalence calculations (number of non-zero values per feature)
  if (filter_method %in% c('prevalence', 'both')) {
    frequencies <- colSums(new_feat_table[,-1]>0) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, frequencies > prevalence_cutoff)]
  }
  
  # Average abundance calculations
  if (filter_method %in% c('avg_abundance', 'both')) {
    avg_abundances <- colSums(new_feat_table[,-1]) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, avg_abundances > avg_abundance_cutoff)]
  }
  
  n_taxa_after_filter <- ncol(new_feat_table)-1
  message(n_taxa_after_filter, ' of ',n_taxa_before_filter, ' taxa were filtered.')
  new_feat_table = new_feat_table %>% column_to_rownames("sample_id")
  norm_new_feat_table = (new_feat_table / rowSums(new_feat_table)) * 100
  return(norm_new_feat_table %>% rownames_to_column("sample_id"))
}

#' CLR-transform a feature table.  
#'    
#' @param feat_table Any feature table representing abundances (taxonomic/
#'    functional...). First column is expected to hold sample IDs.
#' @param pseudo_count Added to all values in the table to deal with zeros.
#' 
#' @return A CLR-transformed feature table 
clr_transform <- function(feat_table, pseudo_count = 0.0000001) {
  # Required libraries loading
  require(vegan)
  require(dplyr)
  
  # CLR transformation for rows, with pseudocount
  clr_df <- decostand(feat_table %>% select(-1), method = 'clr', MARGIN = 1, pseudocount = pseudo_count)
  clr_df <- bind_cols(feat_table %>% select(1), clr_df)
  
  return(clr_df)
}