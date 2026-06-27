#' Run MaAsLin2 (differential abundance) in case there are more than 2 discrete 
#' values in the col in metadata_cols_to_analyze.
#' This function will get the differential abundance of each pair and plot all 
#' the results in one boxplot.
#' In case you want to do "one vs. all" comparison you should doe one-hote encoding
#' to the specific column and run the get_differential_abundance function.
#'    
#' @param feat_table Any table with continuous features.
#' @param disct_metadata A table with the metadata features to analyze. A 
#'   'sample id' column is expected, columns (metadata_cols_to_analyze) should 
#'   be discrete.
#' @param sample_id_column The name of the sample id column.
#' @param metadata_cols_to_analyze The discrete columns in the disct_metadata 
#'   to run MaAaLin2.
#' @param metadata_cols_fix_effects A vector of columns to add as fixed 
#'   effects (e.g. age, sex etc.). 
#'   If numeric with less than 5 or categorical, will perform one-hot encoding.
#'   NULL to ignore.
#' @param metadata_cols_rand_effects A vector of columns to add as random 
#'   effects (e.g. subject idnetifier in the case of multiple samples per 
#'   subject, or batches in case of a multi-batch dataset). NULL to ignore.
#' @param fdr_threshold FDR threshold used to define significant findings. 
#'   Plots will only be generated for significant findings.
#' @param output_dir the directory to keep MaAsLin2 output.
#' @param return_maaslin_plots Boolean. FALSE (default) will return our own 
#'   plots instead of MaAsLin2 plots.
#' @param keep_maaslin_plots Boolean. FALSE (default) will not generate 
#'   MaAsLin2's default output files.
#' @param normalization_method A normalization method supported by MaAsLin2. 
#'   Defaults to 'NONE'. Other options are: CLR, CSS, NONE, TMM
#' @param transformation_method A transformation method supported by MaAsLin2. 
#'   Defaults to 'LOG'. Other options are: LOGIT, AST, NONE
#' @param analysis_method The model MaAsLin2 will run. Defaults to 'LM'.
#'   Other options are: CPLM, NEGBIN, ZINB
#' @param category_name_map raw category values will be mapped to names 
#'   according to this named vector, for plotting purposes. Also defines the 
#'   order of categories. NULL to skip.
#' @param logscale_y_for_plot log scale the Y axis in the plots (default: FALSE).
#' @param ord vector of values corresponding to metadata column to analyse unique
#'   values that will be the order of the x axis valus.
#' 
#' @return A list of two elements:
#'   `sig_results` contains statistical test results. Data frame with rows that 
#'    represent significant results. 
#'   `all_results` contains statistical test results. Data frame with rows that 
#'    represent all the results (provide the option to perform additional FDR 
#'    correction if needed,)
#'   `plots` is a list of relevant plots: boxplots that represent the 
#'    differential abundance of feature regarding the discrete metadata feature.
#'    Plots generated for adjusted p-value < fdr_threshold.
#'    The names of the plots list are the unique values (strings of the rows in 
#'    the results dataframe.)
get_differential_abundance <- function(feat_table, 
                                       disct_metadata, 
                                       sample_id_column = 'sample_id', 
                                       metadata_cols_to_analyze = NULL,
                                       metadata_cols_fix_effects = NULL,
                                       fdr_threshold = 0.1,
                                       analysis_method = "LM",
                                       output_dir = NULL) {
  # Required libraries loading
  require(Maaslin2)
  require(ggplot2)
  require(ggsignif)
  require(grid)
  require(png)
  require(dplyr)
  require(ggpubr)
  require(purrr)
  
  # Convert sample_id column to rownames for MaAsLin2 run.
  disct_metadata <- disct_metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  feat_table <- feat_table %>% remove_rownames %>% column_to_rownames(var=sample_id_column)
  
  # Remove samples not in both tables, and order samples identically
  samples_to_include <- intersect(row.names(disct_metadata), row.names(feat_table))
  disct_metadata <- disct_metadata[samples_to_include,]
  feat_table <- feat_table[samples_to_include,,drop=FALSE]
  
  # One hot for categorical fixd cols
  new_fixed_effects = c()
  for (fixed_col in metadata_cols_fix_effects){
    if (length(unique(disct_metadata[,fixed_col])) > 3 & !is.numeric(disct_metadata[,fixed_col])){
      message("Too many unique values in ", fixed_col, " column. It might affect the FDR correction and omit some significant results!")
      tmp <- disct_metadata %>%
                rownames_to_column(sample_id_column) %>%
                dplyr::select(c(sample_id_column, all_of(fixed_col))) %>%
                mutate(!!as.name(fixed_col) := str_replace(!!as.name(fixed_col), "-", "_")) %>%
                tidyr::pivot_longer(cols = !(!!as.name(sample_id_column)), values_to = fixed_col) %>%
                filter(!is.na(fixed_col)) %>%
                mutate(dummy = TRUE) %>%
                tidyr::pivot_wider(id_cols = !!as.name(sample_id_column),
                                   names_from = !!as.name(fixed_col), 
                                   names_prefix = paste0(fixed_col, "_"), 
                                   values_from = dummy) %>%
                replace(is.na(.), FALSE) 
      fixed_col_new_names = tmp %>% 
                              dplyr::select(-sample_id_column) %>% 
                              names()
      disct_metadata <- disct_metadata  %>%
                          rownames_to_column(sample_id_column) %>% 
                          left_join(tmp, by = c(sample_id_column)) %>%
                          column_to_rownames(sample_id_column) 
      new_fixed_effects = append(new_fixed_effects, fixed_col_new_names)
    } else {
      new_fixed_effects = append(new_fixed_effects, fixed_col)
    }
  }
  
  # Rum MaAsLin2
  maaslin2_signif_results <- data.frame() 
  maaslin2_all_results <- data.frame() 
  
  for (col in metadata_cols_to_analyze) {
      message('Working on column: ', col)
      curr_fixed_effects = append(new_fixed_effects, col)
      tmp_metadata <- disct_metadata %>% 
                        filter(!is.na(!!as.name(col))) %>%
                          dplyr::select(c(col, all_of(new_fixed_effects)))
      
      # Get the final list of samples for which DA will be calculated
      samples_to_include <- intersect(row.names(tmp_metadata), row.names(feat_table))
      tmp_metadata <- tmp_metadata[samples_to_include,, drop = F]
      tmp_feat_table <- feat_table[samples_to_include,, drop = F]
      
      invisible(capture.output(tmp_results <- Maaslin2(
        input_data = tmp_feat_table, 
        input_metadata = tmp_metadata, 
        output = paste(output_dir, "/", col, sep=""),
        fixed_effects = curr_fixed_effects,
        analysis_method = analysis_method,
        max_significance = 0.1,
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        save_scatter = FALSE,
        save_models = FALSE)))
      
      # Save all results
      maaslin2_all_results <- bind_rows(maaslin2_all_results, tmp_results$results)
    }
    maaslin2_sig_results = maaslin2_all_results %>%
      arrange(pval) %>%
      mutate(fdr = p.adjust(pval, method = "BH")) %>%
      filter(metadata %in% metadata_cols_to_analyze) %>%
      filter(fdr <= fdr_threshold)

  return(list(sig_results = maaslin2_sig_results, all_results = maaslin2_all_results))
}
