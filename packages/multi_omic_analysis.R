###############################################################################
#' The following script contains various functions for the analysis of 
#' multiple microbiome-related omics / feature tables.
#'  
#' Last updated: 17/09/2024
#' 
#' Implemented functions:
#' - run_mantel_test
#' - run_procrustes_test
###############################################################################



#' Run a Mantel test on two distance matrices
#'
#' @param distancesX A microbiome distance matrix (column names correspond 
#'    to sample id's, and first column also lists sample id's).
#' @param distancesY As above - for the second omic
#' @param sample_id_column The name of the sample_id column
#' @param corr_method Either "spearman", "kendall" or "pearson". Mantel 
#'    supports all 3.
#' @param n_permutations Number of permutations for significance test
#'
#' @return A list including a (1) ggplot object (mantel statistic) and (2) a 
#'    result table with statistic value and p-value.
#' @export
run_mantel_test <- function(distancesX, 
                            distancesY, 
                            sample_id_column = "sample_id",
                            corr_method = "spearman",
                            n_permutations = 499,
                            seed = 1234,
                            quiet = FALSE) {
  
  require(tibble)
  require(vegan)
  
  # Data verification 
  if (!sample_id_column %in% colnames(distancesX))
    stop("sample_id_column not found in distancesX")
  if (!sample_id_column %in% colnames(distancesY))
    stop("sample_id_column not found in distancesY")
  
  distX_samples <- distancesX[sample_id_column] %>%
    pull(sample_id_column)
  distY_samples <- distancesY[sample_id_column] %>%
    pull(sample_id_column)
  shared_samples <- intersect(distX_samples, distY_samples)
  
  if (length(shared_samples) == 0) stop("No shared samples between matrices")
  
  # Filter matrices to only contain shared samples
  if (length(shared_samples) != length(distX_samples) | length(shared_samples) != length(distY_samples)) {
    distancesX <- distancesX%>%
      filter(!!sym(sample_id_column) %in% shared_samples) %>%
      select(all_of(c(sample_id_column, shared_samples)))
    distancesY <- distancesY%>%
      filter(!!sym(sample_id_column) %in% shared_samples) %>%
      select(all_of(c(sample_id_column, shared_samples )))
    if (!quiet) message(str_c(length(shared_samples), " shared samples will be used for analysis\n"))
  }
  
  # Convert to matrix format
  distX_matrix <- distancesX %>%
    column_to_rownames(sample_id_column)
  
  distY_matrix <- distancesY%>%
    column_to_rownames(sample_id_column)
  
  # Make sure identical sample ordering
  distX_matrix <- distX_matrix[shared_samples, shared_samples]
  distY_matrix <- distY_matrix[shared_samples, shared_samples]
  
  # Make sure no NA's in matrices
  if (sum(is.na(distX_matrix)) > 0) stop("Found NAs in distance matrix X. No NAs allowed.")
  if (sum(is.na(distY_matrix)) > 0) stop("Found NAs in distance matrix Y. No NAs allowed.")
  
  # Calculate Mantel correlations
  set.seed(seed)
  mantel_results <- vegan::mantel(distX_matrix, distY_matrix, method = corr_method, permutations = n_permutations)
  
  result_table <- data.frame(p_value = mantel_results$signif, 
                             mantel_corr = mantel_results$statistic,
                             n_shared_samples = length(shared_samples))
  
  # Plot mantel statistic against shuffled versions
  p <- ggplot(data.frame(corr_stat = mantel_results$perm), aes(x = corr_stat)) +
    geom_histogram(color = 'black', fill = 'grey80', alpha = 0.4, bins = round(n_permutations / 10)) +
    geom_vline(xintercept = mantel_results$statistic, color = 'cadetblue3', linewidth = 2, alpha = 0.8) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    xlab(paste0('Mantel statistic (',corr_method,')')) +
    ylab('Count') +
    ggtitle('Mantel statistic - shuffled vs. true') +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  return(list(results = result_table, plot = p))
}


#' Run a procrustes test on two omics, using distance matrices.
#'   The function returns a procrustes plot visualizing the superimposition and
#'   a p-value as returned from a permutation-based test (implemented in the 
#'   'vegan' package).
#' Note that significant p-values do not always imply "pretty" plots...
#'
#' @param distancesX A microbiome distance matrix (column names correspond 
#'    to sample id's, and first column also lists sample id's).
#' @param distancesY As above - for the second
#' @param omic_name_X Name of the first omic (for plot)
#' @param omic_name_Y Name of the second omic (for plot)
#' @param dist_metric Name of the distance metric used
#' @param sample_id_column The name of the sample_id column
#' @param n_permutations Number of permutations for significance test
#' @param quiet Print messages?
#' @param pcoa_n_dimensions Number of PCoA dimensions to use for procrustes 
#' @param seed For reproducability
#'
#' @return A list including a (1) ggplot object and (2) a p-value.
#' @export
#' 
run_procrustes_test <- function(distancesX, 
                                distancesY, 
                                omic_name_X = "omic X",
                                omic_name_Y = "omic Y",
                                dist_metric = "--", 
                                sample_id_column = 'sample_id', 
                                n_permutations = 499, 
                                quiet = FALSE,
                                pcoa_n_dimensions = 10,
                                transparency_reflecting_residuals = FALSE,
                                seed = 1234) {
  
  # Required libraries loading
  require(vegan)
  
  # Data validations
  # TODO
  
  # Extract only common samples
  common_samples <- intersect(distancesX[[sample_id_column]], distancesY[[sample_id_column]])
  if(!quiet) message('Found ', length(common_samples), ' samples common to both omics')
  d1 <- distancesX %>% filter(!!as.symbol(sample_id_column) %in% common_samples) %>% column_to_rownames(var = sample_id_column)
  d2 <- distancesY %>% filter(!!as.symbol(sample_id_column) %in% common_samples) %>% column_to_rownames(var = sample_id_column)
  d1 <- d1[common_samples, common_samples] 
  d2 <- d2[common_samples, common_samples] 
  
  # Calculate PCoA
  pcoa1 <- cmdscale(d1, k = min(pcoa_n_dimensions, nrow(d1)-1), eig = TRUE)
  pcoa2 <- cmdscale(d2, k = min(pcoa_n_dimensions, nrow(d2)-1), eig = TRUE)
  # Check out: ordiplot(pcoa1, type = "text", display = 'sites')
  
  # Procrustes test (test -> symmetric, roatation -> not symmetric)
  set.seed(seed)
  proc_res <- protest(pcoa1, pcoa2, scale = T, scores = "sites", permutations = n_permutations)
  
  # Extract rotated data, for plotting (also use scale factor)
  rotations_df <- data.frame(
    sample_id = common_samples,
    d1_ax1 = proc_res$X[,1],
    d1_ax2 = proc_res$X[,2],
    d2_rotated_ax1 = proc_res$Yrot[,1] / proc_res$scale,
    d2_rotated_ax2 = proc_res$Yrot[,2] / proc_res$scale,
    residual = residuals(proc_res)
  )
  
  # Generate plot
  title <- paste0('Procrustes plot (P-value: ', proc_res$signif, ')')
  subtitle <- paste('Based on',dist_metric,'distances')
  omic_name_Y <- paste(omic_name_Y,'(superimposed)',sep='\n')
  omic_colors <- c('grey20', 'darkred'); names(omic_colors) <- c(omic_name_X, omic_name_Y)
  
  if (transparency_reflecting_residuals) {
    p <- ggplot(rotations_df) +
      geom_point(aes(x = d1_ax1, y = d1_ax2, color = omic_name_X, alpha = 1/residual), size = 2.2) +
      geom_segment(aes(x = d1_ax1, y = d1_ax2, xend = d2_rotated_ax1, yend = d2_rotated_ax2, alpha = 1/residual), color=omic_colors[2], arrow = arrow(length = unit(0.2,"cm"))) +
      geom_point(aes(x = d2_rotated_ax1, y = d2_rotated_ax2, color = omic_name_Y, alpha = 1/residual), size = 2.2) +
      scale_color_manual(values = omic_colors, breaks = names(omic_colors)) +
      theme_classic() +
      ggtitle(title) +
      labs(subtitle = subtitle) +
      xlab('PCoA Axis 1') +
      ylab('PCoA Axis 2') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      theme(legend.title = element_blank())
  } else {
    p <- ggplot(rotations_df) +
      geom_point(aes(x = d1_ax1, y = d1_ax2, color = omic_name_X), alpha = 0.8, size = 1) +
      geom_segment(aes(x = d1_ax1, y = d1_ax2, xend = d2_rotated_ax1, yend = d2_rotated_ax2), alpha = 0.2, color=omic_colors[2], arrow = arrow(length = unit(0.2,"cm"))) +
      geom_point(aes(x = d2_rotated_ax1, y = d2_rotated_ax2, color = omic_name_Y), alpha = 0.8, size = 1) +
      scale_color_manual(values = omic_colors, breaks = names(omic_colors)) +
      ggtitle(title) +
      labs(subtitle = subtitle) +
      xlab('PCoA Axis 1') +
      ylab('PCoA Axis 2') +
      theme_bw(base_family = "Gafata") + 
      theme(plot.title = element_text(hjust = 0.5, size = 20)) +
      theme(plot.subtitle = element_text(hjust = 0.5, size = 20)) +
      theme(axis.ticks = element_blank()) +
      theme(axis.text = element_blank()) +
      theme(legend.title = element_blank()) +
      theme(text = element_text(size = 20))
  }
  
  # Generate a plot of actual residuals compared to residuals from null models
  null_residuals <- data.frame()
  proc_shuf <- list()
  for(i in 0:10) {
    pcoa2_shuffled <- pcoa2
    
    # Shuffle all but 1st iteration
    if (i>0) {
      pcoa2_shuffled$points <-  pcoa2_shuffled$points[sample(1:nrow(pcoa2_shuffled$points)),]
      rownames(pcoa2_shuffled$points) <- sample(rownames(pcoa2_shuffled$points))
    }
    
    # Run procrustes
    proc_shuf[[as.character(i)]] <- procrustes(pcoa1$points, pcoa2_shuffled$points, scale = T, scores = "sites", symmetric = T)
    
    # Record residuals for plotting
    null_residuals <- bind_rows(
      null_residuals,
      data.frame(residual = residuals(proc_shuf[[as.character(i)]]), label = ifelse(i>0,'Shuffled','True'), shuffle_id = i)
    )
  }
  
  # Residuals plot
  p2 <- ggplot(null_residuals, aes(x = as.character(shuffle_id), y = residual, fill = label)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.05, height = 0, width = 0.2) + 
    theme_classic() +
    ggtitle('Residuals in true vs. shuffled models') +
    ylab('Point-wise residuals') +
    xlab('True/shuffled procrustes models') +
    scale_fill_manual(values = c('True' = 'cadetblue3', 'Shuffled' = 'grey80')) +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  # Generate a plot of the procrustes statistic compared to null models 
  # Note: Procrustes statistic = sqrt(1 - sum-of-squares-of-residuals)
  # Interpretation tip: The further away the "true" line is from the null values, the better.
  p3 <- ggplot(data.frame(t_stat = proc_res$t), aes(x = t_stat)) +
    geom_histogram(color = 'black', fill = 'grey80', alpha = 0.4, bins = 30) +
    geom_vline(xintercept = proc_res$t0, color = 'cadetblue3', linewidth = 2, alpha = 0.8) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    xlab('Procrustes statistic (sqrt(1-SSE))') +
    ylab('Count') +
    ggtitle('Procrustes statistic - shuffled vs. true') +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  results_df = data.frame(
    p_value = proc_res$signif, 
    var_explained = proc_res$t0 ^ 2, # to is the protest statistic, equivalent to a correlation 
    n_shared_samples = length(common_samples))
  
  list(results = results_df, plot = p, plot_residuals = p2, plot_procrustes_stat = p3)
}
