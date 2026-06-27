get_abundance_plot = function(df_lst, res_df, sgb_2_gtdb, res_dir, nrow = 3){
  ord_sgb_2_gtdb = sgb_2_gtdb %>%
    mutate(sgb = str_replace(sgb, "t__", "")) %>%
    select(sgb, s)
  tested_sgb = str_remove(unique(res_df$sgb), "^t__")
  ord_vec = res_df %>% arrange(p.adj) %>% pull(species)
  filt_dfs_lst = list()
  for (organ in names(df_lst)){
    df = df_lst[[organ]]
    filt_df = df %>% 
      select(c("sample_id", tested_sgb)) %>%
      pivot_longer(-sample_id) %>%
      mutate(organ = organ)
    filt_dfs_lst[[organ]] = filt_df
  }
  
  conc_dfs = bind_rows(filt_dfs_lst$rectal, filt_dfs_lst$vaginal) %>%
    group_by(name) %>%
    mutate(same_person = n_distinct(organ) == 2 &
             all(c("vaginal", "rectal") %in% organ)) %>%
    ungroup() %>%
    left_join(ord_sgb_2_gtdb, by = c("name"= "sgb")) %>%
    mutate(s = str_replace_all(s, "_", " ")) %>%
    mutate(organ = str_to_title(organ)) %>%
    mutate(s = factor(s, levels = ord_vec)) %>%
    arrange(sample_id, name, s, organ) 
  maaslin_df = conc_dfs %>%
    pivot_wider(names_from = organ, values_from = value) %>%
    select(-name) %>%
    pivot_longer(cols = c(Rectal, Vaginal), names_to = "organ", values_to = "value") %>%
    pivot_wider(names_from = s, values_from = value) %>%
    mutate(sample_id = paste0(sample_id, "_", organ))
  res <- get_differential_abundance(feat_table = maaslin_df %>% select(c(-c("same_person", "organ"))),  
                                    disct_metadata = maaslin_df %>% select(c(c("sample_id", "same_person", "organ"))), 
                                    sample_id_column = 'sample_id', 
                                    metadata_cols_to_analyze = 'organ', # cat_cols,
                                    metadata_cols_fix_effects = NULL, 
                                    # normalization_method = "NONE",
                                    fdr_threshold = 0.1,
                                    output_dir = res_dir)
  y_position_df = conc_dfs %>%
    group_by(s) %>%
    summarise(y_position = max(value) - 5) %>% #  - 0.005
    ungroup()
  plot_res = res$all_results %>%
    mutate(signif = case_when(qval <= 0.001 ~ "***", 
                              qval <= 0.01 ~ "**",
                              qval <= 0.05 ~ "*",
                              TRUE ~ NA_character_)) %>%
    mutate(signif = case_when(is.na(signif) ~ NA_character_,
                              signif == "ns" ~ NA_character_,
                              TRUE ~ paste0(signif, case_when(coef < 0 ~ ",  >",
                                                              coef > 0 ~ ",  <",
                                                              TRUE ~ "")))) %>%
    mutate(feature = str_replace_all(feature, "\\.", " ")) %>%
    full_join(y_position_df, by = c("feature" = "s")) %>%
    mutate(group1 = "Rectal", group2 = "Vaginal", .y. = "value") %>%
    rename("compare_group" = "value", "s" = "feature") %>%
    mutate(s = factor(s, levels = ord_vec))
  # wilcox_results = conc_dfs %>%
  #   group_by(s) %>%
  #   wilcox_test(value ~ organ, paired = TRUE) %>%
  #   ungroup() %>%
  #   mutate(p.adj = p.adjust(p, method = 'fdr')) %>%
  #   mutate(signif = case_when(p.adj <= 0.001 ~ "***", 
  #                             p.adj <= 0.01 ~ "**",
  #                             p.adj <= 0.05 ~ "*",
  #                             TRUE ~ "ns")) %>%
  #   left_join(y_position_df, by = "s") %>%
  #   mutate(s = factor(s, levels = ord_vec))
  
  
  plot = ggplot(data = conc_dfs, aes(x = organ, y = value)) +
    geom_boxplot(aes(fill = organ), color = 'black', alpha = 0.7, outlier.shape = NA) +
    geom_point(size = 1, alpha = 0.3, fill = "black", color = 'black', shape = 21, position = position_jitter(width = 0.2, height = 0, seed = 42)) +
    theme_bw(base_family = "Gafata") +
    scale_fill_manual(values = c("#C1121F", "#669BBC")) +
    # stat_pvalue_manual( plot_res, label = "signif", y.position = "y_position", geom = "text", size = 15, vjust = 0.2, tip.length = 0.01, hide.ns = TRUE) +
    stat_pvalue_manual(plot_res, label = "signif", y.position = "y_position", size = 15, hide.ns = TRUE, vjust = 0) +
    facet_wrap(~ s, nrow = nrow, scale = "free_y") +
    labs(x = "Organ", y = "Relative abundance") +
    theme(legend.position  = "none",
          strip.background = element_rect(fill = "#FEFCFD"),
          axis.title       = element_text(size = 40),
          axis.title.y     = element_text(size = 50),
          axis.text.y      = element_text(size = 42),
          axis.text.x      = element_text(size = 45),
          strip.text       = element_text(size = 50)) 
  return(list("plot" = plot, "res" = plot_res))         
}

get_all_paired_ngd = function(ngd_lst, stats_results, step.inc = 0.08){
  species_levels = stats_results %>% pull(species)
  pair_levels = c("Vag-Rec\nDiff woman", "Vag-Vag", "Rec-Rec")
  lvl_map <- setNames(seq_along(pair_levels), pair_levels)
  
  boxplot_df = bind_rows(ngd_lst) %>%
    as.data.frame() %>%
    mutate(species = str_replace_all(species, "_", " ")) %>%
    # filter(sample_1_organ != sample_2_organ | sample_id_1 == sample_id_2) %>%
    mutate(pair_type = case_when(
      sample_1_organ == "vaginal" & sample_2_organ == "rectal" & sample_id_1 == sample_id_2 ~ "Vag-Rec\nSame woman",
      sample_1_organ == "vaginal" & sample_2_organ == "rectal" & sample_id_1 != sample_id_2 ~ "Vag-Rec\nDiff woman",
      sample_1_organ == sample_2_organ & sample_1_organ == "vaginal"                        ~ "Vag-Vag",
      sample_1_organ == sample_2_organ & sample_1_organ == "rectal"                         ~ "Rec-Rec",
      TRUE ~ NA_character_)) %>%
    filter(!is.na(pair_type)) %>%
    mutate(pair_type = factor(pair_type, levels = pair_levels)) %>%
    mutate(species = factor(species, levels = species_levels))
    
    kw_df = boxplot_df %>%
                filter(pair_type %in% pair_levels) %>%
                group_by(species) %>%
                summarise(pval = kruskal.test(nGD ~ pair_type)$p.value,
                  .groups = "drop") %>%
                mutate(fdr = p.adjust(pval, method = "fdr")) %>%
                mutate(signif = case_when(fdr <= 0.05 ~ "*", fdr <= 0.01 ~ "**", fdr <= 0.001 ~ "***", TRUE ~ "ns"))
    sig_spec = kw_df %>% filter(fdr <= 0.05) %>% pull(species)
    max_spec <- boxplot_df %>%
                  filter(pair_type %in% pair_levels) %>%
                  group_by(species) %>%
                  summarise(max_nGD = max(nGD, na.rm = TRUE), .groups = "drop")
    dunn_df <- boxplot_df %>%
                filter(pair_type %in% pair_levels) %>%
                mutate(pair_type = factor(pair_type, levels = pair_levels)) %>%
                group_by(species) %>%
                dunn_test(nGD ~ pair_type, p.adjust.method = "fdr") %>%
                mutate(z = statistic,
                  group1 = as.character(group1),
                  group2 = as.character(group2),
                  direction = case_when(is.na(z) ~ NA_character_,
                    z > 0 ~ paste0("<"),
                    z < 0 ~ paste0(">"),
                    TRUE ~ NA_character_)) %>%
                ungroup() %>%
                left_join(max_spec, by = "species") %>%
                group_by(species) %>%
                mutate(y.position = max_nGD * (1.2 + 0.1 * (row_number() - 1))) %>%
                ungroup() %>%
                mutate(p.adj.signif = case_when(p.adj <= 0.05 ~ "*", p.adj <= 0.01 ~ "**", p.adj <= 0.001 ~ "***"),
                  direction.signif = case_when(
                    is.na(p.adj.signif) ~ NA_character_,
                    p.adj.signif == "ns" ~ NA_character_,
                    TRUE ~ paste0(direction, ", ", p.adj.signif)))
  
  pair_colors = c("Vag-Rec\nSame woman" = "#C1121F",
                  "Vag-Rec\nDiff woman" = "#669BBC",
                  "Vag-Vag"             = "#EB5E00",
                  "Rec-Rec"             = "#FFB703")
  
  boxplot_p = ggplot(boxplot_df, aes(x = pair_type, y = nGD)) +
    geom_boxplot(aes(fill = pair_type), alpha = 0.6, outlier.shape = NA, color = "black") +
    geom_point(aes(fill = pair_type), shape = 21, color = "black",
               alpha = 0.2, size = 0.8,
               position = position_jitter(width = 0.2, height = 0)) +
    stat_pvalue_manual(dunn_df, label      = "direction.signif",
                       tip.length = 0.01,
                       hide.ns    = TRUE, size = 15,
                       inherit.aes = FALSE) +
    facet_wrap(~species, scales = "free_y", nrow = 3) +
    scale_fill_manual(values = pair_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_bw(base_family = "Gafata") +
    labs(x = "", y = "nGD") +
    theme(legend.position  = "none",
          strip.background = element_rect(fill = "#FEFCFD"),
          axis.title       = element_text(size = 40),
          axis.title.y     = element_text(size = 50),
          axis.text.y      = element_text(size = 42),
          axis.text.x      = element_text(size = 45),
          strip.text       = element_text(size = 50)) +
    theme(axis.text.x = element_text(lineheight = 0.3))
  
  return(boxplot_p)
}

get_percentile = function(ngd_all, percentile = 0.01, n = 2){
  filt_ngd_all = ngd_all %>%
    group_by(sgb) %>%
    filter(sum(same_person == TRUE) >= n) %>%
    ungroup() 
  percent_df = filt_ngd_all %>%
    filter(same_person == FALSE) %>%
    group_by(species) %>%
    summarise(percentile_nGD = quantile(nGD, probs = percentile, na.rm = TRUE),.groups = "drop")
  full_percent_df = filt_ngd_all %>% left_join(percent_df, by = "species") %>%
    mutate(sig = if_else((same_person == TRUE & nGD <= percentile_nGD), "shared", "not_shared"))
}

kw_comp_plot = function(df_long, results_table, measure_col_name = "Unclassified %", measure_col = "unclassified"){
  sig_organ = results_table %>%
                  filter(p_fdr <= 0.05) %>% distinct(organ) %>% pull(organ)
  max_cst <- df_long %>%
    group_by(CST) %>%
    summarise(max_measure_col = max(!!as.name(measure_col) - 1, na.rm = TRUE), .groups = "drop")
  dunn_df <- df_long %>%
    group_by(organ) %>%
    dunn_test(formula = as.formula(paste(measure_col, "~ CST")), p.adjust.method = "fdr") %>%
    mutate(z = statistic,
           group1 = as.character(group1),
           group2 = as.character(group2),
           direction = case_when(is.na(z) ~ NA_character_,
                                 z > 0 ~ paste0("<"),
                                 z < 0 ~ paste0(">"),
                                 TRUE ~ NA_character_)) %>%
    ungroup() %>%
    add_y_position(data = df_long, formula = reformulate(x_col, measure_col)) %>%
    mutate(p.adj.signif = case_when(organ %in% sig_organ & p.adj <= 0.05 ~ "*", 
                                    organ %in% sig_organ & p.adj <= 0.01 ~ "**", 
                                    organ %in% sig_organ & p.adj <= 0.001 ~ "***"),
           direction.signif = case_when(
             is.na(p.adj.signif) ~ NA_character_,
             p.adj.signif == "ns" ~ NA_character_,
             TRUE ~ paste0(direction, ", ", p.adj.signif)))
  
  boxplot_p = ggplot(df_long, aes_string(x = "CST", y = measure_col)) +
    geom_boxplot(aes(fill = CST), alpha = 0.6, outlier.shape = NA, color = "black") +
    geom_point(aes(fill = CST), shape = 21, color = "black",
               alpha = 0.2, size = 0.8,
               position = position_jitter(width = 0.2, height = 0)) +
    stat_pvalue_manual(dunn_df, label      = "direction.signif",
                       tip.length = 0.01,
                       hide.ns    = TRUE, size = 15,
                       inherit.aes = FALSE) +
    scale_fill_manual(values = cst_colors) +
    facet_wrap(~organ, scales = "free_y", nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_bw(base_family = "Gafata") +
    labs(x = "CST", y = measure_col_name) +
    theme(legend.position  = "none",
          strip.background = element_rect(fill = "#FEFCFD"),
          axis.title       = element_text(size = 40),
          axis.title.y     = element_text(size = 50),
          axis.text.y      = element_text(size = 42),
          axis.text.x      = element_text(size = 45),
          strip.text       = element_text(size = 50)) +
    theme(axis.text.x = element_text(lineheight = 0.3))
}