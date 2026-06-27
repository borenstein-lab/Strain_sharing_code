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

