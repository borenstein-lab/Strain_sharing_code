get_ngd_violin = function(ngd_curr_lst, alpha_false = 0.1, size_false = 1.5, save_file_name = "ngd_all", n = 2, box = FALSE){
  # Statistical analysis
  ngd_all = do.call("rbind", ngd_curr_lst)
  ngd_for_plot = ngd_all %>%
    as.data.frame() %>%
    group_by(sgb) %>%
    filter(sum(same_person == TRUE) >= n) %>%
    ungroup() %>%
    mutate(species = str_replace_all(species, "_", " ")) %>%
    mutate(same_person = factor(same_person, levels = c(FALSE, TRUE))) %>%
    select(-c(sample_1, sample_2)) %>%
    filter(sample_1_organ != sample_2_organ) 
  y_position_df = ngd_for_plot %>%
    group_by(species) %>%
    summarise(y_position = max(nGD) - 0.0005) %>% #  - 0.005
    ungroup()
  wilcox_results = ngd_for_plot %>%
    group_by(species) %>%
    wilcox_test(nGD ~ same_person) %>%
    ungroup() %>%
    mutate(p.adj = p.adjust(p, method = 'fdr')) %>%
    mutate(signif = case_when(p.adj <= 0.001 ~ "***", 
                              p.adj <= 0.01 ~ "**",
                              p.adj <= 0.05 ~ "*",
                              TRUE ~ "ns")) %>%
    left_join(y_position_df, by = "species") %>% arrange(p.adj) %>%
    left_join(sgb_2_gtdb %>% select(sgb, s) %>% mutate(species = str_replace_all(s, "_", " ")), by = "species") %>%
    arrange("p.adg") %>%
    mutate(species = factor(species))
  
  ngd_for_plot = ngd_for_plot %>%
    mutate(species = factor(species, levels = wilcox_results$species))
  
  # Plot
  ngd_all_plot = ggplot(data = ngd_for_plot, aes(x = same_person, y = nGD)) 
  
  if (box){
    ngd_all_plot = ngd_all_plot +
      geom_boxplot(color = 'black', fill = '#D9E1EC', alpha = 0.7, outlier.shape = NA) +
      geom_point(aes(fill = same_person, alpha = same_person, size = same_person), color = 'black', shape = 21, position = position_jitter(width = 0.2))
  } else {
    ngd_all_plot = ngd_all_plot +
      geom_violin(data = ngd_for_plot %>% filter(same_person == FALSE),
                  fill = "#0A9396", alpha = 0.6, color = 'black', scale = "count")
      # geom_beeswarm(aes(fill = same_person, alpha = same_person, size = same_person), color = 'black', groupOnX = TRUE, pch = 21)
  }
  # Save results nicely
  save_wilcox_df = wilcox_results %>% select(c(species, sgb, n1, n2, p, p.adj, signif)) %>%
    rename("fdr" = "p.adj")
  write.csv(save_wilcox_df, file = paste0(save_dir, "wilcoxon_", save_file_name, ".csv"))
  message(paste0("Wilcoxon saved to ", "wilcoxon_", save_file_name, ".csv"))
  # Plot
  ngd_all_plot = ngd_all_plot +
    geom_point(data = ngd_for_plot %>% filter(same_person == TRUE),
               size = 3, alpha = 0.9, fill = "orange", color = 'black', shape = 21, position = position_jitter(width = 0.2)) +
    theme_bw(base_family = "Gafata") +
    facet_wrap(~species, nrow = 3, scales = "free_y") +
    stat_pvalue_manual(wilcox_results, label = "signif", y.position = "y_position", size = 15, hide.ns = TRUE) +
    labs(x = "Same woman", y = "nGD") +
    scale_fill_manual(values = c("#0A9396", "orange"), name = "") +
    scale_alpha_manual(values = c(alpha_false, 0.8), name = "") +
    # scale_size_manual(values = c(size_false, 2.5), name = "") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(drop = FALSE) +
    theme(legend.position = "none") +
    theme(strip.background = element_rect(fill = "#FEFCFD")) +
    theme(axis.title = element_text(size = 60),
          axis.text = element_text(size = 50),
          strip.text = element_text(size = 40)) +
    theme(plot.margin = margin(0, 0, -1.5, 0))
  return(list("plot" = ngd_all_plot, "stats_results" = wilcox_results))
}


get_tile_heatmap = function(df, breaks_vec = c(-6, -4, -2, 0, 2)){
  breaks <- seq(round(min(df$zscore, na.rm = TRUE), 1), round(max(df$zscore, na.rm = TRUE), 1), length.out = 6)  
  pal = c("#CA6702","#fcc411","#fcf1d0","#C47EC4","#802780")
  
  tile_plot = ggplot(df, aes(x = sample_id_1, y = species, fill = zscore)) +
    geom_tile() +
    geom_text(aes(label = signif), color = "black", size = 12) +
    scale_fill_gradientn(
      colors = pal,
      limits = c(min(breaks_vec), max(breaks_vec)),
      breaks = breaks_vec) +
    guides(fill = guide_colourbar(barwidth = unit(5, "cm"), barheight = unit(0.5, "cm"))) +
    labs(x = "Sample ID", y = "", title = "") +
    scale_y_discrete(position = "right") + 
    theme_bw(base_family = "Gafata") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = "bottom") +
    theme(text = element_text(size = 55)) +
    theme(legend.title = element_text(vjust = 1)) +
    theme(plot.margin = margin(-1.5, 0, 0, 0))
  
  return(tile_plot)
}

get_dist_plot = function(df, vertical_lines_df){
  density_plot = ggplot(df, aes(x = nGD), fill = "grey90") + 
    geom_density(alpha = .5) +
    scale_fill_viridis_d() +
    geom_vline(data = vertical_lines_df, aes(xintercept = nGD), color = "#461220", size = 0.1) +
    facet_wrap(~species, ncol = 1) + 
    labs(y = "", x = "nGD") +
    theme_bw(base_family = "Gafata") +
    theme(text = element_text(size = 25)) +
    theme(legend.position = "none") +
    theme(axis.title.y = element_blank()) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(plot.margin = margin(-1.5, -0.8, 0, 0))
  return(density_plot)
}

get_pcoa_plot = function(df, segment_df, nrow = 5){
  points_segment_df = df %>%
                      left_join(segment_df, by = c("sgb", "sample_id", "species")) %>%
                      filter(!is.na(x)) %>%
                      select(-c(x, y, xend, yend))
  pcoa_plot = ggplot(data = df) +
    geom_point(aes(x = PCoA1, y = PCoA2, shape = organ), size = 2, color = "#4B4C4C", alpha = 0.3, fill = "grey") + 
    scale_shape_manual(values = c(21, 24), labels = labels, name = "") +
    # scale_fill_manual(values = c("vaginal" = "orange", "rectal" = "#057D7F"), labels = labels, name = "") +
    geom_point(data = points_segment_df, aes(x = PCoA1, y = PCoA2, shape = organ), size = 2.5, color = "#620508", fill = "#620508", alpha = 1) + 
    geom_segment(data = segment_df, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 1.1, color = "#620508", lineend='round') +
    theme_bw(base_family = "Gafata") +
    facet_wrap(~species, nrow = nrow, scales = "free") +
    theme(legend.position = "bottom") +
    theme(strip.background = element_rect(fill = "#F8EBC3")) +
    theme(text = element_text(size = 40)) +
    theme(strip.text = element_text(size = 43)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    theme(legend.text = element_text(size = 30))
  return(pcoa_plot)
}

## Complex heatmap functions
order_wide_matrix = function(zscore_df, color_column, order_vactor){
  wide_df = zscore_df %>% 
    mutate(species = str_replace_all(species, "_", " ")) %>%
    filter(sample_1_organ != sample_2_organ) %>%
    filter(same_person == TRUE) %>%
    select(c(sample_id_1, !!as.name(color_column), species)) %>%
    spread(species, !!as.name(color_column)) %>% column_to_rownames("sample_id_1") 
  heat_df = wide_df[order(rowSums(is.na(wide_df)), decreasing = T),]
  heat_df = heat_df[order_vector] #heat_df[order(colSums(heat_df), decreasing = T)]
  
  heat_mat = heat_df %>% t() %>% as.matrix()
  return(heat_mat)
}

get_row_anno = function(ngd_long_all, order_vactor){
  ngd_df = ngd_long_all %>% select(c(species, nGD, sample_id_1)) %>% mutate(species = str_replace_all(species, "_", " "))
  ngd_df$species = factor(ngd_df$species, levels = order_vector)
  ngd_lst = split(ngd_df$nGD, ngd_df$species)
  lt = lapply(ngd_lst, function(x) data.frame(density(x)[c("x", "y")]))
  return(lt)
}

get_heatmap = function(zscore_df, color_col, order_vactor){
  main_color_mat = order_wide_matrix(zscore_df, color_col, order_vactor)
  pval_mat = order_wide_matrix(zscore_df, "pval")
  ngd_lst = get_row_anno(ngd_long_all)
  costum_cell_fun = function(j, i, x, y, width, height, fill){
    if(pval_mat[i, j] < 0.0001) {
      grid.text("***", x, y)
    } else if(pval_mat[i, j] < 0.001) {
      grid.text("**", x, y)
    } else if(pval_mat[i, j] < 0.01) {
      grid.text("*", x, y)
    }
  }
  
  ## Plot
  col_fun = colorRamp2(seq(min(main_color_mat, na.rm = TRUE), max(main_color_mat, na.rm = TRUE), length = 5), (c("#CA6702","#fcc411","#fcf1d0","#C47EC4","#802780")))
  
  ha = rowAnnotation(nGD = anno_joyplot(ngd_lst, width = unit(5, "cm"), 
                                        gp = gpar(fill = "#778DA9"), 
                                        transparency = 0.6),
                     annotation_name_gp = gpar(fontsize = 18))
  h = Heatmap(main_color_mat, name = color_col, col = col_fun, na_col = "#EBEDF0",
              cell_fun = function(j, i, x, y, width, height, fill){
                if (!is.na(pval_mat[i, j])){
                  if(pval_mat[i, j] < 0.0001) {
                    grid.text("***", x, y)
                  } else if(pval_mat[i, j] < 0.001) {
                    grid.text("**", x, y)
                  } else if(pval_mat[i, j] < 0.01) {
                    grid.text("*", x, y)
                  }
                } 
              },
              cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
              row_names_gp = gpar(fontsize = 18), left_annotation = ha, 
              row_names_max_width = unit(30, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = 18), 
                                          labels_gp = gpar(fontsize = 14)))
  
  save_filname = paste0(save_dir, "sgbs_to_individual_heatmap_diff_organ_mean", color_col, "_times.png")
  png(save_filname, 
      width = 1300, height = 600)
  draw(h, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 15), "mm"))
  dev.off()
  return(h)
}