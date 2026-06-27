get_ngd_violin = function(ngd_curr_lst, alpha_false = 0.1, size_false = 1.5, save_file_name = "ngd_all", n = 2, box = FALSE, nrow = 3, yax_title = "nGD", xax_col = "same_person"){
  # Statistical analysis
  if(is.data.frame(ngd_curr_lst)){
    ngd_all = ngd_curr_lst
  } else {
    ngd_all = do.call("rbind", ngd_curr_lst)
  }
  
  ngd_for_plot = ngd_all %>%
    as.data.frame() %>%
    group_by(sgb) %>%
    filter(sum(.data[[xax_col]] == TRUE, na.rm = TRUE) >= n) %>%
    ungroup() %>%
    mutate(species = str_replace_all(species, "_", " ")) %>%
    mutate(!!xax_col := factor(.data[[xax_col]], levels = c(FALSE, TRUE))) %>%
    select(-c(sample_1, sample_2)) %>%
    filter(sample_1_organ != sample_2_organ)
  
  y_position_df = ngd_for_plot %>%
    group_by(species) %>%
    summarise(
      y_position = max(nGD, na.rm = TRUE) - 0.0005, .groups = "drop")
  
  wilcox_results = ngd_for_plot %>%
    group_by(species) %>%
    wilcox_test(reformulate(xax_col, response = "nGD")) %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p, method = "fdr"),
      signif = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01  ~ "**",
        p.adj <= 0.05  ~ "*",
        TRUE           ~ "ns")) %>%
    left_join(y_position_df, by = "species") %>%
    arrange(p.adj) %>%
    left_join(sgb_2_gtdb %>% select(sgb, s) %>% mutate(species = str_replace_all(s, "_", " ")), by = "species") %>%
    mutate(species = factor(species))
  
  ngd_for_plot = ngd_for_plot %>%
    mutate(species = factor(species,levels = wilcox_results$species))
  
  # Plot
  ngd_all_plot = ggplot(data = ngd_for_plot,aes(x = .data[[xax_col]],y = nGD))
  
  if (box) {
    ngd_all_plot = ngd_all_plot +
      geom_boxplot(color = "black",fill = "#D9E1EC",alpha = 0.7,outlier.shape = NA) +
      geom_point(aes(fill = .data[[xax_col]],alpha = .data[[xax_col]], size = .data[[xax_col]]),color = "black",shape = 21,position = position_jitter(width = 0.2))
  } else {
    ngd_all_plot = ngd_all_plot +
      geom_violin(data = ngd_for_plot %>%
          filter(.data[[xax_col]] == FALSE), fill = "#0A9396",alpha = 0.6,color = "black",scale = "count")
  }
  # Save results
  save_wilcox_df = wilcox_results %>%
                    select(species, sgb, n1, n2, p, p.adj, signif) %>%
                    rename(fdr = p.adj)
  
  write.csv(save_wilcox_df,file = paste0(save_dir,"wilcoxon_",save_file_name,"_",as.character(n),".csv") )
  message(paste0("Wilcoxon saved to ","wilcoxon_",save_file_name,"_",as.character(n),".csv"))
  
  ngd_all_plot = ngd_all_plot +
    geom_point(data = ngd_for_plot %>%
        filter(.data[[xax_col]] == TRUE), size = 3,
      alpha = 0.9,fill = "orange",color = "black",shape = 21,position = position_jitter(width = 0.1, height = 0)) +
    theme_bw(base_family = "Gafata") +
    facet_wrap(~species, nrow = nrow,scales = "free_y") +
    stat_pvalue_manual(wilcox_results,label = "signif", y.position = "y_position",size = 15,hide.ns = TRUE) +
    labs(x = str_to_sentence(str_replace_all(xax_col, c("phlame_" =  "", "_" = " "))),y = yax_title) +
    scale_fill_manual(values = c("#0A9396", "orange"), name = "") +
    scale_alpha_manual(values = c(alpha_false, 0.8), name = "") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(drop = FALSE) +
    theme(legend.position = "none",
      strip.background = element_rect(fill = "#FEFCFD"),
      axis.title = element_text(size = 60),
      axis.text = element_text(size = 50),
      strip.text = element_text(size = 40),
      plot.margin = margin(0, 0, -1.5, 0) )
  
  return(list(plot = ngd_all_plot, stats_results = wilcox_results))
}


get_tile_heatmap = function(df, breaks_vec = c(-7.5, -4, -2, 0, 1)){
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
    # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.text.x = element_blank()) +
    theme(legend.position = "bottom") +
    theme(text = element_text(size = 55)) +
    theme(legend.title = element_text(vjust = 1)) +
    theme(plot.margin = margin(-1.5, 0, 0, 0))
  
  return(tile_plot)
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
