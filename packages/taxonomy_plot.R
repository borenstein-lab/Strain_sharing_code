create_long_df = function(df, samp_order, n = 15){
  long_df <- df %>% 
    filter(sample_id %in% samp_order) %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = 'taxon', values_to = 'relab')
  
  # Group rare otu's together for plot simplicity
  all_otus = length(unique(long_df$taxon))
  rare_otus <- long_df %>%
    group_by(taxon) %>%
    summarise(mean_relab = mean(relab), .groups = 'drop') %>%
    arrange(mean_relab) %>%
    slice_head(n = (all_otus - n)) %>%
    pull(taxon)
    # filter(mean_relab < rare_otu_cutoff) %>%

  long_df <- long_df %>%
    mutate(otu_grouped = ifelse(taxon %in% rare_otus, 'Others', taxon)) %>%
    group_by(sample_id, otu_grouped) %>%
    summarise(relab = sum(relab), .groups = 'drop')
  
  # Sort samples by the abundance of upmost taxon
  long_df$sample_id <- factor(long_df$sample_id, levels = samp_order)
  return(long_df)
}

order_df_before_plot = function(df, samp_order, title_organ, manual_colors_named, n = 15){
  long_df = create_long_df(df, samp_order, n)
  long_df = long_df %>%
              mutate(organ = title_organ) %>%
              mutate(otu_grouped = str_replace_all(otu_grouped, "_", " "))
  
  # Slightly reorder to make the 'Others' category last
  taxa_list <- sort(unique(long_df$otu_grouped))
  taxa_list <- c(taxa_list[taxa_list != 'Others'], 'Others')
  
  # Plot
  other_taxa <- setdiff(taxa_list, names(manual_colors_named))
  other_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, name = 'Spectral'))(length(other_taxa))
  if (any(other_colors %in% manual_colors_named)){
    message("***Problem in colors***")
  } 
  other_colors_named = setNames(other_colors, other_taxa)
  all_colors = c(manual_colors_named, other_colors_named)
  all_colors = all_colors[c(setdiff(names(all_colors), "Others"), "Others")]
  
  long_df = long_df %>%
    mutate(order_group = case_when(otu_grouped %in% names(manual_colors_named) & otu_grouped != "Others" ~ 1,                           
                                   organ == "Vaginal" & !(otu_grouped %in% names(manual_colors_named)) ~ 2,   
                                   organ == "Rectal"  & !(otu_grouped %in% names(manual_colors_named)) ~ 3,  
                                   otu_grouped == "Others" ~ 4,                                
                                   TRUE ~ 5))
  ordered_levels <- long_df %>% arrange(order_group) %>% distinct(otu_grouped) %>% pull(otu_grouped)
  long_df <- long_df %>%
    mutate(taxa_grouped = factor(otu_grouped, levels = ordered_levels)) %>%
    select(-order_group)
  return(list("long_df" = long_df, "all_colors" = all_colors))
}

barplot_taxonomy = function(long_both, tmp_palette, legend_name, facet = TRUE){
  barplot = ggplot(long_both, aes(fill = taxa_grouped, y = relab, x = sample_id)) +
    geom_bar(position = "fill", stat = "identity", color = 'black', linewidth = 0.01) +
    scale_fill_manual(values = tmp_palette, name = legend_name) +
    theme_bw(base_family = "Gafata") +
    scale_y_continuous(expand = c(0,0)) +
    ylab(paste0(unique(long_both$organ), '\nrelative abundance')) +
    xlab("Sample") +
    guides(fill = guide_legend(ncol = 2, byrow = FALSE,
                               title.position="top",
                               keywidth = 0.2,
                               keyheight = 0.1,
                               default.unit = "cm")) 
  if (facet){
    barplot = barplot + facet_wrap(~organ, nrow = 2) 
  }
  barplot = barplot +
    # guides(fill = guide_legend(override.aes = list(size = 0.2))) +
    theme(legend.key.height = unit(6, "mm"),
      legend.key.width  = unit(6, "mm"),
      legend.spacing.y = unit(4, "mm"),   
      legend.spacing.x = unit(4, "mm"),
      legend.position = "bottom") +
    theme(axis.text.x = element_blank()) +
    theme(strip.background = element_rect(fill = "#FFFCEB", linewidth = 0.1)) +
    theme(strip.text = element_text(size = 40)) +
    theme(text = element_text(size = 35)) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(lineheight = 0.3)) +
    theme(plot.margin = margin(l = 0.1, b = 0.1, t = 0.1, r = 0.1, unit = "cm"))
  return(barplot)
}

plot_taxonomy <- function(vag_df, rec_df, samp_order, n, manual_colors_named) {
  long_vag = create_long_df(vag_df, samp_order, n) %>% mutate(organ = "Vaginal")
  long_rec = create_long_df(rec_df, samp_order, n) %>% mutate(organ = "Rectal")
  long_both = bind_rows(long_vag, long_rec) %>%
    mutate(organ = factor(organ, levels = c("Vaginal", "Rectal"))) %>%
    mutate(otu_grouped = str_replace_all(otu_grouped, "_", " "))
  
  # Slightly reorder to make the 'Others' category last
  taxa_list <- sort(unique(long_both$otu_grouped))
  taxa_list <- c(taxa_list[taxa_list != 'Others'], 'Others')

  # Plot
  other_taxa <- setdiff(taxa_list, names(manual_colors_named))
  other_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, name = 'Spectral'))(length(other_taxa))
  if (any(other_colors %in% manual_colors_named)){
    message("***Problem in colors***")
  } 
  other_colors_named = setNames(other_colors, other_taxa)
  all_colors = c(manual_colors_named, other_colors_named)
  all_colors = all_colors[c(setdiff(names(all_colors), "Others"), "Others")]
  
  long_both = long_both %>%
    mutate(order_group = case_when(otu_grouped %in% names(manual_colors_named) & otu_grouped != "Others" ~ 1,                           
                          organ == "Vaginal" & !(otu_grouped %in% names(manual_colors_named)) ~ 2,   
                          organ == "Rectal"  & !(otu_grouped %in% names(manual_colors_named)) ~ 3,  
                          otu_grouped == "Others" ~ 4,                                
                          TRUE ~ 5))
  ordered_levels <- long_both %>% arrange(order_group) %>% distinct(otu_grouped) %>% pull(otu_grouped)
  long_both <- long_both %>%
    mutate(taxa_grouped = factor(otu_grouped, levels = ordered_levels)) %>%
    select(-order_group)
  
  barplot = barplot_taxonomy(long_both, all_colors, "Taxa", facet = TRUE)
  return(barplot)
}
