# Metaphlan
# Convert metaphlan table to convenient feature table
# level should be one of: 'p','c','o','f','g','s'
reorganize_metaphlan_table <- function(df, level = 's') {
  # Regex for required taxonomic level
  rgx_in <- paste0('\\|',level,'__')
  
  # Regex to exclude (lower taxonomy levels)
  levels <- c('k'=1,'p'=2,'c'=3,'o'=4,'f'=5,'g'=6,'s'=7,'t'=8)
  level_out <- names(levels[1+levels[level]])
  rgx_out <- paste0('\\|',level_out,'__')
  
  df %>%
    dplyr::filter(clade_name == 'UNCLASSIFIED' | grepl(rgx_in, clade_name)) %>%
    dplyr::filter(! grepl(rgx_out, clade_name)) %>%
    mutate(clade_name = gsub(paste0('^.*',rgx_in), '', clade_name)) %>%
    column_to_rownames('clade_name') %>%
    t() %>%
    data.frame() %>%
    rownames_to_column('sample_id') %>%
    mutate(sample_id = gsub('_metaphlan4_bugs_list', '', sample_id))
}

remove_unclassified = function(df){
  tmp = df %>% select(-UNCLASSIFIED) %>% column_to_rownames("sample_id")
  tmp <- (tmp / rowSums(tmp)) * 100
  tmp = tmp %>% 
    rownames_to_column("sample_id") %>%
    dplyr::filter(rowSums(!is.na(select(., -sample_id))) > 0) 
  return(tmp)
}

sep_one_organ_metaphlan <- function(df, rem = TRUE, remove_rare = TRUE, prevalence_cutoff = 0.1, avg_abundance_cutoff = 0.05, remove_unclassified_column = FALSE, unclassified_threshold = 70){
  lst = list()
  for (organ in c("vaginal", "rectal")){
    if(!rem){
      df = df %>%
        mutate(sample_id_orig = sample_id) 
    } 
    tmp <- df %>%
      mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "CP04600_L00" = "", "_B" = "", "_1" = ""))) %>% # , "_" = ""))) %>%
      mutate(sample_id = sub("_(?=-)", "", sample_id, perl = TRUE)) %>%
      dplyr::filter(grepl(organ, sample_id)) %>% 
      dplyr::filter(rowSums(!is.na(select(., -sample_id))) > 0) %>% 
      mutate(sample_id = str_replace_all(sample_id, paste0("-", organ), ""))
      # mutate(sample_id = str_replace_all(sample_id, "_", ""))
    organ_nrow = tmp %>% nrow()
    tmp = tmp %>% 
            mutate(UNCLASSIFIED = as.numeric(UNCLASSIFIED)) %>%
            dplyr::filter(UNCLASSIFIED <= unclassified_threshold)
            # dplyr::filter(trimws(UNCLASSIFIED) <= unclassified_threshold)
    message(paste0(nrow(tmp), " of ", organ_nrow, " samples were filtered due to high UNCLASSIFIED values."))
    if ((tmp %>% nrow()) != (length(unique(tmp$sample_id)))){
      second_occurrence <- duplicated(tmp$sample_id)
      tmp$sample_id[second_occurrence] <- paste0(tmp$sample_id[second_occurrence], "_2")
    }
    if (remove_unclassified_column){
      tmp = remove_unclassified(tmp)
    }
    if (remove_rare){
      tmp = remove_rare_features(tmp, prevalence_cutoff = prevalence_cutoff, avg_abundance_cutoff = avg_abundance_cutoff)
      # tmp = remove_unclassified_features(tmp)
      tmp = tmp %>% column_to_rownames("sample_id")
      tmp <- (tmp / rowSums(tmp)) * 100
      tmp = tmp %>% rownames_to_column("sample_id") %>% filter(rowSums(!is.na(select(., -sample_id))) > 0)
      if (organ == "rectal"){
        tmp = tmp %>% rowwise() %>%
          filter(sum(c_across(-sample_id) != 0) > 3) %>%
          ungroup()
      }
    }
    lst[[organ]] = tmp
  }
  return(lst)
}