# Separate to organs
sep_one_organ_dissimilarity <- function(df){
  lst = list()
  df = df %>%
    rename(sample_id = 1) 
  for (organ in c('rectal', 'vaginal')){
    tmp <- df %>%
           mutate_at(c('sample_id'), as.factor) %>%
           filter(grepl(organ, sample_id)) %>%
           dplyr::select(matches(paste0(organ, "|sample_id"))) %>%
           mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>% 
           mutate(sample_id = str_replace_all(sample_id, paste0("-", organ), "")) %>% 
           rename_with(~ str_replace_all(., c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>%
           rename_with(~ str_replace_all(., paste0("-", organ), ""))
    
    lst[[organ]] = tmp
  }
  return(lst)
}

order_organ_diversity = function(file_path, organ, metric){
  if (organ == "rectal"){
    suffix = "rec"
  } else {
    suffix = "vag"
  }
  file_path = paste0(feat_table_dir, "filtered_diversity/", metric, "_", suffix, ".tsv")
  if (metric %in% c("bc", "wuf")){
    df = read_tsv2(file.path(file_path)) %>% rename("sample_id" = "...1") 
  } else {
    df <- read_delim(file.path(file_path), delim = "\t", escape_double = FALSE, col_names = c('sample_id', metric), skip = 1, trim_ws = TRUE, show_col_types = FALSE)
  }
  ord_df = df %>%
            mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>% 
            mutate(sample_id = str_replace_all(sample_id, paste0("-", organ), ""))
  if (ord_df %>% ncol() > 2){
    ord_df = ord_df %>%
              rename_with(~ str_replace_all(., c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>%
              rename_with(~ str_replace_all(., paste0("-", organ), ""))
  }
  return(ord_df)          
}

sep_one_organ_tax <- function(df, feat_table_dir, rem = TRUE, remove_rare = TRUE, prevalence_cutoff = 0.1, avg_abundance_cutoff = 0.05, pathway = FALSE){
  lst = list()
  for (organ in c('rectal', 'vaginal')){
    tmp <- df %>%
           mutate(sample_id = as.character(sample_id)) %>%
           filter(grepl(organ, sample_id)) %>%
           mutate(sample_id = str_replace_all(sample_id, paste0("-", organ), ""))
    if(rem){
      tmp = tmp %>%
        dplyr::select(-(sample_id_orig))
    } 
    if ((tmp %>% nrow()) != (length(unique(tmp$sample_id)))){
      second_occurrence <- duplicated(tmp$sample_id)
      tmp$sample_id[second_occurrence] <- paste0(tmp$sample_id[second_occurrence], "_2")
    }
    if (pathway){
      tmp = order_pathway_df(tmp)
    } 
    if ((rowSums(tmp[sapply(tmp, is.numeric)], na.rm = TRUE)[1] > 100) & (grepl("humann", feat_table_dir))) {
      divisor = rowSums(tmp %>% column_to_rownames("sample_id"))[1] / 100
      tmp = tmp %>%
        column_to_rownames("sample_id") %>%
        mutate(across(everything(), ~ . / divisor)) %>%
        rownames_to_column("sample_id")
    }
    if (remove_rare){
      tmp = tmp[rowSums(is.na(tmp)) == 0, ]
      tmp = remove_rare_features(tmp, prevalence_cutoff = prevalence_cutoff, avg_abundance_cutoff = avg_abundance_cutoff)
      tmp = remove_unclassified_features(tmp)
      remove_species_lst = c("Chromatium weissei", "Azoarcus sp012911025", 
                             "JACQX01 sp016208715", "CAG-873 sp004553485")
    }
    lst[[organ]] = tmp
  }
  return(lst)
}

# Alpha diversity
order_alpha_diversity <- function(df){
  lst = list()
  for (organ in c('rectal', 'vaginal')){
    tmp <- df %>%
      rename(sample_id = 1) %>%
      filter(grepl(organ, sample_id)) %>%
      mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>% 
      mutate(sample_id = str_replace_all(sample_id, paste0("-", organ), "")) 
    lst[[organ]] = tmp
  }
  return(lst)
}

# Metadata 
order_meta = function(dir, valencia_dir, valencia_filename = "valencia_res", remove_na_cols = TRUE, remove_samples = c("HM-78", "SM-202", "DS-73", "SV-49")){
  # load
  meta = read_csv(dir, na = c("", "NA", "#DIV/0!", "#VALUE!"), show_col_types = FALSE)
  meta = as.data.frame(meta)
  # Duplicate rows for duplicated samples
  duplicated_samples = c('NH-295', 'CP-25', 'OR-226')
  dup_df = meta %>% dplyr::filter(sample_id %in% duplicated_samples) %>% mutate(sample_id = paste(sample_id, "_2", sep = ""))
  meta = meta %>% bind_rows(dup_df)
  # Add classified counts
  proc <- read_delim(paste0(feat_table_dir, "qc_stats.tsv"),
                     delim = "\t", escape_double = FALSE, col_names = TRUE,
                     show_col_types = FALSE)
  save_cols = c('total_reads_before_filtering', 'total_reads_after_qc', 'total_reads_after_bowtie', 'kraken_input_reads')
  proc <- proc %>%
    mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "CP04600_L00" = "", "\\." = "-", "_B" = ""))) %>%
    mutate(sample_id_dup = sample_id) %>%
    mutate(sample_id = str_replace_all(sample_id, '-vaginal|-rectal', "")) %>%
    mutate(sample_id = sub("_(?=-)", "", sample_id, perl = TRUE)) %>%
    extract(sample_id_dup, into = c("sample_id_dup", "organ"), "(.*)-([^-]+)$") %>%
    dplyr::select(-sample_id_dup) %>%
    mutate(across('organ', str_replace, '_2', '')) %>%
    dplyr::select(c(save_cols, 'sample_id', 'organ', 'total_reads_after_bowtie'))
  proc_wide <- pivot_wider(proc,
                           names_from = organ,
                           names_glue = "{organ}_{.value}",
                           values_from = append(save_cols, "total_reads_after_bowtie")) %>%
               mutate_at(c('sample_id'), ~ make.unique(.))
  meta <- meta %>% left_join(proc_wide, by = 'sample_id')
  ## Add CST
  valencia_res = read_csv(paste0(valencia_dir, valencia_filename, ".csv"), show_col_types = FALSE)
  if (!"sample_id" %in% colnames(valencia_res)){
    valencia_res <- valencia_res %>%
      rename("sample_id" = "sampleID")
  }
  meta <- meta %>% left_join(valencia_res %>%
                               mutate_at(c('sample_id'), ~ make.unique(.)) %>%
                               dplyr::select(c(sample_id, subCST, CST, score)), by = 'sample_id')
  # Convert infertility cause to column
  tmp <- meta %>%
    dplyr::select(c(sample_id, starts_with('infertility_cause'))) %>%
    tidyr::pivot_longer(cols = -sample_id, values_to = 'infer_cause') %>%
    dplyr::filter(!is.na(infer_cause)) %>%
    mutate(infer_cause = tolower(infer_cause)) %>%
    mutate(dummy = TRUE) %>%
    tidyr::pivot_wider(id_cols = sample_id,
                       names_from = infer_cause, 
                       names_prefix = 'infer_cause_', 
                       values_from = dummy) %>%
    replace(is.na(.), FALSE) %>%
    # mutate(across(starts_with("infer_cause"), as.numeric)) %>%
    mutate(across(starts_with("infer_cause_"), as.factor))
  meta <- meta %>% left_join(tmp, by = c('sample_id')) 
  # Change unwanted characters
  meta <- meta %>%
    dplyr::select(-(c(serial_number, OPU_date, percent_divided_IVF))) %>%
    mutate_at(vars(contains(c("percent", "total"))), as.numeric) %>%
    mutate(frozen_day_5_6  = ifelse(is.na(frozen_day_5) & is.na(frozen_day_6), NA,
                             ifelse(is.na(frozen_day_5), frozen_day_6,
                             ifelse(is.na(frozen_day_6), frozen_day_5, frozen_day_5 + frozen_day_6)))) %>%
    mutate(pregnancy_test = ifelse(meta$pregnancy_test == "positive", 1, ifelse(meta$pregnancy_test == "negative", 0, NA))) 
  # Complete missing values where needed
  meta$E2_on_hCG[is.na(meta$E2_on_hCG)] <- (meta$'E2_on_hCG_1'[is.na(meta$E2_on_hCG)]) * 1.5
  meta$Prog_on_hCG[is.na(meta$Prog_on_hCG)] <- (meta$'Prog_on_hCG_1'[is.na(meta$Prog_on_hCG)]) * 1.5
  meta <- meta %>%
    mutate_at(c('patient_age'), ~na_if(., 0)) %>%
    mutate(day_of_transfer = case_when(day_of_transfer %in% c("1", "2", "3", "4") ~ "1_4",
                                       day_of_transfer %in% c("5") ~ "5")) %>%
    mutate(no_transfer = case_when(number_of_transferred == 0 ~ 1,
                                   is.na(number_of_transferred) ~ no_transfer,
                                   TRUE ~ no_transfer)) %>%
    mutate(pregnancy_test = case_when(is.na(pregnancy_test) & is.na(no_transfer) ~ 0,
                              no_transfer == 1 ~ NA_real_,
                              TRUE ~ pregnancy_test)) %>%
    mutate(full_clinical_pregnancy = case_when(is.na(clinical_pregnancy) & is.na(no_transfer) ~ 0,
                              no_transfer == 1 ~ NA_real_,
                              TRUE ~ clinical_pregnancy)) %>%
    mutate(full_outcome = case_when(is.na(outcome) & is.na(no_transfer) ~ "chemical",
                              no_transfer == 1 ~ NA_character_,
                              TRUE ~ as.character(outcome))) %>%
    mutate(cat_outcome = ifelse(is.na(outcome), outcome, ifelse(outcome == "LB", 1, 0))) %>% 
    mutate(full_cat_outcome = ifelse(is.na(full_outcome), full_outcome, ifelse(full_outcome == "LB", 1, 0))) %>%
    mutate(lacto_bin = case_when(CST %in% c("I", "II", "III", "V") ~ "lactobacillus",
                                 CST %in% c("IV-A", "IV-B", "IV-C") ~ "diverse")) %>%
    mutate_at(c('probiotics', 'smoking', 'smoke', 'food_additives', 'medications', 'vegan', 'kosher', 'veg', 'clinical_pregnancy', 'progeterone_po_pv', 'progesterone_IM', 'antibiotics_in_the_past_3_months', 'no_transfer'), as.factor) %>%
    mutate_at(c('sperm_volume', 'final_concentration', 'first_beta', 'second_beta', 'frozen_day_5_6'), as.numeric) %>%
    mutate(sample_id = str_replace_all(sample_id, "_", ""))
  ## Remvoe columns with NA
  if (remove_na_cols){
    important_meta = meta[, c('sample_id', 'pregnancy_test', 'clinical_pregnancy', 'outcome', 'cat_outcome', 'full_clinical_pregnancy', 'full_outcome', 'full_cat_outcome')]
    meta <- meta %>%
      select_if(~ sum(is.na(.)) <= (dim(meta)[1] * 0.15)) %>%
      left_join(important_meta, by = c('sample_id'))
  }
  ## Remove samples
  n_samples = meta %>% nrow()
  rectal_thresh = 300000
  vaginal_thresh = 300000
  meta = meta %>%
          dplyr::filter(is.na(rectal_total_reads_after_bowtie) | rectal_total_reads_after_bowtie >= rectal_thresh) %>%
          dplyr::filter(is.na(vaginal_total_reads_after_bowtie) | vaginal_total_reads_after_bowtie >= vaginal_thresh) %>%
          dplyr::filter(!sample_id %in% remove_samples)
  filt_n_samples = meta %>% nrow()
  message(paste0(filt_n_samples, " samples were filtered from ", n_samples, " due to low reads after bowtie"))
  # filter(sample_id != "SV-49") %>%
  return(meta)
}

# Beta diversity 
order_beta_dissimilarity <- function(df) {
  require(dplyr)
  require(tibble)
  df <- df %>%
    rename(sample_id = '...1') %>%
    # remove duplicated samples
    filter(!grepl('CP04600_L002', sample_id)) %>%
    select(-contains("CP04600_L002")) %>%
    # order samples names
    mutate(sample_id = str_replace_all(sample_id, c("_CP04600_L001" = "", "\\." = "-", "_B" = "-"))) %>%
    rename_with(~ sub("_CP04600_L001", '', .) %>%
                  sub("[.]", '-', .)%>%
                  sub("_B", '', .))
  lst = sep_one_organ_dissimilarity(df)
  return(lst)
}
