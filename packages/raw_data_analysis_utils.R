# Read tab-delimited files
read_tsv2 <- function(p, cmnt = "") { 
  return(read_delim(p, comment = cmnt, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)) 
}

