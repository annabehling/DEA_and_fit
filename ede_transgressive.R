## load functions

# to get all EDEs
get_EDEs <- function(sub_classes_df){
  #sub_classes_df : table of each non-NA gene with column of expression category for each gene (df)
  p_sig_genes <- sub_classes_df[(sub_classes_df$padj_p < 0.05), ]
  p_sig <- p_sig_genes[(p_sig_genes$log2FC_p > 5.64 | p_sig_genes$log2FC_p < -5.64), ]
  
  h_sig_genes <- sub_classes_df[(sub_classes_df$padj_h < 0.05), ]
  h_sig <- h_sig_genes[(h_sig_genes$log2FC_h > 5.64 | h_sig_genes$log2FC_h < -5.64), ]
  
  rbind(p_sig, h_sig)
}

# to get only parental EDEs
get_par_EDEs <- function(sub_classes_df){
  #sub_classes_df : table of each non-NA gene with column of expression category for each gene (df)
  sig_genes <- sub_classes_df[(sub_classes_df$padj_p < 0.05), ]
  sig_genes[(sig_genes$log2FC_p > 5.64 | sig_genes$log2FC_p < -5.64), ]
}

# to get only hybrid EDEs
get_hyb_EDEs <- function(sub_classes_df){
  #sub_classes_df : table of each non-NA gene with column of expression category for each gene (df)
  sig_genes <- sub_classes_df[(sub_classes_df$padj_h < 0.05), ]
  sig_genes[(sig_genes$log2FC_h > 5.64 | sig_genes$log2FC_h < -5.64), ]
}

# to get transgressive genes
get_transgressives <- function(exp_file, sub_classes_df){
  #exp_file : HyLiTE expression.txt file read in (df)
  #sub_classes_df : table of each non-NA gene with column of expression category for each gene (df)
  
  #calculate mean read counts for p1 in col 2:3, p2 in col 6:7, hybrid in col 4:5
  exp_file$mean_p1 <- rowMeans(exp_file[, c(2:3)], na.rm = TRUE)
  exp_file$mean_p2 <- rowMeans(exp_file[, c(6:7)], na.rm = TRUE)
  exp_file$mean_hy <- rowMeans(exp_file[, c(4:5)], na.rm = TRUE)
  
  #determine genes with transgressive expression
  #transgressive = mean_hybrid > 2x highest mean_parent, or < 0.5x lowest mean_parent
  exp_file$transgressive <- ifelse((exp_file$mean_hy > 2*(exp_file$mean_p1) & 
                                      exp_file$mean_hy > 2*(exp_file$mean_p2)) |
                                     (exp_file$mean_hy < 0.5*(exp_file$mean_p1) &
                                        exp_file$mean_hy < 0.5*(exp_file$mean_p2)), "yes", "no")
  
  #make a table of only transgressively expressed genes
  trans_genes <- exp_file[exp_file$transgressive == "yes", ]
  
  #merge with genes that have a classification from the two differential expression analyses
  trans_cats <- merge(sub_classes_df, trans_genes, by.x = "gene_id", by.y = "GENE")
  trans_cats
}
