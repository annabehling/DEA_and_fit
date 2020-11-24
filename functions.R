## load libraries

library(DESeq2)

## load functions

# hybrid read count data processing 

# read in hybrid expression data (hyliter)
read_hylite <- function(file_name){
  #file_name : path to HyLiTE read.summary.txt file
  read.table(file_name,sep="\t", stringsAsFactors=FALSE, header=1) #read in the read.summary.txt file
}

# process one of the hybrid expression files (modified from hyliter)
process_one_h_file <- function(file_name, include_N){
  #file_name : path to HyLiTE read.summary.txt file
  #include_N : if true, include the +N columns in the read count totals (logical)
  expr_data <- read_hylite(file_name) #read in the read.summary.txt file
  if(include_N){ #if we want the +N read counts included in the parental totals,
    with_N <- cbind(as.integer(rowSums(expr_data[,2:3])), as.integer(rowSums(expr_data[,4:5]))) #sum up the parent and parent+N columns for each
    colnames(with_N) <- c(names(expr_data)[2], names(expr_data)[4]) #set the column names
    return(with_N) #matrix with parent (including +N) readcount data
  }
  as.matrix(expr_data[,c(2,4)]) #matrix with parent (no +N) readcount data
}

# process hybrid read count data (modified from hyliter)
hybrid_DE <- function(results_dir, species, include_N, parent_1, parent_2){
  #results_dir : path to directory containing HyLiTE read.summary.txt files
  #species : specific file identifier in case the directory contains output files from multiple HyLiTE runs
  #include_N : if true, include the +N columns in the read count totals (logical)
  #parent_1 : name of one parental species - keep constant across functions
  #parent_2 : name of second parental species - keep constant across functions
  fnames <- list.files(results_dir, 
                       paste0(species, ".*.", "read.summary.txt"), #looking in directory for read.summary.txt files
                       full.names=TRUE)
  read_counts <- do.call(cbind, lapply(fnames, process_one_h_file, include_N)) #cbind the process_one_h_file outputs for each replicate
  col_data <- data.frame(parent =colnames(read_counts)) #makes data frame with read_counts column names as rows
  colnames(read_counts) <- paste(colnames(read_counts), rep(1:(ncol(read_counts)/2), each=2), sep="_") #pastes _1 _2 after the parent names to distinguish replicates
  mod <- DESeqDataSetFromMatrix( read_counts, 
                          colData = col_data, 
                          design = ~ parent) #make DESeq dataset
  mod$parent <- factor(mod$parent, levels = c(parent_1, parent_2)) #defining the level order of parent_1 and parent_2 so that these can be arbitrarily defined
  mod
}

# read in parental expression data
read_exp_file <- function(results_dir, species){
  #results_dir : path to directory containing HyLiTE expression.txt file
  #species : specific file identifier in case the directory contains output files from multiple HyLiTE runs
  expression_file <- list.files(results_dir, 
                                paste0(species, ".*.", "expression.txt"), #looking in directory for expression.txt file
                                full.names=TRUE)
  read.table(expression_file, sep="\t", stringsAsFactors=FALSE, header=1) #read in the expression.txt file
}

# process parental read count data (modified from hyliter)
parental_DE <- function(expr_data, parent_1 ,parent_2){
  #expr_data : dataframe output of read_exp_file()
  #parent_1 : name of one parental species - keep constant across functions
  #parent_2 : name of second parental species - keep constant across functions
  spp <- c("GENE", sub("\\..+", "", colnames(expr_data)[-1])) #get the species names from the column headers without the HyLiTE formatting
  parents <- c(parent_1, parent_2)
  read_counts <- as.matrix(expr_data[, spp %in% parents]) #get the read count values in the parental columns 
  col_data <- data.frame(species=spp[ spp %in% parents ]) #get the corresponding names of the parental columns
  mod <- DESeqDataSetFromMatrix( read_counts, 
                          colData = col_data, 
                          design = ~ species) #make DESeq dataset
  mod$species <- factor(mod$species, levels = c(parent_1, parent_2)) #defining the level order of parent_1 and parent_2 so that these can be arbitrarily defined
  mod
}

# function to fit parental/hybrid DESeq2 results (modified from hyliter), classify where there is a bias in parental origin
fit_and_classify <- function(mod, parent_1, parent_2, cutoff, min_p){ 
  #mod : output of parental_DE() or hybrid_DE()
  #parent_1 : name of one parental species - keep constant across functions
  #parent_2 : name of second parental species - keep constant across functions
  #cutoff : log2 fold change threshold
  #min_p : adjusted p value threshold
  fit <- results(DESeq(mod), independentFiltering = FALSE) #differential expression analysis without independent filtering
  fit$padj <- p.adjust(fit$pvalue, method="BH") #calculate adjusted p value with the Benjamini-Hochberg procedure
  parent_class <- cut(fit$log2FoldChange, breaks=c(-Inf, -cutoff,cutoff, Inf), labels=c(parent_1, "equal", parent_2)) #cut fit by log2 fold change threshold
  parent_class[fit$padj > min_p] <- "equal" #DE with adjusted p value higher than threshold is not-significant, therefore equal
  data.frame(fit, parent_class)
}

# classify each gene to one of four gene expression categories
gene_cats <- function(parent_fit, hybrid_fit, parent_1, parent_2, results_dir, species){
  #parent_fit : output of fit_and_classify() that was run on parental_DE()
  #hybrid_fit : output of fit_and_classify() that was run on hybrid_DE()
  #parent_1 : name of one parental species - keep constant across functions
  #parent_2 : name of second parental species - keep constant across functions
  #results_dir : path to directory containing HyLiTE expression.txt file
  #species : specific file identifier in case the directory contains output files from multiple HyLiTE runs
  classes_df <- cbind(as.data.frame(parent_fit$parent_class), as.data.frame(parent_fit$log2FoldChange), as.data.frame(parent_fit$padj),
                      as.data.frame(hybrid_fit$parent_class), as.data.frame(hybrid_fit$log2FoldChange), as.data.frame(hybrid_fit$padj)) #make a df with only relevant columns
  names(classes_df)[1] <- "parent_class_p"
  names(classes_df)[2] <- "log2FC_p"
  names(classes_df)[3] <- "padj_p"
  names(classes_df)[4] <- "parent_class_h"
  names(classes_df)[5] <- "log2FC_h"
  names(classes_df)[6] <- "padj_h"
  
  exp_file <- read_exp_file(results_dir, species) #read in expression file to get gene ids (these are in same order as the classes_df table) 
  classes_df$gene_id <- exp_file$GENE
  
  classes_df$classification <- ""
  for(i in 1:nrow(classes_df)){
    message(i)
    if( any(is.na(classes_df[i,c(1, 4)]))){
      classes_df$classification[i] <- NA #not possible to compare NA (no data?) to classification, so pass on NA to final result
    } 
    else if(classes_df$parent_class_p[i] == classes_df$parent_class_h[i]){
      classes_df$classification[i] <- "PEI" #covers 3 conditions where expression is Parental Expression Inheritance
    }
    else if((classes_df$parent_class_p[i] == parent_1 | classes_df$parent_class_p[i] == parent_2) & classes_df$parent_class_h[i] == "equal"){
      classes_df$classification[i] <- "HEBl" #covers 2 conditions where expression is Homeolog Expression Blending
    }
    else if(classes_df$parent_class_p[i] == "equal" & (classes_df$parent_class_h[i] == parent_1 | classes_df$parent_class_h[i] == parent_2)){
      classes_df$classification[i] <- "HEBi" #covers 2 conditions where expression is Homeolog Expression Bias
    }
    else if((classes_df$parent_class_p[i] == parent_1 & classes_df$parent_class_h[i] == parent_2) | (classes_df$parent_class_p[i] == parent_2 & classes_df$parent_class_h[i] == parent_1)){
      classes_df$classification[i] <- "HER" #covers 2 conditions where expression is Homeolog Expression Reversal
    }
  }
  
  classes_df
}

#make a dataframe with only classified genes (i.e. no 'NA' genes)
get_nonNA <- function(classes_df){
  #classes_df : output dataframe from gene_cats()
  classes_df[(!is.na(classes_df$parent_class_h) & !is.na(classes_df$parent_class_p)), ] #only rows that had non-NA result from both DEAs
}