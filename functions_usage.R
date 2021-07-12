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
    else if(classes_df$parent_class_p[i] == "equal" & classes_df$parent_class_h[i] == "equal"){
      classes_df$classification[i] <- "PEI eq" #covers 1 condition where expression is Parental Equal Expression Inheritance
    }
    else if ((classes_df$parent_class_p[i] == parent_1 & classes_df$parent_class_h[i] == parent_1) | 
             (classes_df$parent_class_p[i] == parent_2 & classes_df$parent_class_h[i] == parent_2)){
      classes_df$classification[i] <- "PEI de" #covers 2 conditions where expression is Parental Differential Expression Inheritance
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

# make a dataframe with only classified genes (i.e. no 'NA' genes)
get_nonNA <- function(classes_df){
  #classes_df : output dataframe from gene_cats()
  classes_df[(!is.na(classes_df$parent_class_h) & !is.na(classes_df$parent_class_p)), ] #only rows that had non-NA result from both DEAs
}

## usage

## allopolyploid fungi

# reading in read count data
hybrid_mod_allo_f <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "allo_f", include_N = TRUE, 
                               parent_1 = "E_elymi", parent_2 = "E_amarillans") #hybrid

exp_file_allo_f <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "allo_f") #parent
parent_mod_allo_f <- parental_DE(exp_file_allo_f, parent_1 = "E_elymi", parent_2 = "E_amarillans")

# identify differentially expressed genes
parent_classes_allo_f <- fit_and_classify(parent_mod_allo_f, parent_1 = "E_elymi", parent_2 = "E_amarillans", cutoff=1, min_p = 0.05)
hybrid_classes_allo_f <- fit_and_classify(hybrid_mod_allo_f, parent_1 = "E_elymi", parent_2 = "E_amarillans", cutoff=1, min_p = 0.05)

classes_df_allo_f <- gene_cats(parent_classes_allo_f, hybrid_classes_allo_f, parent_1 = "E_elymi", parent_2 = "E_amarillans", 
                               results_dir = "~/Desktop/research/frozen_data/hylite/", species = "allo_f")

# exclude genes with very low coverage
sub_classes_df_allo_f <- get_nonNA(classes_df_allo_f)


## homoploid hybrid fungi

# reading in read count data
hybrid_mod_hh_f <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "HH_f", include_N = TRUE, 
                             parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae") #hybrid

exp_file_hh_f <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "HH_f") #parent
parent_mod_hh_f <- parental_DE(exp_file_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae")

# identify differentially expressed genes
parent_classes_hh_f <- fit_and_classify(parent_mod_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", cutoff=1, min_p = 0.05)
hybrid_classes_hh_f <- fit_and_classify(hybrid_mod_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", cutoff=1, min_p = 0.05)

classes_df_hh_f <- gene_cats(parent_classes_hh_f, hybrid_classes_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", 
                             results_dir = "~/Desktop/research/frozen_data/hylite/", species = "HH_f")

# exclude genes with very low coverage
sub_classes_df_hh_f <- get_nonNA(classes_df_hh_f)


## allopolyploid plants

# reading in read count data
hybrid_mod_allo_p <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "allo_p", include_N = TRUE, 
                               parent_1 = "G_raimondii", parent_2 = "G_arboreum") #hybrid

exp_file_allo_p <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "allo_p") #parent
parent_mod_allo_p <- parental_DE(exp_file_allo_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum")

# identify differentially expressed genes
parent_classes_allo_p <- fit_and_classify(parent_mod_allo_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)
hybrid_classes_allo_p <- fit_and_classify(hybrid_mod_allo_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)

classes_df_allo_p <- gene_cats(parent_classes_allo_p, hybrid_classes_allo_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                               results_dir = "~/Desktop/research/frozen_data/hylite/", species = "allo_p")

# exclude genes with very low coverage
sub_classes_df_allo_p <- get_nonNA(classes_df_allo_p)


## homoploid hybrid plants

# reading in read count data
hybrid_mod_hh_p <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "HH_p", include_N = TRUE, 
                             parent_1 = "G_raimondii", parent_2 = "G_arboreum") #hybrid

exp_file_hh_p <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "HH_p") #parent
parent_mod_hh_p <- parental_DE(exp_file_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum")

# identify differentially expressed genes
parent_classes_hh_p <- fit_and_classify(parent_mod_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)
hybrid_classes_hh_p <- fit_and_classify(hybrid_mod_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)

classes_df_hh_p <- gene_cats(parent_classes_hh_p, hybrid_classes_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                             results_dir = "~/Desktop/research/frozen_data/hylite/", species = "HH_p")

# exclude genes with very low coverage
sub_classes_df_hh_p <- get_nonNA(classes_df_hh_p)


## allopolyploid animals

# reading in read count data
hybrid_mod_allo_a <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "allo_a", include_N = TRUE, 
                               parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides") #hybrid

exp_file_allo_a <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "allo_a") #parent
parent_mod_allo_a <- parental_DE(exp_file_allo_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides")

# identify differentially expressed genes
parent_classes_allo_a <- fit_and_classify(parent_mod_allo_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", cutoff=1, min_p = 0.05)
hybrid_classes_allo_a <- fit_and_classify(hybrid_mod_allo_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", cutoff=1, min_p = 0.05)

classes_df_allo_a <- gene_cats(parent_classes_allo_a, hybrid_classes_allo_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", 
                               results_dir = "~/Desktop/research/frozen_data/hylite/", species = "allo_a")

# exclude genes with very low coverage
sub_classes_df_allo_a <- get_nonNA(classes_df_allo_a)


## homoploid hybrid animals

# reading in read count data
hybrid_mod_hh_a <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "HH_a", include_N = TRUE, 
                             parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides") #hybrid

exp_file_hh_a <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "HH_a") #parent
parent_mod_hh_a <- parental_DE(exp_file_hh_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides")

# identify differentially expressed genes
parent_classes_hh_a <- fit_and_classify(parent_mod_hh_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", cutoff=1, min_p = 0.05)
hybrid_classes_hh_a <- fit_and_classify(hybrid_mod_hh_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", cutoff=1, min_p = 0.05)

classes_df_hh_a <- gene_cats(parent_classes_hh_a, hybrid_classes_hh_a, parent_1 = "S_pyrenaicus", parent_2 = "S_alburnoides", 
                             results_dir = "~/Desktop/research/frozen_data/hylite/", species = "HH_a")

# exclude genes with very low coverage
sub_classes_df_hh_a <- get_nonNA(classes_df_hh_a)
