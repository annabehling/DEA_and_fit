## load functions

fit_and_classify_IF <- function(mod, parent_1, parent_2, cutoff, min_p){ 
  #mod : output of parental_DE() or hybrid_DE()
  #parent_1 : name of one parental species - keep constant across functions
  #parent_2 : name of second parental species - keep constant across functions
  #cutoff : log2 fold change threshold
  #min_p : adjusted p value threshold
  fit <- results(DESeq(mod)) #differential expression analysis with independent filtering
  parent_class <- cut(fit$log2FoldChange, breaks=c(-Inf, -cutoff,cutoff, Inf), labels=c(parent_1, "equal", parent_2)) #cut fit by log2 fold change threshold
  parent_class[fit$padj > min_p] <- "equal" #DE with adjusted p value higher than threshold is not-significant, therefore equal
  data.frame(fit, parent_class)
}

## testing the arbitrariness of the parental definition

# where G. raimondii == parent 1 and G. arboreum == parent 2
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

table(classes_df_hh_p$classification)
#HEBi   HEBl    HER PEI de PEI eq 
#721    920     28    159   7714


#where G. arboreum == parent 1 and G. raimondii == parent 2
# reading in read count data
hybrid_mod_hh_p_2 <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "HH_p", include_N = TRUE, 
                             parent_1 = "G_arboreum", parent_2 = "G_raimondii") #hybrid

exp_file_hh_p <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "HH_p") #parent
parent_mod_hh_p_2 <- parental_DE(exp_file_hh_p, parent_1 = "G_arboreum", parent_2 = "G_raimondii")

# identify differentially expressed genes
parent_classes_hh_p_2 <- fit_and_classify(parent_mod_hh_p_2, parent_1 = "G_arboreum", parent_2 = "G_raimondii", cutoff=1, min_p = 0.05)
hybrid_classes_hh_p_2 <- fit_and_classify(hybrid_mod_hh_p_2, parent_1 = "G_arboreum", parent_2 = "G_raimondii", cutoff=1, min_p = 0.05)

classes_df_hh_p_2 <- gene_cats(parent_classes_hh_p_2, hybrid_classes_hh_p_2, parent_1 = "G_arboreum", parent_2 = "G_raimondii", 
                             results_dir = "~/Desktop/research/frozen_data/hylite/", species = "HH_p")

table(classes_df_hh_p_2$classification)
#HEBi   HEBl    HER PEI de PEI eq 
#721    920     28    159   7714


## testing the extraneous effects of no independent filtering

# identify differentially expressed genes, using independent filtering
parent_classes_hh_p_IF <- fit_and_classify_IF(parent_mod_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)
hybrid_classes_hh_p_IF <- fit_and_classify_IF(hybrid_mod_hh_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)

# check that of the DESEQ2 differential expression analysis results columns, no values except the adjusted p values have changed
all(parent_classes_hh_p[, c(1:5)] == parent_classes_hh_p_IF[, c(1:5)], na.rm = TRUE) #TRUE
all(hybrid_classes_hh_p[, c(1:5)] == hybrid_classes_hh_p_IF[, c(1:5)], na.rm = TRUE) #TRUE
