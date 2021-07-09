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
hybrid_mod_hh_f <- hybrid_DE("~/Desktop/research/frozen_data/hylite/", "single_end_test", include_N = TRUE, 
                             parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae") #hybrid

exp_file_hh_f <- read_exp_file("~/Desktop/research/frozen_data/hylite/", "single_end_test") #parent
parent_mod_hh_f <- parental_DE(exp_file_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae")

# identify differentially expressed genes
parent_classes_hh_f <- fit_and_classify(parent_mod_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", cutoff=1, min_p = 0.05)
hybrid_classes_hh_f <- fit_and_classify(hybrid_mod_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", cutoff=1, min_p = 0.05)

classes_df_hh_f <- gene_cats(parent_classes_hh_f, hybrid_classes_hh_f, parent_1 = "S_paradoxus", parent_2 = "S_cerevisiae", 
                             results_dir = "~/Desktop/research/frozen_data/hylite/", species = "single_end_test")

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