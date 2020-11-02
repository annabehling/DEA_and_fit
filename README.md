# DEA_and_fit

write about expression classes, intro

## Description

This code takes the read count matrices outputted by [HyLiTE](https://hylite.sourceforge.io/) analyses, performs differential expression analyses, fits the regression models, and classifies each gene into one of the four gene expression categories or NA.

The majority of functions were taken from [hyliter](https://github.com/dwinter/hyliter), sometimes with modifications made specific for this project. Where applicable, this has been commented in the R code file.

## Installation

To install the required scripts, first clone the **DEA_and_fit** repository.
```
git clone https://github.com/annabehling/DEA_and_fit
```

## The HyLiTE output

HyLiTE outputs a number of files. An outline of each output file and what information they contain can be found [here](https://hylite.sourceforge.io/outformat.html#outformat "HyLiTE output formats").
The important output files for the following analyses are the `expression.txt` and `read.summary.txt` files.

## Usage

This example code uses files from a HyLiTE analysis on a homoploid plant hybrid derived from the parental species *Gossypium raimondii* and *G. arboreum*. These can be found in the `files/` folder.

The assignment of **parent_1** and **parent_2** are arbitrary, but should be kept constant at every point throughout the code.

First load the functions.
```{r}
source("functions.R")
```

The next step is to read in the hybrid and parental read count data, and perform the two differential expression analyses.

Assuming you are in the directory containing the HyLiTE read count files, run:
```{r}
#hybrid
hybrid_DEA_res <- hybrid_DE(".", "HH_p", include_N = TRUE, parent_1 = "G_raimondii", parent_2 = "G_arboreum")

#parental
exp_file_HH_p <- read_exp_file(".", "HH_p")
parent_DEA_res <- parental_DE(exp_file_HH_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum")
```

The, to identify differentially expressed genes, run:
```{r}
parent_classes <- fit_and_classify(parent_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)
hybrid_classes <- fit_and_classify(hybrid_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", cutoff=1, min_p = 0.05)
```

Lastly, to classify each gene into one of the four expression categories, run:
```{r}
classes_df <- gene_cats(parent_classes, hybrid_classes, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                        results_dir = ".", species = "HH_p")
```

To remove the genes that did not get a classification (due to an 'NA' result in one or both of the differential expression analyses) from the above table, run:
```{r}
sub_classes_df <- get_nonNA(classes_df)
```

An example output file can be found at `files/sub_classes_df.txt`.

## Checks

## Acknowledgements

David Winter/hyliter
