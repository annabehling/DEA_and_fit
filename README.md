# DEA_and_fit

Given the extent of variability in genetic complements across Eukarya, the cross-kingdom comparison of gene expression patterns cannot be as simple as comparing read count matrices from representative systems. It is possible to instead combine the results of the parental [parent gene vs parent gene] and hybrid [parentally-derived gene copy vs parentally-derived gene copy] differential expression analysis results to form four expression categories that can be compared across kingdoms.

Briefly, these expression categories are:

1. Parental expression inheritance (PEI): any parental bias (or lack thereof) is maintained in the hybrid.
2. Homeolog expression blending (HEBl): a parental expression bias is lost in the hybrid.
3. Homeolog expression bias (HEBi): a hybrid expression bias has arisen from no parental bias.
4. Homeolog expression reversal (HER): an expression bias in the parents is reversed in the hybrid.

More information on these classes can be found in the following publication:

Cox, M.P., T. Dong, G. Shen, Y. Dalvi, D.B. Scott and A.R.D. Ganley. 2014. An interspecific fungal hybrid reveals cross-kingdom rules for allopolyploid gene expression patterns. *PLoS Genetics* 10: e1004180. [https://doi.org/10.1371/journal.pgen.1004180](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004180)

## Description

This code takes the read count matrices outputted by [HyLiTE](https://hylite.sourceforge.io/) analyses, performs differential expression analyses (DEA), fits the regression models, and classifies each gene into one of the four gene expression categories or NA.

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

The next step is to read in the hybrid and parental read count data, and perform the two differential expression analyses. To do so, run:
```{r}
#hybrid
hybrid_DEA_res <- hybrid_DE(results_dir = "./files", species = "HH_p", include_N = TRUE, 
                            parent_1 = "G_raimondii", parent_2 = "G_arboreum")

#parental
exp_file_HH_p <- read_exp_file(results_dir = "./files", species = "HH_p")
parent_DEA_res <- parental_DE(exp_file_HH_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum")
```

The, to identify differentially expressed genes, run:
```{r}
#hybrid
hybrid_classes <- fit_and_classify(hybrid_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                                   cutoff=1, min_p = 0.05)

#parental
parent_classes <- fit_and_classify(parent_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                                   cutoff=1, min_p = 0.05)
```

Lastly, to classify each gene into one of the four expression categories, run:
```{r}
classes_df <- gene_cats(parent_classes, hybrid_classes, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                        results_dir = "./files", species = "HH_p")
```

To remove the genes that did not get a classification (due to an 'NA' result in one or both of the differential expression analyses) in the above table, run:
```{r}
sub_classes_df <- get_nonNA(classes_df)
```

## Output structure

When the above code has been run to completion, it will output a dataframe called `sub_classes_df`.

An example of what this dataframe should look like can be found in this repository under `files/sub_classes_df.txt`.

The dataframe consists of six columns, with the headings:

1. parent_class_p (the gene class assigned from the parental DEA)
2. log2FC_p (the log2(fold change) in gene expression from the parental DEA)
3. parent_class_h (the gene class assigned from the hybrid DEA)
4. log2FC_h (the log2(fold change) in gene expression from the hybrid DEA)
5. gene_id (the genes that had a non-NA result from both DEAs)
6. classification (the gene expression category assigned to each gene)

## Acknowledgements

Thank you to [David Winter](https://github.com/dwinter "github.com/dwinter") for the [hyliter](https://github.com/dwinter/hyliter "github.com/dwinter/hyliter") code and assistance with its implementation.