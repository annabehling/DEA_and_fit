# DEA_and_fit

Given the extent of variability in the genetic complements of species across Eukarya, the cross-kingdom comparison of gene expression patterns cannot be as simple as comparing read count matrices from representative systems. It is possible to instead combine the results of the parental [parent gene vs parent gene] and hybrid [parentally-derived gene copy vs parentally-derived gene copy] differential expression analyses results to form five expression categories that can be compared across kingdoms.

Briefly, these expression categories are:

1. **Parental differential expression inheritance** (PEI de): a parental expression bias is maintained in the hybrid.
2. **Parental equal expression inheritance** (PEI eq) : equal parental expression is inherited by the hybrid.
3. **Homeolog expression blending** (HEBl): a parental expression bias is lost in the hybrid.
4. **Homeolog expression bias** (HEBi): a hybrid expression bias has arisen from no parental bias.
5. **Homeolog expression reversal** (HER): an expression bias in the parents is reversed in the hybrid.

More information on these classes can be found in the following publications:

Cox, M.P., T. Dong, G. Shen, Y. Dalvi, D.B. Scott and A.R.D. Ganley. 2014. An interspecific fungal hybrid reveals cross-kingdom rules for allopolyploid gene expression patterns. *PLoS Genetics* 10: e1004180. [https://doi.org/10.1371/journal.pgen.1004180](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004180)

Yoo, M.J., Szadkowski, E., & Wendel, J.F. 2013. Homeolog expression bias and expression level dominance in allopolyploid cotton. *Heredity* 110: 171-180. [https://doi.org/10.1038/hdy.2012.94](https://www.nature.com/articles/hdy201294)

## Description

This code takes the read count matrices outputted by [HyLiTE](https://hylite.sourceforge.io/) analyses, performs differential expression analyses (DEA), fits the regression models, and classifies each gene into one of the five gene expression categories or NA. The code has been tested to work on R version 4.0.3.

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

This example code uses files from a HyLiTE analysis on a homoploid plant hybrid derived from the parental species *Gossypium raimondii* and *G. arboreum*. These can be found in the `example_files/` folder.

The assignment of **parent_1** and **parent_2** are arbitrary, but should be kept constant at every point throughout the code.

First load the functions:
```{r}
source("functions.R")
```

The next step is to read in the hybrid and parental read count data, and perform the two differential expression analyses. To do so, run:
```{r}
#hybrid
hybrid_DEA_res <- hybrid_DE(results_dir = "./example_files", species = "HH_p", include_N = TRUE, 
                            parent_1 = "G_raimondii", parent_2 = "G_arboreum")

#parental
exp_file_HH_p <- read_exp_file(results_dir = "./example_files", species = "HH_p")
parent_DEA_res <- parental_DE(exp_file_HH_p, parent_1 = "G_raimondii", parent_2 = "G_arboreum")
```

Then, to identify differentially expressed genes, run:
```{r}
#hybrid
hybrid_classes <- fit_and_classify(hybrid_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                                   cutoff=1, min_p = 0.05)

#parental
parent_classes <- fit_and_classify(parent_DEA_res, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                                   cutoff=1, min_p = 0.05)
```

Lastly, to classify each gene into one of the five expression categories, run:
```{r}
classes_df <- gene_cats(parent_classes, hybrid_classes, parent_1 = "G_raimondii", parent_2 = "G_arboreum", 
                        results_dir = "./example_files", species = "HH_p")
```

To remove the genes that did not get a classification (due to an 'NA' result in one or both of the differential expression analyses) in the above table, run:
```{r}
sub_classes_df <- get_nonNA(classes_df)
```

## Output structure

When the above code has been run to completion, it will output a dataframe called `sub_classes_df`.

An example of what this dataframe should look like can be found in this repository under `example_files/sub_classes_df.txt`.

The dataframe consists of eight columns, with the headings:

1. **parent_class_p**: the gene class assigned from the parental DEA
2. **log2FC_p**: the log2 fold change in gene expression from the parental DEA
3. **padj_p** : the adjusted *p* value from the parental DEA
4. **parent_class_h**: the gene class assigned from the hybrid DEA
5. **log2FC_h**: the log2 fold change in gene expression from the hybrid DEA
6. **padj_h** : the adjusted *p* value from the hybrid DEA
7. **gene_id**: the genes that had a non-NA result from both DEAs
8. **classification**: the gene expression category assigned to each gene

## Next steps

Further analyses can be perfomed on the output dataframe, `sub_classes_df`. These include identifying genes with extreme differential expression or transgressive expression.

First load the following additional functions:
```{r}
source("ede_transgressive.R")
```

### Identifying genes with extreme differential expression

A gene with a fold change > 50 in either the parental or hybrid differential expression analysis is considered extremely differentially expressed.

To identify all extremely differentially expressed genes, run:
```{r}
get_EDEs(sub_classes_df)
```

Alternatively, to identify only the extremely differentially expressed parental genes, run:
```{r}
get_par_EDEs(sub_classes_df)
```

Or, to identify only the extremely differentially expressed hybrid genes, run:
```{r}
get_hyb_EDEs(sub_classes_df)
```

### Identifying genes with transgressive expression

Transgressive expression is defined as mean hybrid expression greater than two times the highest (parent 1 or parent 2) mean parental expression, or less than half of the lowest (parent 1 or parent 2) mean parental expression.

To identify all genes with transgressive expression, run:
```{r}
get_transgressives(exp_file_HH_p, sub_classes_df)
```

## Additional data availability

In addition to the read count matrices provided for the above example, all parental and hybrid HyLiTE read count matrices used in the analysis of all representative systems in the associated research project are available [here](https://github.com/annabehling/DEA_and_fit/tree/master/all_count_matrices "all_count_matrices/").

The `all_count_matrices/` folder contains one parental and two replicate hybrid count matrices for a representative system from each of:

* **allopolyploid fungi** (file prefix 'allo_f')
* **homoploid hybrid fungi** (file prefix 'HH_f')
* **allopolyploid plants** (file prefix 'allo_p')
* **homoploid hybrid plants** (file prefix 'HH_p')
* **allopolyploid animals** (file prefix 'allo_a')
* **homoploid hybrid animals** (file prefix 'HH_a')

More information about the raw genomic and RNA-seq data used to generate these read count matrices can be found [here].

There is also a `functions_usage.R` file that contains the usage of the `functions.R` code on all of the additional available data in `all_count_matrices/`.

## Additional data checks

The `additional_checks.R` file contains code and usage for two validations, run on the `example_files/` files:

1. Testing the arbitrariness of the parental definition
2. Testing for extraneous effects of no independent filtering

## Acknowledgements

Thank you to [David Winter](https://github.com/dwinter "github.com/dwinter") for the [hyliter](https://github.com/dwinter/hyliter "github.com/dwinter/hyliter") code and assistance with its implementation.
