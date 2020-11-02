# DEA_and_fit

## Description

This code takes the read count matrices outputted by [HyLiTE](https://hylite.sourceforge.io/) analyses, performs differential expression analyses, fits the regression models, and classifies each gene into one of the four gene expression categories or NA.

The majority of functions were taken from [hyliter](https://github.com/dwinter/hyliter), sometimes with modifications made specific for this project. Where applicable, this has been commented in the R code file.

## Installation

To install the required scripts, first clone the **DEA_and_fit** repository.
```
git clone https://github.com/annabehling/DEA_and_fit
```

## The HyLiTE output

HyLiTE outputs a number of files. An outline of each output file and what information they contain can be found[here](https://hylite.sourceforge.io/outformat.html#outformat "HyLiTE output formats").
The important output files for the following analyses are the `expression.txt` and `read.summary.txt` files.

