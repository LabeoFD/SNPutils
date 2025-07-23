# SNPutils

<!-- badges: start -->
<!-- badges: end -->

## Overview

SNPutils provides tools for converting numerical genotype data (0, 1, 2) to allele pairs (e.g., "A/G"). The package streamlines genetic data processing workflows by:

- Translating genotype codes into biologically meaningful alleles
- Supporting standard validation rules for SNP data
- Handling both standard SNPs and indels (insertions/deletions)
- Generating TYP format output for downstream analysis

## Installation

You can install SNPutils from GitHub:

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("LabeoFD/SNPutils")
```

## Quick Start
```r
library(SNPutils)

# Load example data
data(example_genotypes)
data(example_annotations)

# Basic workflow
results <- genotype2allele(example_genotypes, example_annotations)
typ_output <- allele2typ(results)

# Check the first few results
head(typ_output)
```

## Key Features

- **Efficient conversion**: Transform numerical codes to nucleotide pairs
- **Data validation**: Apply standard genetic rules to ensure data quality
- **Indel support**: Correctly handle insertions and deletions
- **Formatted output**: Generate standardized TYP format files

## Citation
If you use SNPutils in your research, please cite it as:

> LabeoFD (2025). SNPutils: Tools for Converting Numerical Genotype Data to Allele Pairs. R package.


## License
This package is free and open source software, licensed under GPL-3.
