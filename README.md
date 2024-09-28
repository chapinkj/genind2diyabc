# genind2diyabc

## About
A function to convert genind object to DIYABC input format

## Usage
Accepts a genind object, outputs DIYABC format.

Optionally exclude intrapopulation monomorphic loci, and/or missing data with `remove_intrapop_mono` and `remove_intrapop_allmissing`.

Results are exported as a .snp file,

Optionally provide a `file_name`. If left bank, a name will be generated based on the population names in the genind object

Build in R 4.3

## Citations

[DIYABC](https://diyabc.github.io/)

[adegenet](https://cran.r-project.org/web/packages/adegenet/index.html)