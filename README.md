# What is IMPRINTS.PhosphoQP ?

This is an R package in the continuity of the packages IMPRINTS.CETSA and IMPRINTS.CETSA.app.It will help you process phosphoproteomics data and quantitative proteomics data in the IMPRINTS-CETSA format and will allow you to compare them.

## How to install IMPRINTS.PhosphoQP ?  
First go to Rstudio. This package depends on the last version of the IMRPRINTS.CETSA and IMPRINTS.CETSA.app packages. If you havn't installed them yet, please do so.
For installation from github https://github.com/nkdailingyun/IMPRINTS.CETSA and https://github.com/mgerault/IMPRINTS.CETSA.app .

When this first step is done, you can install IMPRINTS.PhosphoQP with:

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools")
} 
devtools::install_github("mgerault/IMPRINTS.PhosphoQP")

library(IMPRINTS.PhosphoQP)
```

You are now ready to use IMPRINTS.PhosphoQP. To access its documentation type the following :

```c
vignette('IMPRINTS.PhosphoQP_doc', package = 'IMPRINTS.PhosphoQP')
```
