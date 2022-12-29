
"Retrofit One sentence (with paper link)"

## Overview

"`retrofit` overview"

![Figure1](https://user-images.githubusercontent.com/90921267/209952993-4a4a5e49-9638-4dee-acdc-ce6a1f18d870.png)

## Installation

Install `retrofit` from github:

``` r
install.packages("devtools") 
devtools::install_github("qunhualilab/retrofit")
```
Install `retrofit` from [Bioconductor](https://bioconductor.org/packages/devel/bioc/html/retrofit.html).

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')
BiocManager::install("retrofit")
```

## Example

``` r
library(retrofit)
## load built in data
data(ReadmeData)
x           = data$extra5_loc_x
sc_ref      = data$sc_ref
marker_ref  = data$marker_ref

## decompose 
res         = retrofit::decompose(x, L=16, iterations=100, verbose=TRUE)
W           = res$w
H           = res$h
TH          = res$th

## annotate with correlations
res         = retrofit::annotateWithCorrelations(sc_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated = res$w
H_annotated = res$h
cells       = res$ranked_cells

## annotate with markers
res         = retrofit::annotateWithMarkers(marker_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated = res$w
H_annotated = res$h
cells       = res$ranked_cells	  
```
