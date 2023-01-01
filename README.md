RETROFIT enables reference-free cell-type deconvolution of spatial transcriptomics data. The overall approach is detailed in the paper script. [Link to paper]

## Overview

Figure below shows the method schematic for RETROFIT. RETROFIT (Reference-freE spatial TRanscriptOmic FactorIzaTion) only takes as input a ST count matrix $X$ of $G$ genes at $S$ spots, and then projects $X$ to a low-dimension space spanned by L non-negative latent components in an unsupervised manner, where the expression of each gene at each spot for each component is further decomposed into (1) expression specific to this gene and (2) background expression shared by all genes. A tuning parameter $\lamda \ge 0$ controls the size of background expression. Like any unsupervised learning, RETROFIT produces unlabeled results. To annotate them, we develop two simple post hoc strategies of matching the latent components inferred by RETROFIT to known cell types. The first annotation strategy requires a single-cell expression reference and the second strategy requires a curated list of cell-type-specific marker genes for all K cell types present in the ST data. Read more about the strategies in Methods section of the paper.

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
x           = ReadmeData$extra5_loc_x
sc_ref      = ReadmeData$sc_ref
marker_ref  = ReadmeData$marker_ref

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

More details can be found in the vignettes.
[Link to vignettes]
