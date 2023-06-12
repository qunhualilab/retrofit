# RETROFIT: Reference-free deconvolution of cell-type mixtures in spatial transcriptomics

Welcome to RETROFIT, [a *Bioconductor* package](https://doi.org/doi:10.18129/B9.bioc.retrofit) for reference-free learning of cell-type composition and cell-type-specific gene expression in spatial transcriptomics (ST).

If you find this R package or any part of this repository useful for your work,
please kindly cite the following research article:

> Roopali Singh, Xi He, Adam Keebum Park, Ross Cameron Hardison, Xiang Zhu, Qunhua Li.
> RETROFIT: Reference-free deconvolution of cell-type mixtures in spatial transcriptomics.
> *bioRxiv* (2023) <https://doi.org/10.1101/2023.06.07.544126>.

Correspondence should be addressed to X.Z. (`xiangzhu[at]psu.edu`) and Q.L. (`qunhua.li[at]psu.edu`).

## Overview

ST profiles gene expression in intact tissues. However, ST data measured at each spatial location may represent gene expression of multiple cell types, making it difficult to identify cell-type-specific transcriptional variation across spatial contexts. Existing cell-type deconvolutions of ST data often require single-cell transcriptomic references, which can be limited by availability, completeness and platform effect of such references. We present RETROFIT, a reference-free Bayesian method that produces sparse and interpretable solutions to deconvolve cell types underlying each location independent of single-cell transcriptomic references.

Built on a Bayesian hierarchical model, RETROFIT deconvolves the ST data matrix ($X$) into two matrices: one reflecting cell-type-specific gene expression ($\widetilde{W}$) and the other reflecting the proportion of each cell type ($\widetilde{H}$). Implemented with a structured stochastic variational inference algorithm, RETROFIT scales well with thousands of genes and spots in a typical ST dataset.

The following figure shows the method schematic of RETROFIT. First, RETROFIT takes a ST data matrix as the only input and decomposes this matrix into latent components in an unsupervised manner, without using any single-cell transcriptomic information. Second, RETROFIT matches these latent components to known cell types using either a cell-type-specific gene expression reference or a list of cell-type-specific marker genes for the cell types present in the ST sample. Lastly, RETROFIT outputs a cell-type-specific gene expression matrix and a cell-type proportion matrix.

![Figure1](https://user-images.githubusercontent.com/90921267/220766755-daea9d4b-4ac0-4dd3-978c-7e71b31bc36e.png)

## Installation

To install `retrofit` from [Bioconductor](https://doi.org/doi:10.18129/B9.bioc.retrofit),
please start R (version "4.3") and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("retrofit")
```

Alternatively, please follow these steps to install `retrofit` from GitHub:

``` r
install.packages("devtools") 
devtools::install_github("qunhualilab/retrofit")
```

## Example

Use the following code to apply RETROFIT to your data.

``` r
library(retrofit)

## load the built-in demo data
utils::data(testSimulationData)
x           <- testSimulationData$extra5_x
sc_ref      <- testSimulationData$sc_ref
marker_ref  <- testSimulationData$marker_ref

## decompose the ST data matrix into latent components 
res         <- retrofit::decompose(x, L=16, iterations=100, verbose=TRUE)
W           <- res$w
H           <- res$h
TH          <- res$th

## match the latent components to known cell types using a cell-type-specific gene expression reference
res         <- retrofit::annotateWithCorrelations(sc_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells

## match the latent components to known cell types using a list of cell-type-specific marker genes
res         <- retrofit::annotateWithMarkers(marker_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells	  
```

## Tutorial
- [Simulation Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/SimulationVignette.Rmd) is designed to get started with RETROFIT and understand its usage reproducing main results showcased in the paper. A simple practice of following the codes will provide an overall picture of the basic but comprehensive scenario. 

## Paper Reproducibility
- [Colon Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/ColonVignette.Rmd) is a slightly more advanced vignette as it not only uses RETROFIT to deconvolve a real ST data but also reproduces certain results from the paper. This vignette utilizes real data from Human Colon tissue generated using the 10x Genomics Visium platform in this [Paper](https://www.sciencedirect.com/science/article/pii/S009286742031686X). Herein, we demonstrate that our method is effective in identifying biologically relevant spatial patterns using ST tissues.
