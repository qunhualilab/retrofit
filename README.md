RETROFIT enables reference-free cell-type deconvolution of spatial transcriptomics data.
##
![Figure1](https://user-images.githubusercontent.com/90921267/220766755-daea9d4b-4ac0-4dd3-978c-7e71b31bc36e.png) <br />

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
utils::data(testSimulationData)
x           <- testSimulationData$extra5_x
sc_ref      <- testSimulationData$sc_ref
marker_ref  <- testSimulationData$marker_ref

## decompose 
res         <- retrofit::decompose(x, L=16, iterations=100, verbose=TRUE)
W           <- res$w
H           <- res$h
TH          <- res$th

## annotate with correlations
res         <- retrofit::annotateWithCorrelations(sc_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells

## annotate with markers
res         <- retrofit::annotateWithMarkers(marker_ref, K=8, decomp_w=W, decomp_h=H)
W_annotated <- res$w
H_annotated <- res$h
cells       <- res$ranked_cells	  
```

More details can be found in the vignettes.
- [Simulation Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/SimulationVignette.Rmd)
- [Colon Vignette](https://github.com/qunhualilab/retrofit/blob/main/vignettes/ColonVignette.Rmd)
