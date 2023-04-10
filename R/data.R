#' SimulationVignette process supporting data
#'
#' A dataset that supports the process of following SimulationVignette.
#' Items include input synthetic ST data to be decomposed, and resulting data requiring expensive computations. The synthetic sample data and reference is created using scRNA-seq data published in [paper](https://pubmed.ncbi.nlm.nih.gov/30923225/).
#'
#' \itemize{
#'   \item vignetteSimulationData
#' }
#'
#' @docType data
#' @name vignetteSimulationData
#' @usage data(vignetteSimulationData)
#' @format ## `vignetteSimulationData`
#' A data frame with sample input for RETROFIT, reference for annotation purposes, decomposed results of the sample data using 4000 iterations and true cell type proportions for the sample input data for analysis.
#' \describe{
#'   \item{n10m3_x}{500 x 1000 GeneExpressions x Spots}
#'   \item{sc_ref_w}{500 x 10 GeneExpressions x Cell types}
#'   \item{sc_ref_h}{10 x 1000 Cell types x Spots}
#'   \item{results_4k_iterations/}{Results of 4k iterations computation}
#'   \item{results_4k_iterations/decompose/h}{20 x 1000 Components x Spots}
#'   \item{results_4k_iterations/decompose/w}{500 x 20 GeneExpressions x Components}
#'   \item{results_4k_iterations/decompose/th}{20 Components}
#' }
#' @source Roopali Singh, Xi He, Adam Keebum Park, Ross Cameron Hardison, Xiang Zhu, Qunhua Li. RETROFIT: Reference-free deconvolution of cell-type mixtures in spatial transcriptomics, Preprint Forthcoming (2023)
#' @source Rodriques, S. G. et al. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463–1467 (2019)
"vignetteSimulationData"


#' ColonVignette process supporting data
#'
#' A dataset that supports the process of following ColonVignette.
#' Items include input ST data tissue to be decomposed, and resulting data requiring expensive computations. The ST data and reference is from this [Paper](https://www.sciencedirect.com/science/article/pii/S009286742031686X). 
#'
#' \itemize{
#'   \item vignetteColonData
#' }
#'
#' @docType data
#' @name vignetteColonData
#' @usage data(vignetteColonData)
#' @format ## `vignetteColonData`
#' A data frame with the input for RETROFIT, reference for annotation purposes, decomposed results of the data using 4000 iterations.
#' \describe{
#'   \item{a3_x}{722 x 1080 GeneExpressions x Spots}
#'   \item{a3_coords}{1080 Spots x 5 Atriibutes}
#'   \item{sc_ref}{722 x 8 GeneExpressions x Cell types}
#'   \item{marker_ref}{Dictionary of 8 Cell Type:Genes}
#'   \item{marker_ref_df}{Rows of Cell Type:Genes}
#'   \item{fetal12_genes}{14 Genes}
#'   \item{fetal12_genes}{27 Genes}
#'   \item{a3_results_4k_iterations}{Results of 4k iterations computation on a3}
#'   \item{a3_results_4k_iterations/decompose/h}{16 x 1080 Components x Spots}
#'   \item{a3_results_4k_iterations/decompose/w}{722 x 16 GeneExpressions x Components}
#'   \item{a3_results_4k_iterations/decompose/th}{16 Components}
#'   \item{a3_results_4k_iterations/annotateWithCorrelations/w}{413 x 10 GeneExpressions x Cell Types}
#'   \item{a1_results_4k_iterations}{Results of 4k iterations computation on a1}
#'   \item{a1_results_4k_iterations/annotateWithCorrelations/w}{829 x 10 GeneExpressions x Cell Types}
#'   \item{intestine_w_12pcw}{33538 Genes x 8 Cell Types}
#' }
#' @source Roopali Singh, Xi He, Adam Keebum Park, Ross Cameron Hardison, Xiang Zhu, Qunhua Li. RETROFIT: Reference-free deconvolution of cell-type mixtures in spatial transcriptomics, Preprint Forthcoming (2023)
#' @source Fawkner-Corbett, D. et al. Spatiotemporal analysis of human intestinal development at single-cell resolution. Cell 184, 810–826 (2021)
"vignetteColonData"

#' Test simulation data
#'
#' A dataset with input and output of retrofit functions for unit&reproducibility testing
#'
#' \itemize{
#'   \item testSimulationData
#' }
#'
#' @docType data
#' @name testSimulationData
#' @usage data(testSimulationData)
#' @format ## `testSimulationData`
#' A data frame with the input for RETROFIT, reference for annotation purposes, decomposed results of the data using 4000 iterations.
#' \describe{
#'   \item{extra5_x}{500 x 1000 GeneExpressions x Spots}
#'   \item{sc_ref}{500 x 10 GeneExpressions x Cell types}
#'   \item{marker_ref}{Dictionary of 10 Cell Type:Genes}
#'   \item{decompose/}{Results of decomposition}
#'   \item{decompose/w}{500 x 16 GeneExpressions x Components}
#'   \item{decompose/h}{16 x 1000 Components x Spots}
#'   \item{decompose/th}{16 Components}
#'   \item{annotateWithCorrelations/}{Results of cell type annotation with correlations}
#'   \item{annotateWithCorrelations/w}{500 x 10 GeneExpressions x Cell Types}
#'   \item{annotateWithCorrelations/h}{10 x 1000 Cell Types x Spots}
#'   \item{annotateWithCorrelations/w_prop}{500 x 10 (normalized) GeneExpressions x Cell Types}
#'   \item{annotateWithCorrelations/h_prop}{10 x 1000 (normalized) Cell Types x Spots}
#'   \item{annotateWithCorrelations/ranked_cells}{10 Cell Types by rank}
#'   \item{annotateWithCorrelations/ranked_correlations}{10 Correlations by rank}
#'   \item{annotateWithCorrelations/sc_ref_prop}{500 x 10 (normalized) GeneExpressions x Cell types}
#'   \item{annotateWithMarkers/}{Results of cell type annotation with markers}
#'   \item{annotateWithMarkers/w}{500 x 10 GeneExpressions x Cell Types}
#'   \item{annotateWithMarkers/h}{10 x 1000 Cell Types x Spots}
#'   \item{annotateWithMarkers/w_prop}{500 x 10 (normalized) GeneExpressions x Cell Types}
#'   \item{annotateWithMarkers/h_prop}{10 x 1000 (normalized) Cell Types x Spots}
#'   \item{annotateWithMarkers/ranked_cells}{10 Cell Types by rank}
#'   \item{annotateWithMarkers/gene_sums}{10 x 16 Cell Types vs Components gene expression sums}
#' }
#' @source Roopali Singh, Xi He, Adam Keebum Park, Ross Cameron Hardison, Xiang Zhu, Qunhua Li. RETROFIT: Reference-free deconvolution of cell-type mixtures in spatial transcriptomics, Preprint Forthcoming (2023)
"testSimulationData"
