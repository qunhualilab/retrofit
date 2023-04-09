#' SimulationVignette process supporting data
#'
#' A dataset based on "mouse brain Slide-seq"[Roopali] that supports the process of following SimulationVignette.
#' Items include input sequencing data to be decomposed, and resulting data requiring expensive computations. 
#'
#' \itemize{
#'   \item vignetteSimulationData
#' }
#'
#' @docType data
#' @name vignetteSimulationData
#' @usage data(vignetteSimulationData)
#' @format ## `vignetteSimulationData`
#' A data frame with input, references, decomposed results of "mouse brain Slide-seq"[Roopali].
#' \describe{
#'   \item{n10m3_x}{500 x 1000 GeneExpressions x Spots}
#'   \item{sc_ref_w}{500 x 10 GeneExpressions x Cell types}
#'   \item{sc_ref_h}{10 x 1000 Cell types x Spots}
#'   \item{results_4k_iterations/}{Results of 4k iterations computation}
#'   \item{results_4k_iterations/decompose/h}{20 x 1000 Components x Spots}
#'   \item{results_4k_iterations/decompose/w}{500 x 20 GeneExpressions x Components}
#'   \item{results_4k_iterations/decompose/th}{20 Components}
#' }
#' @source [Roopali] < include our future preprint as the source, cited as "[Full Author lists], [Title], Preprint Forthcoming (2023)".>
"vignetteSimulationData"


#' ColonVignette process supporting data
#'
#' A dataset based on "human intestine Visium"[Roopali] that supports the process of following ColonVignette.
#' Items include input sequencing data to be decomposed, and resulting data requiring expensive computations. 
#'
#' \itemize{
#'   \item vignetteColonData
#' }
#'
#' @docType data
#' @name vignetteColonData
#' @usage data(vignetteColonData)
#' @format ## `vignetteColonData`
#' A data frame with input, references, decomposed results of "human intestine Visium"[Roopali].
#' \describe{
#'   \item{a3_x}{722 x 1080 GeneExpressions x Spots}
#'   \item{a3_coords}{1080 x 5 [Roopali]}
#'   \item{sc_ref}{722 x 8 GeneExpressions x Cell types}
#'   \item{marker_ref}{Dictionary of 8 Cell Type:Genes}
#'   \item{marker_ref_df}{Rows of Cell Type:Genes}
#'   \item{fetal12_genes}{14 strings [Roopali]}
#'   \item{fetal12_genes}{27 strings [Roopali]}
#'   \item{a3_results_4k_iterations}{Results of 4k iterations computation on a3}
#'   \item{a3_results_4k_iterations/decompose/h}{16 x 1080 Components x Spots}
#'   \item{a3_results_4k_iterations/decompose/w}{722 x 16 GeneExpressions x Components}
#'   \item{a3_results_4k_iterations/decompose/th}{16 Components}
#'   \item{a3_results_4k_iterations/annotateWithCorrelations/w}{413 x 10 GeneExpressions x Cell Types}
#'   \item{a1_results_4k_iterations}{Results of 4k iterations computation on a1}
#'   \item{a1_results_4k_iterations/annotateWithCorrelations/w}{829 x 10 GeneExpressions x Cell Types}
#'   \item{intestine_w_12pcw}{33538 x 8 [Roopali]}
#' }
#' @source [Roopali] < include our future preprint as the source, cited as "[Full Author lists], [Title], Preprint Forthcoming (2023)".>
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
#' A data frame with input, references, decomposed results based on "mouse brain Slide-seq"[Roopali].
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
#' @source [Roopali] < include our future preprint as the source, cited as "[Full Author lists], [Title], Preprint Forthcoming (2023)".>
"testSimulationData"
