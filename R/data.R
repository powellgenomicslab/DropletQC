#' QC metrics from four 10x Genomics single cell gene expression datasets.
#'
#' This dataset contains a collection of summary and QC metrics from four
#' publicly available 10x Genomics single cell gene expression datasets:
#'
#' 1. Human Glioblastoma Multiforme: 3’v3 Whole Transcriptome Analysis
#' 2. Hodgkin's Lymphoma, Dissociated Tumor: Whole Transcriptome Analysis
#' 3. 10k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells,
#'    Dual Indexed
#' 4. PBMCs from a Healthy Donor: Whole Transcriptome Analysis
#'
#' To be included in the dataset, cells were required to pass
#' DropletUtils::emptyDrops using the default FDR threshold of 1 and have a
#' maximum of 15 mitochondrial gene content. The included variables are as
#' follows:
#'
#' @format A data frame with 27,597 rows and 15 variables:
#' \describe{
#'   \item{sample}{sample name, one of four values; GBM (glioblastoma), HL
#'   (Hodgkin's lymphoma), MB (mouse brain) and PBMC (peripheral blood
#'   mononuclear cells)}
#'   \item{cell_barcode}{the 16 nucleotide cell barcode e.g. AAACCCAAGGCGATAC-1}
#'   \item{umap_1}{UMAP coordinates 1, for visualisation}
#'   \item{umap_2}{UMAP coordinates 2, for visualisation}
#'   \item{seurat_clusters}{cell clusters identified with the Louvain algorithm
#'   using default parameters implemented in the Seurat package}
#'   \item{cell_type}{a very rough cell type annotation}
#'   \item{umi_count}{the number of UMIs detected}
#'   \item{log10_umi_count}{log10 of the number of UMIs detected}
#'   \item{percent_mt}{percentage of UMIs mapping to mitochondrial genes}
#'   \item{empty_drops_log_prob}{the log-probability of observing the cell's
#'    count vector under the null model - from DropletUtils::emptyDrops output}
#'   \item{empty_drops_p_value}{the Monte Carlo p-value against the null model -
#'   from DropletUtils::emptyDrops output}
#'   \item{empty_drops_fdr}{FDR values returned from DropletUtils::emptyDrops
#'   output}
#'   \item{nuclear_fraction_droplet_qc}{nuclear fraction metric calculated with
#'   dropletQC}
#'   \item{nuclear_fraction_velocyto}{nuclear fraction metric calculated using
#'   the output from velocyto}
#'   \item{flag}{assigned cell status, taking one of three values; cell,
#'   damaged_cell, empty_droplet}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/}
"qc_examples"
