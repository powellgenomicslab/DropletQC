#' Identify damaged cells
#'
#' @param nf_umi_ed
#'
#' @return
#' @export
#' @importFrom mclust mclustBIC
#'
#' @examples #1

data("qc_examples")
gbm <- qc_examples[qc_examples$sample=="GBM",]
gbm.ed <- data.frame(nf=gbm$nuclear_fraction_droplet_qc, umi=gbm$umi_count)
gbm.ed <- identify_empty_drops(nf_umi = gbm.ed)
gbm.ed$cell_type <- gbm$cell_type

nf_umi_ed_ct <- gbm.ed

identify_damaged_cells <- function(nf_umi_ed_ct,
                                   nf_sep,
                                   umi_sep){

  # Check and parse arguments

  # Extract data for EM
  em.data <- data.frame(nf = unlist(nf_umi_ed_ct[,1], use.names = FALSE),
                        umi = log10(unlist(nf_umi_ed_ct[,2], use.names = FALSE)),
                        ct = unlist(nf_umi_ed_ct[,4], use.names = FALSE))
  row.names(em.data) <- 1:nrow(em.data)

  # Filter out any empty droplets
  em.data <- em.data[nf_umi_ed_ct[,3]=="cell",]

  # Split by cell type
  em.data.ct <-split(em.data, em.data$ct)

  # Run EM
  em_mod <- lapply(em.data.ct, function(x) mclust::Mclust(data = x[,1:2], G = 1:2, modelNames = "EEI"))
  plot(em_mod[[1]], "classification")
  plot(em_mod[[2]], "classification")
  plot(em_mod[[3]], "classification")
  plot(em_mod[[4]], "classification")
  plot(em_mod[[5]], "classification")

  # Check models pass rescue thresholds

  return(em_mod)
}


