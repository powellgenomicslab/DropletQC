#' Identify empty droplets
#'
#' @description This function is used to identify a suitable nuclear fraction
#'   cut-off point to guide the identification of empty droplets. To do this it
#'   calculates the kernel density estimate of the input nuclear fraction scores
#'   and identifies the trough after the first peak, which is assumed to
#'   represent the population of empty droplets.
#'
#' @param nf_umi data frame, containing two columns; the nuclear fraction
#'   estimates in the first column and the total UMI count for each barcode in
#'   the second column
#' @param nf_rescue numeric, a rescue parameter defining a minimum nuclear
#'   fraction score between zero and one. This is used in combination with
#'   `umi_rescue` to identify cells that were misidentified as empty droplets
#' @param umi_rescue integer, a rescue parameter defining a minimum UMI count.
#'   This is used in combination with `nf_rescue` to identify cells that were
#'   misidentified as empty droplets
#' @param include_plot logical, whether or not to produce a plot illustrating
#'   how the nuclear fraction threshold was identified and which barcodes have
#'   been called as empty droplets. In the plot of nuclear fraction vs log10(UMI
#'   counts), empty droplets are expected to occupy the lower left corner of the
#'   plot.
#'
#' @return data frame, the original data frame is returned plus an additional
#'   column identifying each barcode as a "cell" or "empty_droplet"
#' @export
#'
#' @examples
#' data("qc_examples")
#' gbm <- qc_examples[qc_examples$sample=="GBM",]
#' gbm <- data.frame(nf=gbm$nuclear_fraction_droplet_qc, umi=gbm$umi_count)
#' gbm.ed <- identify_empty_drops(nf_umi = gbm)
#' head(gbm.ed)
#' table(gbm.ed$cell_status)
#'
identify_empty_drops <-
  function(nf_umi,
           nf_rescue = 0.05,
           umi_rescue = 1000,
           include_plot = FALSE) {

    ## Check and parse arguments
    if (any(class(nf_umi) == "data.frame")) {

      # Assume nuclear fraction is in the first column
      nf <- unlist(nf_umi[, 1], use.names = FALSE)
      # Assume UMI counts are in the second column
      umi <- unlist(nf_umi[, 2], use.names = FALSE)

      # Check values are reasonable
      if(any(c(max(nf)>1, min(nf)<0))){
        warning(paste0("The nuclear fraction values provided in the first column of 'nf_umi' should be between 0 and 1, but values outside this range were identified : minimum = ",min(nf),", maximum = ",max(nf)), call.=FALSE)
      }

      if(!all(umi == floor(umi))){
        non_integer_examples <- which(umi != floor(umi))
        if(length(non_integer_examples)>5){
          non_integer_examples <- non_integer_examples[1:5]
          }
        non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
        warning(paste0("Non-integer values detected in the second column of 'nf_umi' (e.g. ",non_integer_examples,") where umi counts were expected"), call.=FALSE)
      }

      if(max(umi)<100){
        umi
        warning(paste0("The total umi counts provided in the second column of 'nf_umi' appear to be quite low (max = ",max(umi),"), are these the total UMI counts per cell?"), call.=FALSE)
      }

    } else {
      stop(paste0("A data frame should be supplied to the nf_umi argument, but an object of class ",paste(class(nf_umi), collapse = "/")," was provided"), call.=FALSE)
    }

  # Density estimation (automatically chosen bandwidth)
  kdde_0 <- ks::kdde(x = nf, deriv.order = 0)
  kdde_0 <- data.frame(estimate = kdde_0[["estimate"]],
                       eval.points = kdde_0[["eval.points"]])

  # Density derivative estimation (automatically chosen bandwidth, but different
  # from kdde_0!)
  kdde_1 <- ks::kdde(x = nf, deriv.order = 1)
  kdde_1 <- data.frame(estimate = kdde_1[["estimate"]],
                       eval.points = kdde_1[["eval.points"]])

  # Find point to place cut-off between empty droplets and cells
  gradient_sign <- rle(kdde_1[["estimate"]]>0)
  nf_cutoff <- kdde_1[["eval.points"]][sum(gradient_sign[["lengths"]][1:2])]

  ## Check if there is more than one peak
  if(length(gradient_sign$values)<4){
    warning(paste0("Could not detect more than one peak in the nuclear fraction distribution. There may not be any empty droplets present. We suggest visualising the density estimation (include_plot=TRUE)."), call.=FALSE)
  }

  # Label cells
  nf_umi$cell_status <- "cell"
  nf_umi$cell_status[nf < nf_cutoff] <- "empty_droplet"

  # Rescue miscalled cells
  nf_umi$cell_status[nf > nf_rescue &  umi > umi_rescue] <- "cell"

  # Plots
  if (include_plot) {
    p1 <-
      ggplot2::ggplot(kdde_0, ggplot2::aes(x = eval.points, y = estimate)) +
      ggplot2::geom_line() +
      ggplot2::ggtitle("Density estimation") +
      ggplot2::xlab("Nuclear fraction") +
      ggplot2::ylab("Density function") +
      ggplot2::geom_vline(xintercept = nf_cutoff,
                          col = "dodgerblue",linetype = "dashed")

    p2 <-
      ggplot2::ggplot(kdde_1, ggplot2::aes(x = eval.points, y = estimate)) +
      ggplot2::geom_line() +
      ggplot2::ggtitle("Density derivative estimation") +
      ggplot2::xlab("Nuclear fraction") +
      ggplot2::ylab("Density derivative function") +
      ggplot2::geom_vline(xintercept = nf_cutoff,
                          col = "dodgerblue",linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0, col = "grey")

    p3_data <-
      data.frame(nf = nf,
                 umi = umi,
                 cell_status = nf_umi$cell_status)
    p3 <-
      ggplot2::ggplot(p3_data,
                      ggplot2::aes(
                        x = nf,
                        y = log10(umi),
                        colour = cell_status
                      )) +
      ggplot2::geom_point() +
      ggplot2::xlim(0, 1) +
      ggplot2::xlab("Nuclear fraction") +
      ggplot2::ylab("log10(UMI count)") +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    # Add three dashed lines to indicate threshold + rescue area
    l1 <- data.frame(x_start = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                     x_end = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                     y_start = log10(min(umi)),
                     y_end = max(log10(umi_rescue), log10(min(umi))))

    l2 <- data.frame(x_start = ifelse(umi_rescue<min(umi) & nf_rescue<nf_cutoff, max(min(nf),nf_rescue), nf_cutoff),
                     x_end = min(max(min(nf),nf_rescue), nf_cutoff),
                     y_start = max(log10(umi_rescue), log10(min(umi))),
                     y_end = max(log10(umi_rescue), log10(min(umi))))

    l3 <- data.frame(x_start = min(max(min(nf),nf_rescue), nf_cutoff),
                     x_end = min(max(min(nf),nf_rescue), nf_cutoff),
                     y_start = max(log10(umi_rescue), log10(min(umi))),
                     y_end = max(log10(umi)))

    p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
      x = l1$x_start,
      y = l1$y_start,
      xend = l1$x_end,
      yend = l1$y_end),
      colour = "dodgerblue",
      linetype = "dashed",
      data = p3_data)

    p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
      x = l2$x_start,
      y = l2$y_start,
      xend = l2$x_end,
      yend = l2$y_end),
      colour = "dodgerblue",
      linetype = "dashed",
      data = p3_data)

    p3 <- p3 + ggplot2::geom_segment(ggplot2::aes(
      x = l3$x_start,
      y = l3$y_start,
      xend = l3$x_end,
      yend = l3$y_end),
      colour = "dodgerblue",
      linetype = "dashed",
      data = p3_data)

    p.grid <- ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, ncol = 2),
                                nrow = 2,
                                p3)
    print(p.grid)
  }

  return(nf_umi)
  }
