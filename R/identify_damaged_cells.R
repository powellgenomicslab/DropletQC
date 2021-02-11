#' Assess EM model results and assign cells as damaged or not
#'
#'@description This function accepts the EM results for a particular cell type,
#'  as well as the UMI and nuclear fraction thresholds, and assigns cells as
#'  "cell" or "damaged_cell" This is a helper function called by
#'  `identify_damaged_cells` and was not intended for more general use.
#'
#' @param em Mclust, result of running EM on the log10(UMI counts) and the
#'   nuclear fraction to asses whether there is likely to be two cell
#'   populations; cells and damaged cells, or just cells
#' @param umi_thresh numeric, percentage of the cell disitribution UMI counts
#'   which the damaged_cell distribution must be below in order to be accepted
#'   as damaged cells. For example if the mean UMI count of the distribution fit
#'   to the cell population is 1,000 and `umi_thresh` is set to 30, the
#'   damaged_cell distribution mean must be less than 700 to be classified as
#'   damaged_cells, otherwise all barcodes will be returned as cells
#' @param nf_thresh numeric, the minimum difference in the nuclear fraction
#'   score required between the two distributions (cells and damaged_cells) for
#'   the damaged_cell distribution to be accepted as real
#'
#' @return data frame, with three columns containing the log10(UMI counts),
#'   nuclear fraction score and the assigned cell status; "cell" or
#'   "damaged_cell"
#'
#' @keywords internal
assess_EM <- function(em, umi_thresh, nf_thresh){

  # Separately for each cell type, assign a barcode as "cell" or "damaged_cell"
  # based on the EM results using the following sequential procedure:

  # 1. If the EM model selected contained only one distribution, score all
  # barcodes as "cell".
  check_1 <- em$G==2

  # 2. If two distirbutions were fit, check the distribution with the higher
  # nuclear_fraction mean also has a lower umi mean. If this criteria is
  # satisfied, we consider the population with the lower umi count and higher
  # nuclear_fraction scores the damaged_cell population and move on to step 3.
  # If this is not the case we score all barcodes as "cell".
  if(check_1){
    nf_means <- em$parameters$mean["nf",]
    umi_means <- em$parameters$mean["umi",]
    check_2 <- umi_means[which.max(nf_means)] < umi_means[which.min(nf_means)]
  } else {
    check_2 <- FALSE
  }

  # 3. If 2. was satisfied, we check the damaged_cell nuclear_fraction
  # distribution mean is at least nf_thresh greater than the cell mean. Also
  # check the damaged_cell umi distribution mean is at least umi_thresh
  # percent lower than the cell umi distribution mean. If these two criteria
  # are satisfied we assign the damaged cells, otherwise all barcodes are
  # labelled "cell".
  if(check_2){
    # Check nuclear fraction threshold is satisfied
    nf_check <- nf_means[which.max(nf_means)] - nf_means[which.min(nf_means)] > nf_thresh
    # Check umi threshold is satisfied
    damaged_cell_umi <- 10^umi_means[which.max(nf_means)]
    cell_umi <- 10^umi_means[which.min(nf_means)]
    umi_check <- damaged_cell_umi < (cell_umi - cell_umi*(umi_thresh/100))

    if(all(nf_check, umi_check)){ check_3 <- TRUE } else { check_3 <- FALSE }
  } else {
    check_3 <- FALSE
  }

  # Assign barcodes
  em_classification <- data.frame(nf = em$data[,"nf"],
                                  umi = em$data[,"umi"],
                                  classification= em$classification)
  row.names(em_classification) <- row.names(em$data)

  if(all(check_1, check_2, check_3)){
    damaged_cells <- em_classification$classification==which.max(em$parameters$mean["nf",])
    em_classification$classification[damaged_cells] <- "damaged_cell"
    em_classification$classification[!damaged_cells] <- "cell"
  } else {
    em_classification$classification <- "cell"
  }

  return(em_classification)

}



#' Plot EM results
#'
#' @description This is a helper function called intarnally by the
#'   `identify_damaged_cells` function (if plots are requested). It's not
#'   intended for general use.
#'
#' @param emMclust Mclust, result of running EM on the log10(UMI counts) and the
#'   nuclear fraction
#' @param em_classified data frame, with three columns; log10(UMI counts),
#'   nuclear fraction and a column defining each cell as a "cell" or
#'   "damaged_cell" - this should be the output from running `assess_EM` on the
#'   model results
#' @param input_cell_type character, the name of the cell type
#'
#' @return ggarrange, three ggplots combined using `ggpubr::ggarrange`
#'
#' @keywords internal
plot_EM <- function(em, em_c, input_cell_type, umi_thresh, nf_thresh){

  # Get model parameters for plotting
  em_params <- list(nf_means = em$parameters$mean["nf",],
                    umi_means = em$parameters$mean["umi",],
                    nf_sd = rep(sqrt(em$parameters$variance$Sigma["nf",1]), em$G),
                    umi_sd = rep(sqrt(em$parameters$variance$Sigma["umi",2]), em$G),
                    amplitude = em$parameters$pro)
  # Add colours
  cell_colour <- "#4575b4"
  damaged_cell_colour <- "#d73027"
  if(em$G==2){
    em_params$distribution_colour <- c(damaged_cell_colour, cell_colour)[order(em_params$nf_means, decreasing = TRUE)]
  } else {
    em_params$distribution_colour <- cell_colour
  }

  # Plot

  ### nf vs. umi ###
  p1 <- ggplot2::ggplot(em_c, ggplot2::aes(x=nf, y=umi, colour=classification)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Nuclear fraction") +
    ggplot2::ylab("log10(UMI count)") +
    ggplot2::scale_colour_discrete(stringr::str_wrap(input_cell_type, width = 15))


  ### nf ###
  p2 <- ggplot2::ggplot(em_c, ggplot2::aes(x = nf)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, len = 51), colour = "black", fill = "white") +
    mapply(
      function(mean, sd, amplitude, n, binwidth, dist_colour) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * amplitude
          },
          size=1.5, geom="area", alpha=0.5, fill=dist_colour)
      },
      mean = em_params[["nf_means"]], #mean
      sd = em_params[["nf_sd"]], #standard deviation
      amplitude = em_params[["amplitude"]], #amplitude
      n = nrow(em_c), #sample size
      binwidth = 1/51, #binwidth used for histogram
      dist_colour = em_params[["distribution_colour"]]
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Nuclear fraction")
  # add lines for damaged cell nuclear fraction mean (solid, red) and threshold (dashed, blue)
  if(em$G==2){
    p2 <- p2 + ggplot2::geom_vline(xintercept = em_params[["nf_means"]][which.max(em_params[["nf_means"]])], colour = damaged_cell_colour)
    p2 <- p2 + ggplot2::geom_vline(xintercept = em_params[["nf_means"]][which.min(em_params[["nf_means"]])] + nf_thresh, linetype="dashed", colour = cell_colour)
    }


  ### umi ###
  p3 <- ggplot2::ggplot(em_c, ggplot2::aes(x = umi)) +
    ggplot2::geom_histogram(breaks = seq(min(em_c$umi), max(em_c$umi), len = 51), colour = "black", fill = "white") +
    mapply(
      function(mean, sd, amplitude, n, binwidth, dist_colour) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * amplitude
          },
          size=1.5, geom="area", alpha=0.5, fill=dist_colour)
      },
      mean = em_params[["umi_means"]], #mean
      sd = em_params[["umi_sd"]], #standard deviation
      amplitude = em_params[["amplitude"]], #amplitude
      n = nrow(em_c), #sample size
      binwidth = c(max(em_c$umi) - min(em_c$umi))/51, #binwidth used for histogram
      dist_colour = em_params[["distribution_colour"]]
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("log10(UMI count)")
  # add lines for damaged cell umi mean (solid, red) and threshold (dashed, blue)
  if(em$G==2){
    # add line for damaged cell umi mean (solid, red)
    p3 <- p3 + ggplot2::geom_vline(xintercept = em_params[["umi_means"]][which.max(em_params[["nf_means"]])], colour = damaged_cell_colour)
    # add lines for threshold (dashed, blue)
    damaged_cell_umi <- 10^(em_params[["umi_means"]][which.max(em_params[["nf_means"]])])
    cell_umi <- 10^(em_params[["umi_means"]][which.min(em_params[["nf_means"]])])
    umi_diff <- log10(cell_umi - cell_umi*(umi_thresh/100))
    p3 <- p3 + ggplot2::geom_vline(xintercept = umi_diff, linetype="dashed", colour = cell_colour)
  }

  p.grid <- ggpubr::ggarrange(ggpubr::ggarrange(p1, p2, p3, ncol = 3))
  return(p.grid)

}



#' Identify damaged cells
#'
#' @description This function uses a combination of the cell UMI counts and the
#'   nuclear fraction score to assign each cell one of two values; "cell" or
#'   "damaged_cell". This is based on the idea that damaged cells have a lower
#'   UMI count and higher nuclear fraction than whole cells. The expected input
#'   is a data frame with four columns. The first three columns should contain;
#'   the nuclear fraction score, total UMIs and a character vector describing
#'   each cell as "cell" or "empty_droplet". This is the format output by the
#'   `identify_empty_drops` function. The fourth column should be a character
#'   vector with user-assigned cell types. Internally, the provided data frame
#'   is split by cell type and a Gaussian mixture model with a maximum of two
#'   components is fit to the umi counts and nuclear fraction scores. The
#'   parameters of the model are estimated using expectation maximisation (EM)
#'   with the `mclust` package. The best model is selected using the Bayesian
#'   Information Criterion (BIC). The two populations (cells and damaged cells)
#'   are assumed to have equal variance (mclust model name "EEI").
#'
#' @param nf_umi_ed_ct data frame, with four columns. The first three columns
#'   should match the output from the `identify_empty_drops` function. The
#'   fourth column should contain cell type names.
#' @param nf_sep numeric, the minimum separation of the nuclear fraction score
#'   required between the cell and damage cell populations
#' @param umi_sep_perc numeric, this is the minimum percentage of UMIs which the
#'   damaged cell population is required to have compared to the cell
#'   population. For example, if the mean UMI of the distribution fit to the
#'   whole cell population is 10,000 UMIs, the mean of the distribution fit to
#'   the damaged cell population must be at less than 7,000 UMIs if the umi_sep
#'   parameter is 30 (%)
#' @param output_plots logical, whether or not to return plots
#' @param verbose logical, whether to print updates and progress while fitting
#'   with EM
#'
#' @return list, of length two. The first element in the list contains a data
#'   frame with the same dimensions input to the `nf_umi_ed_ct` argument, with
#'   "damaged_cell" now recorded in the third column. The second element is NULL
#'   unless `output_plots`=TRUE. If requested, three plots are returned for each
#'   cell type in a named list, combined using `ggpubr::ggarrange`. For each
#'   cell type, the first plot illustrates the cell and damaged cell populations
#'   (if any) in a plot of nuclear fraction vs log10(UMI counts). Damaged cells
#'   are expected to be in the lower right portion of the plot(lower UMI counts
#'   and higher nuclear fraction). The second and third plots show the model
#'   fits to the nuclear fraction and UMI count distributions respectively.
#'   Solid lines indicate the distribution mean, while dashed lines indicate the
#'   positions of the thrsholds controlled by the `nf_sep` and `umi_sep`
#'   parameters.
#'
#' @export
#' @importFrom mclust mclustBIC
#'
#' @examples #1
#' data("qc_examples")
#' gbm <- qc_examples[qc_examples$sample=="MB",]
#' gbm.ed <- gbm[,c("nuclear_fraction_droplet_qc","umi_count")]
#' gbm.ed <- identify_empty_drops(nf_umi = gbm.ed)
#' gbm.ed$cell_type <- gbm$cell_type
#' gbm.ed.dc <- identify_damaged_cells(gbm.ed, verbose=FALSE)
#' gbm.ed.dc <- gbm.ed.dc[[1]]
#' head(gbm.ed.dc)
#' table(gbm.ed.dc$cell_status)
#'
identify_damaged_cells <- function(nf_umi_ed_ct,
                                   nf_sep=0.15,
                                   umi_sep_perc=50, # UMI counts percentage less than cell
                                   output_plots=FALSE,
                                   verbose=TRUE){

  # Check nf_umi_ed_ct argument



  if (any(class(nf_umi_ed_ct) == "data.frame")) {

    # Check four columns exist
    if(ncol(nf_umi_ed_ct)!=4){
      stop(paste0("nf_umi_ed_ct should be a data frame with four columns, see function arguments"), call.=FALSE)
    }

    # Assume nuclear fraction is in the first column
    nf <- unlist(nf_umi_ed_ct[, 1], use.names = FALSE)
    # Assume UMI counts are in the second column
    umi <- unlist(nf_umi_ed_ct[, 2], use.names = FALSE)
    # Assume third column contains "cell" or "empty_droplet"
    ed <- unlist(nf_umi_ed_ct[, 3], use.names = FALSE)
    # Assume fourth column contains cell types
    ct <- unlist(nf_umi_ed_ct[, 4], use.names = FALSE)


    # Check values are reasonable
    if(any(c(max(nf)>1, min(nf)<0))){
      warning(paste0("The nuclear fraction values provided in the first column of 'nf_umi_ed_ct' should be between 0 and 1, but values outside this range were identified : minimum = ",min(nf),", maximum = ",max(nf)), call.=FALSE)
    }

    if(!all(umi == floor(umi))){
      non_integer_examples <- which(umi != floor(umi))
      if(length(non_integer_examples)>5){
        non_integer_examples <- non_integer_examples[1:5]
      }
      non_integer_examples <- paste(umi[non_integer_examples], collapse = ",")
      warning(paste0("Non-integer values detected in the second column of 'nf_umi_ed_ct' (e.g. ",non_integer_examples,") where umi counts were expected"), call.=FALSE)
    }

    if(max(umi)<100){
      warning(paste0("The total umi counts provided in the second column of 'nf_umi_ed_ct' appear to be quite low (max = ",max(umi),"), are these the total UMI counts per cell?"), call.=FALSE)
    }

    if(!all(unique(ed)%in%c("cell", "empty_droplet"))){
      ed_output <- unique(ed)
      if(length(ed_output)>5){
        ed_output <- ed_output[1:5]
        ed_output <- paste(ed_output, collapse = ",")
      }
      warning(paste0("The third column of 'nf_umi_ed_ct' was expected to contain either 'cell' or 'empty_droplet' but contains; ",ed_output), call.=FALSE)
    }

    if(verbose){
      ct <- unique(ct)
      ct <- paste(ct, collapse = ",")
      print(paste0("The following cell types were provided; ", ct))
    }

  } else {
    stop(paste0("A data frame should be supplied to the nf_umi_ed_ct argument, but an object of class ",paste(class(nf_umi), collapse = "/")," was provided"), call.=FALSE)
  }

  # Extract data for EM
  em.data <- data.frame(nf = unlist(nf_umi_ed_ct[,1], use.names = FALSE),
                        umi = log10(unlist(nf_umi_ed_ct[,2], use.names = FALSE)),
                        ct = unlist(nf_umi_ed_ct[,4], use.names = FALSE))
  row.names(em.data) <- 1:nrow(em.data)

  # Filter out any empty droplets
  em.data <- em.data[nf_umi_ed_ct[,3]=="cell",]

  # Split by cell type
  em.data.ct <-split(em.data, em.data$ct)

  # Run EM for all cell types
  if(verbose){ print("Fitting models with EM") }
  em_mods <- lapply(em.data.ct, function(x) mclust::Mclust(data = x[,1:2], G = 1:2, modelNames = "EEI", verbose = verbose))

  # Assign barcodes as "cell" or "damaged_cell" using the `assess_EM` function
  em_mods_assessed  <- lapply(em_mods,
                              assess_EM,
                              nf_thresh = nf_sep,
                              umi_thresh = umi_sep_perc)

  # Create plots if requested
  if(output_plots){
    if(verbose){ print("Creating requested plots") }
    em_plots <- lapply(seq_along(em_mods), function(x) plot_EM(em = em_mods[[x]],
                                                             em_c = em_mods_assessed[[x]],
                                                             input_cell_type = names(em_mods)[x],
                                                             umi_thresh = umi_sep_perc,
                                                             nf_thresh = nf_sep))
    names(em_plots) <- names(em_mods)
  }

  # Update the input data frame "nf_umi_ed_ct" (which still contains empty
  # droplets) any with damaged_cell info
  names(em_mods_assessed) <- NULL
  em_mods_assessed <- do.call(rbind, em_mods_assessed)
  nf_umi_ed_ct$cell_status[as.integer(row.names(em_mods_assessed))] <- em_mods_assessed$classification

  # If plots were not requested, return results as a data frame
  if(output_plots){
    return(list(df=nf_umi_ed_ct, plots=em_plots))
  } else {
    return(list(df=nf_umi_ed_ct, plots=NULL))
    }
  }
