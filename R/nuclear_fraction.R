#'Calculate the nuclear fraction statistic
#'
#'@description This function uses the RE tags in the Cell Ranger Barcoded BAM
#'file to calculate for each input cell barcode nuclear the fraction statistic.
#'This is just the fraction of reads that are intronic:
#'
#'nuclear fraction = # intronic reads / (# intronic reads + # of exonic reads)
#'
#'The row names of the returned data frame will match the order and name of the
#'supplied barcodes. As a minimum you can provide as input a directory
#'containing cellranger the output (outs).
#'
#'@param outs character, the path to the 'outs' directory created by Cell
#'  Ranger. We assume outs is structured this way:
#'
#'  ├── filtered_feature_bc_matrix
#'
#'  │   ├── barcodes.tsv.gz
#'
#'  │   ├── features.tsv.gz
#'
#'  │   └── matrix.mtx.gz
#'
#'  ├── possorted_genome_bam.bam
#'
#'  ├── possorted_genome_bam.bam.bai
#'
#'  ├── raw_feature_bc_matrix
#'
#'  │   ├── barcodes.tsv.gz
#'
#'  │   ├── features.tsv.gz
#'
#'  │   └── matrix.mtx.gz
#'
#'  Note that there will probably be other files in the directory as well. We
#'  don't need to worry about those, as the only three files that the function
#'  will require are; possorted_genome_bam.bam, possorted_genome_bam.bam.bai and
#'  filtered_feature_bc_matrix/barcodes.tsv.gz. This is the only required
#'  argument for the function. If your directory structure no longer matches the
#'  one created by Cell Ranger (e.g. you were given the files from a
#'  collaborator or sequencing facility) you can provide the file paths directly
#'  using the bam, bam_index and barcodes arguments.
#'
#'@param bam character, the path to the input bam file. Not required if the
#'  'outs' directory is provided.
#'@param bam_index character, the path to the input bam file index. Not required
#'  if the 'outs' directory is provided.
#'@param barcodes character, the path to the barcodes.tsv.gz file output by Cell
#'  Ranger. Rather than provide the path to a file on disk you can alternatively
#'  provide a vector of barcode names. If providing the cell barocodes as a
#'  vector, make sure that the format matches the one in the BAM file - be
#'  mindful of the "-1" at the end of the barcode sequence. This argument isn't
#'  required if the 'outs' directory is provided - the function will just look
#'  for "barcodes.tsv.gz" in outs/filtered_feature_bc_matrix
#'@param cores numeric, runs the function in parallel using furrr:future_map()
#'  with the requested number of cores. Setting `cores=1` will cause future_map
#'  to run sequentially.
#'@param tiles numeric, to speed up the processing of the BAM file we can split
#'  the genome up into tiles and process reads in chunks
#'@param verbose logical, whether or not to print progress
#'
#'@return data.frame, the function returns a 1-column data frame containing the
#'  calculated nuclear fraction statistic for each input barcode. The order and
#'  names of the rows will match those of the input cell barcodes.
#'
#'@export
#'
#' @examples
#' nf <- nuclear_fraction(outs = system.file("extdata", "outs", package =
#' "dropletQC"), tiles = 10, cores = 1, verbose = FALSE)
#' head(nf)
#'
nuclear_fraction <- function(outs="outs",
                             bam=NULL,
                             bam_index=paste0(bam,".bai"),
                             barcodes=NULL,
                             cores = future::availableCores() - 1,
                             tiles = 100,
                             verbose = TRUE){

  # Check `cores`, `tiles` and `verbose` arguments are valid
  if(!all(is.numeric(cores), cores>0)) { stop("`cores` argument should be a positive numeric scalar specifying the maximum number of parallel futures used for parsing the BAM file (passed to furrr:future_map())", call.=FALSE) }
  if(!all(is.numeric(tiles), tiles>0)) { stop("`tiles` argument should be a positive numeric scalar specifying the number of pieces to split the BAM file into for parsing", call.=FALSE) }
  if(!is.logical(verbose)) { stop("`verbose` argument should be either TRUE or FALSE, to determine whether or not to print progress messages", call.=FALSE) }

  # Check remaining arguments
  if(!is.null(outs)){

    # If `outs` provided, check the directory exists
    if(!is.character(outs)) { stop("`outs` argument should be a character string specifying the path to the cellranger directory e.g. '/path/to/outs'", call.=FALSE) }
    if(!dir.exists(outs)) { stop(paste0("the cellranger `outs` directory provided: ", outs,", does not appear to exist"), call.=FALSE) }

    # Check barcodes file exists - then read in
    if(!dir.exists(paste0(outs,"/filtered_feature_bc_matrix"))) { stop(paste0("Cellranger outs directory provided, but the required internal directory: '", paste0(outs,"/filtered_feature_bc_matrix"),"', does not appear to exist"), call.=FALSE) }
    if(!file.exists(paste0(outs,"/filtered_feature_bc_matrix/barcodes.tsv.gz"))) { stop(paste0("Cellranger outs directory provided, but the file: '", paste0(outs,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),"', does not appear to exist"), call.=FALSE) }
    barcodes <- utils::read.table(paste0(outs, "/filtered_feature_bc_matrix/barcodes.tsv.gz"))[,1]

    # Check and create reference to BAM file
    if(!file.exists(paste0(outs,"/possorted_genome_bam.bam.bai"))) { stop(paste0("Cellranger outs directory provided, but the BAM index: '", paste0(outs,"/possorted_genome_bam.bam.bai"),"', does not appear to exist"), call.=FALSE) }
    bam_check <- check_bam(paste0(outs,"/possorted_genome_bam.bam"))
    if(bam_check[1]=="pass"){
      bam_file <- Rsamtools::BamFile(file = paste0(outs,"/possorted_genome_bam.bam"), index = paste0(outs,"/possorted_genome_bam.bam.bai"))
      if(verbose){ message(bam_check[2]) }
    } else {
      if(bam_check[1]=="warning"){ warning(bam_check[2], call.=FALSE) }
      if(bam_check[1]=="error"){ stop(bam_check[2], call.=FALSE) }
    }

  } else {

    # If `outs` argument is NULL, use the provided three arguments; bam, bam_index and barcodes
    if(!all(!sapply(list(bam, bam_index, barcodes), is.null))){
      missing_args <- paste(c("bam", "bam_index", "barcodes")[sapply(list(bam, bam_index, barcodes), is.null)], collapse = ", ")
      stop(paste0("If `outs` argument not is not provided, all three arguments `bam`, `bam_index` and `barcodes` must be provided. These arguments were still NULL; ", missing_args), call.=FALSE)
    }

    # Check BAM
    if(!file.exists(bam_index)){ stop(paste0("The supplied bam index: '", bam_index, "' does not appear to exist"), call.=FALSE) }
    bam_check <- check_bam(bam)
    if(bam_check[1]=="pass"){
      bam_file <- Rsamtools::BamFile(file = bam, index = bam_index)
      if(verbose){ message(bam_check[2]) }
    } else {
      if(bam_check[1]=="warning"){ warning(bam_check[2], call.=FALSE) }
      if(bam_check[1]=="error"){ stop(bam_check[2], call.=FALSE) }
    }

    # Check barcodes argument
    if(!is.character(barcodes)){ stop("`barcodes` should be either a path to file e.g. barcodes.tsv.gz with cell barcodes in the first column - as produced by cellranger - or a character vector of cell barcodes", call.=FALSE) }
    # If barcodes is a path to a file (length==1), read in
    if(length(barcodes)==1){
      if(!file.exists(barcodes)) { stop(paste0("The provided barcodes file: '", barcodes,"', does not appear to exist"),call.=FALSE) }
      barcodes <- utils::read.table(barcodes)[,1]
    }
  }

  # Print the start of the barcodes brovided
  if(verbose){ message(paste0(length(barcodes)," cell barcodes were provided as input. The first five are; ", paste(barcodes[1:5], collapse = ", "))) }

  # Create tiles across the genome for processing the BAM file in chunks
  bam_info <- Rsamtools::idxstatsBam(bam_file)
  bam_info <- bam_info[bam_info$mapped>0,] # remove seqs with no reads mapped to them
  genome_tiles <- bam_info$seqlength
  names(genome_tiles) <- bam_info$seqnames
  genome_tiles <- unlist(GenomicRanges::tileGenome(genome_tiles, ntile = tiles))

  # Print progress
  if(verbose){ message(paste0("Processing BAM file containing ",round(sum(bam_info$mapped)/1E6,2), " million reads, using ", cores, " cores:")) }

  # Process the BAM file with future_map
  if(cores>1){
    future::plan(future::multisession, workers = cores)
  } else {
    future::plan(future::sequential)
  }

  start_time <- Sys.time()
  tags <- furrr::future_map(.x = seq_along(genome_tiles),
                     .f = function(i) parse_bam(interval = genome_tiles[i], bam = bam_file, bc = barcodes),
                     .progress = verbose)

  # Close multisession workers
  future::plan(future::sequential)

  # Print progress
  time_passed <- round(as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "mins")), 2)
  if(verbose){ message(paste0("Completed parsing BAM file. Time elapsed ~ ", time_passed, " mins")) }

  # Remove any NULL results
  tags <- tags[!sapply(tags, is.null)]

  # Combine results
  E = as.integer(unlist((Reduce(`+`, lapply(tags, `[`, 'E')))))
  I = as.integer(unlist(Reduce(`+`, lapply(tags, `[`, 'I'))))
  N = as.integer(unlist(Reduce(`+`, lapply(tags, `[`, 'N'))))
  results <- data.frame(CB = barcodes, E=E, I=I, N=N)

  # Calculate the nuclear fraction (intronic reads / sum of intronic and exonic reads)
  nuclear_fraction <- results$N/(results$N + results$E)

  # In some rare cases nuclear_fraction my be NaN (0/0+0). For simplicity, change these to zero.
  nuclear_fraction[is.na(nuclear_fraction)] <- 0

  # Return data frame, with cell barcodes as the row names
  nuclear_fraction <- data.frame(nuclear_fraction = nuclear_fraction)
  row.names(nuclear_fraction) <- results$CB

  return(nuclear_fraction)
}
