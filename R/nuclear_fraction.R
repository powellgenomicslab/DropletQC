# Define function `nuclear_fraction()` that calls IEN internally and returns a data frame with one column
# containing the nuclear fraction statistic. The row names of the data frame will match the supplied barcodes.
# As a minimume we can provide as input a directory containing cellranger the output (outs).
# We assume outs is structured this way:
#├── filtered_feature_bc_matrix
#│   ├── barcodes.tsv.gz
#│   ├── features.tsv.gz
#│   └── matrix.mtx.gz
#├── possorted_genome_bam.bam
#├── possorted_genome_bam.bam.bai
#├── raw_feature_bc_matrix
#│   ├── barcodes.tsv.gz
#│   ├── features.tsv.gz
#│   └── matrix.mtx.gz
# Only three files; possorted_genome_bam.bam, possorted_genome_bam.bam.bai &
# filtered_feature_bc_matrix/barcodes.tsv.gz are required. All files must
# exist if the `outs` argument is given and should be named as above.

# If your directory structure no longer matches the one created by
# CellRanger (e.g. you were given the files from a collaborator or
#sequencing facility) you can provide the file paths directly using the
# bam, bam_index and barcodes arguments.

# For the barcodes argument, rather than provide a path to a file on disk
# you can alternatively provide character vector of barocodes. The barcodes must
# match the format of the barcodes in the BAM file.

# To speed up the processing of the BAM file we break the genome up into `tiles` and process
# them in parallel with furrr:future_map() using `cores` cores. Setting `cores=1` will cause
# future_map to run sequentially.


nuclear_fraction <- function(outs=NULL,
                             bam=NULL,
                             bam_index=paste0(bam,".bai"),
                             barcodes=NULL,
                             cores = availableCores() - 1,
                             tiles = 1000,
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
    barcodes <- read.table(paste0(outs, "/filtered_feature_bc_matrix/barcodes.tsv.gz"))[,1]

    # Check and create reference to BAM file
    if(!file.exists(paste0(outs,"/possorted_genome_bam.bam.bai"))) { stop(paste0("Cellranger outs directory provided, but the BAM index: '", paste0(outs,"/possorted_genome_bam.bam.bai"),"', does not appear to exist"), call.=FALSE) }
    bam_check <- check_bam(paste0(outs,"/possorted_genome_bam.bam"))
    if(bam_check[1]=="pass"){
      bam_file <- BamFile(file = paste0(outs,"/possorted_genome_bam.bam"), index = paste0(outs,"/possorted_genome_bam.bam.bai"))
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
      bam_file <- BamFile(file = bam, index = bam_index)
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
      barcodes <- read.table(barcodes)[,1]
    }
  }

  # Print the start of the barcodes brovided
  if(verbose){ message(paste0(length(barcodes)," cell barcodes were provided as input. The first five are; ", paste(barcodes[1:5], collapse = ", "))) }

  # Create tiles across the genome for processing the BAM file in chunks
  bam_info <- idxstatsBam(bam_file)
  bam_info <- bam_info[bam_info$mapped>0,] # remove seqs with no reads mapped to them
  genome_tiles <- bam_info$seqlength
  names(genome_tiles) <- bam_info$seqnames
  genome_tiles <- unlist(tileGenome(genome_tiles, ntile = tiles))

  # Print progress
  if(verbose){ message(paste0("Processing BAM file containing ",round(sum(bam_info$mapped)/1E6,2), " million reads, using ", cores, " cores:")) }

  # Process the BAM file with future_map
  if(cores>1){
    plan(multisession, workers = cores)
  } else {
    plan(sequential)
  }

  start_time <- Sys.time()
  tags <- future_map(.x = seq_along(genome_tiles),
                     .f = function(i) parse_bam(interval = genome_tiles[i], bam = bam_file, bc = barcodes),
                     .progress = verbose)

  # Close multisession workers
  plan(sequential)

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
