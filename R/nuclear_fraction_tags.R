#'Parse cell barcode and region type tags in input region of provided BAM file
#'
#'@description This function accepts three inputs; a GRanges defining a genomic
#'  interval, a reference to a BAM file and a vector of cell barcodes. It then
#'  uses uses `Rsamtools::scanBam()` to import the provided alignment tags which
#'  contain the cell barcode and region type (exonic or intronic). The tags are
#'  then parsed and returned as a data frame with two integer columns; exon and
#'  intron.
#'
#'@param interval a single GRanges defining a genomic interval to import reads
#'  from
#'@param bam a reference to a BAM file to pull reads from (Rsamtools::BamFile())
#'@param bc a character vector of cell barcodes of interest
#'@param CB_tag character, the BAM tag containing the cell barcode sequence
#'@param RE_tag character, the BAM tag containing the region type
#'@param EXON_tag character, the character string that defines a read as exonic
#'@param INTRON_tag character, the character string that defines a read as
#'  intronic
#'
#'@return data.frame returns a data frame with two integer columns; "exon" - the
#'  number of reads that map to exons in the requested genomic interval for each
#'  cell barcode and "intron" - the number of reads that map to introns, where
#'  the rows are in the same order as the input character vector of cell
#'  barcodes. If there are no reads in the interval requested, the function will
#'  return NULL. This function was written as a helper function to be called
#'  internally by `dropletQC::nuclear_fraction_tags` and isn't intended for more
#'  general use.
#'
#'
#'@keywords internal
#'
parse_bam_tags <- function(interval,
                      bam,
                      bc,
                      CB_tag,
                      RE_tag,
                      EXON_tag,
                      INTRON_tag){

  # Import reads
  tags <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(tag=c(CB_tag, RE_tag), which=interval))[[1]]$tag

  # Account for the possibility that there may be no reads in the interval
  if(any(sapply(tags, is.null))){ return(NULL) }

  # Convert from list of length two (CB_tag and RE_tag) to a data frame
  tags <- as.data.frame(tags)

  # Remove cases where there is no cell barcode
  tags <- tags[stats::complete.cases(tags), ]

  # Remove any barcodes that aren't in bc
  tags <- tags[tags[[CB_tag]]%in%bc,]

  # Count
  tags <- data.frame(table(tags))

  # Account for the possibility that there may be no reads remaining
  if(nrow(tags)==0){ return(NULL) }

  # Create a data frame with two columns; exon and intron, with the same row
  # order as bc
  exon = data.frame(CB = tags$CB[tags[[RE_tag]]==EXON_tag], exon = tags$Freq[tags[[RE_tag]]==EXON_tag])
  intron = data.frame(CB = tags$CB[tags[[RE_tag]]==INTRON_tag], intron = tags$Freq[tags[[RE_tag]]==INTRON_tag])
  tags <- Reduce(function(x, y) merge(x, y, by="CB", all=TRUE),
                 list(data.frame(CB=bc), exon, intron))

  # Replace NA's with zero
  tags[is.na(tags)] <- 0

  return(tags)
}


#'Calculate the nuclear fraction statistic using BAM tags
#'
#'@description This function uses the region type tags in a provided BAM file to
#'  calculate for each input cell barcode the nuclear fraction statistic. This
#'  is just the fraction of reads that are intronic:
#'
#'  nuclear fraction = # intronic reads / (# intronic reads + # of exonic reads)
#'
#'  The row names of the returned data frame will match the order and name of
#'  the supplied barcodes. As a minimum you can provide as input a directory
#'  containing cellranger output (outs).
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
#'  argument for the function. If your directory structure doesn't match the one
#'  created by Cell Ranger you can provide the file paths directly using the
#'  bam, bam_index and barcodes arguments.
#'
#'@param bam character, the path to the input bam file. Not required if an
#'  'outs' directory is provided.
#'@param bam_index character, the path to the input bam file index. Not required
#'  if an 'outs' directory is provided.
#'@param barcodes character, either a vector of barcode names or the path to the
#'  barcodes.tsv.gz file output by Cell Ranger. If providing the cell barcodes
#'  as a vector, make sure that the format matches the one in the BAM file -
#'  e.g. be mindful if there are integers appended to the end of the barcode
#'  sequence. This argument isn't required if an 'outs' directory is provided -
#'  the function will just look for "barcodes.tsv.gz" in
#'  outs/filtered_feature_bc_matrix.
#'@param cores numeric, runs the function in parallel using furrr:future_map()
#'  with the requested number of cores. Setting `cores=1` will cause future_map
#'  to run sequentially.
#'@param tiles numeric, to speed up the processing of the BAM file we can split
#'  the genome up into tiles and process reads in chunks
#'@param cell_barcode_tag character, the BAM tag containing the cell barcode
#'  sequence
#'@param region_type_tag character, the BAM tag containing the region type
#'@param exon_tag character, the character string that defines a read as exonic
#'@param intron_tag character, the character string that defines a read as
#'  intronic
#'@param verbose logical, whether or not to print progress
#'
#'@return data.frame, the function returns a 1-column data frame containing the
#'  calculated nuclear fraction statistic for each input barcode. The order and
#'  names of the rows will match those of the input cell barcodes.
#'
#'@export
#'@examples
#' nf1 <- nuclear_fraction_tags(
#'     outs = system.file("extdata", "outs", package =
#'     "dropletQC"),
#'      tiles = 10, cores = 1, verbose = FALSE)
#' head(nf1)
#'
#' nf2 <- nuclear_fraction_tags(
#'    bam = system.file("extdata", "outs","possorted_genome_bam.bam", package =
#'    "dropletQC"),
#'    barcodes = c("AAAAGTCACTTACTTG-1",
#'                 "AAAAGTGGATCTCTAA-1",
#'                 "AAACACGTTCTCATCG-1"),
#'    tiles = 10, cores = 1,
#'    verbose = FALSE)
#' nf2
#'
nuclear_fraction_tags <- function(outs=NULL,
                             bam=NULL,
                             bam_index=paste0(bam,".bai"),
                             barcodes=NULL,
                             cores = future::availableCores() - 1,
                             tiles = 100,
                             cell_barcode_tag = "CB",
                             region_type_tag = "RE",
                             exon_tag = "E",
                             intron_tag = "N",
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
    bam_check <- check_bam(paste0(outs,"/possorted_genome_bam.bam"),
                           cb_tag=cell_barcode_tag,
                           rt_tag=region_type_tag)
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
    bam_check <- check_bam(bam,
                           cb_tag=cell_barcode_tag,
                           rt_tag=region_type_tag)
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
                     .f = function(i) parse_bam_tags(interval = genome_tiles[i],
                                                bam = bam_file,
                                                bc = barcodes,
                                                CB_tag = cell_barcode_tag,
                                                RE_tag = region_type_tag,
                                                EXON_tag = exon_tag,
                                                INTRON_tag = intron_tag),
                     .progress = verbose)

  # Close multisession workers
  future::plan(future::sequential)

  # Print progress
  time_passed <- round(as.numeric(difftime(time1 = Sys.time(), time2 = start_time, units = "mins")), 2)
  if(verbose){ message(paste0("Completed parsing BAM file. Time elapsed ~ ", time_passed, " mins")) }

  # Remove any NULL results
  tags <- tags[!sapply(tags, is.null)]

  # Combine results
  exon = as.integer(unlist((Reduce(`+`, lapply(tags, `[`, 'exon')))))
  intron = as.integer(unlist(Reduce(`+`, lapply(tags, `[`, 'intron'))))
  results <- data.frame(CB = barcodes, exon=exon, intron=intron)

  # Calculate the nuclear fraction (intronic reads / sum of intronic and exonic reads)
  nuclear_fraction <- results$intron/(results$intron + results$exon)

  # In some rare cases nuclear_fraction my be NaN (0/0+0). For simplicity, change these to zero.
  nuclear_fraction[is.na(nuclear_fraction)] <- 0

  # Return data frame, with cell barcodes as the row names
  nuclear_fraction <- data.frame(nuclear_fraction = nuclear_fraction)
  row.names(nuclear_fraction) <- results$CB

  return(nuclear_fraction)
}
