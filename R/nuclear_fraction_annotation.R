#'Check annotation and BAM file are compatible
#'
#'@description This function checks whether the annotation file and the BAM file
#'  provided to `nuclear_fraction_annotation()` are compatible. It pulls all of
#'  the unique values in the seqid column of the annotation and all of the
#'  seqname values from the BAM index. It then checks to see if all of the
#'  sequence names in the BAM file can be found in the annotation. It returns a
#'  character vector with two elements. The first element contains either
#'  "pass", "warning" or error ; depending on whether all, some or nonw of the
#'  sequences were found. The second element contains a message with more
#'  details. Note that the function was written as a simple helper function to
#'  be called `dropletQC::nuclear_fraction_annotation()` and isn't intended for
#'  more general use.
#'
#'@param bam_file character, should be a character vector pointing to the BAM
#'  file to be checked.
#'@param annotation character, annotation should be a character vector pointing
#'  to the annotation file to be checked.
#'
#'@return character. The function returns a character vector with two elements.
#'  The first element may take one of three possible values; "pass", "warning"
#'  or "error". The second element is a message with provides extra details e.g.
#'  if some of the sequences in the BAM file are missing from the annotation, a
#'  message is provided that contains a sample of those sequences.
#'
#'@keywords internal
check_bam_anno <- function(bam_file, annotation){

  # Check annotation exists
  if (!file.exists(annotation)) {
    return(c(
      "error",
      paste0("The provided annotation: ", annotation, " cannot be found.")
    ))
  }

  # Get sequence IDs from the annotation file
  anno_seqids <- rtracklayer::readGFF(filepath = annotation,
                                      filter=list(type=c("gene")),
                                      columns=c("seqid"),
                                      tags=character(0))
  anno_seqids <- unique(as.character(anno_seqids$seqid))

  # Get sequence IDs from the BAM file
  bam_seqids <- unique(as.character(Rsamtools::idxstatsBam(bam_file)$seqnames))

  # Check if any seqids from the bam file are not in the annotation
  missing_from_anno <- bam_seqids[!bam_seqids%in%anno_seqids]

  if(!any(bam_seqids%in%anno_seqids)){
    return(c("error",paste0("None of the seqids in the BAM file can be found in the provided annotation.")))
  }

  # Collapse vector of seqids for printing output
  if (length(missing_from_anno) >0 & length(missing_from_anno) < 12) {
    missing_from_anno <- paste(missing_from_anno, collapse = ", ")
  }

  if (length(missing_from_anno) >12) {
    missing_from_anno <-
      paste0(paste(utils::head(missing_from_anno), collapse = ','),
             "...",
             paste(utils::tail(missing_from_anno), collapse = ','))
  }

  # Return warning if some seqids are missing
  if(length(missing_from_anno)>0){
    return(c("warning",paste0("The following seqids; ",missing_from_anno," are present in the provided BAM file but cannot be found in the provided annotation.")))
  } else{
    return(c("pass","All seqids in the BAM file were found in the provided annotation."))
  }
}


#'Get transcript ranges
#'
#'@description This function accepts an annotation as input and returns a
#'  GrangesList containing; GRanges defining merged exons and introns. Note that
#'  the function was written as a simple helper function to be called
#'  `dropletQC::nuclear_fraction_annotation()` and isn't intended for more
#'  general use.
#'
#'@param input_annotation character, should be a character vector pointing to
#'  the annotation file
#'@param input_annotation_format character. Can be one of "auto", "gff3" or
#'  "gtf". This is passed to the 'format' argument of
#'  GenomicFeatures::makeTxDbFromGFF(). This should generally just be left as
#'  "auto"
#'@param ntiles integer, the number of blocks that exons and introns should be
#'  grouped into
#'@param verbose_output logical, whether to print messages
#'
#'@return GRangesList. Contains GRanges defining merged exons and introns
#'  grouped into ~ntiles blocks. When defining intronic regions, any regions
#'  that overlap exons are removed.  Strand is ignored when merging all
#'  intervals.
#'
#'@keywords internal
get_transcript_ranges <- function(input_annotation,
                                  input_annotation_format,
                                  ntiles,
                                  verbose_output){

  print("Extracting exon and intron ranges from provided annotation file:")
  print(input_annotation)
  txdb <- suppressWarnings(suppressMessages(
    GenomicFeatures::makeTxDbFromGFF(input_annotation,
                                     format = input_annotation_format)
  ))
  exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names=TRUE)
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)

  # Merge any overlapping exon or intron ranges - ignore strand
  exons <- GenomicRanges::reduce(unlist(exons), ignore.strand=TRUE)
  introns <- GenomicRanges::reduce(unlist(introns), ignore.strand=TRUE)

  # Remove from introns any regions that overlap exons
  introns <- BiocGenerics::setdiff(introns, exons)

  # Name exons and introns then combine them
  names(introns) <- rep("intron", length(introns))
  names(exons) <- rep("exon", length(exons))
  exons_introns <- c(exons, introns)
  exons_introns <- sort(exons_introns)

  # Reduce to ~ntiles ranges or less
  exons_introns <- S4Vectors::split(exons_introns,
                                    ceiling(seq_along(exons_introns) / round(length(exons_introns) /
                                                                               ntiles)))

  return(exons_introns)
}


#' Determine read overlap with exon and intron intervals
#'
#' @description A function that accepts four inputs;
#'
#'   A GRanges defining exon and intron ranges
#'
#'   A reference to a BAM file
#'
#'   A vector of cell barcodes
#'
#'   A character vector describing the BAM tag which contains the cell barcode
#'
#'   and returns an integer twice the length of the input cell barcode vector
#'   containing; c(intron counts, exon counts), in the same order as the input
#'   barcodes. Note that the function was written as a simple helper function to
#'   be called `dropletQC::nuclear_fraction_annotation()` and isn't intended for
#'   more general use.
#'
#' @param block GRanges, contains merged exon and intron intervals
#' @param bam_file_ref character, a reference to the BAM file to import reads
#'   from
#' @param cell_barcodes character, a vector of cell barcodes matching the format
#'   in bam_file_ref
#' @param cb_tag character, defines the BAM tags which contains the cell barcode
#'   e.g. "CB"
#'
#' @return integer. This function returns a vector of integers twice the length
#'   of the input cell barcode vector containing; c(intron counts, exon counts),
#'   in the same order as the input barcodes. If a cell barcode was not detected
#'   in the range, `NA_integer_` is returned.
#'
#' @keywords internal
intron_exon_overlap <- function(block,
                                bam_file_ref,
                                cell_barcodes,
                                cb_tag){

  # Import all reads that fall within the genomic interval containing the exon
  # and intron ranges
  reads <- GenomicAlignments::readGAlignments(
    file = bam_file_ref,
    param = Rsamtools::ScanBamParam(tag = cb_tag,
                                    which = range(block))
  )

  # Account for cases where interval contains no reads - return a vector of NA's
  if(length(reads)==0){ return(rep(NA, length(cell_barcodes)*2)) }

  # Find reads that overlap exons and introns
  all_hits <-
    GenomicAlignments::findOverlaps(
      query = block,
      subject = reads,
      minoverlap = 5
    )


  # Count cell barcodes and exon/intron hits
  region_hits <-
    table(cell_barcode=reads@elementMetadata[[cb_tag]][S4Vectors::to(all_hits)],
  region_type = names(block)[S4Vectors::from(all_hits)])
  region_hits <- data.frame(region_hits)

  # Get named vectors of integers
  # exons
  exon_hits <- region_hits$Freq[region_hits$region_type=="exon"]
  names(exon_hits) <- region_hits$cell_barcode[region_hits$region_type=="exon"]
  exon_hits <- exon_hits[cell_barcodes]
  # introns
  intron_hits <- region_hits$Freq[region_hits$region_type=="intron"]
  names(intron_hits) <- region_hits$cell_barcode[region_hits$region_type=="intron"]
  intron_hits <- intron_hits[cell_barcodes]

  return(c(intron_hits, exon_hits))
}


#'Calculate the nuclear fraction statistic by quantifying intron/exon overlap
#'
#'@description This function calculates the nuclear fraction score by parsing
#'  all of the reads in the provided BAM file and determining their overlap with
#'  intron and exon intervals obtained from the provided annotation file. It is
#'  a more flexible function than `nuclear_fraction_tags()` because it doesn't
#'  rely on the BAM file already containing tags describing whether reads are
#'  aligned to intronic or exonic regions, but it is also slower, because this
#'  information needs to be calculated.
#'
#'@param annotation_path character, should be a character vector pointing to the
#'  annotation file
#'@param annotation_format character. Can be one of "auto", "gff3" or "gtf".
#'  This is passed to the 'format' argument of
#'  GenomicFeatures::makeTxDbFromGFF(). THis should generally just be left as
#'  "auto"
#'@param bam character, should be a character vector pointing to the BAM file
#'@param bam_index character, the path to the input bam file index
#'@param barcodes character, either a vector of barcode names or the path to a
#'  file barcodes.tsv.gz containging the cell barcodes. If providing the cell
#'  barcodes as a vector, make sure that the format matches the one in the BAM
#'  file - e.g. be mindful if there are integers appended to the end of the
#'  barcode sequence.
#'@param cell_barcode_tag character, defines the BAM tage which contains the
#'  cell barcode e.g. "CB"
#'@param cores numeric, parsing of the BAM file can be run in parallel using
#'  furrr:future_map() with the requested number of cores. Setting `cores=1`
#'  will cause future_map to run sequentially.
#' @param tiles integer, to speed up the processing of the BAM file we can split
#'  transcripts up into tiles and process reads in chunks, default=1000.
#'@param verbose logical, whether or not to print progress
#'
#'@return data.frame. Returns a data frame containing the nuclear fraction
#'  score. THis is just the fraction of reads that are intronic:
#'
#'  nuclear fraction = # intronic reads / (# intronic reads + # of exonic reads)
#'
#'  The row names of the returned data frame will match the order and name of
#'  the supplied barcodes.
#'
#'@export
#'
#' @examples
#' nf3 <- nuclear_fraction_annotation(
#'  annotation_path = system.file("extdata/outs/chr1.gff3",
#'   package = "dropletQC"),
#'  bam = system.file("extdata/outs/possorted_genome_bam.bam",
#'   package = "dropletQC"),
#'  barcodes = system.file(
#'  "extdata/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#'   package = "dropletQC"),
#'  tiles = 1, cores = 1, verbose = FALSE)
#' head(nf3)
#'
nuclear_fraction_annotation <- function(
  annotation_path,
  annotation_format = "auto",
  bam,
  bam_index=paste0(bam,".bai"),
  barcodes,
  cell_barcode_tag ="CB",
  cores=future::availableCores() - 1,
  tiles=1000,
  verbose=TRUE){

  ## Check "simple" arguments are valid
  if(!all(is.character(annotation_path), length(annotation_path)==1)) { stop("`annotation_path` argument should be a single character string specifying the path to the annotation file", call.=FALSE) }
  if(!(annotation_format%in%c("auto","gff3","gtf"))) { stop("`annotation_format` should be one of 'auto', 'gff3' or 'gtf'", call.=FALSE) }
  if(!all(is.character(cell_barcode_tag), length(cell_barcode_tag)==1)) { stop("`cell_barcode_tag` argument should be a single character string specifying the BAM file tag which contains the cellular barcode sequence", call.=FALSE) }
  if(!all(is.numeric(cores), cores>0)) { stop("`cores` argument should be a positive numeric scalar specifying the maximum number of parallel futures used for parsing the BAM file (passed to furrr:future_map())", call.=FALSE) }
  if(!all(is.numeric(tiles), tiles>0)) { stop("`tiles` argument should be a positive numeric scalar specifying the number of merged gene intervals to use for parsing the BAM file", call.=FALSE) }
  if(!is.logical(verbose)) { stop("`verbose` argument should be either TRUE or FALSE, to determine whether or not to print progress messages", call.=FALSE) }


  ## Get barcodes
  if (!is.character(barcodes)) {
    stop(
      "`barcodes` should be either a path to file e.g. barcodes.tsv.gz with cell barcodes in the first column (no header - gets passed to `utils::read.table(barcodes)[, 1]`), or a character vector defining cell barcodes",
      call. = FALSE
    )
  }

  # If barcodes is a path to a file (length==1), read in
  if (length(barcodes) == 1) {
    if (!file.exists(barcodes)) {
      stop(
        paste0(
          "The provided barcodes file: '",
          barcodes,
          "', does not appear to exist"
        ),
        call. = FALSE
      )
    }
    barcodes <- utils::read.table(barcodes)[, 1]
  }

  ## Check BAM
  if(!file.exists(bam_index)){ stop(paste0("The supplied bam index: '", bam_index, "' does not appear to exist"), call.=FALSE) }
  bam_check <- check_bam(bam, cb_tag=cell_barcode_tag)
  if(bam_check[1]=="pass"){
    if(verbose){ message(bam_check[2]) }
  } else {
    if(bam_check[1]=="warning"){ warning(bam_check[2], call.=FALSE) }
    if(bam_check[1]=="error"){ stop(bam_check[2], call.=FALSE) }
  }

  ## Check the annotation is compatible with the BAM file
  anno_check <- check_bam_anno(bam_file = bam,
                               annotation = annotation_path)
  if(anno_check[1]=="pass"){
    if(verbose){ message(anno_check[2]) }
  } else {
    if(anno_check[1]=="warning"){ warning(anno_check[2], call.=FALSE) }
    if(anno_check[1]=="error"){ stop(anno_check[2], call.=FALSE) }
  }

  ## Get transcript, exon and intron ranges
  tx_ranges <- get_transcript_ranges(input_annotation = annotation_path,
                                     input_annotation_format = annotation_format,
                                     ntiles=tiles,
                                     verbose_output=verbose)

  # Remove ranges from sequences not in the BAM file
  bam_info <- Rsamtools::idxstatsBam(bam)
  bam_seqs <- bam_info$seqnames
  tx_seqs <- GenomicRanges::seqnames(GenomicRanges::seqinfo(tx_ranges))
  tx_ranges <- GenomeInfoDb::dropSeqlevels(x = tx_ranges,
                                   value = tx_seqs[!tx_seqs%in%bam_seqs],
                                   pruning.mode = "tidy")


  ## Loop through transcript blocks to get intron and exon counts
  # Print progress
  if (verbose) {
    message(
      paste0(
        "Processing BAM file containing ",
        round(sum(bam_info$mapped) / 1E6, 2),
        " million reads, using ",
        cores,
        " cores:"
      )
    )
    start_time <- Sys.time()
  }

  # Process the BAM file with future_map
  future::plan(future::multisession, workers = cores)
  ie <-
    furrr::future_map(
      .x = seq_along(tx_ranges),
      .f = function(i)
        intron_exon_overlap(
          block = tx_ranges[[i]],
          bam_file_ref = Rsamtools::BamFile(file = bam, index = bam_index),
          cell_barcodes = barcodes,
          cb_tag = cell_barcode_tag
        ),
      .progress = TRUE
    )

  # Close multisession workers
  future::plan(future::sequential)

  # Print progress
  if (verbose) {
    time_passed <-
      round(as.numeric(difftime(
        time1 = Sys.time(),
        time2 = start_time,
        units = "mins"
      )), 2)
    message(paste0(
      "Completed parsing BAM file. Time elapsed ~ ",
      time_passed,
      " mins"
    ))
  }

  # Combine exon and intron counts for each cell barcode into a data frame
  ie <- colSums(do.call(rbind, ie), na.rm=TRUE)
  ie <- data.frame(I=ie[1:length(barcodes)],
                   E=ie[(length(barcodes)+1):length(ie)])

  # Calculate nuclear fraction
  nf <- data.frame(nuclear_fraction = ie$I/(ie$E + ie$I))
  row.names(nf) <- barcodes

  return(nf)
}
