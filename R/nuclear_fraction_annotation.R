#' Check annotation and BAM file are compatible
#'
#'@description This function checks whether the annotation file and the BAM file
#'  provided to `nuclear_fraction_annotation()` are compatible. It pulls all of
#'  the unique values in the seqid column of the annotation and all of the
#'  seqname values from the BAM index. It then checks to see if all of the
#'  sequence names in the BAM file can be found in the annotation. It returns a
#'  character vector with two elements. The first element contains either "pass"
#'  or "warning" ; depending whether all of the sequences were found or not. The
#'  second element contains a message with more details. Note that the function
#'  was written as a simple helper function to be called
#'  `dropletQC::nuclear_fraction_annotation()` and isn't intended for more
#'  general use.
#'
#' @param bam_file character, should be a character vector pointing to
#'   the BAM file to be checked.
#' @param annotation character, annotation should be a character vector pointing
#'   to the annotation file to be checked.
#'
#' @return character. The function returns a character vector with two elements.
#'   The first element may take one of two possible values; "pass" or "warning".
#'   The second element is a message with provides extra details e.g. if some of
#'   the sequences in the BAM file are missing from the annotation, a message is
#'   provided that contains a sample of those sequences.
#'
#' @keywords internal
check_bam_anno <- function(bam_file, annotation){

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
  if (length(missing_from_anno) < 12) {
    missing_from_anno <- paste(missing_from_anno, collapse = ", ")
  } else{
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


#' Get transcript ranges
#'
#' @description This function accepts an annotation as input and returns three
#'   sets of genomic intervals;
#'
#'   transcript_blocks - GRanges that groups any overlapping genes into a single
#'   GRanges
#'
#'   exons - GRangesList containing merged exons grouped by transcript_blocks
#'
#'   introns - GRangesList containing merged introns grouped by
#'   transcript_blocks
#'
#' Note that the function was written as a simple helper function to be called
#' `dropletQC::nuclear_fraction_annotation()` and isn't intended for more
#' general use.
#'
#' @param input_annotation character, should be a character vector pointing to
#'   the annotation file
#' @param input_annotation_format character. Can be one of "auto", "gff3" or
#'   "gtf". This is passed to the 'format' argument of
#'   GenomicFeatures::makeTxDbFromGFF(). THis should generally just be left as
#'   "auto"
#'
#' @return list. A list is returned containing three elements;
#'
#'   transcript_blocks - GRanges defining intervals that contain any overlapping
#'   genes
#'
#'   exons - GRangesList containing exons grouped by transcript_blocks.
#'   Overlapping exons are merged into a single GRanges. Strand is ignored.
#'
#'   introns - GRangesList containing introns grouped by transcript_blocks.
#'   Overlapping introns are merged into a single GRanges. All intronic regions
#'   that overlap exons are removed.  Strand is ignored.
#'
#' @keywords internal
get_transcript_ranges <- function(input_annotation, input_annotation_format){

  print("Extracting exon and intron ranges from provided annotation file:")
  print(input_annotation)
  txdb <- suppressWarnings(suppressMessages(
    GenomicFeatures::makeTxDbFromGFF(input_annotation,
                                     format = input_annotation_format)
  ))
  exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names=TRUE)
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)

  # Get transcript blocks: GRanges that stretch from the most 5' to most 3'
  # coordinate of overlapping transcripts (ignores strand)
  tx_blocks <- GenomicRanges::reduce(c(unlist(range(exons)),
                                       unlist(range(introns))),
                                     ignore.strand=TRUE)

  # Merge any overlapping exon or intron ranges - ignore strand
  exons <- GenomicRanges::reduce(unlist(exons), ignore.strand=TRUE)
  introns <- GenomicRanges::reduce(unlist(introns), ignore.strand=TRUE)

  # Remove from introns any regions that overlap exons
  introns <- BiocGenerics::setdiff(introns, exons)

  # Group introns and exons into transcript blocks - creates a GRangesList,
  # where the names of each element are integers corresponding to the tx_blocks
  #the introns/exons fall into

  exons <- split(exons,
                 as.factor(S4Vectors::to(
                   GenomicAlignments::findOverlaps(query = exons,
                                                   subject = tx_blocks)
                 )))
  introns <- split(introns,
                   as.factor(S4Vectors::to(
                     GenomicAlignments::findOverlaps(query = introns,
                                                     subject = tx_blocks)
                   )))

  return(list(
    transcript_blocks = tx_blocks,
    exons = exons,
    introns = introns
  ))
}


#' Determine read overlap with exon and intron intervals
#'
#' @description A function that accepts seven inputs;
#'
#' GRanges transcript blocks
#'
#' An integer defining which transcript block to process
#'
#' A GrangesList containing exons in the input transcript block
#'
#' A GrangesList containing introns in the input transcript block
#'
#' A reference to a BAM file
#'
#' A vector of cell barcodes
#'
#' A character vector describing the BAM tag which contains the cell barcode
#'
#' and returns an integer twice the length of the input cell barcode vector
#' containing; c(intron counts, exon counts), in the same order as the input
#' barcodes. Note that the function was written as a simple helper function to
#' be called `dropletQC::nuclear_fraction_annotation()` and isn't intended for
#' more general use.
#'
#' @param transcript_block Granges, transcript_blocks contains GRanges grouping
#'   overlapping genes
#' @param current_block integer, defining which transcript block to process
#' @param exon_ranges GRangesList, contains merged exon intervals grouped by
#'   transcript_blocks
#' @param intron_ranges GRangesList, contains merged intron intervals grouped by
#'   transcript_blocks
#' @param bam_file_ref character, a reference to the BAM file to import reads
#'   from
#' @param cell_barcodes character, a vector of cell barcodes matching the format
#'   in bam_file_ref
#' @param cb_tag character, defines the BAM tage which contains the cell barcode
#'   e.g. "CB"
#'
#' @return integer. This function returns a vector of integers twice the length
#'   of the input cell barcode vector containing; c(intron counts, exon counts),
#'   in the same order as the input barcodes. If a cell barcode was not detected
#'   in the range, NA_integer_ is returned.
#'
#' @keywords internal
intron_exon_overlap <- function(transcript_block,
                                current_block,
                                exon_ranges,
                                intron_ranges,
                                bam_file_ref,
                                cell_barcodes,
                                cb_tag="CB"){

  # Import all reads that fall within the transcript block as GAlignmentPairs
  # objects
  reads <- GenomicAlignments::readGAlignments(
    file = bam_file_ref,
    param = Rsamtools::ScanBamParam(tag = cb_tag,
                                    which = transcript_block[current_block])
  )

  # Account for cases where interval contains no reads - return a vector of NA's
  if(length(reads)==0){ return(rep(NA, length(cell_barcodes)*2)) }

  # Find reads that overlap introns within the input transcript block

  if(!(as.character(current_block) %in% names(intron_ranges))){
    # Account for situation where there are no introns in the transcript block
    intron_hits <- rep(NA, length(cell_barcodes))
  } else{
    intron_hits <-
      GenomicAlignments::findOverlaps(
        query = unlist(intron_ranges[as.character(current_block)]),
        subject = reads,
        minoverlap = 5
      )

    # Count cell barcodes from hits - will remove any NAs for us
    intron_hits <-
      table(reads@elementMetadata[[cb_tag]][S4Vectors::to(intron_hits)])

    # Change from table to a named vector of integers
    intron_hits <- rep(intron_hits)

    # Remove any cell barcodes not in our list. Also changes vector to the same
    # length as our list of barcodes, inserting NA for missing barcodes
    intron_hits <- intron_hits[cell_barcodes]
  }


  # Repeat for exons
  if(!(as.character(current_block) %in% names(exon_ranges))){
    exon_hits <- rep(NA, length(cell_barcodes))
  } else{
    exon_hits <-
      GenomicAlignments::findOverlaps(
        query = unlist(exon_ranges[as.character(current_block)]),
        subject = reads,
        minoverlap = 5
      )
    exon_hits <- table(reads@elementMetadata[[cb_tag]][S4Vectors::to(exon_hits)])
    exon_hits <- rep(exon_hits)
    exon_hits <- exon_hits[cell_barcodes]
  }

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
#'@param cell_barcode_tag character, defines the BAM tage which contains the cell barcode
#'   e.g. "CB"
#'@param cores numeric, parsing of the BAM file can be run in parallel using
#'  furrr:future_map() with the requested number of cores. Setting `cores=1`
#'  will cause future_map to run sequentially.
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
#' #gbm_nf <- nuclear_fraction_annotation(
#' #  annotation_path = "data/human.gtf",
#' #  annotation_format = "auto",
#' #  bam = "data/GBM/outs/possorted_genome_bam.bam",
#' #  bam_index =  "data/GBM/outs/possorted_genome_bam.bam.bai",
#' #  barcodes = "data/GBM/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#' #  cell_barcode_tag = "CB",
#' #  cores = 16,
#' #  verbose = TRUE)
#'
nuclear_fraction_annotation <- function(
  annotation_path,
  annotation_format = "auto",
  bam,
  bam_index=paste0(bam,".bai"),
  barcodes,
  cell_barcode_tag ="CB",
  cores=future::availableCores() - 1,
  verbose=TRUE){

  ## Get barcodes
  if (!is.character(barcodes)) {
    stop(
      "`barcodes` should be either a path to file e.g. barcodes.tsv.gz with cell barcodes in the first column, or a character vector of cell barcodes",
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

  ## Check the annotation is compatible with the BAM file
  anno_check <- check_bam_anno(bam_file = bam,
                               annotation = annotation_path)
  if(anno_check[1]=="pass"){
    if(verbose){ message(anno_check[2]) }
  } else {
    if(anno_check[1]=="warning"){ warning(anno_check[2], call.=FALSE) }
  }


  ## Get transcript, exon and intron ranges
  tx_ranges <- get_transcript_ranges(input_annotation = annotation_path,
                                     input_annotation_format = annotation_format)


  ## Loop through transcript blocks to get intron and exon counts
  # Print progress
  if (verbose) {
    message(
      paste0(
        "Processing BAM file containing ",
        round(sum(Rsamtools::idxstatsBam(bam)$mapped) / 1E6, 2),
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
      .x = seq_along(tx_ranges[["transcript_blocks"]]),
      .f = function(i)
        intron_exon_overlap(
          transcript_block = tx_ranges[["transcript_blocks"]],
          current_block = i,
          exon_ranges = tx_ranges[["exons"]],
          intron_ranges = tx_ranges[["introns"]],
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
