#' Parse CB and RE tags in input region of barcoded BAM file
#'
#' @description
#' This function accepts three inputs; a GRanges defining a genomic interval, a reference to a Cell Ranger Barcoded BAM file and a vector of cell barcodes.
#' It then uses uses `Rsamtools::scanBam()` to import the Cell Ranger Barcoded BAM Alignment Tags:
#' CB - error-corrected cell barcode
#' RE - region type (E = exonic, N = intronic, I = intergenic)
#' parse them and return a data frame with three integer columns; E, I and N.
#'
#' @param interval
#' a single GRanges defining a genomic interval to import reads from
#' @param bam
#' a reference to a Cell Ranger Barcoded BAM file to pull reads from (Rsamtools::BamFile())
#' @param bc
#'a character vector of cell barcodes of interest
#'
#' @return data.frame
#' returns a data frame with three integer columns;
#' E - the number of reads that map to exons in the requested genomic interval for each cell barcode
#' I - the number of reads that map to intergenic regions
#' N - the number of reads that map to introns
#' where the rows are in the same order as the input character vector of cell barcodes. If there are no reads in the interval requested, the function will return NULL. This function was written as a helper function to be called internally by `dropletQC::nuclear_fraction` and isn't intended for more general use.
#'
#' @export
#'
#' @examples
#' #parse_bam(interval = genome_tiles[1], bam = bam_file, bc = barcodes)
parse_bam <- function(interval, bam, bc){

  # Import reads
  tags <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(tag=c("CB", "RE"), which=interval))[[1]]$tag

  # Account for the possibility that there may be no reads in the interval
  if(any(sapply(tags, is.null))){ return(NULL) }

  # Convert from list of length two (CB and RE) to a data frame
  tags <- as.data.frame(tags)

  # Remove cases where there is no cell barcode
  tags <- tags[stats::complete.cases(tags), ]

  # Remove any barcodes that aren't in bc
  tags <- tags[tags$CB%in%bc,]

  # Count
  tags <- data.frame(table(tags))

  # Account for the possibility that there may be no reads remaining
  if(nrow(tags)==0){ return(NULL) }

  # Create a data frame with three columns; E, I & N, with the same row order as bc
  E = data.frame(CB = tags$CB[tags$RE=="E"], E = tags$Freq[tags$RE=="E"])
  I = data.frame(CB = tags$CB[tags$RE=="I"], I = tags$Freq[tags$RE=="I"])
  N = data.frame(CB = tags$CB[tags$RE=="N"], N = tags$Freq[tags$RE=="N"])
  tags <- Reduce(function(x, y) merge(x, y, by="CB", all=TRUE),
                 list(data.frame(CB=bc), E, I, N))

  # Replace NA's with zero
  tags[is.na(tags)] <- 0

  return(tags)
}
