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
#'  internally by `dropletQC::nuclear_fraction` and isn't intended for more
#'  general use.
#'
#'
#'@keywords internal
#'
parse_bam <- function(interval,
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
