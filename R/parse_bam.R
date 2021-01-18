# Define a function that takes as input:
# a single GRanges defining a genomic interval to import reads from
# a reference to the BAM file to pull reads from (Rsamtools::BamFile())
# a character vector of cell barcodes of interest
# and uses the 10x BAM Alignment Tags:
# CB - error-corrected cell barcode
# RE - region type (E = exonic, N = intronic, I = intergenic)
# to return a tibble with three integer columns;
# E - the number of reads that map to exons in the requested genomic interval for each cell barcode
# I - the number of reads that map to intergenic regions
# N - the number of reads that map to introns
# where the rows are in the same order as the input character vector of cell barcodes
parse_bam <- function(interval, bam, bc){

  # Import reads
  tags <- scanBam(bam, param = ScanBamParam(tag=c("CB", "RE"), which=interval))[[1]]$tag

  # Account for the possibility that there may be no reads in the interval
  if(any(sapply(tags, is.null))){ return(NULL) }

  # Convert from list of length two (CB and RE) to a data frame
  tags <- as.data.frame(tags)

  # Remove cases where there is no cell barcode
  tags <- tags[complete.cases(tags), ]

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
