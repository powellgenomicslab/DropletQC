#'Check file is a valid Cell Ranger Barcoded BAM file
#'
#'@description
#'
#'Accepts the path to a BAM file. The function executes two checks on the file.
#'The first check passes the file through `gunzip` and tests whether the bytes
#'at the beginning of the BAM file are correct. This ide was taken from a
#'helpful suggestion by Konrad Rudolph in a StackExchange response to "Is there
#'an efficient way to check an input BAM in R?". The second check reads in the
#'forst 10,000 reads from the BAM file and looks to see if the two required
#'tags; CB (cell barcode) and RE (region type : E = exonic, N = intronic, I =
#'intergenic) are present. The function returns a character vector with two
#'elements. The first element may take one of three possible values; "error",
#'"warning", "pass". The second element is a message with provides extra details
#'about any warnings/errors. Note that if the file passes the first file check,
#'an "error" is returned as the first value. If the file fails the second check
#'a "warning" is returned as the forst value. Note that the function was written
#'as a simple helper function to be called internally by
#'`doubletQC::nuclear_fraction()` and isn't intended for more general use.
#'
#'@param bam_path character, bam_path should be a character vector pointing to
#'  the BAM file you want to check.
#'
#'@return character. The function returns a character vector with two elements.
#'  The first element may take one of three possible values; "error", "warning",
#'  "pass". The second element is a message with provides extra details about
#'  any warnings/errors.
#'
#'@keywords internal
#'
check_bam <- function(bam_path){

  # Check the file exists
  if(!file.exists(bam_path)){ return(c("error", paste0("The BAM file '", bam_path, "', cannot be found")))}

  # Check if it is a BAM file
  bgzf_stream = gzfile(bam_path, 'r')
  bam_peek = suppressWarnings(readChar(bgzf_stream, 4))
  close(bgzf_stream)
  if(!identical(bam_peek, 'BAM\1')){ return(c("error", paste0("The file provided: ,", bam_path,", may not be a valid BAM file. Failed BAM check.")))}

  # Scan the first 10,000 reads to see if the BAM file contains the required CB and RE tags
  bamFile <- Rsamtools::BamFile(bam_path)
  Rsamtools::yieldSize(bamFile) <- 10
  open(bamFile)
  tags <- Rsamtools::scanBam(bamFile, param = Rsamtools::ScanBamParam(tag=c("CB", "RE")))[[1]]$tag
  close(bamFile)
  Rsamtools::yieldSize(bamFile) <- NA
  if(all(is.null(tags$CB), is.null(tags$RE))){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ", bam_path,". Could not find either of the required 'CB' or 'RE' tags. Is this a valid Cell Ranger BAM file?"))) }
  if(is.null(tags$CB)){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ,", bam_path,",. Could not find the required 'CB' tag. Is this a valid Cell Ranger BAM file?"))) }
  if(is.null(tags$RE)){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ,", bam_path,",. Could not find the required 'RE' tag. Is this a valid Cell Ranger BAM file?"))) }

  # Passed all checks
  return(c("pass",paste0("The provided BAM file: '",bam_path, "' passed all of the BAM QC checks")))
}
