# Define a function to check if a bam file is valid for our purposes.
# Accepts the path to the BAM file and returns a character vector with two elements:
# The first element may take one of three possible values; "error", "warning", "pass"
# The second element is a message with provides details about any warnings/errors
check_bam <- function(bam_path){

  # Check the file exists
  if(!file.exists(bam_path)){ return(c("error", paste0("The BAM file '", bam_path, "', cannot be found")))}

  # Check if it is a BAM file
  bgzf_stream = gzfile(bam_path, 'r')
  bam_peek = suppressWarnings(readChar(bgzf_stream, 4))
  close(bgzf_stream)
  if(!identical(bam_peek, 'BAM\1')){ return(c("error", paste0("The file provided: ,", bam_path,", may not be a valid BAM file. Failed BAM check.")))}

  # Scan the first 10,000 reads to see if the BAM file contains the required CB and RE tags
  bamFile <- BamFile(bam_path)
  yieldSize(bamFile) <- 10
  open(bamFile)
  tags <- scanBam(bamFile, param = ScanBamParam(tag=c("CB", "RE")))[[1]]$tag
  close(bamFile)
  yieldSize(bamFile) <- NA
  if(all(is.null(tags$CB), is.null(tags$RE))){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ", bam_path,". Could not find either of the required 'CB' or 'RE' tags. Is this a valid Cell Ranger BAM file?"))) }
  if(is.null(tags$CB)){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ,", bam_path,",. Could not find the required 'CB' tag. Is this a valid Cell Ranger BAM file?"))) }
  if(is.null(tags$RE)){ return(c("warning", paste0("Scanned the first 10,000 reads in the supplied BAM file: ,", bam_path,",. Could not find the required 'RE' tag. Is this a valid Cell Ranger BAM file?"))) }

  # Passed all checks
  return(c("pass",paste0("The provided BAM file: '",bam_path, "' passed all of the BAM QC checks")))
}
