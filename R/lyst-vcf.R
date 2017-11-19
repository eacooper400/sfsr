#' Polar Bear LYST gene SNP data
#'
#' Data from 7 polar bear, 2 grizzly bears, and 1 black bear
#' in VCF format, selected to contain the LYST gene.  
#'
#' @docType data
#'
#' @usage data(lyst)
#'
#' @format A list of 3 data sets: the vcf is in a data.frame, while
#' the fasta and gff3 are stored as vectors of lines
#'
#' @keywords datasets
#'
#' @references Liu et al. (2014) Cell 157(4): 785-794
#'
#' @examples
#' read.vcf(lyst.vcf.txt)
