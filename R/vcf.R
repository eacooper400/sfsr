## vcf.R
## Read VCF file and calculate SFS


#' Calculate either a 1 or 2 dimensional SFS from VCF format data
#'
#' addition to orig sfsr pkg by L.Cooper 11-2017
#'
#' @param vcf.list a list of length>=1 & length<3 VCF format data frames
#' use read.vcf to read VCF file into a table/data frame before using this fxn
#'
#' @return either a named Vector (1D) or Matrix (2D) that is an SFS
#'
#' @export
sfs.fromVcf <- function(vcf.list) {

	## calculate num. chrom. (i.e. max derived count); assumes diploid
	nchr1 = (ncol(vcf.list[[1]]) - 9) * 2

	## get the count of derived alleles at every SNP
	## assumes ALT allele is derived!!
	y=lapply(vcf.list, function(x) apply(x, 1, derivedCount))
	y=lapply(y, unname)

	## If only 1 VCF table provided, create a 1d sfs
	if (length(y)==1) {
        	z=table(factor(y[[1]], levels=0:nchr1))
        	return(z)
	} 

	## If 2 VCF tables, get the 2d joint sfs between them
	else {
		### allow for diff. number of samples in vcf 2
		nchr2 = (ncol(vcf.list[[2]]) - 9) * 2
		
		### set up the empty matrix (allow for 0 counts)
		sfs = matrix(0, ncol=(nchr2+1), nrow=(nchr1+1))

		### for every possible combo of counts, fill in the matrix
		for (i in 0:nchr1) {
            		s1 = which(y[[1]]==i)
            		for (j in 1:nchr2) {
                		s2 = which(y[[2]]==j)
                		sfs[i,j] = sfs[i,j] + length(intersect(s1,s2))
            		}
        	}
		return(sfs)
	}
}

#' Read in VCF format file into R as a table
#'
#' addition to orig sfsr pkg by L.Cooper 11-2017
#'
#' @param file a file in uncompressed VCF format
#'
#' @return a Data Frame with VCF data
#'
#' @export
read.vcf <- function(file, special.char="##", ...) {

	## Search specifically for the ##HEADER lines,
	## and remove them   	
    	my.search.term=paste0(special.char, ".*")
    	all.lines=readLines(file)
    	clean.lines=gsub(my.search.term, "",  all.lines)
	
	## Fix the column headers so that R won't skip them
    	clean.lines=gsub("#CHROM", "CHROM", clean.lines)
    
	## Read the cleaned lines in with read.table
	read.table(..., text=paste(clean.lines, collapse="\n"))
	
}

#' Get the count of derived (ALT) alleles from a row in a vcf
#'
#' addition to orig sfsr pkg by L.Cooper 11-2017
#'
#' @param row A vector of data in VCF format corresponding to 1 SNP
#'
#' @return a single integer of the number of derived alleles in a
#' 	sample at this SNP
#'
#' @export

derivedCount <- function(row) {

	## coerce row to a character vector
	row=as.vector(row, mode="character")

	## extract just the GT fields for the sample columns
	gt=gsub(":.*", "", row[10:length(row)])
	
	## count the "1" allele
	dc=length(which(unlist(strsplit(gt, "/"))=="1"))
	return(dc)

}

