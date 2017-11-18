## io.R
## read and write SFS from disk

#' Read 'flattened' site frequency spectrum (SFS) from a file
#' 
#' @param ff filename
#' @param dims sample sizes (number of chromosomes) along each dimension of SFS
#' @param dtype what sort of frequences to expect; best let the procedure auto-detect
#' @param bootstraps logical; if \code{TRUE}, read bootstrap replicates (one per line) if they are present
#' @param repolarize logical; if \code{TRUE}, swap ancestral and derived states (see Details)
#' @param ... ignored
#' 
#' @return a site frequency spectrum (SFS): a k-dimensional array representing the joint frequencies of derived alleles in each of k populations
#' 
#' @details The SFS is expected to be provided in a text file as space-separated numbers (integers or floating-point), in row-major order.  
#' An \emph{unfolded} SFS (ie. polarized against the derived allele) is expected, although this may be relaxed in future
#'   
#' NB: Recent versions of \code{ANGSD} apparently get the ancestral and derived alleles backwards.  Use \code{repolarize = TRUE} to correct
#' this issue at runtime.
#' 
#' @export
read_sfs <- function(ff, dims, dtype = double, bootstraps = TRUE, repolarize = FALSE, ...) {
	
	#hack
	angsd <- FALSE
	
	## read comment lines with metadata
	message("Reading SFS metadata...")
	infile <- file(ff)
	open(infile)
	meta <- vector("list")
	mn <- character()
	j <- 1
	while (length(line <- readLines(infile, n = 1, warn = FALSE))) {
		
		if (!grepl("^#", line)) {
			break
		}
		
		vals <- stringr::str_match(line, "#+(\\w+)=(.+)")[1,2:3]
		vv <- unlist(stringr::str_split(vals[2], ","))
		mn <- c(mn, vals[1])
		vv.num <- suppressWarnings(as.numeric(vv))
		vv.bool <- suppressWarnings(as.logical(vv))
		if (!any(is.na(vv.num)))
			meta[[j]] <- vv.num
		else if (!any(is.na(vv.bool)))
			meta[[j]] <- vv.bool
		else
			meta[[j]] <- vv
		
		j <- j+1
		
	}
	close(infile)
	names(meta) <- mn
	
	if ("pops" %in% names(meta))
		names(dims) <- meta[["pops"]]
	
	## convenience wrapper for base::scan()
	read_one <- function(ff, ...) {
		line <- scan(ff, what = dtype(), comment.char = "#", quiet = TRUE,
					 nlines = 1, strip.white = TRUE, ...)
		return(line)
	}
	
	## now read spectrum itself
	message("Reading SFS: expecting ", prod(dims+1), " entries...")
	boots <- list()
	infile <- file(ff, "r")
	ii <- 1
	if (!angsd) {
		while(length(x <- read_one(infile)) || ii == 1) {
			if (length(x) == 0)
				next
			else if (ii > 1 && !bootstraps)
				break
			sfs <- aperm( array(x, dim = rev(dims+1)) )
			dimnames(sfs) <- lapply(dims, function(z) c(0,seq_len(z)))
			boots[[ii]] <- sfs
			ii <- ii+1
		}
	}
	else {
		while(length(x <- read_one(infile))) {
			if (ii > 1 && !bootstraps)
				break
			sfs <- array(x, dim = dims+1)
			dimnames(sfs) <- lapply(dims, function(z) c(0,seq_len(z)))
			boots[[ii]] <- sfs
			ii <- ii+1
		}
	}
	close(infile)
	
	## add metadata
	rez <- boots[[1]]
	message("Read ", length(rez), " items.")
	if (length(boots) > 1) {
		message("\t[ ", length(boots), " bootstrap replicates ]")
		attr(rez, "bootstraps") <- boots
	}
	if (length(meta)) {
		for (m in names(meta))
			attr(rez, m) <- meta[[m]]
	}
	
	if (repolarize)
		rez <- repolarize_sfs(rez)
	
	message("First SFS has ", sum(rez), " sites.")
	if (!inherits(rez, "sfs"))
		class(rez) <- c("sfs", class(rez))
	
	
	return(rez)
	
}

#' Convert SFS to vector representation for saving in a file
#' 
#' @param sfs an SFS
#' 
#' @return a vector containing the 'flattened' representation of the SFS, in row-major order
#' 
#' @export
flatten_sfs <- function(sfs) {
	
	## because numpy flattens row-wise and R column-wise
	as.vector(aperm(sfs))
	
}

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
