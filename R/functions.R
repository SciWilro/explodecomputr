################
# READING DATA #
################


#' Read binary GRM files into R
#'
#' @param rootname
#' @export
#' @return List of GRM and id data frames
readGRM <- function(rootname)
{
	bin.file.name <- paste(rootname, ".grm.bin", sep="")
	n.file.name <- paste(rootname, ".grm.N.bin", sep="")
	id.file.name <- paste(rootname, ".grm.id", sep="")

	cat("Reading IDs\n")
	id <- read.table(id.file.name)
	n <- dim(id)[1]
	cat("Reading GRM\n")
	bin.file <- file(bin.file.name, "rb")
	grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(bin.file)
	cat("Reading N\n")
	n.file <- file(n.file.name, "rb")
	N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(n.file)

	cat("Creating data frame\n")
	l <- list()
	for(i in 1:n)
	{
		l[[i]] <- 1:i
	}
	col1 <- rep(1:n, 1:n)
	col2 <- unlist(l)
	grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	

	ret <- list()
	ret$grm <- grm
	ret$id <- id
	return(ret)
}


#' Check if files exist
#'
#' Either checks of all files for a rootname exist, or returns a specific suffix
#'
#' @param rootname
#' @param software Either "plink" or "gcta"
#' @param suffix Either "all" or a specific suffic
#'
#' @export
#' @return Character string
checkRootname <- function(rootname, software="plink", suffix="all")
{
	sfs <- list(
		gcta = c("all", "grm.bin", "grm.N.bin", "grm.id"),
		plink = c("all", "bim", "fam", "bed")
	)

	stopifnot(software %in% names(sfs))
	stopifnot(suffix %in% sfs[[software]])

	if(suffix=="all")
	{
		nom <- paste(rootname, sfs[[software]][-1], sep=".")
		a <- file.exists(nom)
		if(!any(a))
		{
			stop(paste(nom[!a], collapse=" "))
		} else {
			return(rootname)
		}
	} else {
		if(file.exists(rootname))
		{
			return(rootname)
		} else if (file.exists(paste(rootname, suffix, sep="."))) {
			return(paste(rootname, suffix, sep="."))
		} else {
			stop(paste(suffix, "file does not exist"))
		}
	}
}



#' Read bim file
#' 
#' @param rootname
#' @export
#' @return Data frame
readBim <- function(rootname)
{
	nom <- checkRootname(rootname, "plink", "bim")
	bim <- read.table(nom, colClasses=c("character", "character", "numeric", "numeric", "character", "character"))
	names(bim) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
	return(bim)	
}


#' Read fam file
#' 
#' @param rootname
#' @export
#' @return Data frame
readFam <- function(rootname)
{
	nom <- checkRootname(rootname, "plink", "fam")
	fam <- read.table(nom, colClasses=c("character", "character", "numeric", "numeric", "character", "character"))
	names(fam) <- c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHEN")
	return(fam)	
}


#' Read output from plink --linear
#'
#' @param filename
#' @param h Is there a header
#'
#' @export
#' @return Data.frame
readPlinkLinear <- function(filename, h=TRUE)
{
	a <- read.table(filename, header=h, colClass=c("character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric"))
	if(!h)
	{
		names(a) <- c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P")
	}
	return(a)
}



##########
# GRAPHS #
##########



#' <brief desc>
#'
#' <full description>
#'
#' @param gwas <what param does>
#' @param  filename=NULL <what param does>
#'
#' @export
#' @return
qqplotpval <- function(P, filename=NULL)
{
	require(GenABEL)
	l <- estlambda(P, method="median")
	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
	if(!is.null(filename))
	{
		png(filename)
	}
	estlambda(P, method="median", plot=TRUE, main=nom)
	if(!is.null(filename))
	{
		dev.off()
	}
}

