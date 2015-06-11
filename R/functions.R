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


#' Convert long GRM format to matrix format
#'
#' @param grm Result from readGRM
#' @export
#' @return Matrix of n x n
makeGRMmatrix <- function(grm)
{
	mat <- diag(nrow(grm$id))
	mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
	mat <- t(mat)
	nsnpvec <- subset(grm$grm, id1 != id2)$N
	mat[upper.tri(mat, diag=FALSE)] <- nsnpvec
	return(mat)
}



#' Write readGRM style output back to binary GRM for use with GCTA
#'
#' @param grm Output from \link{readGRM}
#' @param rootname
#' @export
writeGRM <- function(grm, rootname)
{
	bin.file.name <- paste(rootname, ".grm.bin", sep="")
	n.file.name <- paste(rootname, ".grm.N.bin", sep="")
	id.file.name <- paste(rootname, ".grm.id", sep="")
	write.table(grm$id, id.file.name, row=F, col=F, qu=F)
	n <- dim(grm$id)[1]
	bin.file <- file(bin.file.name, "wb")
	writeBin(grm$grm$grm, bin.file, size=4)
	close(bin.file)
	n.file <- file(n.file.name, "wb")
	writeBin(grm$grm$N, n.file, size=4)
	close(n.file)
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
	fam <- read.table(nom, colClasses=c("character", "character", "character", "character", "character", "character"))
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


#' Read GCTA xmat file
#'
#' @param filename
#'
#' @export
#' @return list
readGctaXmat <- function(filename)
{
        a <- read.table(filename, colClass="character")
        snps <- as.character(a[1,-c(1:2)])
        ids <- a[-c(1:2), 1:2]
        alleles <- as.character(a[2,-c(1:2)])
        xmat <- matrix(as.numeric(as.matrix(a[-c(1,2), -c(1,2)])), nrow(a)-2, ncol(a)-2)
        ids <- data.frame(a)
        names(ids) <- c("FID", "IID")
        snps <- data.frame(snp=snps, allele=alleles)
        return(list(xmat=xmat, snps=snps, ids=ids))
}



##########
# GRAPHS #
##########



#' QQ plot pvalues
#'
#' Creates a qq plot using GenABEL's function. Calculates lambda and prints to title
#'
#' @param P array of pvalues
#' @param filename=NULL If specified the plot will be printed to png
#'
#' @export
#' @return NULL
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



#' Manhattan Plot
#'
#' @param p P values
#' @param  chr Chromosomes (won't handle X, must not be factors)
#' @param  pos Physical position
#' @param  filename=NULL If specified will print to file
#' @param  width=15
#' @param  height=7
#' @param  threshold=-log10(0.05/1000000) For thresholds
#'
#' @export
#' @return NULL
manhattanPlot <- function(p, chr, pos, filename=NULL, width=15, height=7, threshold=-log10(0.05/1000000))
{
	require(ggplot2)
	dat <- data.frame(chrom=as.numeric(chr), bp=pos, pval=-log10(p))
	dat <- dat[order(dat$chrom, dat$bp), ]
	dat$col <- dat$chr %% 2 + 1
	dat <- subset(dat, !is.na(pval))
	
	pl <- ggplot(dat, aes(x=bp, y=pval)) +
	geom_point(aes(colour=factor(col))) +
	facet_grid(. ~ chrom, scale="free_x", space="free_x") +
	theme(legend.position="none") +
	scale_colour_manual(values=c("#404040", "#ca0020")) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	ylim(0, max(c(threshold, dat$pval, na.rm=TRUE))) +
	labs(y=expression(-log[10]*p), x="Position") +
	geom_hline(yintercept=threshold)

	if(!is.null(filename))
	{
		ggsave(filename, pl, width=width, height=height)		
	} else {
		print(pl)
	}
}







#################
# DATA CLEANING #
#################


removeOutliersRecursive <- function(x, niter=5, d=3)
{
	for(i in 1:niter)
	{
		s <- sd(x, na.rm=T)
		m <- mean(x, na.rm=T)
		index <- x > m + d*s | x < m - d*s
		x[index] <- NA
	}
	return(x)
}

scaleRankTransform <- function(x) 
{
	require(GenABEL)
	s <- sd(x, na.rm=TRUE)
	m <- mean(x, na.rm=TRUE)
	x <- rntransform(x) * s + m
	return(x)
}

ztransform <- function(x)
{
	(x-mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

adjustMeanVariance <- function(y, x, keep.scale=TRUE)
{
	require(plyr)
	s <- sd(y, na.rm=TRUE)
	m <- mean(y, na.rm=TRUE)
	d <- data.frame(y=y, x=x, index=1:length(y))

	d <- ddply(d, .(x), mutate, y1 = ztransform(y))
	if(keep.scale)
	{
		d$y1 <- d$y1 * s + m
	}
	d <- d[order(d$index), ]
	return(d$y1)
}

adjustCov <- function(y, x, keep.scale=TRUE)
{
	m <- mean(y, na.rm=TRUE)
	s <- sd(y, na.rm=TRUE)
	res <- residuals(lm(y ~ x, na.action=na.exclude))
	if(keep.scale)
	{
		res <- res * s + m
	}
	return(res)
}

