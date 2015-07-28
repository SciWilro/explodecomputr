#' Create a genotype in HWE
#'
#' @param n number of individuals
#' @param p allele frequency
#' @export
#' @return array of 0s,1s,2s of length n
make_geno <- function(n, p)
{
	x <- rbinom(n, 2, p)
	return(x)
}

make_geno <- function(n, p, m = 1)
{
	if(length(p) == 1) p <- rep(p, m)
	stopifnot(length(p) == length(m))
	geno <- drop(sapply(p, function(x) rbinom(10, 2, x)))
	return(geno)
}


#' Create a phenotype using arbitrary number of known causal inputs
#'
#' @param effs array of variances for each input
#' @param indep matrix of independent variables corresponding to effs
#' @return simulated phenotype
#' @examples \dontrun{
#' g1 <- make_geno(1000, 0.5)
#' g2 <- make_geno(1000, 0.3)
#' p <- make_phen(cors=c(0.2, 0.1, 0.15, 0.4, 0.15), g1, g2, rnorm(1000), rnorm(1000))
#'}
make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)))
{
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors) <= 1)
	cors <- c(cors, 1-sum(cors))
	n <- nrow(indep)
	indep <- cbind(indep, rnorm(n))
	l <- ncol(indep)

	for(i in 1:l)
	{
		indep[,i] <- (indep[,i] - mean(indep[,i])) / sd(indep[,i]) * cors[i]
	}

	y <- apply(indep, 1, sum)
	return(y)
}


g1 <- make_geno(1000, 0.5)
g1 <- make_geno(1000, 0.5)

conf <- 

