#' Create a genotype or many genotypes in HWE
#'
#' Under Hardy Weinberg a genotype is binomially distributed, with p = allele frequency
#'
#' @param n Sample size
#' @param p Array of allele frequencies to simulate. Can be a single value.
#'
#' @export
#' @return Array or matrix of genotypes
#' @examples \dontrun{
#' g <- make_geno(100, 0.5) # To make a single genotype
#' g <- make_geno(100, c(0.5, 0.3, 0.6)) # To make 3 genotypes in a matrix
#' }
make_geno <- function(n, p)
{
	stopifnot(all(p >= 0 & p <= 1))
	geno <- drop(sapply(p, function(x) rbinom(n, 2, x)))
	return(geno)
}


#' Create a phenotype using arbitrary number of known causal inputs
#'
#' For a set of x variables and effect sizes for each variable (b) 
#' y is constructed such that 
#' y = Xb + e
#' Given that the variance explained in y by X, is
#' r^2 = sum(b * var(x) / (sqrt(x)*sqrt(y)))
#' we can model e ~ N(0, 1 - r^2)
#'
#' @param effs Array of beta values for each input. Leave the vx and vy values to default (=1) to allow effs to be equal to the correlation between y and each x
#' @param indep Matrix of independent variables corresponding to effs
#' @param vy=1 The output variance of y
#' @param vx=rep(1, length(effs)) The desired scaled variance of x
#'
#' @export
#' @return Numeric array, simulated phenotype
#'
#' @examples \dontrun{
#' g1 <- make_geno(1000, 0.5)
#' g2 <- make_geno(1000, 0.3)
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' y <- make_phen(effs=c(0.2, 0.1, 0.15, 0.4), cbind(g1, g2, x1, x2))
#' 
#'}
make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)))
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy))
	return(y)
}


b = cov(x,y) / var(x)

cor(x, y) = cov(x, y) / (sd(x) * sd(y))





x <- rnorm(10000)
y <- make_phen(0.2, x, vy = 2)
cor(x, y)
coef(lm(y ~ x))[2]

var(x)
var(y)



g <- make_geno(100000, 0.5, 30)
x <- make_phen(effs = rnorm(30)/30, g)

y <- make_phen(effs = 0.3, x)
gx <- coef(lm(x ~ g))[-1]
gy <- coef(lm(y ~ g))[-1]
cor(y, x)
plot(gy ~ gx)

y <- make_phen(effs = c(0.3, rnorm(30)/40), cbind(x, g))
gx <- coef(lm(x ~ g))[-1]
gy <- coef(lm(y ~ g))[-1]
cor(y, x)
plot(gy ~ gx)

y <- make_phen(effs = c(0.3, abs(rnorm(30)/40)), cbind(x, g))
gx <- coef(lm(x ~ g))[-1]
gy <- coef(lm(y ~ g))[-1]
cor(y, x)
plot(gy ~ gx)

y <- make_phen(effs = c(0.3, rep(0.01,30)), cbind(x, g))
gx <- coef(lm(x ~ g))[-1]
gy <- coef(lm(y ~ g))[-1]
cor(y, x)
plot(gy ~ gx)

y <- make_phen(effs = c(0.3, rep(-0.01,30)), cbind(x, g))
gx <- coef(lm(x ~ g))[-1]
gy <- coef(lm(y ~ g))[-1]
cor(y, x)
plot(gy ~ gx)

