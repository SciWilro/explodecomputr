#' Estimate standard error of regression coefficient
#'
#' @param b Effect size
#' @param n Sample size
#' @param vx Variance of x
#' @param vy variance of y
#'
#' @export
#' @return Expected standard error of b
expected_se_of_beta <- function(b, n, vx, vy)
{
	r <- b * vx / (sqrt(vx) * sqrt(vy))
	adjusted_r2 <- 1 - (n-1)/(n-2) * (1 - r^2)
	se_reg <- sqrt((1 - adjusted_r2) * vy)
	se_beta <- sqrt(se_reg^2 * 1 / (vx * (n-1)))
	return(se_beta)
}


#' Get the summary stats for use in 2 sample MR
#'
#' @param y vector of dependant variable
#' @param x matrix of independent variables
#'
#' @export
#' @return data frame of effects and standard errors
get_summary_stats <- function(y, x)
{
	mod <- as.data.frame(coefficients(summary(lm(y ~ x)))[-1,])
	names(mod) <- c("b", "se", "tval", "pval")
	as.data.frame(mod)
}
