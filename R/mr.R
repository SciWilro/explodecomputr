#' Egger's regression for Mendelian randomisation
#'
#' 
#'
#' @param b_exp <what param does>
#' @param  b_out <what param does>
#' @param  se_exp <what param does>
#' @param  se_out <what param does>
#'
#' @export
#' @return
eggers_regression <- function(b_exp, b_out, se_exp, se_out)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))
	b_out = b_out*sign(b_exp)
	b_exp = abs(b_exp) 
	dat <- data.frame(b_out=b_out, b_exp=b_exp, se_exp=se_exp, se_out=se_out)	
	mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
	smod <- summary(mod)
	tab <- coefficients(smod)
	rownames(tab)[2] <- paste(xlab, "on", ylab)
	return(list(mod = smod, dat = dat))
}

plot_eggers_regression <- function(b_exp, b_out, se_exp, se_out, ylab = "Gene-outcome", xlab = "Gene-exposure")
{
	require(ggplot2)
	er <- eggers_regression(b_exp, b_out, se_exp, se_out)
	reg <- data.frame(a = coefficients(er$smod)[1,1], b = coefficients(er$mod)[2,1])

	p <- ggplot(er$dat, aes(y = b_out, x = b_exp)) +
		geom_errorbar(aes(ymax = b_out + se_out, ymin = b_out - se_out), width=0, colour="grey") +
		geom_errorbarh(aes(xmax = b_exp + se_exp, xmin = b_exp - se_exp), height=0, colour="grey") +
		geom_point() +
		geom_abline(data=reg, aes(intercept = a, slope = b)) +
		labs(y = ylab, x = xlab)
	return(p)
}

eggers_regression_bootstrap <- function(x,y,xse,yse,n){
	# Do bootstraps
	res <- array(0, c(n+1, 4))
	for (i in 1:n)
	{
		#sample from distributions of SNP betas
		xs <- rnorm(length(x),x,xse)
		ys <- rnorm(length(y),y,yse)

		# Use absolute values for Egger reg
		ys <- ys*sign(xs)
		xs <- abs(xs)

		#weighted regression with given formula
		r <- summary(lm(ys ~ xs, weights=1/yse^2))

		#collect coefficient from given line.
		res[i, 1] <- r$coefficients[1,1]
		res[i, 2] <- r$coefficients[1,2]
		res[i, 3] <- r$coefficients[2,1]
		res[i, 4] <- r$coefficients[2,2]
	}

	# Run original analysis
	y <- y*sign(x)
	x <- abs(x)
	dat <- data.frame(y, x, yse, xse)
	orig <- coefficients(summary(lm(y ~ x, weights=1/yse^2)))
	res[n+1, ] <- c(orig[1,1], orig[1,2], orig[2,1], orig[2,2])
	res <- as.data.frame(res)
	res$what <- "bootstrap"
	res$what[n+1] <- "original"
	return(list(res, dat))
}

two_sample_iv_ml <- function(b_exp, b_out, se_exp, se_out, n=10000)
{
	loglikelihood <- function(param) {
		return(1/2*sum((b_exp-param[1:length(b_exp)])^2/se_exp^2)+1/2*sum((b_out-param[length(b_exp)+1]*param[1:length(b_exp)])^2/se_out^2))
	}
	opt <- optim(
		c(b_exp, sum(b_exp*b_out/se_out^2)/sum(b_exp^2/se_out^2)),
		loglikelihood, 
		hessian=TRUE, 
		control = list(maxit=25000))

	b <- opt$par[length(b_exp)+1]
	se <- sqrt(solve(opt$hessian)[length(b_exp)+1,length(b_exp)+1])
	pval <- pt(abs(b) / se, df = n, low=FALSE)
	return(list(b=b, se=se, pval=pval))
}
