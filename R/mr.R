#' Egger's regression for Mendelian randomisation
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param bootstrap=NULL Number of bootstraps to estimate standard error. If NULL then don't use bootstrap
#'
#' @export
#' @return List of results from MR Egger
eggers_regression <- function(b_exp, b_out, se_exp, se_out, bootstrap=NULL)
{
	stopifnot(length(b_exp) == length(b_out))
	stopifnot(length(se_exp) == length(se_out))
	stopifnot(length(b_exp) == length(se_out))

	sign0 <- function(x)
	{
		x[x==0] <- 1
		return(sign(x))
	}

	to_flip <- sign0(b_exp) == -1
	b_out = b_out*sign0(b_exp)
	b_exp = abs(b_exp) 
	dat <- data.frame(b_out=b_out, b_exp=b_exp, se_exp=se_exp, se_out=se_out, flipped=to_flip)
	mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
	smod <- summary(mod)

	b <- coefficients(smod)[2,1]
	se <- coefficients(smod)[2,2]
	pval <- coefficients(smod)[2,4]
	b_i <- coefficients(smod)[1,1]
	se_i <- coefficients(smod)[1,2]
	pval_i <- coefficients(smod)[1,4]

	if(!is.null(bootstrap))
	{
		boots <- eggers_regression_bootstrap(b_exp, b_out, se_exp, se_out, bootstrap)
		se <- boots$boots$se[boots$boots$param == "slope" & boots$boots$stat == "b" & boots$boots$what == "bootstrap"]
		se_i <- boots$boots$se[boots$boots$param == "intercept" & boots$boots$stat == "b" & boots$boots$what == "bootstrap"]
	}
	return(list(b = b, se = se, pval = pval, b_i = b_i, se_i = se_i, pval_i = pval_i, mod = smod, dat = dat))
}


plot_eggers_regression <- function(b_exp, b_out, se_exp, se_out, ylab = "Gene-outcome", xlab = "Gene-exposure")
{
	require(ggplot2)
	er <- eggers_regression(b_exp, b_out, se_exp, se_out)
	reg <- data.frame(a = coefficients(er$mod)[1,1], b = coefficients(er$mod)[2,1])

	p <- ggplot(er$dat, aes(y = b_out, x = b_exp)) +
		geom_errorbar(aes(ymax = b_out + se_out, ymin = b_out - se_out), width=0, colour="grey") +
		geom_errorbarh(aes(xmax = b_exp + se_exp, xmin = b_exp - se_exp), height=0, colour="grey") +
		geom_point(aes(colour=flipped)) +
		geom_abline(data=reg, aes(intercept = a, slope = b))
		labs(y = ylab, x = xlab, colour = "Exposure\nsign flipped")
	return(p)
}

eggers_regression_bootstrap <- function(b_exp, b_out, se_exp, se_out, nboot)
{
	require(reshape2)
	require(plyr)
	# Do bootstraps
	res <- array(0, c(n+1, 4))
	pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style=3) 
	for (i in 1:nboot)
	{
		setTxtProgressBar(pb, i)
		#sample from distributions of SNP betas
		xs <- rnorm(length(b_exp),b_exp,se_exp)
		ys <- rnorm(length(b_out),b_out,se_out)

		# Use absolute values for Egger reg
		ys <- ys*sign(xs)
		xs <- abs(xs)

		#weighted regression with given formula
		r <- summary(lm(ys ~ xs, weights=1/se_out^2))

		#collect coefficient from given line.
		res[i, 1] <- r$coefficients[1,1]
		res[i, 2] <- r$coefficients[1,2]
		res[i, 3] <- r$coefficients[2,1]
		res[i, 4] <- r$coefficients[2,2]
	}
	cat("\n")

	# Run original analysis
	b_out <- b_out*sign(b_exp)
	b_exp <- abs(b_exp)
	dat <- data.frame(b_out, b_exp, se_out, se_exp)
	orig <- coefficients(summary(lm(b_out ~ b_exp, weights=1/se_out^2)))
	res[n+1, ] <- c(orig[1,1], orig[1,2], orig[2,1], orig[2,2])
	res <- as.data.frame(res)
	res$what <- "bootstrap"
	res$what[n+1] <- "original"

	datl <- melt(res, measure.vars=c("V1", "V2", "V3", "V4"))
	datl$param <- "slope"
	datl$param[datl$variable %in% c("V1", "V2")] <- "intercept"
	datl$stat <- "b"
	datl$stat[datl$variable %in% c("V2", "V4")] <- "se"

	qu <- ddply(datl, .(param, stat, what), summarise, 
		m=mean(value),
		se=sd(value),
		q05=quantile(value, 0.05),
		q95=quantile(value, 0.95),
		pval=sum(value < 0)/length(value))

	res <- as.data.frame(res)
	names(res) <- c("b_i", "se_i", "b", "se")

	return(list(boots=qu, res=res, data=dat))
}


#' Perform 2 sample IV 
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Standard errors of genetic effects on exposure
#' @param se_out Standard errors of genetic effects on outcome
#' @param n=10000 Sample size
#' @param method Use "standard" method (default) or "ml"
#'
#' @export
#' @return List of results from 2 sample IV
two_sample_iv <- function(b_exp, b_out, se_exp, se_out, n=10000, method="standard")
{
	if(method == "standard")
	{
		b <- sum(b_exp*b_out / se_out^2) / sum(b_exp^2/se_out^2)
		se <- sqrt(1 / sum(b_exp^2/se_out^2))
		pval <- pt(abs(b) / se, df = n, low=FALSE)		
	} else if (method == "ml") {
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
	} else {
		stop("Must specify standard or ml for method")
	}
	return(list(b=b, se=se, pval=pval))
}

