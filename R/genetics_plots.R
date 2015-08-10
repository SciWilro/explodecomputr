#' Estimate the inflation factor for a distribution of P-values (from GenABEL)
#'
#' Estimate the inflation factor for a distribution of P-values or 1df chi-square test. The major use of this procedure is the Genomic Control, but can also be used to visualise the distribution of P-values coming from other tests. Methods implemented include 'median' (median(chi2)/0.455...) and regression (of observed onto expected)
#'
#' @param data A vector of reals. If all are <=1, it is assumed that this is a vector of P-values, else it is treated as a vector of chi-squares
#' @param plot Whether the plot should be shown or not (default).
#' @param proportion The proportion of lowest P (or chi^2) values to be used when estimating the inflation factor lambda. Default = 1
#' @param method "regression" (default) or "median"
#' @param filter if the test statistics with 0-value of chi^2 should be excluded prior to estimation of lambda. Default = TRUE
#' @param df Number of degrees of freedom. Default = 1
#' @param ... Arguments passed to plot function
#'
#' @export
#' @return A list with elements estimate and se
est_lambda <- function (data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1, ...)
{
	data <- data[which(!is.na(data))]
	if (proportion > 1 || proportion <= 0) 
		stop("proportion argument should be greater then zero and less than or equal to one")
	ntp <- round(proportion * length(data))
	if (ntp < 1) 
		stop("no valid measurements")
	if (ntp == 1) {
		warning(paste("One measurement, lambda = 1 returned"))
		return(list(estimate = 1, se = 999.99))
	}
	if (ntp < 10) 
		warning(paste("number of points is too small:", ntp))
	if (min(data) < 0) 
		stop("data argument has values <0")
	if (max(data) <= 1) {
		data <- qchisq(data, 1, lower.tail = FALSE)
	}
	if (filter) {
		data[which(abs(data) < 1e-08)] <- NA
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
	out <- list()
	if (method == "regression") {
		s <- summary(lm(data ~ 0 + ppoi))$coeff
		out$estimate <- s[1, 1]
		out$se <- s[1, 2]
	}
	else if (method == "median") {
		out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
			df)
		out$se <- NA
	}
	else {
		stop("'method' should be either 'regression' or 'median'!")
	}
	if (plot) {
		lim <- c(0, max(data, ppoi, na.rm = TRUE))
		oldmargins <- par()$mar
		par(mar = oldmargins + 0.2)
		plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
			ylab = expression("Observed " ~ chi^2), ...)
		abline(a = 0, b = 1)
		abline(a = 0, b = out$estimate, col = "red")
		par(mar = oldmargins)
	}
	out
}

#' QQ plot pvalues
#'
#' Creates a qq plot using GenABEL's function. Calculates lambda and prints to title
#'
#' @param P array of pvalues
#' @param filename=NULL If specified the plot will be printed to png
#'
#' @export
#' @return NULL
qqplot_pval <- function(P, filename=NULL)
{
	l <- est_lambda(P, method="median")
	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
	if(!is.null(filename))
	{
		png(filename)
	}
	est_lambda(P, method="median", plot=TRUE, main=nom)
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
manhattan_plot <- function(p, chr, pos, filename=NULL, width=15, height=7, threshold=-log10(0.05/1000000))
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
