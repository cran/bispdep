spcorrelogram.bi <- function (neighbours, varX, varY, order = 1, method = "corr",
                              style = "W", randomisation = TRUE, zero.policy = NULL,
                              spChk=NULL,alternative="greater",drop.EI2=FALSE)
{
    if (!inherits(neighbours, "nb"))
        stop("not a neighbours list")
    stopifnot(is.vector(varX))
    if (any(is.na(varX)))
        stop("no NAs permitted in variable")
    stopifnot(is.vector(varY))
    if (any(is.na(varY)))
        stop("no NAs permitted in variable")
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = get('.spdepOptions', envir = asNamespace('spdep'), inherits = FALSE))
    stopifnot(is.logical(zero.policy))
    if (is.null(spChk))
        spChk <- spdep::get.spChkOption()
    if (spChk && !chkIDs(varX, nb2listw(neighbours, zero.policy = zero.policy)))
        stop("Check of data of variable X and weights ID integrity failed")
    if (spChk && !chkIDs(varY, nb2listw(neighbours, zero.policy = zero.policy)))
        stop("Check of data of variable Y and weights ID integrity failed")
    if (order < 1)
        stop("order less than 1")
    nblags <- nblag(neighbours, maxlag = order)
    cardnos <- vector(mode = "list", length = order)
    for (i in 1:order) cardnos[[i]] <- table(card(nblags[[i]]))
    nobs <- sapply(cardnos, function(x) sum(x[names(x) > "0"]))
    if (any(nobs < 3))
        stop("spcorrelogram.bi: too few included observations in higher lags:\n\treduce order.")
    if (method == "corr") {
        lags.y <- matrix(0, nrow = length(varY), ncol = order)
        for (i in 1:order) lags.y[, i] <- lag.listw(nb2listw(nblags[[i]],
            style = style, zero.policy = zero.policy), varY, zero.policy = zero.policy)
        res <- cor(cbind(varX, lags.y))[1, -1]
        names(res) <- 1:order
    }
    else if ((method == "I") || (method == "C")) {
        res <- matrix(NA, nrow = order, ncol = 3)
        for (i in 1:order) {
            listw <- nb2listw(nblags[[i]], style = style, zero.policy = zero.policy)
            if (method == "I") {
                res[i, ] <-  moranbi.test(varX,varY,listw,randomisation = randomisation,
                             zero.policy = zero.policy,alternative = alternative, drop.EI2=drop.EI2)$estimate
            }
            else {
                res[i, ] <- gearybi.test(varX, varY, listw=listw, randomisation = randomisation,
                  zero.policy = zero.policy,alternative = alternative)$estimate
            }
        }
        rownames(res) <- 1:order
    }
    else stop("method unknown")
    obj <- list(res = res, method = method, cardnos = cardnos,
                varX = deparse(substitute(varX)),varY = deparse(substitute(varY)),
                alternative=alternative)
    class(obj) <- "spcorbi"
    obj
}

print.spcorbi <- function (x, p.adj.method="none", ...)
{
    cat("Bivariate spatial correlogram for", x$varX, "-",x$varY,"\nmethod: ")
    if (x$method == "I") {
        cat( "Bivariate Moran's I_{xy}\n")       # Modif. PL
    } else if(x$method == "C") {
        cat( "Bivariate Geary's C_{xy}\n")   # Modif. PL
    } else {
        cat("Bivariate spatial autocorrelation\n")          # Modif. PL
    }
    if ((x$method == "I") || (x$method == "C")) {
        res <- as.matrix(x$res)
	ZIxy <- (res[,1]-res[,2])/sqrt(res[,3])
#	pv <- p.adjust(2*pnorm(abs(ZI), lower.tail=FALSE), method=p.adj.method)
	if (x$alternative == "two.sided")
	  PrIxy <- 2 * pnorm(abs(ZIxy), lower.tail=FALSE)
	else if (x$alternative == "greater")
	  PrIxy <- pnorm(ZIxy, lower.tail=FALSE)
	else PrIxy <- pnorm(ZIxy)
	res <- cbind(res, ZIxy, PrIxy)
        rownames(res) <- paste(rownames(x$res), " (", sapply(x$cardnos,
          function(x) sum(x[names(x) > "0"])), ")", sep="")
        colnames(res) <- c("estimate", "expectation", "variance",
	    "Z", paste("Pr(Z) ",x$alternative))
        printCoefmat(res, ...)
    } else {
        res <- as.vector(x$res)
        names(res) <- names(x$res)
	print(res)
    }
    invisible(res)
}


plot.spcorbi <- function (x, main, ylab, ylim, ...)
{
    if (missing(main))
        main <- paste(x$varX,"-",x$varY)
    if ((x$method == "I") || (x$method == "C")) {
        lags <- as.integer(rownames(x$res))
        to.plot <- which((x$res[,3] > 0) & (x$res[,3] != Inf))  # Addition PL
        sd2 <- rep(0, nrow(x$res))                              # Addition PL
        sd2[to.plot] <- 2 * sqrt(x$res[to.plot, 3])             # Modif. PL
#        sd2 <- 2*sqrt(x$res[,3])
        if (missing(ylim)) {
            ylim <- range(c(x$res[,1]+sd2, x$res[,1]-sd2))
	}
        if (missing(ylab))
            if(x$method == "I") ylab <- "Bivariate Moran's  I_{xy}"            # Addition PL
            if(x$method == "C") ylab <- "Bivariate Geary's  C_{xy}"            # Addition PL
        plot(lags, x$res[,1], type="p", pch=18, ylim = ylim, main = main, ylab = ylab, xaxt = "n")
#        segments(lags, x$res[,1], lags, x$res[,2], lwd=4, col="grey")
        arrows(lags, x$res[,1]+sd2, lags, x$res[,1]-sd2, length=0.1, angle=90)
        arrows(lags, x$res[,1]-sd2, lags, x$res[,1]+sd2, length=0.1, angle=90)
        axis(1, at = lags)
        abline(h = x$res[1,2])
    }
    else {
        res <- as.vector(x$res)
        lags <- as.integer(names(x$res))
        if (missing(ylim))
            ylim <- c(-1, 1)
        if (missing(ylab))
            ylab <- "Bivariate spatial autocorrelation"
        plot(lags, res, type = "h", ylim = ylim, main = main, ylab = ylab,
            lwd = 4, xaxt = "n")
        axis(1, at = lags)
        abline(h = 0)
    }
}


