#' Initialize new classes
#'
#' Creates new classes that were needed for this package
#'
#' @param NA
#' @return NA
#' @example
#' init_classes()
init_classes <- function(){
  setClass("rateTable", representation(x="matrix"))
  setClass("riskTable", representation(x="matrix"))
  setClass("ccTable", representation(x="matrix"))
  setClass("mccTable", representation(x="matrix"))
}

#' Converts data to rate table
#'
#' Creates new matrix represented as a table with person time
#'
#' @param eCases number of exposed cases
#' @param ePT amount of exposed person-time
#' @param neCases number of not exposed Cases
#' @param nePT amount of unexposed person-time
#' @return NA
#' @example
#' as.rateTable(...)
as.rateTable <- function(eCases, ePT, neCases, nePT, ...) {
  init_classes()
	if (nargs() == 4) {
		obj <- matrix(c(eCases, ePT, neCases, nePT), nrow=2)
		class(obj) <- "rateTable"
	}
	else {
		if (nargs() %% 4) {stop("Incorrect number of table entries (need multiple of 4)")}
		arglist <- as.list(match.call())[-1]
		ntable <- nargs()/4
		obj <- vector("list", length=ntable)
		class(obj) <- c("rateTable", "list")
		for (i in 1:ntable) {
			a <- 4*(i-1)
			obj[[i]] <- as.rateTable(arglist[[1+a]], arglist[[2+a]], arglist[[3+a]], arglist[[4+a]])
		}
	}
	return(obj)
}

#' Converts data to rate table (updated)
#'
#' Creates new matrix represented as a table with person time
#'
#' @param cases number of cases
#' @param person.time amount of person time
#' @return NA
#' @example
#' as.rateTable.new(...)
as.rateTable.new <- function(cases,person.time, ...) {
  init_classes()
  if (length(cases) == 2) {
    obj <- matrix(c(as.numeric(cases[2]), as.numeric(person.time[2]), as.numeric(cases[1]), as.numeric(person.time[1])), nrow=2)
    class(obj) <- "rateTable"
  }
  else {
    if (length(cases)>2 & length(cases) %% 2) {stop("Incorrect number of table entries (need multiple of 2)")}
    ntable <- length(cases)/2
    obj <- vector("list", length=ntable)
    class(obj) <- c("rateTable", "list")
    for (i in 1:ntable) {
      a <- 2*(i-1)
      obj[[i]] <- as.rateTable(as.numeric(cases[2+a]), as.numeric(person.time[2+a]), as.numeric(cases[1+a]), as.numeric(person.time[1+a]))
    }
  }
  return(obj)
}

#' Prints rate table
#'
#' Creates output of rate table
#'
#' @param x rate table object
#' @return NA
#' @example
#' print.rateTable(...)
print.rateTable <- function(x,...) {
  init_classes()
	if (!("list" %in% class(x))) {
		cat("\t\t\t\t\tExposed\t\tNot exposed\n")
		cat("\t", paste("Cases", format(x[1,1], width=7), format(x[1,2], width=7), sep="\t\t\t"), "\n")
		cat("\t", paste("Person-time", format(x[2,1], width=7), "\t", format(x[2,2], width=7), sep="\t"))
	} else {
		for (i in 1:length(x)) {
			cat(paste("Stratum ", i, ":\n", sep=""))
			print.rateTable(x[[i]])
			cat("\n\n")
		}
		cat("Crude data:\n")
		print.rateTable(Reduce("+", x))
	}
}

#' Prints summary rate table
#'
#' Creates summary output of rate table
#'
#' @param x rate table object
#' @param alpha alpha threshold
#' @return NA
#' @example
#' summary.rateTable(...)
summary.rateTable <- function(x, alpha=.05, ...) {
  init_classes()
	mult <- qnorm(1-alpha/2)
	if (!("list" %in% class(x))) {
		rat <- (x[1,1]/x[2,1])/(x[1,2]/x[2,2])
		dif <- (x[1,1]/x[2,1]) - (x[1,2]/x[2,2])
		cir <- exp(log(rat) + c(-mult, mult)*sqrt(1/x[1,1]+1/x[1,2]))
		cid <- dif + c(-mult, mult)*sqrt(x[1,1]/x[2,1]^2+x[1,2]/x[2,2]^2)

		ctot <- x[1,1]+x[1,2]
		ttot <- x[2,1]+x[2,2]
		teststat <- (x[1,1]-x[2,1]*ctot/ttot)^2/(x[2,1]*x[2,2]*ctot/ttot^2)
		pval <- 1-pchisq(teststat, 1)

		cat("Point estimate of IRR: ", format(rat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for IRR: ", format(cir[1], digits=4, nsmall=4),
			" - ", format(cir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Point estimate of IRD: ", format(dif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for IRD: ", format(cid[1], digits=4, nsmall=4),
			" - ", format(cid[2], digits=4, nsmall=4), sep=""), "\n")
		if (pval < .0001) {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p =",
				format(pval, digits=4, nsmall=4)))
		}
	}
	else {
		collapseMat <- Reduce("+", x)
		a <- collapseMat[1,1]
		b <- collapseMat[1,2]
		c <- collapseMat[2,1]
		d <- collapseMat[2,2]
		cruderat <- (a*d)/(b*c)
		crudedif <- (a/c) - (b/d)
		crudecir <- exp(log(cruderat) + c(-mult, mult)*sqrt(1/a+1/b))
		crudecid <- crudedif + c(-mult, mult)*sqrt(a/c^2+b/d^2)
		crudestat <- (a-c*(a+b)/(c+d))^2/(c*d*(a+b)/(c+d)^2)
		crudepval <- 1-pchisq(crudestat, 1)

		ornum <- sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]/(y[2,1]+y[2,2])))))
		orden <- sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]/(y[2,1]+y[2,2])))))
		ormh <- ornum/orden
		A <- sum(unlist(lapply(x, function(y) return((y[1,1]+y[1,2])*y[2,1]*y[2,2]/(y[2,1]+y[2,2])^2))))
		B <- sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]/(y[2,1]+y[2,2])))))
		CC <- sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]/(y[2,1]+y[2,2])))))
		varln <- A/(B*CC)
		mhcir <- exp(log(ormh) + c(-mult, mult)*sqrt(varln))

		wgts <- unlist(lapply(x, function(y) return(y[2,1]^2*y[2,2]^2/(y[2,2]^2*y[1,1]+y[2,1]^2*y[1,2]))))
		irds <- unlist(lapply(x, function(y) return(y[1,1]/y[2,1]-y[1,2]/y[2,2])))
		ivdif <- sum(wgts*irds)/sum(wgts)
		ivcid <- ivdif + c(-mult, mult)*sqrt(1/sum(wgts))
		stratnum <- (a - sum(unlist(lapply(x, function(y) return(y[2,1]*(y[1,1]+y[1,2])/(y[2,1]+y[2,2]))))))^2
		stratden <- sum(unlist(lapply(x, function(y) return((y[1,1]+y[1,2])*y[2,1]*y[2,2]/(y[2,1]+y[2,2])^2))))
		stratstat <- stratnum/stratden
		stratpval <- 1-pchisq(stratstat, 1)

		rhomstat <- sum(unlist(lapply(x, function(y) return((log(y[1,1]*y[2,2]/y[1,2]/y[2,1])-log(ormh))^2/(1/y[1,1]+1/y[1,2])))))
		rhompval <- 1-pchisq(rhomstat, length(x)-1)
		dhomstat <- sum(unlist(lapply(x, function(y) return((y[1,1]/y[2,1]-y[1,2]/y[2,2]-ivdif)^2/(y[1,1]/y[2,1]^2+y[1,2]/y[2,2]^2)))))
		dhompval <- 1-pchisq(dhomstat, length(x)-1)

		cat("Crude estimate of IRR: ", format(cruderat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for crude IRR: ", format(crudecir[1], digits=4,
			nsmall=4), " - ", format(crudecir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Crude estimate of IRD: ", format(crudedif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for crude IRD: ", format(crudecid[1], digits=4,
			nsmall=4), " - ", format(crudecid[2], digits=4, nsmall=4), sep=""), "\n")
		if (crudepval < .0001) {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p =",
				format(crudepval, digits=4, nsmall=4)))
		}
		cat("\n\nMantel-Haenszel estimate of IRR: ", format(ormh, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for MH IRR: ", format(mhcir[1], digits=4,
			nsmall=4), " - ", format(mhcir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Inverse-variance estimate of IRD: ", format(ivdif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for IV IRD: ", format(ivcid[1], digits=4,
			nsmall=4), " - ", format(ivcid[2], digits=4, nsmall=4), sep=""), "\n")
		if (stratpval < .0001) {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4),
				"with p < 0.0001"))
		} else {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4), "with p =",
				format(stratpval, digits=4, nsmall=4)))
		}

		cat(paste("\n\nTest of IRR homogeneity: X^2 =", format(rhomstat, digits=4, nsmall=4), "with p =",
				format(rhompval, digits=4, nsmall=4)))
		cat(paste("\nTest of IRD homogeneity: X^2 =", format(dhomstat, digits=4, nsmall=4), "with p =",
				format(dhompval, digits=4, nsmall=4)))
	}
}

#' Converts data to risk table
#'
#' Creates new matrix represented as a table for a closed cohort (or cross-sectional)
#'
#' @param eCases number of exposed cases
#' @param eNonCases number of exposes non-cases
#' @param neCases number of non-exposed cases
#' @param neNonCases number of non-exposed non-cases
#' @return NA
#' @example
#' as.riskTable(...)
as.riskTable <- function(eCases, eNonCases, neCases, neNonCases, ...) {
	init_classes()
  if (nargs() == 4) {
		obj <- matrix(c(eCases, eNonCases, neCases, neNonCases), nrow=2)
		class(obj) <- "riskTable"
	}
	else {
		if (nargs() %% 4) {stop("Incorrect number of table entries (need multiple of 4)")}
		arglist <- as.list(match.call())[-1]
		ntable <- nargs()/4
		obj <- vector("list", length=ntable)
		class(obj) <- c("riskTable", "list")
		for (i in 1:ntable) {
			a <- 4*(i-1)
			obj[[i]] <- as.riskTable(arglist[[1+a]], arglist[[2+a]], arglist[[3+a]], arglist[[4+a]])
		}
	}
	return(obj)
}

#' Converts data to new risk table
#'
#' Creates new matrix represented as a table for a closed cohort (or cross-sectional)
#'
#' @param Cases number of cases
#' @param NonCases number of non-cases
#' @return NA
#' @example
#' as.riskTable.new(...)
as.riskTable.new <- function(cases,noncases, ...) {
  init_classes()
  if (length(cases) == 2) {
    obj <- matrix(c(as.numeric(cases[2]), as.numeric(noncases[2]), as.numeric(cases[1]), as.numeric(noncases[1])), nrow=2)
    class(obj) <- "riskTable"
  }
  else {
    if (length(cases)>2 & length(cases) %% 2) {stop("Incorrect number of table entries (need multiple of 2)")}
    ntable <- length(cases)/2
    obj <- vector("list", length=ntable)
    class(obj) <- c("riskTable", "list")
    for (i in 1:ntable) {
      a <- 2*(i-1)
      obj[[i]] <- as.riskTable(as.numeric(cases[2+a]), as.numeric(noncases[2+a]), as.numeric(cases[1+a]), as.numeric(noncases[1+a]))
    }
  }
  return(obj)
}

#' Prints risk table
#'
#' Prints new matrix represented as a table for a closed cohort (or cross-sectional)
#'
#' @param x riskTable object
#' @return NA
#' @example
#' print.riskTable(...)
print.riskTable <- function(x,...) {
  init_classes()
  if (!("list" %in% class(x))) {
		cat("\t\t\t\t\tExposed\t\tNot exposed\n")
		cat("\t", paste("Cases", format(x[1,1], width=7), format(x[1,2], width=7), sep="\t\t\t"), "\n")
		cat("\t", paste("Non-cases", format(x[2,1], width=11), "\t", format(x[2,2], width=7), sep="\t"))
	} else {
		for (i in 1:length(x)) {
			cat(paste("Stratum ", i, ":\n", sep=""))
			print.riskTable(x[[i]])
			cat("\n\n")
		}
		cat("Crude data:\n")
		print.riskTable(Reduce("+", x))
	}
}

#' Summarizes risk table
#'
#' Summarizes new matrix represented as a table for a closed cohort (or cross-sectional)
#'
#' @param x riskTable object
#' @param alpha alpha threshold
#' @return NA
#' @example
#' summary.riskTable(...)
summary.riskTable <- function(x, alpha=.05, ...) {
  init_classes()
  mult <- qnorm(1-alpha/2)
	if (!("list" %in% class(x))) {
		N1 <- as.numeric(x[1,1]+x[2,1])
		N0 <- as.numeric(x[1,2]+x[2,2])
		rat <- (x[1,1]/N1)/(x[1,2]/N0)
		dif <- x[1,1]/N1 - x[1,2]/N0
		cir <- exp(log(rat) + c(-mult, mult)*sqrt(x[2,1]/(x[1,1]*N1)+x[2,2]/(x[1,2]*N0)))
		cid <- dif + c(-mult, mult)*sqrt(x[1,1]*x[2,1]/N1^3+x[1,2]*x[2,2]/N0^3)

		M1 <- as.numeric(x[1,1]+x[1,2])
		M0 <- as.numeric(x[2,1]+x[2,2])
		teststat <- (x[1,1]-N1*M1/(M1+M0))^2/(as.numeric(M1*M0*N1*N0)/(M1+M0)^3)
		pval <- 1-pchisq(teststat, 1)

		cat("Point estimate of CIR: ", format(rat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for CIR: ", format(cir[1], digits=4, nsmall=4),
			" - ", format(cir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Point estimate of CID: ", format(dif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for CID: ", format(cid[1], digits=4, nsmall=4),
			" - ", format(cid[2], digits=4, nsmall=4), sep=""), "\n")
		if (pval < .0001) {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p =",
				format(pval, digits=4, nsmall=4)))
		}
	}
	else {
		collapseMat <- Reduce("+", x)
		a <- collapseMat[1,1]
		b <- collapseMat[1,2]
		c <- collapseMat[2,1]
		d <- collapseMat[2,2]
		N1 <- as.numeric(a+c)
		N0 <- as.numeric(b+d)
		M1 <- as.numeric(a+b)
		M0 <- as.numeric(c+d)
		cruderat <- (a/N1)/(b/N0)
		crudedif <- a/N1 - b/N0
		crudecir <- exp(log(cruderat) + c(-mult, mult)*sqrt(c/(a*N1)+d/(b*N0)))
		crudecid <- crudedif + c(-mult, mult)*sqrt(a*c/N1^3+b*d/N0^3)
		crudestat <- (a-N1*M1/(M1+M0))^2/(M1*M0*N1*N0/(M1+M0)^3)
		crudepval <- 1-pchisq(crudestat, 1)

		cirnum <- sum(unlist(lapply(x, function(y) return(y[1,1]*(y[1,2]+y[2,2])/sum(y)))))
		cirden <- sum(unlist(lapply(x, function(y) return(y[1,2]*(y[1,1]+y[2,1])/sum(y)))))
		cirmh <- cirnum/cirden
		A <- sum(unlist(lapply(x, function(y) return(((y[1,1]+y[1,2])*(y[1,1]+y[2,1])*(y[1,2]+y[2,2])-y[1,1]*y[1,2]*sum(y))/sum(y)^2))))
		B <- sum(unlist(lapply(x, function(y) return(y[1,1]*(y[1,2]+y[2,2])/sum(y)))))
		CC <- sum(unlist(lapply(x, function(y) return(y[1,2]*(y[1,1]+y[2,1])/sum(y)))))
		varln <- A/(B*CC)
		mhcir <- exp(log(cirmh) + c(-mult, mult)*sqrt(varln))

		wgts <- unlist(lapply(x, function(y) return((y[1,1]+y[2,1])^3*(y[1,2]+y[2,2])^3/((y[1,2]+y[2,2])^3*y[1,1]*y[2,1]+(y[1,1]+y[2,1])^3*y[1,2]*y[2,2]))))
		irds <- unlist(lapply(x, function(y) return(y[1,1]/(y[1,1]+y[2,1])-y[1,2]/(y[1,2]+y[2,2]))))
		ivdif <- sum(wgts*irds)/sum(wgts)
		ivcid <- ivdif + c(-mult, mult)*sqrt(1/sum(wgts))
		stratnum <- (a - sum(unlist(lapply(x, function(y) return((y[1,1]+y[2,1])*(y[1,1]+y[1,2])/sum(y))))))^2
		stratden <- sum(unlist(lapply(x, function(y) return((y[1,1]+y[1,2])*(y[2,1]+y[2,2])*(y[1,1]+y[2,1])*(y[1,2]+y[2,2])/sum(y)^3))))
		stratstat <- stratnum/stratden
		stratpval <- 1-pchisq(stratstat, 1)

		rhomstat <- sum(unlist(lapply(x, function(y) return((log(y[1,1]*(y[1,2]+y[2,2])/y[1,2]/(y[1,1]+y[2,1]))-log(cirmh))^2/((y[2,1]/y[1,1]/(y[1,1]+y[2,1]))+y[2,2]/y[1,2]/(y[1,2]+y[2,2]))))))
		rhompval <- 1-pchisq(rhomstat, length(x)-1)
		dhomstat <- sum(unlist(lapply(x, function(y) return((y[1,1]/(y[1,1]+y[2,1])-y[1,2]/(y[1,2]+y[2,2])-ivdif)^2/(y[1,1]*y[2,1]/(y[1,1]+y[2,1])^3+y[1,2]*y[2,2]/(y[1,2]+y[2,2])^3)))))
		dhompval <- 1-pchisq(dhomstat, length(x)-1)

		cat("Crude estimate of CIR: ", format(cruderat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for crude CIR: ", format(crudecir[1], digits=4,
			nsmall=4), " - ", format(crudecir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Crude estimate of CID: ", format(crudedif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for crude CID: ", format(crudecid[1], digits=4,
			nsmall=4), " - ", format(crudecid[2], digits=4, nsmall=4), sep=""), "\n")
		if (crudepval < .0001) {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p =",
				format(crudepval, digits=4, nsmall=4)))
		}
		cat("\n\nMantel-Haenszel estimate of CIR: ", format(cirmh, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for MH CIR: ", format(mhcir[1], digits=4,
			nsmall=4), " - ", format(mhcir[2], digits=4, nsmall=4), sep=""), "\n")
		cat("Inverse-variance estimate of CID: ", format(ivdif, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for IV CID: ", format(ivcid[1], digits=4,
			nsmall=4), " - ", format(ivcid[2], digits=4, nsmall=4), sep=""), "\n")
		if (stratpval < .0001) {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4),
				"with p < 0.0001"))
		} else {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4), "with p =",
				format(stratpval, digits=4, nsmall=4)))
		}

		cat(paste("\n\nTest of CIR homogeneity: X^2 =", format(rhomstat, digits=4, nsmall=4), "with p =",
				format(rhompval, digits=4, nsmall=4)))
		cat(paste("\nTest of CID homogeneity: X^2 =", format(dhomstat, digits=4, nsmall=4), "with p =",
				format(dhompval, digits=4, nsmall=4)))
	}
}

###########################################################################################
# methods for case-control data

#' Converts data into ca-co format
#'
#' Creates new matrix represented as a table for a case control study
#'
#' @param eCases number of exposed cases
#' @param eControls number of exposes controls
#' @param neCases number of non-exposed cases
#' @param neControls number of non-exposed controls
#' @return NA
#' @example
#' as.ccTable(...)
as.ccTable <- function(eCases, eControls, neCases, neControls, ...) {
	init_classes()
  if (nargs() == 4) {
		obj <- matrix(c(eCases, eControls, neCases, neControls), nrow=2)
		class(obj) <- "ccTable"
	}
	else {
		if (nargs() %% 4) {stop("Incorrect number of table entries (need multiple of 4)")}
		arglist <- as.list(match.call())[-1]
		ntable <- nargs()/4
		obj <- vector("list", length=ntable)
		class(obj) <- c("ccTable", "list")
		for (i in 1:ntable) {
			a <- 4*(i-1)
			obj[[i]] <- as.ccTable(arglist[[1+a]], arglist[[2+a]], arglist[[3+a]], arglist[[4+a]])
		}
	}
	return(obj)
}

#' Converts data into ca-co format (updated)
#'
#' Creates new matrix represented as a table for a case control study
#'
#' @param Cases number of cases
#' @param Controls number of controls
#' @return NA
#' @example
#' as.ccTable.new(...)
as.ccTable.new <- function(cases, controls, ...) {
  init_classes()
  if (length(cases) == 2) {
    obj <- matrix(c(as.numeric(cases[2]), as.numeric(controls[2]), as.numeric(cases[1]), as.numeric(controls[1])), nrow=2)
    class(obj) <- "ccTable"
  }
  else {
    if (length(cases)>2 & length(cases) %% 2) {stop("Incorrect number of table entries (need multiple of 2)")}
    ntable <-length(cases)/2
    obj <- vector("list", length=ntable)
    class(obj) <- c("ccTable", "list")
    for (i in 1:ntable) {
      a <- 2*(i-1)
      obj[[i]] <- as.ccTable(as.numeric(cases[2+a]), as.numeric(controls[2+a]), as.numeric(cases[1+a]), as.numeric(controls[1+a]))
    }
  }
  return(obj)
}

#' Prints data from ca-co format
#'
#' Creates new matrix represented as a table for a case control study
#'
#' @param x ccTable object
#' @return NA
#' @example
#' as.ccTable(...)
print.ccTable <- function(x,...) {
  init_classes()
	if (!("list" %in% class(x))) {
		cat("\t\t\t\t\tExposed\t\tNot exposed\n")
		cat("\t", paste("Cases", format(x[1,1], width=7), format(x[1,2], width=7), sep="\t\t\t"), "\n")
		cat("\t", paste("Controls", format(x[2,1], width=11), "\t", format(x[2,2], width=7), sep="\t"))
	} else {
		for (i in 1:length(x)) {
			cat(paste("Stratum ", i, ":\n", sep=""))
			print.ccTable(x[[i]])
			cat("\n\n")
		}
		cat("Crude data:\n")
		print.ccTable(Reduce("+", x))
	}
}

#' Summarizes data from ca-co format
#'
#' Creates new matrix represented as a table for a case control study
#'
#' @param x ccTable object
#' @param alpha alpha threshold
#' @return NA
#' @example
#' as.ccTable(...)
summary.ccTable <- function(x, alpha=.05, ...) {
  init_classes()
	mult <- qnorm(1-alpha/2)
	if (!("list" %in% class(x))) {
		rat <- (x[1,1]/x[2,1])/(x[1,2]/x[2,2])
		cir <- exp(log(rat) + c(-mult, mult)*sqrt(1/x[1,1]+1/x[1,2]+1/x[2,1]+1/x[2,2]))

		N1 <- x[1,1]+x[2,1]
		N0 <- x[1,2]+x[2,2]
		M1 <- x[1,1]+x[1,2]
		M0 <- x[2,1]+x[2,2]
		tot <- sum(x)
		teststat <- (x[1,1]-N1*M1/tot)^2/(M1*M0*N1*N0/(tot^2*(tot-1)))
		pval <- 1-pchisq(teststat, 1)

		cat("Point estimate of OR: ", format(rat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for OR: ", format(cir[1], digits=4, nsmall=4),
			" - ", format(cir[2], digits=4, nsmall=4), sep=""), "\n")
		if (pval < .0001) {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("H0 test: X^2 =", format(teststat, digits=4, nsmall=4), "with p =",
				format(pval, digits=4, nsmall=4)))
		}
	}
	else {
		collapseMat <- Reduce("+", x)
		a <- collapseMat[1,1]
		b <- collapseMat[1,2]
		c <- collapseMat[2,1]
		d <- collapseMat[2,2]
		cruderat <- (a*d)/(b*c)
		crudecir <- exp(log(cruderat) + c(-mult, mult)*sqrt(1/a+1/b+1/c+1/d))
		N1 <- a+c
		N0 <- b+d
		M1 <- a+b
		M0 <- c+d
		tot <- sum(collapseMat)
		crudestat <- (a-N1*M1/tot)^2/(M1*M0*N1*N0/(tot^2*(tot-1)))
		crudepval <- 1-pchisq(crudestat, 1)

		ornum <- sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]/sum(y)))))
		orden <- sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]/sum(y)))))
		ormh <- ornum/orden
		num1 <- sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]*(y[1,1]+y[2,2])/sum(y)^2))))
		num2 <- sum(unlist(lapply(x, function(y) return((y[1,1]*y[2,2]*(y[2,1]+y[1,2])+y[1,2]*y[2,1]*(y[1,1]+y[2,2]))/sum(y)^2))))
		num3 <- sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]*(y[1,2]+y[2,1])/sum(y)^2))))
		den1 <- (sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]/(sum(y)))))))^2
		den2 <- sum(unlist(lapply(x, function(y) return(y[1,1]*y[2,2]/sum(y)))))*sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]/(sum(y))))))
		den3 <- (sum(unlist(lapply(x, function(y) return(y[1,2]*y[2,1]/sum(y))))))^2
		varln <- (num1/den1+num2/den2+num3/den3)/2
		mhcir <- exp(log(ormh) + c(-mult, mult)*sqrt(varln))

		stratnum <- (a - sum(unlist(lapply(x, function(y) return((y[1,1]+y[2,1])*(y[1,1]+y[1,2])/sum(y))))))^2
		stratden <- sum(unlist(lapply(x, function(y) return((y[1,1]+y[1,2])*(y[2,1]+y[2,2])*(y[1,1]+y[2,1])*(y[1,2]+y[2,2])/sum(y)^2/(sum(y)-1)))))
		stratstat <- stratnum/stratden
		stratpval <- 1-pchisq(stratstat, 1)

		rhomstat <- sum(unlist(lapply(x, function(y) return((log(y[1,1]*y[2,2]/y[1,2]/y[2,1])-log(ormh))^2/(1/y[1,1]+1/y[1,2]+1/y[2,1]+1/y[2,2])))))
		rhompval <- 1-pchisq(rhomstat, length(x)-1)

		cat("Crude estimate of OR: ", format(cruderat, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for crude OR: ", format(crudecir[1], digits=4,
			nsmall=4), " - ", format(crudecir[2], digits=4, nsmall=4), sep=""), "\n")
		if (crudepval < .0001) {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p < 0.0001"))
		} else {
			cat(paste("Crude H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p =",
				format(crudepval, digits=4, nsmall=4)))
		}
		cat("\n\nMantel-Haenszel estimate of OR: ", format(ormh, digits=4, nsmall=4), "\n")
		cat(paste((1-alpha)*100, "% confidence interval for MH OR: ", format(mhcir[1], digits=4,
			nsmall=4), " - ", format(mhcir[2], digits=4, nsmall=4), sep=""), "\n")
		if (stratpval < .0001) {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4),
				"with p < 0.0001"))
		} else {
			cat(paste("Stratified H0 test: X^2 =", format(stratstat, digits=4, nsmall=4), "with p =",
				format(stratpval, digits=4, nsmall=4)))
		}

		cat(paste("\n\nTest of OR homogeneity: X^2 =", format(rhomstat, digits=4, nsmall=4), "with p =",
				format(rhompval, digits=4, nsmall=4)))
	}
}

#' Converts data into matched ca-co format
#'
#' Creates new matrix represented as a table for a case control study
#'
#' @param matchRatio ratio used to calculate number of controls chosen for each case
#' @return NA
#' @example
#' as.mccTable(...)
as.mccTable <- function(matchRatio, ...) {
  init_classes()
	arglist <- as.list(match.call())[-c(1,2)]
	if (length(arglist) != 2*(matchRatio+1)) {stop("Incorrect number of table entries")}
	obj <- matrix(unlist(arglist), nrow=2, byrow=TRUE)
	class(obj) <- "mccTable"

	return(obj)
}

#' Prints data into matched ca-co format
#'
#' Prints new matrix represented as a table for a case control study
#'
#' @param x mccTable object
#' @return NA
#' @example
#' print.mccTable(...)
print.mccTable <- function(x,...) {
  init_classes()
	cat("Matched data:\n")
	cat("\t\t\t\t  Exposed controls\n\t\t\t\t\t")
	for (i in (ncol(x)-1):0) {
		cat(paste(i, "\t\t"))
	}
	cat("\nCase exposed\t")
	for (i in 1:ncol(x)) {
		cat(paste(format(x[1,i], width=5), "\t"))
	}
	cat("\nCase unexposed\t")
	for (i in 1:ncol(x)) {
		cat(paste(format(x[2,i], width=5), "\t"))
	}

	crudeData <- matrix(0, nrow=2, ncol=2)
	crudeData[1,1] <- apply(x, 1, sum)[1]
	crudeData[1,2] <- apply(x, 1, sum)[2]
	for (i in 1:ncol(x)) {
		crudeData[2,1] <- crudeData[2,1] + (ncol(x)-i)*(x[1,i]+x[2,i])
		crudeData[2,2] <- crudeData[2,2] + (i-1)*(x[1,i]+x[2,i])
	}
	cat("\n\nCrude data:\n")
	print.ccTable(crudeData)
}

#' Summarizes data into matched ca-co format
#'
#' Summarizes new matrix represented as a table for a case control study
#'
#' @param x mccTable object
#' @param alpha alpha threshold
#' @return NA
#' @example
#' print.mccTable(...)
summary.mccTable <- function(x, alpha=.05, ...) {
  init_classes()
	mult <- qnorm(1-alpha/2)

	crudeData <- matrix(0, nrow=2, ncol=2)
	crudeData[1,1] <- apply(x, 1, sum)[1]
	crudeData[1,2] <- apply(x, 1, sum)[2]
	a <- rowSums(x)[1] - x[1,1]
	G <- H <- GP <- GQHP <- HQ <- Ea <- Vara <- 0
	f10 <- f01 <- 0
	for (i in 1:ncol(x)) {
		crudeData[2,1] <- crudeData[2,1] + (ncol(x)-i)*(x[1,i]+x[2,i])
		crudeData[2,2] <- crudeData[2,2] + (i-1)*(x[1,i]+x[2,i])
		f10 <- f10 + (i-1)*x[1,i]
		f01 <- f01 + (ncol(x)-i)*x[2,i]
		G <- G + (i-1)*x[1,i]/ncol(x)
		H <- H + (ncol(x)-i)*x[2,i]/ncol(x)
		GP <- GP + (i-1)*i*x[1,i]/(ncol(x))^2
		GQHP <- GQHP + (i-1)*(ncol(x)-i)*(x[1,i]+x[2,i])/(ncol(x))^2
		HQ <- HQ + (ncol(x)-i)*(ncol(x)-i+1)*x[2,i]/(ncol(x))^2
		if(i < ncol(x)) {
			Ea <- Ea + ((x[1,i+1]+x[2,i])*(ncol(x)-i))/ncol(x)
			Vara <- Vara + ((x[1,i+1]+x[2,i])*(ncol(x)-i)*i)/(ncol(x))^2
		}
	}

	crudeOR <- (crudeData[1,1]/crudeData[2,1])/(crudeData[1,2]/crudeData[2,2])
	crudeORci <- exp(log(crudeOR) + c(-mult, mult)*sqrt(1/crudeData[1,1]+1/crudeData[1,2]+1/crudeData[2,1]+1/crudeData[2,2]))
	N1 <- crudeData[1,1]+crudeData[2,1]
	N0 <- crudeData[1,2]+crudeData[2,2]
	M1 <- crudeData[1,1]+crudeData[1,2]
	M0 <- crudeData[2,1]+crudeData[2,2]
	tot <- sum(crudeData)
	crudestat <- (crudeData[1,1]-N1*M1/tot)^2/(M1*M0*N1*N0/(tot^2*(tot-1)))
	crudepval <- 1-pchisq(crudestat, 1)

	mhOR <- f10/f01
	mhORvar <- (GP/G^2+GQHP/G/H+HQ/H^2)/2
	mhORci <- exp(log(mhOR) + c(-mult, mult)*sqrt(mhORvar))
	mhstat <- (a-Ea)^2/Vara
	mhpval <- 1-pchisq(mhstat, 1)

	cat("Crude estimate of OR: ", format(crudeOR, digits=4, nsmall=4), "\n")
	cat(paste((1-alpha)*100, "% confidence interval for crude OR: ", format(crudeORci[1], digits=4, nsmall=4), " - ", format(crudeORci[2], digits=4, nsmall=4), sep=""), "\n")
	if (crudepval < .0001) {
		cat(paste("H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p < 0.0001"))
	} else {
	cat(paste("H0 test: X^2 =", format(crudestat, digits=4, nsmall=4), "with p =",
		format(crudepval, digits=4, nsmall=4)))
	}
	cat("\n\nMatched MH OR estimate: ", format(mhOR, digits=4, nsmall=4), "\n")
	cat(paste((1-alpha)*100, "% confidence interval for MH OR: ", format(mhORci[1], digits=4, nsmall=4), " - ", format(mhORci[2], digits=4, nsmall=4), sep=""), "\n")
	if (mhpval < .0001) {
		cat(paste("H0 test: X^2 =", format(mhstat, digits=4, nsmall=4), "with p < 0.0001"))
	} else {
	cat(paste("H0 test: X^2 =", format(mhstat, digits=4, nsmall=4), "with p =",
		format(mhpval, digits=4, nsmall=4)))
	}
}

#' Convert wide data format into long data format
#'
#' Takes in wide data format for closed cohort and case control studies and converts it into long format
#'
#' @param x wide data format data frame
#' @return long data format data frame
#' @example
#' as.data.frame.riskTable(x)
as.data.frame.riskTable <- function(x, ...) {
	nstrat <- length(x)
	n <- sum(unlist(lapply(x, sum)))
	out <- data.frame(Y=rep(0, n), X=rep(0, n), C=rep(LETTERS[1:nstrat], length.out=n))
	temp <- 1
	for (i in 1:nstrat) {
		ncases <- x[[i]][1,1]+x[[i]][1,2]
		out$Y[temp:(temp+ncases-1)] <- 1
		out$X[temp:(temp+x[[i]][1,1]-1)] <- 1
		out$X[(temp+ncases):(temp+ncases+x[[i]][2,1]-1)] <- 1
		out$C[temp:(temp+sum(x[[i]])-1)] <- LETTERS[i]
		temp <- temp + sum(x[[i]])
	}
	return(out)
}

#' Convert wide data format into long data format
#'
#' Takes in wide data format for closed cohort and case control studies and converts it into long format
#'
#' @param x wide data format data frame
#' @return long data format data frame
#' @example
#' as.data.frame.ccTable(x)
as.data.frame.ccTable <- function(x, ...) {
  nstrat <- length(x)
  n <- sum(unlist(lapply(x, sum)))
  out <- data.frame(Y=rep(0, n), X=rep(0, n), C=rep(LETTERS[1:nstrat], length.out=n))
  temp <- 1
  for (i in 1:nstrat) {
    ncases <- x[[i]][1,1]+x[[i]][1,2]
    out$Y[temp:(temp+ncases-1)] <- 1
    out$X[temp:(temp+x[[i]][1,1]-1)] <- 1
    out$X[(temp+ncases):(temp+ncases+x[[i]][2,1]-1)] <- 1
    out$C[temp:(temp+sum(x[[i]])-1)] <- LETTERS[i]
    temp <- temp + sum(x[[i]])
  }
  return(out)
}
