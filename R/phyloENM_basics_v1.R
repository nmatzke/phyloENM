################################################# 
# Utility functions
################################################# 

#' logit: transform probabilities to log-odds-ratio (wraps R base stats::qlogis)
#' 
#' Logit converts probability to log-odds ratio
#' logit = function(prob) { log(prob) - log(1-prob) }
#' 
#' This function does the logit transform, converting a probability into a 
#' log-odds ratio.  The inverse operation is sometimes called
#' "inverse logit" or "expit". 
#'
#' This is just a wrapper for \code{stats::qlogis}, which calls a C function \code{C_qlogis}, which is 
#' presumably the fastest way to do it.
#' 
#' The logit operations is useful for building linear models (logistic regression, 
#' or other models) that attempt to 
#' predict probabilities, because predicting logit(probability) avoids many of 
#' the issues that occur when trying to predict just probability. I.e., probabilities
#' are bounded at 0 and 1. Logit(prob), on the other hand, ranges from -Inf to +Inf,
#' a probability of 0.5 equates to logit(prob) of 0, etc., so a model of logit(prob)
#' can look at how various predictors increase or decrease the probability.
#' 
#' Note from ?stats:qlogis: "qlogis(p) is the same as the well-known 'logit' function, 
#' \code{logit(p) = log(p/(1-p))}." 
#'
#' The inputs, below, are copied from \code{stats::plogis} (the parameter descriptions are copied, but 
#' somewhat clarified here). 
#'
#' @param prob A vector of probabilities.
#' @param location This is the mean, sometimes written as \emph{m} or \emph{mu}. I.e., the center of the distribution.
#' @param scale	The scale parameter (often written as s), is "proportional to the standard deviation" (\url{https://en.wikipedia.org/wiki/Logistic_distribution}).
#' @param log.p Logical. If \code{TRUE}, the input probabilities \emph{prob} are being input as log(\emph{prob}). Default is \code{FALSE}.
#' @param lower.tail	Logical; if \code{TRUE} (default), probabilities are \emph{P[X ≤ x]}, otherwise, <i>P[X > x]</i>.
#' @return On defaults, this returns the \code{log_odds_ratio} for the input \emph{prob}. Defaults or any other inputs give the quantile function of the specified logistic distribution.
#' @export
#' @seealso \code{\link{stats::dlogis}}, \code{\link{stats::plogis}}, \code{\link{stats::qlogis}}, \code{\link{stats::rlogis}}, all part of \code{\link{stats::Logistic}}.
#' @note The inverse operation (log-odds ratio to probability) is known as "inverse logit", "expit", or the \code{stats::plogis} function.
#' @author Nicholas J. Matzke \email{nickmatzke.ncse@@gmail.com} 
#' @references
#' \url{https://en.wikipedia.org/wiki/Logistic_distribution}
#' \url{https://github.com/nmatzke/phyloENM}
#' @examples
#' test=1
#' 
#' 
#' Titanic survivor data
#' # Original Source: http://amunategui.github.io/binary-outcome-modeling/
#' remotedir = "http://math.ucdenver.edu/RTutorial/"
#' # Local source:
#' tmpdir = system.file("extdata/", package="phyloENM")
#' if (tmpdir == "")
#' 	{
#' 	# Matzke development source
#' 	tmpdir = "/GitHub/phyloENM/inst/extdata/data/"
#' 	} else {
#' 	if (tmpdir == "")
#' 		{
#' 		tmpdir = remotedir
#' 		}
#' 	}
#' fn = paste0(tmpdir, "titanic.txt")
#' titanic = read.table(file=fn, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#' head(titanic)
#' 
#' # Fraction survived total
#' sum(titanic$Survived) # 450
#' nrow(titanic)					# 1313
#' sum(titanic$Survived) / nrow(titanic)
#' # 0.3427266
#' 
#' # Convert to odds ratio
#' prob = 0.3427266
#' oddsratio = prob / (1-prob)
#' oddsratio
#' # 0.5214369
#' 
#' # Convert to log odds ratio
#' log(oddsratio)
#' # -0.651167
#' 
#' # Compare to logit()
#' logit(prob=prob)
#' 
#' 
#' # Predict survival by Intercept-only
#' fit = glm(Survived ~1, family=binomial(link="logit"), data=titanic)
#' coef(fit)
#' # (Intercept) 
#' #  -0.6511671 
#' exp(coef(fit))
#' # (Intercept) 
#' #   0.5214368 
#' oddsratio = exp(coef(fit))
#' oddsratio / (1+oddsratio)
#' # (Intercept) 
#' #   0.3427266 
#' 
#' 
#' # Predict survival by age
#' fit = glm(Survived ~1+Age, family=binomial(link="logit"), data=titanic)
#' 
#' # (natural) Log odds ratios
#' coef(fit)
#' # (Intercept)         Age 
#' # -0.08142783 -0.00879462 
#' 
#' # Relative odds ratios
#' exp(coef(fit))
#' # (Intercept)         Age 
#' #   0.9217992   0.9912439 
#' oddsratio = exp(coef(fit))
#' # (Intercept)         Age 
#' #   0.9217992   0.9912439 
#' 
#' # Plot predicted log(oddsratios)
#' plot(x=titanic$Age[!is.na(titanic$Age)], y=predict(fit))
#' 
#' # Plot predicted probabilities
#' oddsratio = exp(predict(fit))
#' predicted_probs = oddsratio / (1+oddsratio)
#' plot(x=titanic$Age[!is.na(titanic$Age)], y=predicted_probs, ylim=c(-0.1,1.1))
#' points(x=titanic$Age[!is.na(titanic$Age)], y=jitter(titanic$Survived[!is.na(titanic$Age)], factor=0.2))
#' 
#' # Not really significant:
#' summary(fit)
#' 
logit <- function(prob, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
	{
	log_odds_ratio = stats::qlogis(p=prob, location=location, scale=scale, lower.tail=lower.tail, log.p=log.p)
	return(log_odds_ratio)
	} # END logit


# qlogis(p) is the same as the well known ‘logit’ function, 
# logit(p) = log(p/(1-p)), and plogis(x) has consequently 
# been called the ‘inverse logit’.

# Expit = inverse logit
# Converts log-odds ratio to probability
# aka inverse logit

#' invlogit: transform log-odds-ratio to probabilities (wraps R base stats::plogis)
#' 
#' \code{invlogit} converts log-odds ratio to probability. Also known as \code{expit} or \code{stats::plogis}.
#' logit converts probability to log-odds ratio
#' logit = function(prob) { log(prob) - log(1-prob) }
#' expit = function(logoddsratio) function(eta) {1/(1+exp(-logoddsratio))}
#' 
#' This function does inverse of the logit transform, converting a  
#' log-odds ratio into a probability.  The inverse operation is sometimes called
#' "inverse logit" or "expit". 
#'
#' This is just a wrapper for \code{stats::plogis}, which calls a C function \code{C_plogis}, which is 
#' presumably the fastest way to do it.
#' 
#' The inputs, below, are copied from \code{stats::plogis} (the parameter descriptions are copied, but 
#' somewhat clarified here). 
#'
#' @param logoddsratio A vector of log(odds ratios) (or quantiles from the logistic distribution generally).
#' @param location This is the mean, sometimes written as \emph{m} or \emph{mu}. I.e., the center of the distribution.
#' @param scale	The scale parameter (often written as s), is "proportional to the standard deviation" (\url{https://en.wikipedia.org/wiki/Logistic_distribution}).
#' @param log.p Logical. If \code{TRUE}, the input probabilities \emph{prob} are being input as log(\emph{prob}). Default is \code{FALSE}.
#' @param lower.tail	Logical; if \code{TRUE} (default), probabilities are \emph{P[X ≤ x]}, otherwise, <i>P[X > x]</i>.
#' @return On defaults, this returns the \code{log_odds_ratio} for the input \emph{prob}. Defaults or any other inputs give the quantile function of the specified logistic distribution.
#' @export
#' @seealso \code{\link{stats::dlogis}}, \code{\link{stats::plogis}}, \code{\link{stats::qlogis}}, \code{\link{stats::rlogis}}, all part of \code{\link{stats::Logistic}}.
#' @note The inverse operation (log-odds ratio to probability) is known as "inverse logit", "expit", or the \code{stats::plogis} function.
#' @author Nicholas J. Matzke \email{nickmatzke.ncse@@gmail.com} 
#' @references
#' \url{https://en.wikipedia.org/wiki/Logistic_distribution}
#' \url{https://github.com/nmatzke/phyloENM}
#' @examples
#' test=1
#' 
#' # Run some inverse logits
#' invlogit(logoddsratio=-10)
#' invlogit(logoddsratio=-1)
#' invlogit(logoddsratio=0)
#' invlogit(logoddsratio=1)
#' invlogit(logoddsratio=10)
#' 
invlogit <- function(logoddsratio, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
	{
	prob = stats::plogis(q=logoddsratio, location=location, scale=scale, lower.tail=lower.tail, log.p=log.p)
	return(prob)
	} # END invlogit


#' expit: transform log-odds-ratio to probabilities (wraps R base stats::plogis)
#' 
#' \code{expit} converts log-odds ratio to probability. Also known as \code{expit} or \code{stats::plogis}.
#' logit converts probability to log-odds ratio
#' logit = function(prob) { log(prob) - log(1-prob) }
#' expit = function(logoddsratio) function(eta) {1/(1+exp(-logoddsratio))}
#' 
#' This function does inverse of the logit transform, converting a  
#' log-odds ratio into a probability.  The inverse operation is sometimes called
#' "inverse logit" or "expit". 
#'
#' This is just a wrapper for \code{stats::plogis}, which calls a C function \code{C_plogis}, which is 
#' presumably the fastest way to do it.
#' 
#' The inputs, below, are copied from \code{stats::plogis} (the parameter descriptions are copied, but 
#' somewhat clarified here). 
#'
#' @param logoddsratio A vector of log(odds ratios) (or quantiles from the logistic distribution generally).
#' @param location This is the mean, sometimes written as \emph{m} or \emph{mu}. I.e., the center of the distribution.
#' @param scale	The scale parameter (often written as s), is "proportional to the standard deviation" (\url{https://en.wikipedia.org/wiki/Logistic_distribution}).
#' @param log.p Logical. If \code{TRUE}, the input probabilities \emph{prob} are being input as log(\emph{prob}). Default is \code{FALSE}.
#' @param lower.tail	Logical; if \code{TRUE} (default), probabilities are \emph{P[X ≤ x]}, otherwise, <i>P[X > x]</i>.
#' @return On defaults, this returns the \code{log_odds_ratio} for the input \emph{prob}. Defaults or any other inputs give the quantile function of the specified logistic distribution.
#' @export
#' @seealso \code{\link{stats::dlogis}}, \code{\link{stats::plogis}}, \code{\link{stats::qlogis}}, \code{\link{stats::rlogis}}, all part of \code{\link{stats::Logistic}}.
#' @note The inverse operation (log-odds ratio to probability) is known as "inverse logit", "expit", or the \code{stats::plogis} function.
#' @author Nicholas J. Matzke \email{nickmatzke.ncse@@gmail.com} 
#' @references
#' \url{https://en.wikipedia.org/wiki/Logistic_distribution}
#' \url{https://github.com/nmatzke/phyloENM}
#' @examples
#' test=1
#' 
#' # Run some inverse logits
#' expit(logoddsratio=-10)
#' expit(logoddsratio=-1)
#' expit(logoddsratio=0)
#' expit(logoddsratio=1)
#' expit(logoddsratio=10)
#' 
expit <- function(logoddsratio, location=0, scale=1, lower.tail=TRUE, log.p=FALSE)
	{
	prob = stats::plogis(q=logoddsratio, location=location, scale=scale, lower.tail=lower.tail, log.p=log.p)
	return(prob)
	} # END expit




#######################################################
# Calculate the base occurrence rate
#######################################################


