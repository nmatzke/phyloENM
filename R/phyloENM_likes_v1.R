# phyloENM_likes: likelihoods of presence/absence data

#' 
#' Likelihood of occurrence data assuming a spline model with 1 breakpoint, 
#' at the optimum value
#' @param predictor The values of the niche value being used for prediction
#' @param response The occurence data, 0=absent, 1=present
#' @param opt    The position of the optimum point on the x-axis (e.g. temperature). 
#'               typically inferred.
#' @param slope1 The slope, in logit-space, leading up to the optimum point. Must be >=0.
#' @param slope2 The slope, in logit-space, leading down from the optimum point. Must be <0.
#' @param intercept logit(prob(occ)) = average occurence probability in the data, in 
#'               log-odds form = intercept of a logistic regression. This is 
#'               auto-calculated by default, making the spline relative to the 
#'               mean occurrence probability.
#' @param maxadd The log-odds to add to the intercept, to form the maximum value in 
#'               the spline curve. Default is logit(0.99)=4.59512; but could be inferred.
#' 

like_1spline <- function(predictor, response, opt, intercept=qlogis(mean(response)), maxadd=qlogis(0.99), slope1=1, slope2=-1)
	{
	# The optimum (highest y-value), in logit space
	intensity_at_opt = intercept + maxadd
	
	# Calculate the predicted probability of occurrence
	logit_prob_occurrence = rep(NA, times=length(response))
	
	# Batch filling in of logit(prob)
	TF = predictor <= opt
	dist_below_x = opt - predictor[TF]
	logit_prob_occurrence[TF] = intensity_at_opt - (dist_below_x * slope1)

	TF = predictor > opt
	dist_below_x = predictor[TF] - opt
	logit_prob_occurrence[TF] = intensity_at_opt + (dist_below_x * slope2)
	
	
	# The log-likelihood of the presences
	prob_of_1s = plogis(q=logit_prob_occurrence[response==1])
	lnL_presences = sum(log(prob_of_1s))

	# The log-likelihood of the absences
	prob_of_0s = 1-plogis(q=logit_prob_occurrence[response==0])
	lnL_absences = sum(log(prob_of_0s))
	
	lnL = lnL_presences + lnL_absences

	if (is.finite(lnL) == FALSE)
		{
		lnL = -1e80
		} else {
		printwait = FALSE
		
		if (printwait == TRUE)
			{
			numsecs = 20
			remainder = round((as.numeric(Sys.time())-1495000000), 0) %% numsecs
			if (remainder == 0)
				{
				cat(paste0("LnL: ", round(lnL,1), ", opt: ", round(par[1],1), ", slope1: ", round(par[2],1), ", slope2: ", round(par[3],1), "\n"))
				}
			} # END if (printresult == TRUE)
		}
		
	return(lnL)
	}


like_1spline_for_GenSA <- function(par)
	{
	lnL = like_1spline(predictor, response, opt=par[1], intercept=intercept, maxadd=maxadd, slope1=par[2], slope2=par[3])
	

	
	cat(paste0("LnL: ", lnL, "	opt: ", round(par[1],1), "	slope1: ", round(par[2],1), "	slope2: ", round(par[3],1), "\n"))
	
	if (lnL > 1e-50)
		{
		# Loads to last_lnl
		load(file="last_lnl.Rdata")
		if (lnL > last_lnl)
			{
			plot_1split_against_data(predictor, response, opt=par[1], intercept=intercept, maxadd=maxadd, slope1=par[2], slope2=par[3])
			title(paste0("lnL = ", lnL))
			}
		} # END if (lnL > 1e-50)
	
	return(-lnL)
	}


plot_1split_against_data <- function(predictor, response, opt, intercept=qlogis(mean(response)), maxadd=qlogis(0.99), slope1=1, slope2=-1)
	{
	hist_object1 = hist(predictor, breaks=200, plot=FALSE)
	presenceTF = this_species$response == 1
	hist_object2 = hist(predictor[presenceTF], breaks=hist_object1$breaks, plot=FALSE)

	# Make a plot of proportions
	#plot(hist_object1, xlab="Temperature (C)", main="Histogram of all pixels, and occurrences")
	#plot(hist_object2, col="blue", add=TRUE)

	# Presence counts vs. total counts
	hist_object3 = hist_object2
	proportion_presences = hist_object2$counts / hist_object1$counts
	hist_object3$counts = proportion_presences
	plot(hist_object3, xlab="Temperature (C)", main=NULL, ylim=c(0,1), xlim=c(-50,450))
	
	
	# Plot the model on top of this
	# The optimum (highest y-value), in logit space
	intensity_at_opt = intercept + maxadd
	
	
	# Batch filling in of logit(prob)
	xvals = seq(min(predictor), max(predictor), by=1)

	# Calculate the predicted probability of occurrence
	logit_prob_occurrence = rep(NA, times=length(xvals))
	
	TF = xvals <= opt
	dist_below_x = opt - xvals[TF]
	logit_prob_occurrence[TF] = intensity_at_opt - (dist_below_x * slope1)

	TF = xvals > opt
	dist_below_x = xvals[TF] - opt
	logit_prob_occurrence[TF] = intensity_at_opt + (dist_below_x * slope2)
	
	
	# The log-likelihood of the presences
	prob_of_1s = plogis(q=logit_prob_occurrence)
	lines(x=xvals, y=prob_of_1s, col="blue", lwd=2)
	
	}

