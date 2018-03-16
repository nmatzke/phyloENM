

# Plots the PNOs (predicted niche occupancies)
#' @param clim A data.frame, where the first column is the environmental variable.
#'             Second and later columns are the individual species.
#' @param varname The name of the climate variable. Default is "climate variable"
#' 

plot_PNOs <- function(clim, varname="climate variable")
	{
	xvals = clim[,1]
	
	colors = rainbow(n=(ncol(clim)-1))
	
	for (i in 2:ncol(clim))
		{
		yvals = clim[,i]
		titletxt = paste0("PNO for ", names(clim)[i])
		plot(xvals, yvals, xlab=varname, ylab="Predicted Niche Occupancy (PNO)", pch=".", col="white")
		title(titletxt)
		lines(xvals, yvals, col=colors[i-1])
		points(xvals, yvals, cex=0.5, col=colors[i-1])
		} # END for (i in 2:ncol(clim))
	} # END plot_PNOs <- function(clim, varname="climate variable")

extract_id_from_list <- function(x, id)
	{
	x <- x[[id]]
	}

# For a particular phylogeny and e.g. 100 samples from a PNO, calculate
# ancestral niche estimates under Brownian Motion
# For application to a list of phylogenies
anc_clim_core <- function(phy, x, returnwhat="original", CIs=c(0.025,0.975))
	{
	defaults='
	returnwhat="original"
	'
	
	x <- x[, match(phy$tip.label, colnames(x))]
	if (method == "GLS")
		{
		# Perform ACE 100 times on the 100 samples
		out <- apply(x, 1, FUN=ace, phy=phy, type="continuous", method="GLS", CI=TRUE, model="BM", kappa=1, corStruct=corBrownian(1,phy))
		}
		
	if (method == "ML")
		{
		out <- apply(x, 1, FUN=ace, phy=phy, type="continuous", method="ML", CI=FALSE, model="BM", kappa=1)
		}
	
	id <- which(c("GLS", "ML") %in% method)
	
	# The 100 estimates of the MEAN ancestral niche value
	ancnode_means = lapply(X=out, FUN=extract_id_from_list, id=id)
	ancnode_means = matrix(unlist(ancnode_means), nrow=length(ancnode_means), ncol=phy$Nnode, byrow=TRUE)
	
	# Take the MEAN of the 100 estimates
	means_of_ancnode_means = apply(ancnode_means, 2, mean)
	
	# The original output was just the means of mean ancestral node values
	# for the 100 samples from the PNOs
	if (returnwhat == "original")
		{
		return(means_of_ancnode_means)
		}
	
	# Otherwise, make a list of various interesting outputs!
	if (returnwhat == "samples")
		{
		outlist = list()
		outlist$ancnode_means = ancnode_means
		outlist$means_of_ancnode_means = means_of_ancnode_means
		outlist$sds_of_ancnode_means = apply(ancnode_means, 2, sd)
		outlist$ci025s_of_ancnode_means = apply(ancnode_means, 2, quantile, probs=CIs[1])
		outlist$ci975s_of_ancnode_means = apply(ancnode_means, 2, quantile, probs=CIs[2])
		outlist$ci025s_of_ancnode_means
		outlist$ci975s_of_ancnode_means
		return(outlist)
		}

	return(stop("Shouldn't get here; error in 'returnwhat' input."))
	} # END anc_clim_core <- function(phy, x, returnwhat="original")



anc_clim_wUncertainty <- function (target, posterior=NULL, pno, n=100, method="GLS", returnwhat="samples", CIs=c(0.025,0.975))
	{
	defaults='
	target=tree; posterior=NULL; pno=clim; n=100; method="GLS"; returnwhat="original"; CIs=c(0.025,0.975)
	'
	
	# Make 2 trees in the posterior, so that lapply() works later;
	# kind of a hack!
	if (is.null(posterior))
		{
		posterior <- list(target, target)
		}
	
	# Make a samples matrix, n=100 replicates for 5 taxa
	ntax <- dim(pno)[2] - 1
	x <- matrix(nrow=n, ncol=ntax)
	colnames(x) <- colnames(pno)[-1]
	
	# Go through the PNOs of each species, collect 100 samples from each species' PNO
	for (i in 2:dim(pno)[2])
		{
		# Sample from the PNO, according to the density (kind of weird sampling from 
		# the empirical distribution of the PNO, rather than the smoothed curve!)
		x[, i-1] = sample(pno[, 1], size=n, replace=TRUE, prob=pno[, i])
		} # END for (i in 2:dim(pno)[2])
	
	# For each tree, get the ancestral climate estimates
	if (returnwhat == "original")
		{
		list_of_means_of_ancnode_means <- lapply(X=posterior, FUN=anc_clim_core, x=x, returnwhat="original", CIs=CIs)
		matrix_of_means_of_ancnode_means <- matrix(unlist(list_of_means_of_ancnode_means), nrow=length(list_of_means_of_ancnode_means), ncol=posterior[[1]]$Nnode, byrow=TRUE)
		} # END if (returnwhat == "original")

	if (returnwhat == "samples")
		{
		outlist <- lapply(X=posterior, FUN=anc_clim_core, x=x, returnwhat="samples", CIs=CIs)
		list_of_means_of_ancnode_means = lapply(X=outlist, FUN=function(x) {return(x$means_of_ancnode_means)})
		matrix_of_means_of_ancnode_means <- matrix(unlist(list_of_means_of_ancnode_means), nrow=length(list_of_means_of_ancnode_means), ncol=posterior[[1]]$Nnode, byrow=TRUE)
		
		list_of_ci025s_of_ancnode_means = lapply(X=outlist, FUN=function(x) {return(x$ci025s_of_ancnode_means)})
		matrix_of_ci025s_of_ancnode_means <- matrix(unlist(list_of_ci025s_of_ancnode_means), nrow=length(list_of_ci025s_of_ancnode_means), ncol=posterior[[1]]$Nnode, byrow=TRUE)
		list_of_ci975s_of_ancnode_means = lapply(X=outlist, FUN=function(x) {return(x$ci975s_of_ancnode_means)})
		matrix_of_ci975s_of_ancnode_means <- matrix(unlist(list_of_ci975s_of_ancnode_means), nrow=length(list_of_ci975s_of_ancnode_means), ncol=posterior[[1]]$Nnode, byrow=TRUE)
		} # END if (returnwhat == "samples")
		
		
	for (i in seq(along=posterior))
		{
		# Get the indexes that match the nodes of the posterior trees to the nodes of the target tree
		id <- AC.node.trans(new=target, old=posterior[[i]], index=TRUE)
		nomatch <- which(!seq(along=id) %in% id)
		id[is.na(id)] <- nomatch
		matrix_of_means_of_ancnode_means[i, ] <- matrix_of_means_of_ancnode_means[i, id] # insert the values for the nodes that matched
		matrix_of_means_of_ancnode_means[i, nomatch] <- NA  # insert NAs for the nodes that didn't match
	
		if (returnwhat == "samples")
			{
			matrix_of_ci025s_of_ancnode_means[i, ] <- matrix_of_ci025s_of_ancnode_means[i, id] # insert the values for the nodes that matched
			matrix_of_ci025s_of_ancnode_means[i, nomatch] <- NA  # insert NAs for the nodes that didn't match
			matrix_of_ci975s_of_ancnode_means[i, ] <- matrix_of_ci975s_of_ancnode_means[i, id] # insert the values for the nodes that matched
			matrix_of_ci975s_of_ancnode_means[i, nomatch] <- NA  # insert NAs for the nodes that didn't match
			}
		} # END for (i in seq(along=posterior))
	
	# Means across tree posterior
	
	mean_of_ancnode_means <- apply(matrix_of_means_of_ancnode_means, 2, mean, na.rm=TRUE)
	if (returnwhat == "samples")
		{
		mean_of_ancnode_ci025s_means <- apply(matrix_of_ci025s_of_ancnode_means, 2, mean, na.rm=TRUE)
		mean_of_ancnode_ci975s_means <- apply(matrix_of_ci975s_of_ancnode_means, 2, mean, na.rm=TRUE)
		} # END if (returnwhat == "samples")

	tips <- apply(x, 2, mean)
	tips <- tips[match(target$tip.label, names(tips))]
	central_density_tips <- apply(x, 2, quantile, probs=CIs)
	central_density_tips <- central_density_tips[, match(target$tip.label, colnames(central_density_tips))]
	out <- c(tips, mean_of_ancnode_means)
	if (any(is.na(out)))
		{
		nas <- which(is.na(out))
		nas <- paste(nas, collapse=", ")
		stop("Reconstruction failed for nodes: ", nas)
		} # END if (any(is.na(out)))
	if (returnwhat == "original")
		{
		out <- list(tree=target, means=out, central_density_tips=central_density_tips)
		}
	if (returnwhat == "samples")
		{
		out <- list(tree=target, means=out, central_density_tips=central_density_tips, mean_of_ancnode_ci025s_means=mean_of_ancnode_ci025s_means, mean_of_ancnode_ci975s_means=mean_of_ancnode_ci975s_means)
		}
		
	return(out)
	} # END anc_clim_wUncertainty 


abb <- function(tips)
	{
	tmp_tips = unlist(strsplit(tips, ""))
	paste(head(tmp_tips, nchar), collapse="")
	} # END abb


plotAncClim_wUncertainty <- function (x, clades=NULL, col, density=TRUE, tipmode=1, nchar=3, cex=1, tipspace=NULL, cladespace=1, lwd=1, ylab="", plot_CIs=TRUE) 
	{
	tr <- x$tree
	tips <- tr$tip.label
	nbtips <- length(tips)
	if (is.null(tipspace))
		{
		tipspace = 1 - (4/nbtips)
		}
	
	clim <- x$means
	if (density) 
			central_density_tips <- x$central_density_tips
	
	if (plot_CIs == TRUE)
		{
		if (is.null(x$mean_of_ancnode_ci025s_means) == TRUE)
			{
			stop("STOP ERROR in plotAncClim_wUncertainty(): the input x, which comes from anc_clim_wUncertainty(), must have x$mean_of_ancnode_ci025s_means and x$mean_of_ancnode_ci975s_means. These are produced by running anc_clim_wUncertainty() with returnwhat='samples' instead of returnwhat='original'.")
			}
		if (is.null(x$mean_of_ancnode_ci975s_means) == TRUE)
			{
			stop("STOP ERROR in plotAncClim_wUncertainty(): the input x, which comes from anc_clim_wUncertainty(), must have x$mean_of_ancnode_ci025s_means and x$mean_of_ancnode_ci975s_means. These are produced by running anc_clim_wUncertainty() with returnwhat='samples' instead of returnwhat='original'.")
			}
		mean_of_ancnode_ci025s_means = x$mean_of_ancnode_ci025s_means
		mean_of_ancnode_ci975s_means = x$mean_of_ancnode_ci975s_means
		}
	
	if (is.null(clades))
		{
		ts <- terminal.sisters(tr)
		nts <- tr$tip.label[!tr$tip.label %in% ts]
		for (i in seq(along=nts))
			{
			ts <- rbind(ts, rep(nts[i],2))
			}
		lts <- vector(mode="list", length=dim(ts)[1])
		for (i in seq(along=lts)) lts[[i]] <- ts[i, ]
		clades <- lts
		}
	tips <- tr$tip
	
	
	tips <- sapply(tips, abb)
	nodeages <- c(rep(0, nbtips), -branching.times(tr))
	max_age <- -max(branching.times(tr))
	xy <- cbind(nodeages, clim)
	
	if (!density)
		{
		yrange <- range(xy[, 2])
		} else {
		yrange <- c(min(central_density_tips[1, ]), max(central_density_tips[2, ]))
		}
	if (tipmode == 3) 
			yrange[2] <- yrange[2] + 0.05 * diff(yrange)
	if (tipmode == 0) 
			tipspace <- 0
	tipspace <- -tipspace * max_age
	plot(x=c(max_age, tipspace), y=yrange, type="n", xlab="Time (Ma)", 
			ylab=ylab, bty="c", xaxp=c(floor(max_age), 0, 5))
	tipxpos <- max(strwidth(tips, units="f", cex=cex)) * 1:nbtips
	tipxpos <- tipxpos * cladespace
	tex <- cbind(tipxpos, xy[1:nbtips, 2])
	rownames(tex) <- tr$tip.label
	for (i in 2:length(clades))
		{
		id <- which(rownames(tex) %in% clades[[i]])
		extraspace <- (-tex[1, 1] + tex[1, 1] * 1) * (i - 1)
		tex[id, 1] <- tex[id, 1] + extraspace
		}
	rownames(tex) <- tips
	maxtip <- max(tex[, 1])
	tex[, 1] <- tex[, 1]/maxtip * tipspace
	totalspace <- max(tex[, 1])
	
	n <- noi(tr, clades) 	# noi=Number of internal?
	
	# Descendants of internal nodes
	n <- lapply(n, phyloclim::descendants, tree=tr, internal=TRUE)
	lincol <- rep("grey", dim(xy)[1])
	if (missing(col)) 
			col <- rainbow(length(n))
	for (i in seq(along=n))
		{
		lincol[n[[i]]] <- col[i]
		}
	colind <- vector()
	for (i in seq(along=tr$edge[, 1]))
		{
		ind <- tr$edge[i, ]
		xvals = xy[ind, 1]
		yvals = xy[ind, 2]
		lines(x=xvals, y=yvals, lwd=lwd, col=lincol[ind[2]])
		colind <- c(colind, ind[2])
		
		anc_nodenums_done = NULL
		if (plot_CIs == TRUE)
			{
			# The older node is the one to plot CIs on
			anc_nodenum = ind[1]
			if ((anc_nodenum %in% anc_nodenums_done) == FALSE)
				{
				colnum = anc_nodenum - nbtips
				xval = xy[anc_nodenum,1]
				yval = xy[anc_nodenum,2]
			
				CIs = c(mean_of_ancnode_ci025s_means[colnum], mean_of_ancnode_ci975s_means[colnum])
				lines(x=c(xval,xval), y=CIs, col=lincol[ind[2]], lty="dotted")
			
				# Plot tickmarks
				tick_xmin = xval - 0.01*abs(max_age)
				tick_xmax = xval + 0.01*abs(max_age)
				segments(x0=tick_xmin, x1=tick_xmax, y0=CIs[1], y1=CIs[1], col=lincol[ind[2]], lty="solid")
				segments(x0=tick_xmin, x1=tick_xmax, y0=CIs[2], y1=CIs[2], col=lincol[ind[2]], lty="solid")
			
				anc_nodenums_done = c(anc_nodenums_done, anc_nodenum)
				} # END if ((anc_nodenum %in% anc_nodenums_done) == FALSE)
			} # END if (plot_CIs == TRUE)
		} # END for (i in seq(along=tr$edge[, 1]))

	if (min(xy[, 2]) < 0 & max(xy[, 2] > 0)) 
			lines(x=c(min(xy[, 1]), 0), y=rep(0, 2), lwd=0.8 * 
					lwd, lty="12", col="gray50")
	
	# Add the CIs for the tips
	if (density & tipmode != 0)
		{
		for (i in seq(along=central_density_tips[1, ]))
			{
			xvals = rep(tex[i, 1], 2)
			yvals = central_density_tips[, i]
			lines(x=xvals, y=yvals, col=lincol[1:nbtips][i], lwd=lwd, lty="12")

			# Plot tickmarks
			tick_xmin = xvals[1] - 0.01*abs(max_age)
			tick_xmax = xvals[1] + 0.01*abs(max_age)
			segments(x0=tick_xmin, x1=tick_xmax, y0=yvals[1], y1=yvals[1], col=lincol[1:nbtips][i], lty="solid")
			segments(x0=tick_xmin, x1=tick_xmax, y0=yvals[2], y1=yvals[2], col=lincol[1:nbtips][i], lty="solid")
			}
		} # END if (density & tipmode != 0)
	
	if (tipmode != 0) 
			if (tipmode == 1) 
					text(tex[, 1], tex[, 2], rownames(tex), cex=cex, 
							adj=0.5, col=lincol[1:nbtips], font=4)
	
	if (tipmode == 2)
		{
		points(tex, col=lincol[1:nbtips], pch=18)
		id <- mean(yrange)
		id <- tex[, 2] < id
		offset <- strheight(".", cex=cex, font=4)
		text(tex[id, 1], central_density_tips[2, id] + offset, rownames(tex)[id], 
			cex=cex, adj=c(0, 0.5), font=4, srt=90, col=lincol[1:nbtips][id])
		text(tex[!id, 1], central_density_tips[1, !id] - offset, rownames(tex)[!id], 
			cex=cex, adj=c(1, 0.5), srt=90, col=lincol[1:nbtips][!id], 
			font=4)
		}

	if (tipmode == 3)
		{
		points(tex, col=lincol[1:nbtips], pch=18)
		text(tex[, 1], max(central_density_tips[2, ]), rownames(tex), cex=cex, 
			adj=c(0, 0.5), col=lincol[1:nbtips], font=4, 
			srt=45)
		}
	} # END plotAncClim_wUncertainty

