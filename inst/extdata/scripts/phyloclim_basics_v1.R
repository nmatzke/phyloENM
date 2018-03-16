
#######################################################
# Use PNOs (predicted niche occupancies) and
# phyloclim::anc.clim to estimate 
# ancestral niche profiles, using the 
# Smith et al. (2009) method
#######################################################

library(phyloclim)
source("/GitHub/phyloENM/R/")


# load phylogeny and PNOs of Oxalis sect. Palmatifoliae
data(tree)
data(PNO)

# choose summer precipitation for analysis
clim <- PNO$PrecipitationWarmestQuarter

plot_PNOs(clim, varname="PrecipitationWarmestQuarter")

# estimate ancestral tolerances
target=tree; posterior=NULL; pno=clim; n=100; method="GLS"; returnwhat="original"; CIs=c(0.025,0.975)
ac = anc_clim_wUncertainty(target=tree, posterior=NULL, pno=clim, n=100, method="GLS", returnwhat="samples", CIs=c(0.025,0.975))


# visualize results
plotAncClim_wUncertainty(x=ac, clades=NULL, density=TRUE, tipmode=1, nchar=3, ylab="Precipitation of warmest quarter (mm)", plot_CIs=FALSE)


#######################################################
# What is the statistical uncertainty of these 
# ancestral niche estimates?
#######################################################
x=ac; clades=NULL; density=TRUE; tipmode=1; nchar=3; ylab="Precipitation of warmest quarter (mm)"
cex=NULL; tipspace=NULL; cladespace=1; lwd=NULL; plot_CIs=TRUE

plotAncClim_wUncertainty(x=ac, clades=NULL, density=TRUE, tipmode=1, nchar=3, ylab="Precipitation of warmest quarter (mm)", plot_CIs=TRUE)

# Not bad actually...



#######################################################
# NEXT:

# ? PNOs from Dan's fit models?
# 1. Model features from Dan's fit models
#######################################################
