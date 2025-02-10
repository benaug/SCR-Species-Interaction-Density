# SCR-Species-Interaction-Density

This is a 2 species SCR interaction model where species 2 either avoids or is attracted to species 1. Species 2 activity centers are a function of the local density* of species 1.

By "density" I mean the expected density of the realized individuals at a snapshot in time. This needs to be unpacked! There is a within home range resource selection process that combines an individual availability distribution with an RSF function to produce individual use distributions, or "utilization distributions".
Assuming individual site use is proportional to the utilization distribution, the expected density of realized individuals in each cell at any snapshot in time
is just the sum of the individual utilization distributions in each cell.

You many want to convert this density to different units, e.g. if we interpret density in units of km^2, converting to units of 100km^2 leads to plausible beta.interact values that are closer to 0. Generally, I've found plausible values for repulsion between -5 and 0 and for attraction, between 0 and 1 or so.
Attraction is much more likely to produce unidentifiable process model parameters. You need species 2 individuals living far away enough from species 1 individuals such that they are not attracted to them to be able to estimate the species 2 density intercept well.

The species interaction can be added to the species 2 linear predictor like this (I have the density intercept separated out so the activity center likelihoods don't need to be computed when updating the intercept):

lambda2.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta12*D.cov2[1:n.cells] + beta.interact*(100*use.dist1.sum[1:n.cells]/cellArea))

Note, that we can also write it in terms of an "interaction function" for easier comparison to the soft-core point process approach where the interaction function is a function of distance instead of density.

h[1:n.cells] <- exp(beta.interact*(100*use.dist1.sum[1:n.cells]/cellArea))
 
lambda2.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta12*D.cov2[1:n.cells])*h[1:n.cells]

As with the soft-core point process approach, this model seems to require heroic SCR data sets. The testscripts are set up with 144 traps and detection turned up higher than is realistic in most cases.
Identifiability is better if no RSF covariate is used. There is a test script that uses the RSF covariate and another that does not.