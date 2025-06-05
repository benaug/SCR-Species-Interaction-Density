NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Density covariates
  D01 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D02 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D.beta11 ~ dnorm(0,sd=10) #sps 1 response to D.cov1
  D.beta12 ~ dnorm(0,sd=10) #sps 2 response to D.cov2
  beta.interact ~ dunif(-10,10) #sps1 effect on sps2 expected D
  #RSF coefficients
  rsf.beta1 ~ dnorm(0,sd=10) 
  rsf.beta2 ~ dnorm(0,sd=10)
  #availability distribution spatial scale
  sigma1 ~ dunif(0,20)
  sigma2 ~ dunif(0,20)
  #detection intensity
  lambda.detect1 ~ dunif(0,10)
  lambda.detect2 ~ dunif(0,10)
  #--------------------------------------------------------------
  #Density model
  D1.intercept <- D01*cellArea
  D2.intercept <- D02*cellArea
  lambda1.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta11*D.cov1[1:n.cells])
  use.dist1.sum[1:n.cells] <- sumUse(use.dist=use.dist1[1:M1,1:n.cells],z=z1[1:M1])
  # 100*use.dist1.sum[1:n.cells]/cellArea this is expected density of realized individuals in N/100km^2 units
  #given that realized individuals are moving around their home range, what is expected number per unit area in each cell
  #at any snapshot in time? Assuming individual site use choices are categorical(use.dist[i,1:n.cells]) and independent.
  h[1:n.cells] <- exp(beta.interact*(100*use.dist1.sum[1:n.cells]/cellArea))
  lambda2.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta12*D.cov2[1:n.cells])*h[1:n.cells]
  pi1.cell[1:n.cells] <- lambda1.cell[1:n.cells]/pi1.denom #expected proportion of total N1 in cell c
  pi2.cell[1:n.cells] <- lambda2.cell[1:n.cells]/pi2.denom #expected proportion of total N2 in cell c
  pi1.denom <- sum(lambda1.cell[1:n.cells])
  pi2.denom <- sum(lambda2.cell[1:n.cells])
  lambda1 <- D1.intercept*pi1.denom #Expected N1
  lambda2 <- D2.intercept*pi2.denom #Expected N2
  N1 ~ dpois(lambda1) #Realized N1
  N2 ~ dpois(lambda2) #Realized N2
  #Resource selection function evaluated across all cells
  rsf1[1:n.cells] <- exp(rsf.beta1*rsf.cov1[1:n.cells])
  rsf2[1:n.cells] <- exp(rsf.beta2*rsf.cov2[1:n.cells])
  for(i in 1:M1){#sps 1
    s1.cell[i] ~ dcat(pi1.cell[1:n.cells])
    s1[i,1:2] <- dSS[s1.cell[i],1:2] #pull out X,Y for this cell
    #Individual availability distributions conditioned on cells, bivariate Normal centered on activity center
    avail.dist1[i,1:n.cells] <- getAvail(s=s1[i,1:2],sigma=sigma1,res=res,x.vals=x.vals[1:n.cells.x],
                                        y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    use.dist1[i,1:n.cells] <- getUse(rsf=rsf1[1:n.cells],avail=avail.dist1[i,1:n.cells])
    #extract ind by trap use probability from cells with traps, multiply by detection intensity
    for(j in 1:J){
      # expected detections at trap is proportional to use of the cell containing trap
      #convert counts to detections
      p1[i,j] <- 1 - exp(-lambda.detect1*use.dist1[i,trap.to.cell[j]]) 
    }
    y1[i,1:J] ~ dBinomialVector(size=K1D[1:J],prob=p1[i,1:J],z=z1[i]) #vectorized obsmod, skips z[i]=0 calculations
  }
  for(i in 1:M2){#sps 2
    s2.cell[i] ~ dcat(pi2.cell[1:n.cells])
    s2[i,1:2] <- dSS[s2.cell[i],1:2] #pull out X,Y for this cell
    #Individual availability distributions conditioned on cells, bivariate Normal centered on activity center
    avail.dist2[i,1:n.cells] <- getAvail(s=s2[i,1:2],sigma=sigma2,res=res,x.vals=x.vals[1:n.cells.x],
                                         y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    use.dist2[i,1:n.cells] <- getUse(rsf=rsf2[1:n.cells],avail=avail.dist2[i,1:n.cells])
    #extract ind by trap use probability from cells with traps, multiply by detection intensity
    for(j in 1:J){
      #expected detections at trap is proportional to use of the cell containing trap
      #convert counts to detections
      p2[i,j] <- 1 - exp(-lambda.detect2*use.dist2[i,trap.to.cell[j]])
    }
    y2[i,1:J] ~ dBinomialVector(size=K1D[1:J],prob=p2[i,1:J],z=z2[i]) #vectorized obsmod, skips z[i]=0 calculations
  }
})# end model