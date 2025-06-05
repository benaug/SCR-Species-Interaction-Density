sumUse <- nimbleFunction(
  run = function(use.dist = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(use.dist)[1]
    n.cells <- nimDim(use.dist)[2]
    use.dist.sum <- rep(0,n.cells)
    for(i in 1:M){
      if(z[i]==1){
        use.dist.sum <- use.dist.sum + use.dist[i,]
      }
    }
    return(use.dist.sum)
  }
)

dBinomialVector <- nimbleFunction(
  run = function(x = double(1), size = double(1), prob = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual (never occurs with all known IDs)
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size=size, prob=prob, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),size = double(1), prob = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(prob)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)


getCell <- nimbleFunction(#cell 0 not allowed in this model, but leaving in as an error check
  run = function(u = double(1),res=double(0),cells=integer(2),xlim=double(1),ylim=double(1)) {
    returnType(double(0))
    inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
    if(inout==1){
      u.cell <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
    }else{
      u.cell <- 0
    }
    return(u.cell)
  }
)

getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),n.cells.x=integer(0),n.cells.y=integer(0)) {
    returnType(double(1))
    avail.dist.x <- rep(0,n.cells.x)
    avail.dist.y <- rep(0,n.cells.y)
    delta <- 1e-8 #this sets the degree of trimming used to get individual availability distributions
    x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
    y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
    #convert to grid edges instead of centroids
    x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x]+0.5*res)
    y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y]+0.5*res)
    #trim in x and y direction
    if(x.vals.edges[1]<x.limits[1]){
      x.start <- sum(x.vals.edges<x.limits[1])
    }else{
      x.start <- 1
    }
    if(x.vals.edges[n.cells.x]>x.limits[2]){
      x.stop <- which(x.vals.edges>x.limits[2])[1]
    }else{
      x.stop <- n.cells.x
    }
    if(y.vals.edges[1]<y.limits[1]){
      y.start <- sum(y.vals.edges<y.limits[1])
    }else{
      y.start <- 1
    }
    if(y.vals.edges[n.cells.y]>y.limits[2]){
      y.stop <- which(y.vals.edges>y.limits[2])[1]
    }else{
      y.stop <- n.cells.y
    }
    #get pnorms
    pnorm.x <- rep(0,n.cells.x+1)
    pnorm.y <- rep(0,n.cells.y+1)
    pnorm.x[x.start:(x.stop+1)] <- pnorm(x.vals.edges[x.start:(x.stop+1)], mean=s[1], sd=sigma)
    pnorm.y[y.start:(y.stop+1)] <- pnorm(y.vals.edges[y.start:(y.stop+1)], mean=s[2], sd=sigma)
    # Compute availability distributions
    avail.dist.x[x.start:x.stop] <- pnorm.x[(x.start+1):(x.stop+1)] - pnorm.x[x.start:x.stop]
    avail.dist.y[y.start:y.stop] <- pnorm.y[(y.start+1):(y.stop+1)] - pnorm.y[y.start:y.stop]
    avail.dist.tmp <- matrix(0,n.cells.x,n.cells.y)
    sum.dist <- 0
    for(i in x.start:x.stop){
      for(j in y.start:y.stop){
        avail.dist.tmp[i,j] <- avail.dist.x[i]*avail.dist.y[j]
        sum.dist <- sum.dist + avail.dist.tmp[i,j]
      }
    }
    avail.dist <- c(avail.dist.tmp)
    #if any probability mass is outside state space, normalize
    if(sum.dist<1){
      avail.dist <- avail.dist/sum.dist
    }
    return(avail.dist)
  }
)

getUse <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1)) {
    returnType(double(1))
    use.dist <- rsf*avail.dist
    use.dist <- use.dist/sum(use.dist)
    return(use.dist)
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#Required custom update for N/z
zSampler1 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    n.cells <- control$n.cells
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    sps2.nodes <- control$sps2.nodes
    sps2.lik.nodes <- control$sps2.lik.nodes
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    #track these "manually" so computations faster than nimble will do them
    use.dist1.sum.initial <- model$use.dist1.sum
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z1==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.sps2 <- model$getLogProb(sps2.lik.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N1[1] <<-  model$N1[1] - 1
          model$z1[pick] <<- 0

          #sps2 stuff
          use.dist1.sum.proposed <- use.dist1.sum.initial - model$use.dist1[pick,] #subtract this ind out
          model$use.dist1.sum <<- use.dist1.sum.proposed #put into model to skip resumming over all individuals
          model$calculate(sps2.nodes) #update sps2 nodes

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.sps2 <- model$calculate(sps2.lik.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y + lp.proposed.sps2) -
                                    (lp.initial.N + lp.initial.y + lp.initial.sps2)
          accept <- decide(log_MH_ratio)

          if(accept) {
            mvSaved["N1",1][1] <<- model[["N1"]]
            mvSaved["use.dist1.sum",1][1:n.cells] <<- model[["use.dist1.sum"]][1:n.cells]
            mvSaved["lambda2",1][1] <<- model[["lambda2"]]
            mvSaved["lambda2.cell",1][1:n.cells] <<- model[["lambda2.cell"]][1:n.cells]
            mvSaved["pi2.cell",1][1:n.cells] <<- model[["pi2.cell"]][1:n.cells]
            mvSaved["pi2.denom",1][1] <<- model[["pi2.denom"]]
            mvSaved["z1",1][pick] <<- model[["z1"]][pick]
            use.dist1.sum.initial <- use.dist1.sum.proposed
          }else{
            model[["N1"]] <<- mvSaved["N1",1][1]
            model[["use.dist1.sum"]][1:n.cells] <<- mvSaved["use.dist1.sum",1][1:n.cells]
            model[["lambda2"]] <<- mvSaved["lambda2",1][1]
            model[["lambda2.cell"]][1:n.cells] <<- mvSaved["lambda2.cell",1][1:n.cells]
            model[["pi2.cell"]][1:n.cells] <<- mvSaved["pi2.cell",1][1:n.cells]
            model[["pi2.denom"]] <<- mvSaved["pi2.denom",1][1]
            model[["z1"]][pick] <<- mvSaved["z1",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
            model$calculate(sps2.lik.nodes)
          }
        }
      }else{#add
        if(model$N1[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z1==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.sps2 <- model$getLogProb(sps2.lik.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0

          #propose new N/z
          model$N1[1] <<-  model$N1[1] + 1
          model$z1[pick] <<- 1

          #sps2 stuff
          use.dist1.sum.proposed <- use.dist1.sum.initial + model$use.dist1[pick,] #add this ind in
          model$use.dist1.sum <<- use.dist1.sum.proposed #put into model to skip resumming over all individuals
          model$calculate(sps2.nodes) #update sps2 nodes

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.sps2 <- model$calculate(sps2.lik.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick])

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y + lp.proposed.sps2) -
                                      (lp.initial.N + lp.initial.y + lp.initial.sps2)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N1",1][1] <<- model[["N1"]]
            mvSaved["use.dist1.sum",1][1:n.cells] <<- model[["use.dist1.sum"]][1:n.cells]
            mvSaved["lambda2",1][1] <<- model[["lambda2"]]
            mvSaved["lambda2.cell",1][1:n.cells] <<- model[["lambda2.cell"]][1:n.cells]
            mvSaved["pi2.cell",1][1:n.cells] <<- model[["pi2.cell"]][1:n.cells]
            mvSaved["pi2.denom",1][1] <<- model[["pi2.denom"]]
            mvSaved["z1",1][pick] <<- model[["z1"]][pick]
            use.dist1.sum.initial <- use.dist1.sum.proposed
          }else{
            model[["N1"]] <<- mvSaved["N1",1][1]
            model[["use.dist1.sum"]][1:n.cells] <<- mvSaved["use.dist1.sum",1][1:n.cells]
            model[["lambda2"]] <<- mvSaved["lambda2",1][1]
            model[["lambda2.cell"]][1:n.cells] <<- mvSaved["lambda2.cell",1][1:n.cells]
            model[["pi2.cell"]][1:n.cells] <<- mvSaved["pi2.cell",1][1:n.cells]
            model[["pi2.denom"]] <<- mvSaved["pi2.denom",1][1]
            model[["z1"]][pick] <<- mvSaved["z1",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
            model$calculate(sps2.lik.nodes)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z
zSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z2==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N2[1] <<-  model$N2[1] - 1
          model$z2[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N2",1][1] <<- model[["N2"]]
            mvSaved["z2",1][pick] <<- model[["z2"]][pick]
          }else{
            model[["N2"]] <<- mvSaved["N2",1][1]
            model[["z2"]][pick] <<- mvSaved["z2",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N2[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z2==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N2[1] <<-  model$N2[1] + 1
          model$z2[pick] <<- 1
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N2",1][1] <<- model[["N2"]]
            mvSaved["z2",1][pick] <<- model[["z2"]][pick]
          }else{
            model[["N2"]] <<- mvSaved["N2",1][1]
            model[["z2"]][pick] <<- mvSaved["z2",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)