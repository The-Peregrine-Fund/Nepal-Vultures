model{ 
  # JAGS code
  ####################################
  # Indices and parameter descriptions
  ####################################
  # lmu.N log-transformed abundance index after effort (effort.effect), circle, and distance corrections (dist offset)
  # r= population growth rate on log scale, 
  # lambda= backtransformed population growth rate
  # Effort correction: parameters B and p
  # sig.obs= observation error of counts, sig.proc= process error of abundance
  # indices
  # k= count in ncounts, 
  # y= year (also yr), 
  # s= bird conservation region (also str) in nstrata ,
  # c= route
  
  # Priors and constants
  p <- 0.001
  B ~ dnorm(0, 0.001)
  mean.r ~ dnorm(0, 0.001)
  mean.lambda<- exp(mean.r)
  tau.proc <-  1/(sig.proc*sig.proc)
  sig.proc ~ dnorm(0,1)T(0,)
  tau.obs <-  1/(sig.obs*sig.obs)
  sig.obs ~ dnorm(0,1)T(0,)
  
  #### Data model ###### 
  for(k in 1 : ncounts) { 	
    lN.est[k] <- lmu.N[yr[k],str[k]] + effort.effect[k]  
    logtheta[k] ~ dnorm (lN.est[k], tau.obs)
    log(theta[k]) <- logtheta[k]
    count[k] ~ dpois(theta[k]) 
    # effort effects
    ScaledEffort[k] <- effort[k]/xbar.effort
    trans.effort[k] <- (pow(ScaledEffort[k],p)-1)/p   # Common p # 
    effort.effect[k] <-  B*trans.effort[k] 
  }  # k
  
  for(s in 1 : nstratas) { 
    # Baseline year and strata expectations ######
    lmn1st[s] <- log(mn1st[s]+0.01) # add a small constant to prevent log(zero) which will result in failure to run
    lmu.N[1,s] ~ dnorm(lmn1st[s], 0.1) # weak priors for 1st yr abundance are observed values
    # Dynamics
    for( y in startyear : (endyear-1)) { 
      lmu.N[y+1,s] <- lmu.N[y,s] + r[y,s] 
      r[y,s] ~ dnorm(mean.r, tau.proc)
      lambda[y,s] <- exp(r[y,s])
    } # y
} # s

  # Derived parameters 
  for( y in startyear : (endyear-1 )) {
    mn.r[y] <- mean ( r[y,] )
    mn.lam[y] <- exp( mn.r[y] )
  } # y

  for(s in 1 : nstratas) { 
    mn.N.site[s] <- exp( mean(lmu.N[,s]) )
  } # s
  mn.N <-  exp( mean(log(mn.N.site[]) ) )

  } #model end

