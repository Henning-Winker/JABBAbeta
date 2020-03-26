##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## Stock Assessment execution File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><


jabba_libs <- function(){
  # Install required packages if missing
  list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales",'snpar')
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  # Load Packages
  library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);library('snpar')
  library(mvtnorm);library(scales)}  

#---------------------------------------------------
# jabba2jags() to write jabba model JAGS file
#--------------------------------------------------


jabba2jags = function(assessment,scenario,model,add.catch.CV,psi.dist,sigma.proc,sigma.est,projection,output.dir){
  cat("\n","><> Write JABBA model for",assessment,scenario,"in JAGS to",output.dir," <><","\n","\n")
  # JAGS MODEL Standard
  sink(paste0(output.dir,"/JABBA.jags"))
  cat("
    
    model {
    
    # Prior specifications  
    eps <- 0.0000000000000000000000000000000001 # small constant    
    
    #Catchability coefficients
    for(i in 1:nq)
    {   
    q[i] ~ dunif(q_bounds[1],q_bounds[2])    
    }  
    
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])
    
    # Carrying Capacity SB0
    K ~ dlnorm(log(K.pr[1]),pow(K.pr[2], -2))
    
    # informative priors for Hmsy as a function of r
    r ~ dlnorm(log(r.pr[1]),pow(r.pr[2],-2))  
    
    ")
  
  if(model==4){
    cat("
      # Shape m prior
      m ~ dlnorm(log(mu.m),pow(m.CV,-2))
      ",append=TRUE)  
  }else{ cat(" 
      m <- mu.m     
    ",append=TRUE)}
  
  if(psi.dist =="beta"){
    cat("
      # Beta Prior for Biomass depletion at the start (deteministic)
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
  } else {
    cat("
      # Lognormal for Biomass depletion at the start (deteministic)
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)    
      ",append=TRUE)  
  }
  
  if(sigma.proc==TRUE){
    cat("
      # Process variance
      isigma2 <- isigma2.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
      ",append=TRUE)  
  }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)
           
           ",append=TRUE)}
  
  if(sigma.est==TRUE){
    cat("
      # Obsevation variance
      for(i in 1:nvar)
      {
      # Observation error
      itau2[i]~ dgamma(0.001,0.001)
      tau2[i] <- 1/itau2[i]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance
      
      }}
      ",append=TRUE)  
  }else{ cat(" 
      # Obsevation variance
           for(i in 1:nvar)
           {
           # Observation error
           itau2[i]~ dgamma(4,0.01)
           tau2[i] <- 1/itau2[i]
           }
           
           for(i in 1:nI)
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}
  
  # Run standard JABBA
  if(add.catch.CV==FALSE){
    cat("  
for(t in 1:N){    
estC[t] <- TC[t]
}
",append=TRUE)} else {
  cat("  
for(t in 1:N){    
      estC[t] ~ dlnorm(log(TC[t]),pow(CV.C[t],-2))
}
",append=TRUE)}
cat("  
  
    #Process equation
    Pmean[1] <- log(psi)
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    P[1] ~ dlnorm(Pmean[1],iPV[1]) # set to small noise instead of isigma2
    penB[1]  <- ifelse(P[1]<P_bound[1],log(K*P[1])-log(K*P_bound[1]),ifelse(P[1]>P_bound[2],log(K*P[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    penBK[1] <- 0
    
    # Process equation
    for (t in 2:N) 
    {
    Pmean[t] <- ifelse(P[t-1] > Plim,
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1)) - estC[t-1]/K,0.001)),
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1))*P[t-1]*slope.HS - estC[t-1]/K,0.001)))
    iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    P[t] ~ dlnorm(Pmean[t],iPV[t])
    penB[t]  <- ifelse(P[t]<(P_bound[1]),log(K*P[t])-log(K*(P_bound[1])),ifelse(P[t]>P_bound[2],log(K*P[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    # Depletion prior
    penBK[t] <- ifelse(b.yr[t] < 1,0,log(P[t])-log(b.pr[1]))
    }
    
    
    
    # Process error deviation 
    for(t in 1:N){
    Proc.Dev[t] <- log(P[t]*K)-log(exp(Pmean[t])*K)} 
    
    # Enforce soft penalties on bounds for P
    for(t in 1:N){
    pen.P[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    pen.bk[t] ~ dnorm(penBK[t],pow(b.pr[2],-2))
    }
    
    Hmsy <- r*pow(m-1,-1)*(1-1/m) 
    
    for (t in 1:N) 
    { 
    SB[t] <- K*P[t]    
    H[t] <- TC[t]/SB[t] 
    }
    
    # Observation equation in related to EB
    
    for(i in 1:nI)
    {
    for (t in 1:N) 
    { 
    Imean[t,i] <- log(q[sets.q[i]]*P[t]*K);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))#q[[i]]*P[t]*SB0*EBtoSB[t,i]
    Ihat[t,i]  <- exp(Imean[t,i])
  
}}
    
  
    #Management quantaties
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    
    MSY <- SBmsy*Hmsy
    for (t in 1:N)
    {
    # use x y to put them towards the end of the alphabetically sorted  mcmc object
    #SP[t] <- pow(r.pella,-(m-1))*SB[t]*(1-pow(P[t],m-1))
    BtoBmsy[t] <- SB[t]/SBmsy
    HtoHmsy[t] <- H[t]/(Hmsy) 
    }
    
    
    # Enforce soft penalty on K if < K_bounds >  
    K.pen ~ dnorm(penK,1000) # enforce penalty 
    penK  <- ifelse(K<(K_bounds[1]),log(K)-log(K_bounds[1]),ifelse(K>K_bounds[2],log(K)-log(K_bounds[2]),0)) # penalty if Pmean is outside viable biomass
    

    # Enforce soft penalty on process deviance if sigma.proc > 0.2 
    proc.pen ~ dnorm(penProc,1000) # enforce penalty 
    penProc  <- ifelse(sigma>sigmaproc_bound,log(sigma)-log(sigmaproc_bound),0) 
     
    
    # Enforce soft penalty on observation error if sigma.obs > sigma_bound 
    for(i in 1:nvar){
    obs.pen[i] ~ dnorm(penObs[i],1000) # enforce penalty 
    penObs[i]  <- ifelse(pow(tau2[i],0.5)>sigmaobs_bound,log(pow(tau2[i],0.5))-log(sigmaobs_bound),0) 
    }
    
    ", append=TRUE)

# PROJECTION
if(projection==TRUE){
  cat("
      for(i in 1:nTAC){
      # Project first year into the future
      prPmean[1,i] <- ifelse(P[N] > Plim,
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TACint/K,0.005)),
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))*4*P[N] - TACint/K,0.001)))
      prP[1,i] ~ dlnorm(prPmean[1,i],isigma2) 
      # Project all following years
      for(t in 2:pyrs){
      prPmean[t,i] <- ifelse(prP[t-1,i] > Plim,
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/K,0.001)),
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))*slope.HS*prP[t-1,i] - TAC[t-1,i]/K,0.001)))
      # process error (as monte-carlo simular)
      prP[t,i] ~ dlnorm(prPmean[t,i],isigma2)}
      for(t in 1:pyrs){
      prB[t,i] <- prP[t,i]*K
      prH[t,i] <- TAC[t,i]/prB[t,i]
      prHtoHmsy[t,i] <- prH[t,i]/Hmsy
      prBtoBmsy[t,i] <- prB[t,i]/SBmsy
      }}  
      ",append=TRUE)} else {
        cat("
            #Prevent error for unused input    
            fakeTAC <-  TAC
            fakepyrs <- pyrs 
            fakenTAC <- nTAC    
            fakeTACint <- TACint
            prHtoHmsy <- 1
            prP <- 1 
            prBtoBmsy <- 1    
            ", append=TRUE)}  

cat("
} # END OF MODEL
    ",append=TRUE,fill = TRUE)
sink()
} # #END jabba JAGS  


build_jabba <- function(
  catch = NULL,
  cpue = NULL,
  se =  NULL,
  assessment = "example",
  scenario = "s1",
  output.dir = NULL,
  model.type = c("Schaefer","Fox","Pella","Pella_m"),
  add.catch.CV = TRUE, # to match original assessment
  catch.cv = 0.15,
  Plim = 0, # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  r.dist = c("lnorm","range"), 
  r.prior = c(0.2,0.5), # mu, lod.sd
  K.dist = c("lnorm","range"), # K prior distribution
  K.prior = NULL,
  psi.dist = c("lnorm","beta"),
  psi.prior = c(0.9,0.25), 
  b.prior = c(FALSE,0.3,NA,c("bk","bbmy")[1]), # alternativ set as b.prior = c(mean,cv,yr,type=c("bk","bbmsy"))   
  BmsyK = 0.4, # Required specification for Pella-Tomlinson (model = 3 | model 4)
  shape.CV = 0.3, # Required specification for Pella-Tomlinson (Model 4)
  sets.q = 1:(ncol(cpue)-1), 
  sigma.est = FALSE, # Estimate additional observation variance
  sets.var = 1:(ncol(cpue)-1), # estimate individual additional variace
  fixed.obsE = c(0.01), # Minimum fixed observation erro 
  sigma.proc =  TRUE, # TRUE: Estimate observation error, else set to value
  proc.dev.all = TRUE, # TRUE: All year, year = starting year 
  igamma = c(0.3,0.01), # informative mean 0.07, CV 0.4
  projection = FALSE, # Switch on by Projection = TRUE 
  TACs = NULL,
  TACint =  NULL, # default avg last 3 years
  imp.yr = NULL, # default last year plus ONE
  pyrs = NULL, # Set number of projections years
  # Penalties
  P_bound = c(0.02,1.3),  # Soft penalty bounds for b/k
  sigmaobs_bound = 1, # Adds an upper bound to the observation variance  
  sigmaproc_bound = 0.2, # Adds an upper bound to the process variance  
  q_bounds= c(10^-30,1000), # Defines lower and upper bounds for q 
  K_bounds= c(0.01,10^10), # Defines lower and upper bounds for q 
  # Settings
  KOBE.plot = TRUE, # Produces JABBA Kobe plot 
  KOBE.type = c("ICCAT","IOTC")[2], # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
  Biplot= FALSE, # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
  save.trajectories =FALSE, # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
  harvest.label = c("Hmsy","Fmsy")[2], # choose label preference H/Hmsy versus Fmsy
  CPUE.plot= TRUE, # Runs state-tool to produce "alligned" multi-CPUE plot  
  meanCPUE = FALSE, # Uses averaged CPUE from state-space tool instead of individual indices  
  save.projections = FALSE, # saves projection posteriors as .RData object 
  catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
){
  
  
  # define model typue
  mod.names = model.type[1]
  model = which(c("Schaefer","Fox","Pella","Pella_m")%in%mod.names)
  # Set up output structure
  if(is.null(output.dir)) output.dir=paste0(getwd(),"/",assessment,"/",scenario)
  dir.create(output.dir,showWarnings = FALSE)
  
  #-------------------------
  # Prepare input data
  #-------------------------
  cat("\n","><> Prepare JABBA input data <><","\n","\n")
  
  
  if(is.null(catch)) stop("\n","\n","><> Catch Time series not provided <><","\n","\n")
  years = catch[,1]
  styr = min(years)
  endyr = max(years)
  n.years = length(years)
  if(length(which(styr:endyr%in%catch[,1]==FALSE))>0) stop("\n","\n","><> catch time series must have years in sequential order <><","\n","\n")
  
  if(is.null(se)) SE.I = FALSE
  if(is.null(cpue)){ 
    CatchOnly = TRUE
    cat(paste0("\n","><> Running Catch-Only mode: CatchOnly = TRUE <><","\n","\n"))
    cpue= catch[,1:2]
    colnames(cpue) = c("Year","Dummy Index")
    cpue[,2] = NA  
    cpue[1,2] = 1  
    SE.I = FALSE
    add.catch.CV =FALSE
    sigma.est = FALSE
    sets.q =1
    sets.var=1
    n.indices = 1
  } else {
    CatchOnly = FALSE
    
  }
  if(nrow(catch)!=nrow(cpue)) stop("\n","\n","><> cpue and catch differ in number of year <><","\n","\n") 
  catches = names(catch)[2:ncol(catch)]
  n.catches = length(catches)
  if(nrow(catch)!=nrow(cpue)) stop("\n","\n","><> cpue and catch differ in number of year <><","\n","\n") 
  styr.catch = min(catch[,1])
  styr.C = styr.catch-styr+1 
  conv.catch = as.numeric(rbind(matrix(rep(NA,(styr.C-1)*n.catches),styr.C-1,n.catches),as.matrix(catch[,-1])))
  
  if(length(which(conv.catch<0.0001)>0)) {
    cat(paste0("\n","><> Warnning: Replacing 0 Catch by small constant 0.0001 <><","\n","\n"))
    conv.catch[conv.catch< 0.001]=0
  }
  if(length(which(is.na(conv.catch)))>0) stop("\n","\n","><> Missing Catch values currently not permitted (NAs should be 0) <><","\n","\n")
  
  # Build Catch input
  Catch=matrix(conv.catch,nrow=n.years,ncol=n.catches)
  if(ncol(Catch)>1){
    cat(paste0("\n","><> Aggrigating multiple catch colums (assumed by fleet) to a single total catch column <><","\n","\n"))
    Catch = apply(Catch,1,sum)  
  }
  TC = as.numeric(Catch) # Total Catch
  
  # Catch CV option.
  if(add.catch.CV==TRUE){
    if(length(catch.cv)>1) {CV.C =catch.cv[,2]} else {CV.C = rep(catch.cv,length(TC))}   
    cat("\n","><> Assume Catch with error CV = ",mean(CV.C)," <><","\n","\n")
  } else {
    cat("\n","><> Assume Catch to be known without error <><","\n","\n")
  }
  
  
  
  # Build CPUE input
  indices = names(cpue)[2:ncol(cpue)]
  n.indices = max(length(indices),1)
  styr.cpue = min(cpue[,1])
  styr.I = styr.cpue-styr+1 
  # Convert input data to matrices for JAGS input
  conv.cpue = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
  CPUE=matrix(conv.cpue,nrow=n.years,ncol=n.indices)
  
  # Build Standard Error matrix
  if(is.null(se)){SE.I=FALSE} else {SE.I=TRUE} 
  if(SE.I==FALSE){
    cat("\n","><> SE.I=FALSE: Creatinng SE dummy matrix <><","\n","\n")
    se = cpue  
    conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
    se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
  } else{
    if(ncol(cpue)!=ncol(se) | nrow(cpue)!=nrow(se)) stop("\n","\n","><> SE and CPUE matrix must match <><","\n","\n")
    conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(se[,-1])))
    se2 = matrix(ifelse(is.na(conv.se),0.3^2,conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
  }
  
  
  # Add depletion prior
  if(b.prior[1]==FALSE){
    b.yr = rep(0,n.years) # activate by setting one year to one
    b.pr = c(0.5,0.1) # not active    
  } else {
    b.yr = rep(0,n.years)
    if(length(b.prior)>2 & b.prior[3] %in% years){
      b.pr = b.prior[1:2]     
      b.yr[which(years %in% b.prior[3])] =1
    } else {
      b.yr[n.years] = 1
      b.pr = b.prior[1:2]
      
    }}
  
  
  #---------------------
  # Index color palette
  #---------------------
  jabba.colors = as.character(c('#e6194b', "#3cb44b", "#ffe119",
                                "#0082c8","#f58231", "#911eb4",
                                "#46f0f0", "#f032e6", "#d2f53c",
                                "#fabebe", "#008080","#e6beff", "#aa6e28",rainbow(10)))
  
  
  #----------------------------------------------------------
  # Get JABBA parameterization and suplus production function
  #----------------------------------------------------------
  
  # For Pella-Tomlinson
  getshape = function(model,BmsyK){
    if(model==3 | model==4){ 
      #-----------------------------------------------
      # find inflection point
      ishape = NULL
      # Find shape for  SBmsytoK 
      ishape = seq(0.1,10,0.001)
      
      check.shape =((ishape)^(-1/(ishape-1))-BmsyK)^2
      #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
      shape =  ishape[check.shape==min(check.shape)] 
    } else {shape=FALSE}
    return(shape) 
  }
  shape = getshape(model=model,BmsyK=BmsyK)
  
  # Set shape m for Fox and Schaefer: Fox m ~1; Schaefer m =2
  if(shape==FALSE){
    if(model == 1){m=2} else {m = 1.001}}else{m=shape}
  
  #------------------------------------------------
  
  
  #----------------------------------------------------
  #----------------------------------------------------
  if(r.dist[1]=="range"){
    # initial range of r based on resilience (FishBase.org)
    if(length(r.prior)>1){ start.r = r.prior} else
      if(r.prior == "High") {
        start.r <- c(0.6,1.5)} else if(r.prior == "Medium") {
          start.r <- c(0.2,0.8)}    else if(r.prior == "Low") {
            start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
              start.r <- c(0.015,0.1)}  
    
    log.r = mean(log(start.r))
    sd.r = abs(log.r - log(start.r[1]))/2
    r.prior = c(exp(log.r),sd.r)  
    CV.r = sqrt(exp(sd.r^2)-1) 
  } else {
    log.r = log(r.prior[1])
    sd.r = r.prior[2]
    CV.r = sqrt(exp(sd.r^2)-1) 
  }
  
  #----------------------------------------------------
  # Prepare K prior
  #----------------------------------------------------
  if(is.null(K.prior)) K.prior = c(4*max(Catch),1) # mean and CV
  if(length(K.prior)<2) K.prior[2] = 1
  
  
  if(K.dist[1]=="range"){
    log.K =  mean(log(K.prior)) 
    sd.K= abs(log.K - log(K.prior[1]))/2
    log.K = mean(log(K.prior))+0.5*sd.K^2 # Will be back-bias corrected in lognorm function
    CV.K = sqrt(exp(sd.K^2)-1)
  } else {
    
    CV.K = K.prior[2]
    sd.K=sqrt(log(CV.K^2+1))
    log.K = log(K.prior[1])#-0.5*sd.K^2
  }
  
  
  
  
  # Get input priors
  K.pr = plot_lnorm(exp(log.K),CV.K,Prior="K")
  r.pr = plot_lnorm(mu=exp(log.r),CV=CV.r,Prior="r")
  
  psi.dist = psi.dist[1]
  if(psi.dist=="beta"){
    psi.pr = get_beta(mu=psi.prior[1],CV=psi.prior[2],Min=0,Prior=paste0("Prior B(",years[1],")/K"))} else {
      psi.pr = plot_lnorm(mu=psi.prior[1],CV=psi.prior[2],Prior=paste0("Prior B(",years[1],")/K"))  
    }
  if(model==4){
    dm = dlnorm((seq(0.001,5,0.1)),log(m),shape.CV)
    dm = dm/max(dm)
    bmsyk  = (seq(0.001,5,0.1))^(-1/(seq(0.001,5,0.1)-1))
  }
  
  
  # Note PRIORS and save input subfolder
  Priors =rbind(c(K.pr[1],CV.K),psi.prior,c(r.pr[1],CV.r))
  row.names(Priors) = c("K","Psi","r")
  colnames(Priors) = c("Mean","CV")                          
  Priors = data.frame(Priors)
  Priors$log.sd = sqrt(log(Priors[,2]^2+1))
  
  #----------------------------------------------------------
  # Set up JABBA 
  #----------------------------------------------------------
  cat("\n","><> Set up JAGS input data structure <><","\n")
  # Plot MSY
  # remove scientific numbers
  options(scipen=999)
  
  #----------------------------------------------------------
  # Setup TAC projection
  #---------------------------------------------------------
  if(is.null(TACs)) TACs = seq(0.5,1.3,0.1)*TC[n.years]
  if(is.null(TACint)) TACint =  mean(c(TC[length(TC)-1],TC[length(TC)])) # avg last 3 years
  # default avg last 3 years
  if(is.null(imp.yr))imp.yr = years[n.years]+1 # default last year plus ONE
  if(is.null(pyrs)) pyrs = 10 # Set number of projections years
  
  if(projection==TRUE){
    nTAC = length(TACs)
    TAC = mat.or.vec(pyrs,nTAC)
    yr.last = max(years) # assessment year  
    
    for(i in 1:nTAC){
      TAC[,i] = c(rep(TACint,imp.yr-yr.last-1),rep(TACs[i],pyrs-(imp.yr-yr.last-1)))  
    }
    
  }else{
    nTAC = 1  
    TAC = TC[n.years] #  
    pyrs = 1
  }
  
    
  
  
  #---------------------------------------------------------------
  # JABBA Schaefer/Fox Models 1-2, Pella 3, 4 estimate shape
  #---------------------------------------------------------------
  
  # Slope of hockey-stick
  slope.HS = ifelse(Plim==0,1/10^-10,1/Plim)
  nSel = 1 # setup for JABBA-SELECT version (in prep)
  nI = ncol(CPUE) # number of CPUE series
  stI = ifelse(proc.dev.all==TRUE,1, c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1]) #first year with CPUE
  
  # starting values
  nq = length(unique(sets.q))
  nvar = length(unique(sets.var))
  
  
  # JABBA input data 
  surplus.dat = list(N=n.years, TC = TC,I=CPUE,SE2=se2,r.pr=r.pr,psi.pr=psi.pr,K.pr = K.pr,
                     nq=nq,nI = nI,nvar=nvar,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),
                     sets.var=sets.var, sets.q=sets.q,Plim=Plim,slope.HS=slope.HS,
                     nTAC=nTAC,pyrs=pyrs,TAC=TAC,igamma = igamma,stI=stI,TACint =TACint,pen.P = rep(0,n.years) ,pen.bk = rep(0,n.years),proc.pen=0,K.pen = 0,
                     obs.pen = rep(0,nvar),P_bound=P_bound,q_bounds=q_bounds,sigmaobs_bound=sigmaobs_bound,sigmaproc_bound=sigmaproc_bound,K_bounds=K_bounds,mu.m=m,b.yr=b.yr, b.pr = b.pr)
  
  
  # PARAMETERS TO MONITOR
  params <- c("K","r", "q", "psi","sigma2", "tau2","m","Hmsy","SBmsy", "MSY", "BtoBmsy","HtoHmsy","CPUE","Ihat","Proc.Dev","P","SB","H","prP","prBtoBmsy","prHtoHmsy","TOE")
  
  
  #-----------------------------------------------
  # If shape parameter is estimated (Model =4)
  if(model==4){
    surplus.dat$m.CV = shape.CV }
  #-----------------------------------------------
  # If Catch Estimation with CV is used
  if(add.catch.CV==TRUE){
    surplus.dat$CV.C = CV.C  
    params = c(params,"estC")
    
  }
  
  
  # JAGS model file
  JABBA = "JABBA.jags"
  
  #--------------------------
  # Capture Settings
  #--------------------------
  jbinput = list()
  jbinput$data = list() 
  jbinput$jagsdata = list()   
  jbinput$settings = list()
  jbinput$data$yr = years
  jbinput$data$catch = catch
  jbinput$data$cpue = cpue 
  jbinput$data$se = se
  jbinput$jagsdata = surplus.dat
  jbinput$settings$params = params
  jbinput$settings$psi.dist = psi.dist
  jbinput$settings$psi.prior.raw = psi.prior
  jbinput$settings$SE.I = SE.I
  jbinput$settings$sigma.proc = sigma.proc
  jbinput$settings$model.id = model
  jbinput$settings$model.type = mod.names
  jbinput$settings$add.catch.CV = add.catch.CV
  jbinput$settings$catch.cv = catch.cv
  jbinput$settings$CatchOnly = CatchOnly
  jbinput$settings$proc.dev.all = proc.dev.all
  jbinput$settings$do.projections = projection
  jbinput$settings$TAC.implementation = imp.yr
  jbinput$settings$catch.metric = catch.metric
  jbinput$settings$harvest.label = harvest.label
  jbinput$settings$assessment = assessment
  jbinput$settings$scenario = scenario
  jbinput$settings$output.dir = output.dir
  jbinput$settings$cols = jabba.colors
  capture.output(jbinput, file=paste0(output.dir,"/Settings.txt"))
  #-------------------------------------------------------------------
  # write JAGS MODEL
  jabba2jags(assessment=assessment,scenario=scenario,model=model,add.catch.CV=add.catch.CV,psi.dist=psi.dist,sigma.proc=sigma.proc,sigma.est=sigma.est,projection=projection,output.dir=output.dir)
  
  return(jbinput)

  } # end of build_jabba()

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# fit_jabba()
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

fit_jabba = function(jbinput,
                     # MCMC settings
                     ni = 30000, # Number of iterations
                     nt = 5, # Steps saved
                     nb = 5000, # Burn-in
                     nc = 2, # number of chains
                     # init values
                     init.values = FALSE,
                     K.init = NULL,
                     r.init = NULL,
                     q.init = NULL,# vector
                     peels = NULL, # retro peel option 
                     save.all = FALSE,
                     save.jabba = TRUE,
                     save.csvs = FALSE,
                     save.projections = NULL
                     ){
  # mcmc saved
  nsaved = (ni-nb)/nt*nc
  # jabba model data
  jbd = jbinput$jagsdata
  # Initial starting values (new Eq)
  if(init.values==FALSE){
    inits = function(){list(K= rlnorm(1,log(jbd$K.pr[1])-0.5*0.3^2,0.3),r = rlnorm(1,log(jbd$r.pr[1]),jbd$r.pr[2]) ,q = runif(jbd$nq,min(jbd$I,na.rm=T)/max(jbd$TC,na.rm=T),min(jbd$I,2,na.rm=T)/max(jbd$TC,na.rm=T)), isigma2.est=runif(1,20,100), itau2=runif(jbd$nvar,80,200))}
  }else {
    if(is.null(init.K)) 
      stop("\n","\n","><> Provide init.K guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.r)) 
      stop("\n","\n","><> Provide init.r guess for option init.values=TRUE  <><","\n","\n")
    if(is.null(init.q)) 
      stop("\n","\n","><> Provide init.q vector guess for option init.values=TRUE  <><","\n","\n")
    if(length(init.q)!=nq)
      stop("\n","\n","><> init.q vector must match length of estimable q's, length(unique(sets.q))   <><","\n","\n")
    inits = function(){ list(K= init.K,r=init.r,q = init.q, isigma2.est=runif(1,20,100), itau2=runif(nvar,80,200))}
  }
  
  # retrospecitive
  years = jbinput$settings$yr
  if(is.null(peels)) peels = 0  
  if(peels > 0){
  jbinput$jagsdata$I[(length(years)-peels+1) : length(years),]  = NA
  }
    
  
  
  output.dir = jbinput$settings$output.dir
  params = jbinput$settings$params
  # jabba model building conditions
  
  
  ptm <- proc.time()
  
  mod <- jags(jbd, inits,params,paste0(output.dir,"/JABBA.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)  # adapt is burn-in
  
  proc.time() - ptm
  save.time = proc.time() - ptm
  
  # unpack
  settings= c(jbinput$data,jbinput$jagsdata,jbinput$settings)
  catch = settings$catch
  cpue= settings$cpue
  se = settings$se
  n.years = settings$N
  years = settings$yr
  assessment = settings$assessment
  scenario = settings$scenario
  CPUE = settings$I
  n.indices = settings$nI
  
  # if run with library(rjags)
  posteriors = mod$BUGSoutput$sims.list
  cat(paste0("\n","><> Produce results output of ",settings$model.type," model for ",settings$assessment," ",settings$scenario," <><","\n"))
  
  #-----------------------------------------------------------
  # <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
  #-----------------------------------------------------------
  
  # run some mcmc convergence tests
  par.dat= data.frame(posteriors[params[c(1:7)]])
  geweke = geweke.diag(data.frame(par.dat))
  pvalues <- 2*pnorm(-abs(geweke$z))
  pvalues
  
  heidle = heidel.diag(data.frame(par.dat))
  
  # postrior means + 95% BCIs
  #Model  parameter
  apply(par.dat,2,quantile,c(0.025,0.5,0.975))
  
  man.dat = data.frame(posteriors[params[8:10]])
  #Management quantaties
  apply(man.dat,2,quantile,c(0.025,0.5,0.975))
  
  # Depletion
  Depletion = posteriors$P[,c(1,n.years)]
  colnames(Depletion) = c(paste0("P",years[1]),paste0("P",years[n.years]))
  
  # Current stock status (Kobe posterior)
  H_Hmsy.cur = posteriors$HtoHmsy[,c(n.years)]
  B_Bmsy.cur = posteriors$BtoBmsy[,c(n.years)]
  
  
  # Prepare posterior quantaties
  man.dat = data.frame(man.dat,Depletion,B_Bmsy.cur,H_Hmsy.cur)
  
  results = round(t(cbind(apply(par.dat,2,quantile,c(0.025,0.5,0.975)))),3)
  
  results = data.frame(Median = results[,2],LCI=results[,1],UCI=results[,3],Geweke.p=round(pvalues,3),Heidel.p = round(heidle[,3],3))
  
  ref.points = round(t(cbind(apply(man.dat,2,quantile,c(0.025,0.5,0.975)))),3)
  
  ref.points = data.frame(Median = ref.points[,2],LCI=ref.points[,1],UCI=ref.points[,3])
  # get number of parameters
  npar = length(par.dat)
  # number of years
  N=n.years
  
  #-------------------------------------------------------------------------
  # Save parameters, results table and current status posterior in csv files
  #-------------------------------------------------------------------------
  
  # Safe posteriors (Produces large object!)
  if(save.all==TRUE) save(posteriors,file=paste0(output.dir,"/",assessment,"_",settings$model.type,"_",settings$scenario,"_posteriors.rdata"))
  
  # Make standard results table with parameter estimates and reference points
  Table = rbind(data.frame(results)[c("K","r","psi","sigma2","m"),1:3],data.frame(ref.points))  
  Table[4,] = round(sqrt((Table[4,])),3) 
  rownames(Table)[4] = "sigma.proc"
  colnames(Table) <- c("mu","lci","uci")
  
  #-----------------------------------------------
  # Stock trajectories
  #-----------------------------------------------
  #Bt = posteriors$SB
  #Ht = posteriors$H
  #Bt_Bmsy = posteriors$BtoBmsy
  #Ht_Hmsy = posteriors$HtoHmsy
  #Bt_K = posteriors$P
  
  Stock_trj = array(data=NA,dim=c(ncol(posteriors$SB),3,6),dimnames = list(years,c("mu","lci","uci"),c("B","F","BBmsy","FFmsy","BB0","procB")))
  for(i in 1:3){
    Stock_trj[,i,] =  cbind(t(apply(posteriors$SB,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$H,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$BtoBmsy,2,quantile,c(0.5,0.025,0.975)))[,i],
                            t(apply(posteriors$HtoHmsy,2,quantile,c(0.5,0.025,0.975)))[,i],t(apply(posteriors$P,2,quantile,c(0.5,0.025,0.975)))[,i],t(apply(posteriors$Proc.Dev,2,quantile,c(0.5,0.025,0.975)))[,i])
    
  }
  
  
  #------------------------
  # Production function
  #------------------------
  m.sp = median(posteriors$m)
  Bit = seq(1,median(posteriors$K),median(posteriors$K)/500)
  Cmsy = Bit*median(posteriors$Hmsy)
  B.sp = apply(posteriors$SB,2,mean)
  Hmsy.sp = median(posteriors$Hmsy) 
  SB0.sp = median(posteriors$K)
  SP = Hmsy.sp/(1-1/m.sp)*Bit*(1-(Bit/SB0.sp)^(m.sp-1))
  Bmsy.sp = median(posteriors$SBmsy)
  MSY.sp = quantile(posteriors$SBmsy*posteriors$Hmsy,c(0.025,0.5,0.975))
  
  #------------------
  # Goodness-of-Fit
  #------------------
  DIC =round(mod$BUGSoutput$DIC,1)
  
  if(settings$CatchOnly==FALSE){
    # get residuals
    Resids = NULL
    for(i in 1:n.indices){
      Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))   
    }
    
    # Standardized Residuals
    StResid = NULL
    for(i in 1:n.indices){
      StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                       apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))        
    }
    
    Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
    DF = Nobs-npar
    RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)
    SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
    Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
    # Produce statistice describing the Goodness of the Fit
  } else {
    
    Nobs =DF = RMSE = SDNR = Crit.value = NA
    
  }
  
  # Save Obs,Fit,Residuals
  jabba.res = NULL
  if(settings$CatchOnly==FALSE){
    for(i in 1:n.indices){
      
      Yr = years
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      exp.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.5))[is.na(cpue[,i+1])==F]
      expLCI.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.025))[is.na(cpue[,i+1])==F]
      expUCI.i = apply(posteriors$CPUE[,,i],2,quantile,c(0.975))[is.na(cpue[,i+1])==F]
      
      obs.i = cpue[is.na(cpue[,i+1])==F,i+1]
      sigma.obs.i = (apply(posteriors$TOE[,,i],2,quantile,c(0.5)))[is.na(cpue[,i+1])==F]
      
      yr.i = Yr[is.na(cpue[,i+1])==F]
      jabba.res = rbind(jabba.res,data.frame(scenario=settings$scenario,name=names(cpue)[i+1],year=yr.i,obs=obs.i,obs.err=sigma.obs.i,hat=exp.i,hat.lci=expLCI.i,hat.uci=expUCI.i,residual=log(obs.i)-log(exp.i),retro.peels=peels))  
    }
  }
  
  #----------------------------------
  # Predicted CPUE
  #----------------------------------
  cpue.hat = array(data=NA,dim=c(N,5,n.indices),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
  for(i in 1:n.indices){
    cpue.hat[,,i] = cbind(t(apply(posteriors$Ihat[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$Ihat[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
  }  
  #------------------------------------
  # Posterior Predictive Distribution
  #------------------------------------  
  
  cpue.ppd = array(data=NA,dim=c(N,5,n.indices),list(years,c("mu","lci","uci","se","obserror"),names(cpue[,-1])))
  for(i in 1:n.indices){
    cpue.ppd[,,i] = cbind(t(apply(posteriors$CPUE[,,i],2,quantile,c(0.5,0.025,0.975))),apply(log(posteriors$CPUE[,,i]),2,sd),(apply(posteriors$TOE[,,i],2,quantile,c(0.5))))
  }  
  
  
  
  #-----------------------------------
  # Note posteriors of key parameters
  #-----------------------------------
  sel.par = c(1,2,7,4,3,5)
  out=data.frame(posteriors[params[sel.par]])
  
  cat(paste0("\n","\n",paste0("><> Scenario ", jbinput$settings$scenario,"_",jbinput$settings$model.type," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))
  
  
  
  #-------------------------------
  # summarize results in jabba list
  #-------------------------------
  jabba = list()
  jabba$assessment  = assessment 
  jabba$scenario = scenario
  jabba$settings = c(jbinput$jagsdata,jbinput$settings)
  jabba$inputseries = list(cpue=cpue,se=se,catch=catch)
  jabba$pars=results
  jabba$estimates=Table
  jabba$yr = years 
  jabba$catch = settings$TC
  if(settings$add.catch.CV ==TRUE){jabba$est.catch = data.frame(yr=years,mu=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[1,],
                                                       lci=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[2,],
                                                       uci=apply(posteriors$estC,2,quantile,c(0.5,0.025,0.975))[3,])} else {
                                                         jabba$est.catch = "Require option: add.catch.cv = TRUE"
                                                       }
  jabba$cpue.hat=cpue.hat
  jabba$cpue.ppd=cpue.ppd
  jabba$timeseries=Stock_trj
  jabba$refpts = data.frame(factor=assessment,level=settings$scenario,quant = c("hat","logse"), k=c(median(posteriors$K),sd(log(posteriors$K))),bmsy=c(median(posteriors$SBmsy),sd(log(posteriors$SBmsy))),
                            fmsy=c(median(posteriors$Hmsy),sd(log(posteriors$Hmsy))),msy=c(median(posteriors$SBmsy*posteriors$Hmsy),sd(log(posteriors$SBmsy*posteriors$Hmsy)))) 
  jabba$pfunc = data.frame(factor=assessment,level=scenario,SB_i=round(Bit,3),SP=round(SP,3),Hmsy=round(Hmsy.sp,4),r=round(Hmsy.sp*(m.sp-1)/(1-1/m.sp),4),m=round(m.sp,3),MSY=round(as.numeric(MSY.sp[2]),3),SB0=round(SB0.sp,3),Cmsy=round(Cmsy,3))
  
  if(settings$CatchOnly==FALSE){ 
    jabba$diags = data.frame(factor=assessment,level=settings$scenario,name=jabba.res$name,year=jabba.res$year,season=1,obs=jabba.res$obs,hat=jabba.res$hat,hat.lci=jabba.res$hat.lci,hat.uci=jabba.res$hat.uci,residual=jabba.res$residual,retro.peels=jabba.res$retro.peels)
    jabba$residuals = array(Resids,dim=c(nrow(Resids),ncol(Resids)),dimnames = list(names(cpue)[-1],years))
    jabba$std.residuals = array(StResid,dim=c(nrow(StResid),ncol(StResid)),dimnames = list(names(cpue)[-1],years))
    
    
  } else {
    jabba$diags = "Not Available for Catch-Only option"
    jabba$residuals = "Not Available for Catch-Only option"
  }
  jabba$stats = data.frame(Stastistic = c("N","p","DF","SDNR","RMSE","DIC"),Value = c(Nobs,npar,DF,SDNR,RMSE,DIC))
  jabba$pars_posterior = out
  jabba$kobe = data.frame(factor=assessment,level=scenario,yr=years[N],stock=posteriors$BtoBmsy[,N],harvest=posteriors$HtoHmsy[,N])
  
  save(jabba,file=paste0(output.dir,"/jabba.rdata"))
  
  
  return(jabba)
  
  } # end of fit_jabba()
  
  





#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Generic jabba functions
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#--------------------------------------------------
# Function to get beta prior parameters
#--------------------------------------------------
get_beta <- function(mu,CV,Min=0,Prior="x",Plot=FALSE){
  a = seq(0.0001,1000,0.001)
  b= (a-mu*a)/mu
  s2 = a*b/((a+b)^2*(a+b+1))
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = (a-mu*a)/mu
  x = seq(Min,1,0.001)  
  pdf = dbeta(x,a,b)  
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}

#--------------------------------------------------
# Function to get gamma prior parameters
#--------------------------------------------------

get_gamma <- function(mu,CV,Prior="x", Plot=FALSE){
  a = seq(0.00001,10000,0.0001)
  b = a/mu
  s2 = (a/b^2)
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = a/mu
  x = sort(rgamma(1000,a,b))  
  pdf = dgamma(x,a,b)  
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}


#--------------------------------------------------
# Function to get lognormal prior parameters
#--------------------------------------------------
plot_lnorm <- function(mu,CV,Prior="x",Plot=FALSE){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))  
  pdf = dlnorm(x,log(mu),sdev)  
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(mu,sdev))
}


#------------------------------------
# Function kobeJabba for FLR
#------------------------------------
kobeJabba<-function(x,minyear=1){
  
  out=cbind(melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year+minyear-1
  out}

#-------------------------------------------------
# Function kobeJabbaProj for projections with FLR
#-------------------------------------------------
kobeJabbaProj<-function(x,minyear=1,tac=NULL){
  
  out=cbind(melt(x[,,,2]),c(x[,,,3]))
  names(out)=c("iter","year","tac","stock","harvest")
  out$year=out$year+minyear-1
  
  out}



#---------------------------------------------------------------------------------
# JABBA plotting functions
#---------------------------------------------------------------------------------- 

#----------------
# Total Landings
#----------------
jbplot_catch <- function(jabba,output.dir=NULL){
  
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  
  cat(paste0("\n","><> jbplot_catch()  <><","\n"))
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Landings_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  cord.x <- c(jabba$yr,rev(jabba$yr))
  y<-rep(0,length(jabba$yr))
  plot(jabba$yr,(jabba$catch),type="l",ylim=c(0,max(jabba$catch,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
  polygon(cord.x,c(jabba$catch,rev(y)),col="gray",border=1,lty=1)
  dev.off()
}

#------------------------------
# Catch estimated with CV
#------------------------------
jbplot_catcherror <- function(jabba,output.dir=NULL){  
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  if(jabba$settings$add.catch.CV==TRUE){
    cat(paste0("\n","><> jbplot_catcherror()  <><","\n"))  
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.7,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/Catch.fit_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
        res = 200, units = "in")
    par(Par)
    # estimated Catch
    predC = jabba$est.catch 
    years = jabba$yr
    cord.x <- c(jabba$yr,rev(jabba$yr))
    cord.y<-c(predC[,3],rev(predC[,4]))
    plot(years,(jabba$catch),type="n",ylim=c(0,max(predC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
    polygon(cord.x,cord.y,col="gray",border=0,lty=1)
    lines(years,predC[,2],lwd=2,col=4)
    points(years,(jabba$catch),pch=21,bg=0,cex=1.5)
    legend("topright",c("Observed","Predicted"),pch=c(21,-1),bg=0,lwd=c(-1,2),col=c(1,4),bty="n")
    dev.off() 
  } else {
    cat(paste0("\n","><> jbplot_catcherror() only available if add.catch.CV=TRUE <><","\n"))
  }
}

#------------------------------
# Plot Posteriors
#------------------------------

jbplot_ppdist <- function(jabba, output.dir=NULL){  
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_ppist() - prior and posterior distributions  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)
  #informative priors
  Prs = as.matrix(cbind(jabba$settings$K.pr,jabba$settings$r.pr,c(0,0),jabba$settings$psi.pr))
  
  #Posteriors
  Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Posteriors_",jabba$assessment,"_",jabba$scenario,".png"),width  = 8, height = 2.5*round(length(node_id)/3,0), 
      res = 200, units = "in")
  par(Par)
  
  for(i in 1:length(node_id))
  {
    
    post.par = as.numeric(unlist(out[paste(node_id[i])]))
    
    if(i==1){
      
      rpr =  rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
      plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab="K",ylab="",xaxs="i",yaxs="i",main="")
      
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
      
    }  
    if(i==2){
      rpr = rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
      pdf = stats::density(post.par,adjust=2) 
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
      plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
      
      
    }
    
    if(i==3){
      if(jabba$settings$model.id<4){
        plot(1,1,type="n",xlim=range(0.5,2.5),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
        abline(v=jabba$pars["m",1],lwd=2)}
      if(jabba$settings$model.id==4){
        mpr = rlnorm(10000,log(jabba$jagsdata$mu.m),jabba$jagsdata$m.CV) 
        pdf = stats::density(post.par,adjust=2) 
        prior = dlnorm(sort(mpr),log(m),shape.CV)   
        plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
        
        polygon(c(sort(mpr),rev(sort(mpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
        polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
        PPVM = round(mean(post.par)/mean(rpr),3)
        legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
        
      }
    }
    
    
    if(i==4){
      if(jabba$settings$psi.dist=="beta"){
        parm = fitdist(post.par[post.par<1 & post.par>0.01], "beta")$estimate
        rpr = rbeta(10000,(Prs[1,4]),Prs[2,4]) 
        pdf = stats::density(post.par,adjust=2)  
        prior = dbeta(sort(rpr),psi.pr[1],psi.pr[2])   
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
        PPVM = round(mean(post.par)/mean(rpr),3)
        legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
        
      } else {
        rpr = rlnorm(10000,log(Prs[1,4]),Prs[2,4]) 
        pdf = stats::density(post.par,adjust=2)  
        prior = dlnorm(sort(rpr),log(Prs[1,4]),Prs[2,4])}
      plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,quantile(rpr,c(0.001,0.999)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
      
      #legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
    }        
    
    if(i>4){
      if(jabba$settings$sigma.proc!=TRUE & i==length(node_id)) {
        plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
        abline(v=sigma.proc^2,lwd=2)} else {
          
          pdf = stats::density(post.par,adjust=2)  
          plot(pdf,type="l",xlim=range(0,post.par),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
          if(i==length(node_id)& jabba$settings$igamma[1]>0.9){
            rpr = 1/rgamma(10000,jabba$settings$igamma[1],jabba$settings$igamma[2])
            prior = stats::density(rpr,adjust=2)
            polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
            PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
            PPVM = round(mean(post.par)/mean(rpr),3)
            legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
            
          }
          
          polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
          #legend('topright',c("Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.8,0.6)),bty="n")
        } }         
    
  }
  mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  dev.off()   
} # End of ppdist plot


#-----------------------------
# MCMC chains of posteriors
#-----------------------------

jbplot_mcmc <- function(jabba, output.dir=NULL){
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_mcmc() - mcmc chains  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)
  
  Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/MCMC_",jabba$assessment,"_",jabba$scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0), 
      res = 200, units = "in")
  par(Par)
  for(i in 1:length(node_id)){
    
    post.par = as.numeric(unlist(out[paste(node_id[i])]))
    plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=4)
    lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)   
  }
  dev.off()
}

#---------------------------------------------------------------------------
# Plot CPUE with expectected mean CIs and Posterior Predictive Distributions
#----------------------------------------------------------------------------

jbplot_cpue <- function(jabba, output.dir=NULL){  
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_cpue() - fits to CPUE <><","\n"))
    
    N = jabba$settings$N 
    years = jabba$yr
    n.indices = jabba$settings$nI
    series = 1:jabba$settings$nI
    CPUE = jabba$settings$I 
    indices = unique(jabba$diags$name)
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    #CPUE FITS
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/Fits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
        res = 200, units = "in")
    par(Par)
    
    
    for(i in 1:n.indices){
      
      # set observed vs predicted CPUE
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      fit =  t(jabba$cpue.ppd[,c(2,1,3),i])   
      fit.hat = t(jabba$cpue.hat[,c(2,1,3),i])
      mufit = mean(fit[2,])
      fit = fit/mufit
      fit.hat = fit.hat/mufit
      
      cpue.i = CPUE[is.na(CPUE[,i])==F,i]
      yr.i = Yr[is.na(CPUE[,i])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,i])==F,(i)])
      
      ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
      
      cord.x <- c(Yr,rev(Yr))
      cord.y <- c(fit[1,yr],rev(fit[3,yr]))
      cord.yhat <- c(fit.hat[1,yr],rev(fit.hat[3,yr]))
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(jabba$yr),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      polygon(cord.x,cord.yhat,col=grey(0.3,0.5),border=grey(0.3,0.5),lty=2)
      
      lines(Yr,fit[2,yr],lwd=2,col=1)
      if(jabba$settings$SE.I  ==TRUE | max(jabba$settings$SE2)>0.01){ plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
        points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")}
      
      legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
    }
    
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  } else {
    cat(paste0("\n","><> jbplot_cpue() not available CatchOnly=TRUE <><","\n"))
  }
} # End of CPUE plotting function

#-------------------------
# Plot logfits CPUE
#-------------------------

jbplot_logfits <- function(jabba, output.dir= NULL){  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_logfits()  <><","\n"))
    
    N = jabba$settings$N 
    years= jabba$yr
    n.indices = jabba$settings$nI
    series = 1:jabba$settings$nI
    CPUE = jabba$settings$I 
    indices = unique(jabba$diags$name)
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    
    #log CPUE FITS
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/logFits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
        res = 200, units = "in")
    par(Par)
    
    
    for(i in 1:n.indices){
      
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      fit = t(jabba$cpue.hat[,c(2,1,3),i])
      mufit = mean(fit[2,])
      fit = fit/mufit
      cpue.i = CPUE[is.na(CPUE[,i])==F,i]
      yr.i = Yr[is.na(CPUE[,i])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,i])==F,(i)])
      
      ylim = log(c(min(fit[,yr[is.na(CPUE[,i])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,i])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
      
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      
      
      lines(Yr,log(fit[2,yr]),lwd=2,col=4)
      if(jabba$settings$SE.I ==TRUE | max(jabba$settings$SE2)>0.01){ plotCI(yr.i,log(cpue.i/mufit),ui=log(exp(log(cpue.i)+1.96*se.i)/mufit),li=log(exp(log(cpue.i)-1.96*se.i)/mufit),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
        points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")}
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
    }
    
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_logfit() not available CatchOnly=TRUE <><","\n"))
  }
} # End of logfit


#-------------------------------------------------
# JABBA Residual Plot
#-------------------------------------------------
jbplot_residuals <- function(jabba, output.dir=NULL){
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_residuals() - JABBA residual plot  <><","\n"))
    years = jabba$yr  
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jabba$residuals
    Yr = jabba$yr
    n.years = length(Yr) 
    n.indices = jabba$settings$nI
    indices = unique(jabba$diags$name)
    series = 1:jabba$settings$nI
    
    # JABBA-residual plot
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/Residuals_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
        res = 200, units = "in")
    par(Par)
    
    
    plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
    boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
    abline(h=0,lty=2)
    positions=runif(n.indices,-0.2,0.2)
    
    for(i in 1:n.indices){
      for(t in 1:n.years){
        lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=jabba$settings$cols[i])}
      points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=jabba$settings$cols[i])}
    mean.res = apply(Resids[,as.numeric(colnames(Resids))%in%cpue.yrs],2,mean,na.rm =TRUE)
    smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
    lines(cpue.yrs,smooth.res,lwd=2)
    # get degree of freedom
    Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
    
    RMSE = round(jabba$stats[5,2],1)
    
    legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
    legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$col[series],1),lwd=c(rep(-1,n.indices),2))
    dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_residuals() not available CatchOnly=TRUE <><","\n"))
  }
} # End of functions

#---------------------------------------
# Plot Stadardized Residuals
#--------------------------------------

jbplot_stdresiduals <- function(jabba, output.dir=NULL){
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_staresiduals() - standardized residuals  <><","\n"))
    years = jabba$yr  
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jabba$residuals
    Yr = jabba$yr
    n.years = length(Yr) 
    n.indices = jabba$settings$nI
    indices = unique(jabba$diags$name)
    series = 1:jabba$settings$nI
    StResid = jabba$std.residuals
    
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/StandardizedResids_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
        res = 200, units = "in")
    par(Par)
    # Standardized Residuals
    plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
    boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
    abline(h=0,lty=2)
    positions=runif(n.indices,-0.2,0.2)
    
    for(i in 1:n.indices){
      for(t in 1:n.years){
        lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=jabba$settings$cols[i])}
      points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=jabba$settings$cols[i])}
    mean.res = apply(StResid[,as.numeric(colnames(StResid))%in%cpue.yrs],2,mean,na.rm =TRUE)
    smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
    lines(cpue.yrs,smooth.res,lwd=2)
    SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(jabba$stats[1,2]-1)),2)
    Crit.value = (qchisq(.95, df=(jabba$stats[1,2]-1))/(jabba$stats[1,2]-1))^0.5
    legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
    legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$cols[series],1),lwd=c(rep(-1,n.indices),2))
    dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_stdresiduals() not available CatchOnly=TRUE <><","\n"))
  }
  
} # end of function



#-------------------------------------------------
# Function to do runs.test and 3 x sigma limits  
#------------------------------------------------
runs.sig3 <- function(x,type=NULL) {
  if(is.null(type)) type="resid"
  if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)} 
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){ 
    runstest = snpar::runs.test(x) 
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001  
    }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}

#----------------------------------------------------
# Runs test plots
#---------------------------------------------------


jbplot_runstest <- function(jabba, output.dir=NULL){
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_runstest()   <><","\n"))
    
    
    years = jabba$yr  
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    Resids = jabba$residuals
    n.years = length(years) 
    n.indices = jabba$settings$nI
    indices = unique(jabba$diags$name)
    series = 1:jabba$settings$nI
    
    
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/Residual_RunsTests_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
        res = 200, units = "in")
    par(Par)
    
    
    for(i in 1:n.indices){
      
      resid = (Resids[i,is.na(Resids[i,])==F])  
      res.yr = years[is.na(Resids[i,])==F]
      runstest = runs.sig3(x=as.numeric(resid),type="resid")
      # CPUE Residuals with runs test
      plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab=expression(log(cpue[obs])-log(cpue[pred])))
      abline(h=0,lty=2)
      lims = runstest$sig3lim
      cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
      rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
      for(j in 1:length(resid)){
        lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))  
      }
      points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
    }  
    
    
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(expression(log(cpue[obs])-log(cpue[pred])), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_stdresiduals() not available CatchOnly=TRUE <><","\n"))
  }
  
} # end of function




#------------------------------
# Plot process error deviation
#------------------------------

jbplot_procdev <- function(jabba, output.dir=NULL){ 
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_procdev() - Process error diviations on log(biomass)  <><","\n"))
  
  years=jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/ProcDev_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  ylim = c(min(-0.22,jabba$timeseries[,,"procB"]),max(0.22,jabba$timeseries[,,"procB"]))#range(proc.dev)*1.1
  cord.x <- c(years,rev(years))
  cord.y <- c(jabba$timeseries[,2,"procB"],rev(jabba$timeseries[,3,"procB"]))
  # Process Error
  plot(years,jabba$timeseries[,1,"procB"],ylab="Process Error Deviates",xlab="Year",ylim=ylim,type="n")
  polygon(cord.x,cord.y,col='grey',border=0,lty=2)
  lines(years,jabba$timeseries[,1,"procB"],lwd=2)
  lines(years,rep(0,length(years)),lty=5)
  dev.off()
} # end of plot function


#-----------------------------------------------------------
# <><<><<><<><<>< JABBA Management Plots ><>><>><>><>><>><>
#-----------------------------------------------------------

#-----------------------
# Plot trajectores: B,F,Bmsy,FFmsy,BB0 
#-----------------------

jbplot_trj <-  function(jabba, type = c("B","F","BBmsy","FFmsy","BB0") ,output.dir=NULL){ 
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  for(i in 1:length(type)){  
    
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/",type[i],"_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
        res = 200, units = "in")
    par(Par)
    cat(paste0("\n","><> jbplot_trj() - ", type[i]," trajectory  <><","\n"))
    
    ylabs = c(paste("Biomass",jabba$settings$catch.metric),ifelse(jabba$settings$harvest.label=="Fmsy","Fishing mortality F","Harvest rate H"),expression(B/B[MSY]),ifelse(jabba$settings$harvest.label=="Fmsy",expression(F/F[MSY]),expression(H/h[MSY])),expression(B/B[0])) 
    trj = jabba$timeseries[,,paste(type[i])] 
    years = jabba$yr
    ylim = c(0, max(trj[,3]))
    cord.x <- c(years,rev(years))
    cord.y <- c(trj[,2],rev(trj[,3]))
    plot(years,trj[,1],ylab=ylabs[i],xlab="Year",ylim=ylim,type="n")
    polygon(cord.x,cord.y,col='grey',border=0,lty=2)
    lines(years,trj[,1],lwd=2,col=1)
    if(i==1) lines(years,rep(jabba$refpts$bmsy[1],length(years)),lty=5)
    if(i==2) lines(years,rep(jabba$refpts$fmsy[1],length(years)),lty=5)
    if(i>2 & i <5) lines(years,rep(1,length(years)),lty=5)
    if(i==1) text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]*1.11,expression(paste(B[MSY])))
    if(i==2) text((max(years)-min(years))/30+years[1],jabba$refpts$fmsy[1]*1.11,ifelse(jabba$settings$harvest.label=="Fmsy",expression(F[MSY]),expression(H[MSY])))
    if(i==5){
      lines(years,rep(jabba$refpts$bmsy[1]/jabba$refpts$k[1],length(years)),lty=5)
      text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]/jabba$refpts$k[1]*1.11,expression(paste(B[MSY])))
    }
    dev.off()
  }} # end of plot function



#-----------------------------------------
# Produce JABBA SP-phase plot
#-----------------------------------------
jbplot_spphase <-  function(jabba ,output.dir=NULL){ 
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_spphase() - JABBA Surplus Production Phase Plot  <><","\n"))
  
  # extract pars  
  m = jabba$pfun$m[1]
  Bit = jabba$pfun$SB_i
  Cmsy = jabba$pfunc$Cmsy
  B = jabba$timeseries[,"mu","B"]
  Hmsy.sp = jabba$pfunc$Hmsy[1] 
  SB0.sp =jabba$pfunc$SB0[1]
  SP = jabba$pfunc$SP
  Bmsy.sp = jabba$estimates["SBmsy",1]
  MSY.sp = jabba$estimates["MSY",1]
  N = jabba$settings$N
  years = jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SPphase_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(SB0.sp,0,max(SP),SB0.sp,SB0.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(jabba$catch,na.rm=T)*1.05,max(MSY.sp*1.1)))),xlim=c(0,max(Bit,B)),ylab="Surplus Production",xlab="Spawning Biomass",xaxs="i",yaxs="i")
  rect(0,0,SB0.sp*1.1,SB0.sp*1.1,col="green",border=0)
  rect(0,0,SB0.sp,SB0.sp,col="yellow",border=0)
  rect(0,max(SP),SB0.sp,SB0.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)  
    #i = i+1
  }
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY.sp[1],2),rep(MSY.sp[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,jabba$catch,lty=1,lwd=1)
  points(B,jabba$catch,cex=0.8,pch=16)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],jabba$catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(Bmsy.sp,2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright', 
         c(expression(B[MSY]),"MSY","SP","Catch",paste(sel.years)), 
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)), 
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
} #end of plotting function


#------------------------------------------------------
# status plot (kobe, biplot)
#------------------------------------------------------

jbplot_kobe <-  function(jabba ,output.dir=NULL){ 
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_kobe() - Stock Status Plot  <><","\n"))
  
  mu.f = jabba$timeseries[,,"FFmsy"]   
  mu.b = jabba$timeseries[,,"BBmsy"]
  f = jabba$kobe$harvest
  b = jabba$kobe$stock
  years=jabba$yr
  N = length(years)
  # fit kernel function
  kernelF <- gplots::ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Kobe_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), ylim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,ylab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  zb2 = c(0,1)
  zf2  = c(1,100)
  zb1 = c(1,100)
  zf1  = c(0,1)
  polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
  polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
  polygon(c(1,100,100,1),c(1,1,100,100),col="orange",border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  
  points(mu.b[,1],mu.f[,1],pch=16,cex=1)
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[,1],mu.f[,1], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[sel.yr,1],mu.f[sel.yr,1],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
  Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
  
  
  
  sel.years = c(years[sel.yr])
  ## Add legend
  
  legend('topright', 
         c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
         lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")  
  
  dev.off()
} # End of Kobe plot 


#---------------------------------------------------------
# Produce 'post-modern' biplot (see Quinn and Collie 2005)
#---------------------------------------------------------

jbplot_biplot <-  function(jabba ,output.dir=NULL){ 
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  cat(paste0("\n","><> jbplot_biplot() - Stock Status Plot  <><","\n"))
  mu.f = jabba$timeseries[,,"FFmsy"]   
  mu.b = jabba$timeseries[,,"BBmsy"]
  f = jabba$kobe$harvest
  b = jabba$kobe$stock
  years=jabba$yr
  N = length(years)
  
  # fit kernel function
  kernelF <- gplots::ci2d(f,b,nbins=201,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Biplot_",jabba$assessment,"_",jabba$cenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", ylim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), xlim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,xlab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")
  
  # and fill areas using the polygon function
  fint = seq(0.001,100,0.01)
  # read ftarget,bthreshold
  ftarget<-0.8
  bthreshold<-0.2
  
  #Zone X
  xb=bthreshold+(1.0-bthreshold)/ftarget*fint
  xf =  ifelse(xb>1,0.8,fint)
  polygon(c(0,0,xf),c(max(xb),bthreshold,xb),col="green")
  zb = bthreshold+(1.0-bthreshold)*fint
  zf  = ifelse(zb>1,1,fint) 
  polygon(c(zf,rep(max(fint),2),rep(0,2)),c(zb,max(zb),0,0,bthreshold),col="red")
  
  polygon(c(xf,rev(zf)),c(xb,rev(zb)),col="yellow")
  
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.f[2,],mu.b[2,],pch=16,cex=1)
  
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.f[,1],mu.b[,1], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.f[sel.yr,1],mu.b[sel.yr,1],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  
  sel.years = years[sel.yr]
  ## Add legend
  legend('topright', 
         c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I."), 
         lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"), 
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,4),1.7,1.7,1.7),bty="n")
  
  Zone  = NULL
  Status = NULL
  X  = 0.15
  Y = 0
  Z = -0.15
  
  for(i  in 1:length(f))
  {
    if(b[i]>1.0){
      if(f[i]<ftarget){
        Zone[i]<-X
      } else if (f[i]>1.0){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }
    } else {
      if(b[i]>bthreshold+(1.0-bthreshold)/ftarget*f[i]){
        Zone[i]<-X
      } else if(b[i]<bthreshold+(1.0-bthreshold)*f[i]){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }
      
      
    }}
  
  perGreen = round(length(Zone[Zone==0.15])/length(Zone)*100,1) 
  perYellow = round(length(Zone[Zone==0])/length(Zone)*100,1) 
  perRed = round(length(Zone[Zone==-0.15])/length(Zone)*100,1)
  
  mtext(expression(paste(B/B[MSY])), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  mtext(ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)
  
  text(0.65,2.4,paste0(perGreen,"%"))
  text(0.9,2.4,paste0(perYellow,"%"))
  text(1.2,2.4,paste0(perRed,"%"))
  
  dev.off()
  
} # End of biplot function

#-------------------------
# wrapper plot function
#------------------------
jbplots = function(jabba,output.dir = NULL,statusplot ="kobe"){
  if(is.null(output.dir)){
    if(file.exists(jabba$settings$output.dir)){output.dir = jabba$settings$output.dir} else {getwd()}}  
  
  jbplot_catch(jabba,output.dir=output.dir) # catch.metric
  jbplot_catcherror(jabba,output.dir=output.dir) # posteriors
  jbplot_cpue(jabba,output.dir=output.dir) # check years
  jbplot_logfits(jabba,output.dir=output.dir) # check n.indices
  jbplot_mcmc(jabba,output.dir=output.dir)
  jbplot_ppdist(jabba,output.dir=output.dir) # check m
  jbplot_procdev(jabba,output.dir=output.dir)
  jbplot_trj(jabba,output.dir=output.dir) 
  jbplot_spphase(jabba,output.dir=output.dir) # check TC
  jbplot_residuals(jabba,output.dir=output.dir) # check years
  jbplot_stdresiduals(jabba,output.dir=output.dir)
  jbplot_runstest(jabba,output.dir=output.dir)
  if(statusplot =="kobe"){
    jbplot_kobe(jabba,output.dir=output.dir)} else {
      jbplot_biplot(jabba,output.dir=output.dir)}
}


