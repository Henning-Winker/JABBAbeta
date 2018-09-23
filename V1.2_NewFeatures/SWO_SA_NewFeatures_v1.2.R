##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## JABBA: Just Another Bayesian Biomass Assessment
## Input File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

rm(list=ls())
# required packages
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library("fitdistrplus")
library(reshape)


#----------------------------------------------------------------
# Setup working directories and output folder labels 
#-----------------------------------------------------------------
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/GitHub/JABBA_testruns" 
setwd(File) # Writes JABBA model in this file
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/GitHub/JABBAbeta"
# JABBA version
version = "v1.2beta"
# Set Assessment file: assement folder within File that includes .csv input files
assessment = "SWO_SA" 
# add specifier for assessment (File names of outputs)


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= FALSE # Runs state-tool to produce "alligned" multi-CPUE plot  
meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
Projection = TRUE # Use Projections: requires to define TACs vectors 
save.projections = TRUE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
Reproduce.seed = FALSE # If FALSE a random seed assigned to each run, if TRUE set.seed(123)
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Optional: Note Scenarios
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Ref: Base-case Model with P_dev for all years 
# Est.m: Estimate m with a prior
# Catch.CV: Add Catch.CV = 0.2, i.e allow Catches to be
# CatchOnly: Run model only on prior and no CPUE 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Specify Scenario name for output file names
Scenarios = c("Ref","Est.m","Catch.CV","CatchOnly") 


#Note for a comparative plots 
Run = "_NewFeatures"
Labels =  Scenarios
B.rs = SP.rs = SP.max = Bx = MSY.rs = Bmsy.rs = Fmsy.rs = BtoBmsy.rs = FtoFmsy.rs = H.rs = P.rs =DICs = RMSEs = NULL

# Execute multiple JABBA runs in loop 
for(s in 1:length(Scenarios)){
  Scenario = Scenarios[s] 
  
  #><>><>><>><>><>><>><>><>><>
  # Add new Feature conditions
  #><>><>><>><>><>><>><>><>><>
  add.catch.CV = FALSE
  CatchOnly = FALSE
  Model = 3
  if(s==2) Model = 4 # shape.CV define in L105
  if(s==3){ 
    add.catch.CV = TRUE
    catch.cv = 0.2 
  }
  if(s==4) CatchOnly = TRUE 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Suplus Production model specifications
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # Choose model type: 
  # 1: Schaefer
  # 2: Fox  
  # 3: Pella-Tomlinsson  
  # 4: Pella-Tomlinsson with m prior
  
  Mod.names = c("Schaefer","Fox","Pella","Pella_m")[Model]
  
  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
  Plim = 0
  
  # Required specification for Pella-Tomlinson (Model = 3)
  BmsyK = 0.4 # Set Surplus Production curve inflection point
  shape.CV = 0.35 # Must be defined if Model = 4!
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>  
  
  #--------------------------------------------------
  # Read csv files
  #--------------------------------------------------
  
  # Use SEs from csv file for abudance indices (TRUE/FALSE)
  SE.I = TRUE
  
  # Load assessment data
  catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
  cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#
  
  if(SE.I ==TRUE){
    se =  read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
  }
  
  names(cpue)
  names(catch)
  
  #--------------------------------------------------
  # option to exclude CPUE time series or catch year
  #--------------------------------------------------
  
  # Set up Base-Case for SWO

  # Remove BrazilI
  cpue = cpue[,-c(2)]
  se = se[,-c(2)]
  
  
  
  names(cpue)
  ncol(catch)
  ncol(cpue)
  
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  meanCPUE = FALSE
  
  #------------------------------------------------
  # Prior for unfished biomass K
  #------------------------------------------------
  # The option are: 
  # a) Specify as a lognormal prior with mean and CV 
  # b) Specify as range to be converted into lognormal prior
  
  K.dist = c("lnorm","range")[1]
  
  # if lnorm use mean and CV; if range use lower,upper bound
  K.prior = c(200000,1) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  
  psi.dist= c("lnorm","beta")[1]
  # specify as mean and CV 
  psi.prior = c(1,0.25) 
  
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  
  
  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are: 
  # a) Specifying a lognormal prior 
  # b) Specifying a resiliance category after Froese et al. (2017; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)
  
  # use [1] lognormal(mean,stdev) or [2] range (min,max) or
  r.dist = c("lnorm","range")[1] 
  
  r.prior = c(0.42,0.37) 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.add = TRUE
  sigma.est = TRUE
  
  # Series
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = c(0.2) # Important if SE.I is not availble
  
  # Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  # Determines if process error deviation are estimated for all years (TRUE)  
  # or only from the point the first abundance index becomes available (FALSE)
  proc.dev.all = TRUE 
  #------------------------------------------
  if(sigma.proc == TRUE){
    igamma = c(4,0.01) #specify inv-gamma parameters
    
    # Process error check
    gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
  }else{
    sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  }
  #--------------------------------------------
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = TRUE # Switch on by Projection = TRUE 
  
  # Check final year catch 
  catch[nrow(catch),]
  
  # Set range for alternative TAC projections
  TACs = seq(10000,18000,1000) #example
  
  # Intermitted TAC to get to TAC implementation year
  #TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # avg last 3 years
  TACint = 10058 # Catch for 2016
  # Set year of first TAC implementation
  imp.yr = 2018
  # Set number of projections years
  pyrs = 10
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # MCMC settings
  ni <- 30000 # Number of iterations
  nt <- 5 # Steps saved
  nb <- 5000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc
  
  # Run model (BSPSPexe file must be in the same working directory)
  source(paste0(JABBA.file,"/JABBA",version,".R")) 
  # Run plot function
  source(paste0(JABBA.file,"/JABBA_plots_",version,".R")) 
  
  
  # Note results
  Bx = cbind(Bx,Bit)
  SP.rs = cbind(SP.rs,SP)
  B.rs = cbind(B.rs,apply(posteriors$SB,2,median))
  H.rs = cbind(H.rs,apply(posteriors$H,2,median))
  P.rs = cbind(P.rs,apply(posteriors$P,2,median))
  MSY.rs = c(MSY.rs,median(posteriors$MSY))
  Bmsy.rs = c(Bmsy.rs,median(posteriors$SBmsy))
  Fmsy.rs = c(Fmsy.rs,median(posteriors$Hmsy))
  BtoBmsy.rs = cbind(BtoBmsy.rs,apply(posteriors$BtoBmsy,2,median))
  FtoFmsy.rs = cbind(FtoFmsy.rs,apply(posteriors$HtoHmsy,2,median))
  DICs = c(DICs,DIC)
  RMSEs = c(RMSEs,RMSE)
}# THE END

folder = "V1.2_NewFeatures"
dir.create(paste0(File,"/",assessment,"/",folder),showWarnings = F)

cols = jabba.colors[c(4:6,11,3)] # Here S1-S3  
# 3x2 
Par = list(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/",folder,"/All_",assessment,Run,".png"), width = 6, height = 8, 
    res = 200, units = "in")
par(Par)
#><> plot B
y = B.rs

plot(years,years,type="n",ylim=c(0,max(y)),ylab="Biomass (t)",xlab="Year")

for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=1,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}

legend("bottomleft",paste(Labels),col=c(cols),bty="n",cex=0.8,pt.cex=0.7,pch=c(seq(14,14+length(Labels),1)),lwd=c(2,rep(1,length(Labels))))

y = H.rs

plot(years,years,type="n",ylim=c(0,max(y)),ylab="Fishing Mortality F",xlab="Year")

for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=1,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}

y = P.rs

plot(years,years,type="n",ylim=c(0,1),ylab="B/K",xlab="Year")

for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=1,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}
abline(h=0.4,lty=2)
#abline(h=1,lty=2)




# Plot SP
plot(years,years,type="n",ylim=c(0,max(SP.rs)),xlim=c(0,max(Bx)),ylab="Surplus Production (t)",xlab="Biomass (t)")


for(i in 1:length(Scenarios)){
  lines(Bx[,i],SP.rs[,i],col=cols[i],lwd=1,lty=5)
  points(Bx[SP.rs[,i]==max(SP.rs[,i]),i],max(SP.rs[,i]),col=cols[i],pch=14+i,cex=1.2)
  
}
#lines(rep(Bx[SP.rs[,1]==max(SP.rs[,1]),1],2),c(0,max(SP.rs[,1])),lty=2)
#abline(h = max(SP.rs[,1]),lty=2)


# Plot BtoBmsy
y = BtoBmsy.rs
plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(B/B[MSY])),xlab="Year")

for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=1,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}

abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)


#><> Plot FtoFmsy
y = FtoFmsy.rs

plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(F/F[MSY])),xlab="Year")

for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=1,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}

#legend("topright",c("Reference",rev(paste(years[(n.years-nback):(n.years-1)]))),col=c(4,grey(cols,1)),bty="n",cex=0.8,pt.cex=0.7,pch=c(-1,seq(15,14+nback,1)),lwd=c(2,rep(1,nback)))
abline(h=1,lty=2)
#mtext(paste("Year"), side=1, outer=T, at=0.5,line=1,cex=0.9)
dev.off()


GOFdat = data.frame(Scenario=Labels, DIC=DICs,RMSE=RMSEs)
write.csv(GOFdat,paste0(File,"/",assessment,"/",folder,"/GOF",assessment,Run,".csv"))

