##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## JABBA: Just Another Bayesian Biomass Assessment
## Input File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

rm(list=ls())
# Install missing packages
list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape)

#----------------------------------------------------------------
# Setup working directories and output folder labels 
#-----------------------------------------------------------------
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/GitHub/JABBA_testruns" 
setwd(File) # Writes JABBA model in this file
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/GitHub/JABBAbeta"
# JABBA version
version = "v1.5beta"
# Set Assessment file: assement folder within File that includes .csv input files
assessment = "SWO_SA" 
# add specifier for assessment (File names of outputs)


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= FALSE # Runs state-tool to produce "alligned" multi-CPUE plot  
meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
Projection = TRUE # Use Projections: requires to define TACs vectors 
save.projections = TRUE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
Reproduce.seed = FALSE # If FALSE a random seed assigned to each run, if TRUE set.seed(123)
jabba2FRL = TRUE
P_bound = c(0.03,1) # Depletion penalty as used in original SWO assessment
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Optional: Note Scenarios
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# S1: Model including Brazil1 
# S2: Model excluding Brazil1
# S3: Base-case Model with time blocks on ESP and JPN 
# S4: Add scenario as example for using average CPUE
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Specify Scenario name for output file names
Scenarios = c(paste0("SWO",c("1.1","1.5beta"))) 

#Note for plotting 
Run = "Versions"
Labels =  Scenarios
B.rs = SP.rs = SP.max = Bx = MSY.rs = Bmsy.rs = Fmsy.rs = BtoBmsy.rs = FtoFmsy.rs = H.rs = P.rs =DICs = RMSEs = NULL


for(s in 1:2){
  Scenario = Scenarios[s] 
  if(s==2) version = "v1.5beta"
  if(s==1) version = "v1.1"
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Suplus Production model specifications
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # Choose model type: 
  # 1: Schaefer
  # 2: Fox  
  # 3: Pella-Tomlinsson  
  # 4: Pella-Tomlinsson with m prior
  
  Model = 3
  
  Mod.names = c("Schaefer","Fox","Pella","Pella_m")[Model]
  
  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
  Plim = 0
  
  # Required specification for Pella-Tomlinson (Model = 3)
  BmsyK = 0.4 # Set Surplus Production curve inflection point
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
  
  # Reverse bias correction to match v1.5
  if(s==2){
    sdev= sqrt(log(K.prior[2]^2+1))
    K.prior[1] = exp(log(K.prior[1])+0.5*sdev^2)  
  } 
  K.prior
  
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  
  psi.dist= c("lnorm","beta")[1]
  # specify as mean and CV 
  psi.prior = c(1,0.25) 
  if(s==2){
    sdev= sqrt(log(psi.prior[2]^2+1))
    psi.prior[1] = exp(log(psi.prior[1])+0.5*sdev^2)  
  } 
  psi.prior
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
  
  # Reverse bias correction to match v1.1
  if(s==2){
  sdev= sqrt(log(r.prior[2]^2+1))
  r.prior[1] = exp(log(r.prior[1])+0.5*sdev^2)  
  } 
  r.prior
  
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
  proc.dev.all = FALSE
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
  Projection = FALSE # Switch on by Projection = TRUE 
  
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

  
  # Note results
  Bx = cbind(Bx,Bit)
  SP.rs = cbind(SP.rs,SP)
  B.rs = cbind(B.rs,apply(posteriors$SB,2,median))
  MSY.rs = c(MSY.rs,median(posteriors$MSY))
  Bmsy.rs = c(Bmsy.rs,median(posteriors$SBmsy))
  Fmsy.rs = c(Fmsy.rs,median(posteriors$Hmsy))
  BtoBmsy.rs = cbind(BtoBmsy.rs,apply(posteriors$BtoBmsy,2,median))
  FtoFmsy.rs = cbind(FtoFmsy.rs,apply(posteriors$HtoHmsy,2,median))
  
}# THE END



folder = "Compare"
dir.create(paste0(File,"/",assessment,"/",folder),showWarnings = F)

cols = c(2,4) # Here S1-S3  
# 2x2 
Par = list(mfrow=c(2,2),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/",folder,"/All_",assessment,Run,".png"), width = 6, height = 6, 
    res = 200, units = "in")
par(Par)
#><> plot B
y = B.rs
plot(years,years,type="n",ylim=c(0,max(y)),ylab="Biomass (t)",xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=2,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}

legend("bottomleft",paste(Labels),col=c(cols),bty="n",cex=0.8,pt.cex=0.7,pch=c(seq(14,14+length(Labels),1)),lwd=c(2,rep(1,length(Labels))))
#abline(h=1,lty=2)


# Plot SP
plot(years,years,type="n",ylim=c(0,max(SP.rs)),xlim=c(0,max(Bx)),ylab="Surplus Production (t)",xlab="Biomass (t)")
for(i in 1:length(Scenarios)){
  lines(Bx[,i],SP.rs[,i],col=cols[i],lwd=2,lty=5)
  points(Bx[SP.rs[,i]==max(SP.rs[,i]),i],max(SP.rs[,i]),col=cols[i],pch=14+i,cex=1.2)
}

# Plot BtoBmsy
y = BtoBmsy.rs
plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(B/B[MSY])),xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=2,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}
abline(h=1,lty=2)

#><> Plot FtoFmsy
y = FtoFmsy.rs
plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(F/F[MSY])),xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=cols[i],lwd=2,lty=1)
  points(years,y[,i],col=cols[i],pch=14+i,cex=0.6)
}
abline(h=1,lty=2)
dev.off()


