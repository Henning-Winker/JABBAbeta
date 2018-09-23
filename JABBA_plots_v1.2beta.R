#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# JABBA Plot Function for version 1.2
# Completed Sep 2018
# Developed by Henning Winker, Felipe Carvalho and Maia Kapur
# NOTE: Best run from "Prime file" as
# source(paste0(JABBA.file,"/JABBA_plots_",version,".R")) 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#----------------
# Total Landings
#----------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/Landings_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

cord.x <- c(years,rev(years))
y<-rep(0,length(years))
plot(years,(TC),type="l",ylim=c(0,max(TC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",catch.metric),main="")
polygon(cord.x,c(TC,rev(y)),col="gray",border=1,lty=1)
dev.off()

# if Catch estimated with CV
if(add.catch.CV==TRUE){
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.7,0), tck = -0.02,cex=0.8)
  
  png(file = paste0(output.dir,"/Catch.fit_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  # estimated Catch
  estC = posteriors$estC
  predC = apply(estC,2,quantile,c(0.5,0.025,0.975)) 
  cord.x <- c(years,rev(years))
  cord.y<-c(predC[2,],rev(predC[3,]))
  plot(years,(TC),type="n",ylim=c(0,max(predC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",catch.metric),main="")
  polygon(cord.x,cord.y,col="gray",border=0,lty=1)
  lines(years,predC[1,],lwd=2,col=4)
  points(years,(TC),pch=21,bg=0,cex=1.5)
  legend("topright",c("Observed","Predicted"),pch=c(21,-1),bg=0,lwd=c(-1,2),col=c(1,4),bty="n")
  dev.off() 
  
}


#------------------------------
# Plot Posteriors
#------------------------------
sel.par = c(1,2,7,4,3,5)

out=data.frame(posteriors[params[sel.par]])
if(nSel>1) out=out[,-c(3:(3+nSel-2))]

node_id = names(out)
#informative priors
Prs = as.matrix(cbind(K.pr,r.pr,c(0,0),psi.pr))

#Posteriors
Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Posteriors_",assessment,"_",Scenario,".png"),width  = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)


node_id = names(out)


#par(mfrow=c(4,2),oma=c(0,1,1,0), mar=c(4,4,1,1))

for(i in 1:length(node_id))
{
  
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  
  if(i==1){
    
    rpr = rlnorm(10000,log(K.pr[1]),K.pr[2]) 
    pdf = stats::density(post.par,adjust=2)  
    prior = dlnorm(sort(rpr),log(K.pr[1]),K.pr[2])   
    plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab="K",ylab="",xaxs="i",yaxs="i",main="")
    
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
  }  
  
  
  if(i==2){
    
    rpr = rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
    pdf = stats::density(post.par,adjust=2) 
    prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
    plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
    
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
  }
  
  if(i==3){
    if(Model<4){
      plot(1,1,type="n",xlim=range(0.5,2.5),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
      abline(v=m,lwd=2)}
    if(Model==4){
      mpr = rlnorm(10000,log(m),shape.CV) 
      pdf = stats::density(post.par,adjust=2) 
      prior = dlnorm(sort(mpr),log(m),shape.CV)   
      plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
      
      polygon(c(sort(mpr),rev(sort(mpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      
      
    }
  }
  
  
  if(i==4){
    if(psi.dist=="beta"){
      parm = fitdist(post.par[post.par<1 & post.par>0.01], "beta")$estimate
      rpr = rbeta(10000,(psi.pr[1]),psi.pr[2]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dbeta(sort(rpr),psi.pr[1],psi.pr[2])   
    } else {
      rpr = rlnorm(10000,log(psi.prior[1]),psi.prior[2]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dlnorm(sort(rpr),log(psi.prior[1]),psi.prior[2])}
    plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,quantile(rpr,c(0.001,0.999)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    #legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
  }        
  
  if(i>4){
    if(sigma.proc!=TRUE & i==length(node_id)) {
      plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
      abline(v=sigma.proc^2,lwd=2)} else {
        
        pdf = stats::density(post.par,adjust=2)  
        plot(pdf,type="l",xlim=range(0,post.par),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
        if(i==length(node_id)& igamma[1]>0.9){
          rpr = 1/rgamma(10000,igamma[1],igamma[2])
          prior = stats::density(rpr,adjust=2)
          polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
        }
        
        polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
        #legend('topright',c("Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.8,0.6)),bty="n")
      } }         
  
}
mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
dev.off()   



#-----------------------------
# MCMC chains of posteriors
#-----------------------------

Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/MCMC_",assessment,"_",Scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)
for(i in 1:length(node_id)){
  
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=4)
  lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)   
}
dev.off()

#-----------------------------------------------------------
# <><<><<><<>< Produce JABBA Model Diagnostics ><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing JABBA Model Fit Diagnostics <><","\n"))


# extract predicted CPUE + CIs

N = n.years
series = 1:n.indices

check.yrs = apply(CPUE,1,sum,na.rm=TRUE)
cpue.yrs = years[check.yrs>0]

#CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Fits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- c(fit[1,yr],rev(fit[3,yr]))
  
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(years),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(Yr,fit[2,yr],lwd=2,col=1)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")}
  
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()

#log CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/logFits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = log(c(min(fit[,yr[is.na(CPUE[,i])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,i])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- log(c(fit[1,yr],rev(fit[3,yr])))
  
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(Yr,log(fit[2,yr]),lwd=2,col=4)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,log(cpue.i/mufit),ui=log(exp(log(cpue.i)+1.96*se.i)/mufit),li=log(exp(log(cpue.i)-1.96*se.i)/mufit),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")}
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()


if(CatchOnly!=TRUE){
# JABBA-residual plot
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Residuals_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

Resids = NULL
for(i in 1:n.indices){
  Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))   
}

plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(Resids,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
# get degree of freedom
Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
DF = Nobs-npar

RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)

legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))

dev.off()

#Save Residuals 
Res.CPUE = data.frame(Resids)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/ResCPUE_",assessment,"_",Scenario,".csv"))

#---------------------------------------
# Stadardized Residuals
#--------------------------------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/StandardizedResids_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)


# Standardized Residuals
StResid = NULL
for(i in 1:n.indices){
  StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                   apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))        
}

plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(StResid,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))


dev.off()

#Save standardized Residuals 
StRes.CPUE = data.frame(StResid)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/StResCPUE_",assessment,"_",Scenario,".csv"))
}

#------------------------------
# Plot process error deviation
#------------------------------

proc.dev = apply(posteriors$Proc.Dev,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/ProcDev_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(-0.22,0.22)#range(proc.dev)*1.1
cord.x <- c(years,rev(years))
cord.y <- c(proc.dev[1,],rev(proc.dev[3,]))
# Process Error
plot(years,proc.dev[2,],ylab="Process deviation P[t]",xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,proc.dev[2,],lwd=2)
lines(years,rep(0,length(years)),lty=5)


dev.off()

#-----------------------------------------------------------
# <><<><<><<><<>< JABBA Management Plots ><>><>><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing Management Plots <><","\n"))

#-----------------------
# Plot Biomass B_t
#-----------------------

B_t = posteriors$SB
mu.B = apply(B_t,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Biomass_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.B ))
cord.x <- c(years,rev(years))
cord.y <- c(mu.B [1,],rev(mu.B [3,]))

# B_t
plot(years,mu.B[2,],ylab=paste0("Biomass ",catch.metric),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.B[2,],lwd=2,col=1)
lines(years,rep(mean(posteriors$SBmsy),length(years)),lty=5)
text((max(years)-min(years))/30+years[1],mean(posteriors$SBmsy)*1.11,expression(paste(B[MSY])))
dev.off()



#-----------------------
# Plot B/Bmsy and H/Hmsy
#-----------------------

HtoHmsy = posteriors$HtoHmsy
BtoBmsy = posteriors$BtoBmsy

mu.f = apply(HtoHmsy,2,quantile,c(0.025,0.5,0.975))
mu.b = apply(BtoBmsy,2,quantile,c(0.025,0.5,0.975))

f = HtoHmsy[,N]
b = BtoBmsy[,N]
Par = list(mfrow=c(1,2),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/TrendMSY_",assessment,"_",Scenario,".png"), width = 7, height = 3, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.f))
cord.x <- c(years,rev(years))
cord.y <- c(mu.f[1,],rev(mu.f[3,]))

# H/Hmsy
plot(years,mu.f[2,],ylab=ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.f[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)

ylim = c(0, max(mu.b,1.1))

cord.x <- c(years,rev(years))
cord.y <- c(mu.b[1,],rev(mu.b[3,]))
plot(years,mu.b[2,],ylab=expression(paste(B/B[MSY])),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.b[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)
dev.off()

if(SP.plot!="phase"){
  #-----------------------------------------
  # Produce simple Production function plot
  #-----------------------------------------
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SP_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  Bit = seq(1,mean(posteriors$K),mean(posteriors$K)/500)
  SP = mean(posteriors$r)/(m-1)*Bit*(1-(Bit/median(posteriors$K))^(m-1))  
  B = apply(posteriors$SB,2,mean)
  MSY = quantile(posteriors$MSY,c(0.025,0.5,0.975)) 
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(Catch,na.rm=T),max(MSY*1.1)))),ylab=paste0("Surplus Production ",catch.metric),xlab="Biomass (t)")
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY[1],2),rep(MSY[3],2)),border = 0,col=grey(0.5,0.4))
  lines(Bit,SP,col=2,lwd=2)
  lines(B,Catch,lty=1)
  points(B,Catch,cex=0.5,pch=4)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],Catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.5)
  abline(h=max(SP),col=4)
  sel.years = c(min(years),years[sel.yr[2]],max(years))
  lines(rep(median(posteriors$SBmsy),2),c(-1000,max(SP)),lty=2,col=2)
  legend('topright', 
         c(expression(B[MSY]),"MSY","Catch",paste(sel.years)), 
         lty=c(2,1,1,1,1,1),pch=c(-1,-1,4,22,21,24),pt.bg=c(0,0,0,rep("white",3)), 
         col=c(2,4,rep(1,4)),lwd=1,cex=0.9,pt.cex=c(-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
}

if(SP.plot=="phase"){
  #-----------------------------------------
  # Produce JABBA SP-phase plot
  #-----------------------------------------
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SPphase_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  Bit = seq(1,median(posteriors$K),median(posteriors$K)/500)
  Cmsy = Bit*median(posteriors$Hmsy)
  SP = median(posteriors$r)/(m-1)*Bit*(1-(Bit/median(posteriors$K))^(m-1))  
  B = apply(posteriors$SB,2,mean)
  MSY = quantile(posteriors$MSY,c(0.025,0.5,0.975)) 
  Bmsy.sp = median(posteriors$SBmsy)
  K.sp = median(posteriors$K)
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(K.sp,0,max(SP),K.sp,K.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(Catch,na.rm=T)*1.05,max(MSY*1.1)))),xlim=c(0,max(Bit,B)),ylab=paste0("Surplus Production ",catch.metric),xlab="Biomass (t)",xaxs="i",yaxs="i")
  rect(0,0,K.sp*1.1,K.sp*1.1,col="green",border=0)
  rect(0,0,K.sp,K.sp,col="yellow",border=0)
  if(KOBE.type!="ICCAT") rect(0,max(SP),K.sp,K.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)  
    #i = i+1
  }
  
  gy.sp = Bit[Bit>Bmsy.sp]
  for(i in (length(ry.sp)+1):length(Bit)){
    
    #lines(rep(Bit[i],2),c(max(SP),Cmsy[i]),col=ifelse(i %% 2== 0,ifelse(KOBE.type=="ICCAT","yellow","orange"),"green"),lty=3)  
    #i = i+1
  }
  
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY[1],2),rep(MSY[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,Catch,lty=1,lwd=1)
  points(B,Catch,cex=0.8,pch=16)
  lines(Bit,Cmsy,col=1,lwd=1,lty=2)
  N=n.years
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],Catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(median(posteriors$SBmsy),2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright', 
         c(expression(B[MSY]),"MSY","SP","Catch",paste(sel.years)), 
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)), 
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
}


#----------------------
# Produce Kobe plot
#----------------------
# extract vectors BtoBmsy and FtoFmsy


if(KOBE.plot==TRUE){
  # prepare 
  
  # fit kernel function
  kernelF <- ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Kobe_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,max(1/BmsyK,mu.b[2,]) +0.05), ylim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,ylab=ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")
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
  polygon(c(1,100,100,1),c(1,1,100,100),col=ifelse(KOBE.type=="ICCAT","yellow","orange"),border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.b[2,],mu.f[2,],pch=16,cex=1)
  
  
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[2,],mu.f[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[2,sel.yr],mu.f[2,sel.yr],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  
  if(KOBE.type=="ICCAT"){               
    Pr.yellow = (sum(ifelse(b<1 & f<1,1,0))+sum(ifelse(b>1 & f>1,1,0)))/length(b)*100} else {
      Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
      Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
    }
  
  
  sel.years = c(years[sel.yr])
  ## Add legend
  if(KOBE.type=="ICCAT"){
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,7)),pch=c(22,21,24,rep(22,7)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.1,3)),bty="n")
  }else{
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")  
    
  }
  dev.off()
}


if(Biplot==TRUE){
  
  #---------------------------------------------------------
  # Produce 'post-modern' biplot (see Quinn and Collie 2005)
  #---------------------------------------------------------
  
  # read ftarget,bthreshold
  ftarget<-0.8
  bthreshold<-0.2
  
  # fit kernel function
  kernelF <- ci2d(f,b,nbins=201,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Biplot_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", ylim=c(0,1/BmsyK+0.05), xlim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,xaxs="i",yaxs="i")
  
  # and fill areas using the polygon function
  fint = seq(0.001,100,0.01)
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
  lines(mu.f[2,],mu.b[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.f[2,sel.yr],mu.b[2,sel.yr],col=
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
  mtext(ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)
  
  text(0.65,2.4,paste0(perGreen,"%"))
  text(0.9,2.4,paste0(perYellow,"%"))
  text(1.2,2.4,paste0(perRed,"%"))
  
  dev.off()
  
}



#--------------------------------------
# Plot projections and safe posteriors
#--------------------------------------
if(Projection ==TRUE){
  cat(paste0("\n","><> Producing Future TAC Projections <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Projections_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  
  proj.yrs =  years[n.years]:(years[n.years]+pyrs)
  # Dims 1: saved MCMC,2: Years, 3:alternatic TACs, 4: P, H/Hmsy, B/Bmsy
  projections = array(NA,c(nsaved,length(proj.yrs),nTAC,3))
  for(i in 1:nTAC){
    projections[,,i,1] = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
  }
  
  for(i in 1:nTAC){
    projections[,,i,2] = cbind(posteriors$BtoBmsy[,(n.years):n.years],posteriors$prBtoBmsy[,,i])
  }
  for(i in 1:nTAC){
    projections[,,i,3] = cbind(posteriors$HtoHmsy[,(n.years):n.years],posteriors$prHtoHmsy[,,i])
  }
  
  kjp = kobeJabbaProj(projections,proj.yrs[1])
  
  save(kjp,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_projections"))
  
  # Change here for ICCAT Bmsy plot
  Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,nTAC])
  
  #plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)+length(proj.yrs)*0.2),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  
  cols = rev(seq(0.4,0.9,0.5/nTAC))
  plot.order = (1:nTAC)
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    polygon(c(proj.yrs,rev(proj.yrs)),c(apply(Traj,2,quantile,0.05),rev(apply(Traj,2,quantile,0.95))),col=grey(cols[i],1),border=NA)
  }
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    lines(proj.yrs,apply(Traj,2,median),col=rev(jabba.colors[j]),lwd=2)
  }
  lines(proj.yrs[1:(imp.yr-yr.last+1)],apply(Traj,2,median)[1:(imp.yr-yr.last+1)],col=1,lwd=2)
  BmsyK=(m)^(-1/(m-1))
  abline(h=BmsyK,lty=2,lwd=2)
  legend("topleft",paste(TACs,"(t)"),col=(jabba.colors[1:nTAC]),lwd=2,cex=0.8)   
  
  dev.off()    
  
}


if(save.trajectories==TRUE){
  cat(paste0("\n","><> Saving Posteriors of FRP trajectories <><","\n"))
  
  # FRP trajectories
  trajectories = array(NA,c(nsaved,n.years,3))
  trajectories[,,1] = posteriors$P 
  trajectories[,,2] = posteriors$BtoBmsy 
  trajectories[,,3] = posteriors$HtoHmsy
  
  kb=kobeJabba(trajectories,years[1])
  save(kb,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_trajectories"))
  
}


cat(paste0("\n","\n","><> Scenario ",Mod.names,"_",Scenario," for ",assessment," - DONE! <><","\n"))





