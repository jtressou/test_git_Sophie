# je fouine dans ton code

nYear=30; 
biomass <- rep(NA,nYear); capt <- rep(NA,nYear)
K=300; 
r = 1.02; 
q=0.5

set.seed(15); 
#set.seed(sample(1:100,1))
sigmaBiomass <- 0.3 
biomass[1] <- rnorm(1,mean = 0.9*K,sd=K/10)
capt[1] <- runif(1,0.1,0.5)*biomass[1]
for (year in 2:nYear){
  Bk <- biomass[year-1]; 
  capt[year-1] <- runif(1,0.1,0.5)*Bk
  innov = rnorm(1,mean=-sigmaBiomass^2/2,sd = sigmaBiomass)
  biomass[year] <- ((Bk+r*Bk*(1-Bk/K))-capt[year-1])*exp(innov)
  capt[year] <- runif(1,0.1,0.5)*biomass[year]
}

obs_error <- rnorm(nYear,mean =-sigmaBiomass^2/2,sd = sigmaBiomass)
abundanceIndex<- q*biomass*exp(obs_error)

plot(1:nYear,abundanceIndex,type='b')



################ FILTRAGE 
M = 1000; 
xi <- matrix(NA,M,nYear)
w <- matrix(NA,M,nYear)
xi[,1] <-abundanceIndex[1]/q*exp(-rnorm(M,mean =-sigmaBiomass^2/2,sd = sigmaBiomass))
ww =dnorm(xi[,1],mean = 0.9*K,sd=K/10,log=TRUE)
w[,1] = exp(ww-mean(ww))/sum(exp(ww-mean(ww)))


#plot(density(xi[,1],weights = w[,1]))
#abline(v=biomass[1])
#lines(density(xi[,1]),col='red')
ESS = rep(0,nYear)
ESS[1] = 1/sum(w[,1]^2)

geneal = matrix(0,M,nYear)
geneal[,1] = 1:M
for (t in 2:nYear){
  
  S = sample(1:M,M,replace=TRUE,prob=w[,t-1])
  geneal = geneal[S,]
  geneal[,t] = 1:M
  
  xi_previous = xi[S,t-1]
     
  omega.proposal.t <- sigmaBiomass^2/2 
  
  mean.prior = (xi_previous+r*xi_previous*(1-xi_previous/K)-capt[t-1])
  cond = which(mean.prior>0)
  if(length(cond)<M){w[-cond,t] = 0}
  mu.0= log(mean.prior[cond])- sigmaBiomass^2/2
  mu.proposal.t <- (mu.0 + (log(abundanceIndex[t])-log(q)+sigmaBiomass^2/2))/2

  
  xi[cond,t] = rlnorm(length(cond),meanlog = mu.proposal.t,sdlog = sqrt(omega.proposal.t))
  xi[-cond,t]=0
  update = dlnorm(abundanceIndex[t]/(q*xi[cond,t]), meanlog=-sigmaBiomass^2/2,sdlog = sigmaBiomass,log=TRUE)
  update2 = dlnorm(xi[cond,t]/mean.prior[cond], meanlog=-sigmaBiomass^2/2,sdlog = sigmaBiomass,log=TRUE)
  update3 = dlnorm(xi[cond,t],meanlog = mu.proposal.t,sdlog = sqrt(omega.proposal.t),log = TRUE)
  log.ww <-  update + update2- update3
  w[cond,t] = exp(log.ww-max(log.ww))/sum( exp(log.ww-max(log.ww)))
  ESS[t]  = 1/sum(w[,t]^2)
  
}
#####
matplot(t(geneal),type='l')
plot(ESS,type='b')
mean_filtr = colMeans(xi[sample(1:M,M,replace=TRUE,prob=w[,30]),])
S = sample(1:M,20,replace=TRUE,prob=w[,30])

echan = xi[S,]
matplot(t(echan),type='l',col='grey')
lines(biomass,lwd=3,col='red',type='b')
lines(mean_filtr,lwd=3,col='blue',type='b')


############### LISSAGE
### calcul des poids de lissage

# 
l=29
xi_l= xi[,l]; xi_lp1 = xi[,l+1]
kernel_m = function(xi_l,xi_lp1,l){
  M  = length(xi_lp1)
  logmean = xi_l+r*xi_l*(1-xi_l/K)-capt[l]
  ld = vapply(1:M,function(i){dlnorm(xi_lp1[i]/logmean,meanlog=-sigmaBiomass^2/2,sdlog = sigmaBiomass)},rep(1,M))
  return(ld)
}

wliss = matrix(0,M,nYear)
wliss[,nYear] = w[,nYear]
for (l in (nYear-1):1){
  m =kernel_m(xi[,l],xi[,l+1],l)
  G = matrix(matrix(w[,l],nrow=1)%*%m,M,M,byrow=TRUE)
  wliss[,l]=w[,l]*(m/G)%*%  wliss[,l+1]
   
}
colSums(wliss)  
# 
#   
#   
  
  
