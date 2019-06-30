###JAGS script
###running the Correlated Chain-Ladder (CCL) model
#estimate parameters alpha and beta
#predicting logloss, mu[i] and logloss[i]

nburn=10000

modelString = "model {
mu[1]<-alpha[w[1]]+beta[d[1]]
logloss[1]~dnorm(mu[1],1/sig2[1])
for (i in 2:length(w)){
mu[i]<-alpha[w[i]]+beta[d[i]]+rho*(logloss[i-1]-mu[i-1])*wne1[i]
logloss[i]~dnorm(mu[i],1/sig2[i])
}

### set up sig2
for (i in 1:length(w)){
sig2[i]<-sigd2[d[i]]
}

for (j in 1:10){
sigd2[j]<-sum(a[j:10])
}
for (k in 1:10){
a[k]~dunif(0.000001,1)
}
#
# specify priors
#
for (i in 1:numlev){
alpha[i]~dnorm(log(premium[i])+logelr,.1)
}
logelr~dnorm(logelr_mean,logelr_sig)
logelr_mean ~ dnorm(0,1)
logelr_sig ~ dunif(0,2)
#
for (i in 1:9){
beta[i]~dnorm(beta_mu[i],beta_sig[i])
beta_mu[i] ~ dnorm(0,1)
beta_sig[i] ~ dunif(0,2)
}
beta[10]<-0 
rho~dunif(-1,1)
# rho~dunif(-.00001,.00001) # Use for LCL model
}" 

###Initialize JAGS model
inits1=list(.RNG.name= "base::Wichmann-Hill",
            .RNG.seed= 12341)
inits2=list(.RNG.name= "base::Marsaglia-Multicarry",
            .RNG.seed= 12342)
inits3=list(.RNG.name= "base::Super-Duper",
            .RNG.seed= 12343)
inits4=list(.RNG.name= "base::Mersenne-Twister",
            .RNG.seed= 12344)

data.for.jags=list(premium= premium[1:10],
                   logloss = log(rloss),
                   numlev  = numw,
                   w       = rdata$w,
                   wne1    = rdata$wne1,
                   d       = rdata$d)

###run the model
nthin=2
maxpsrf=2
while (maxpsrf>1.05){
  nthin=nthin*2
  print(paste("nthin =",nthin))
  jagout=run.jags(model=modelString,monitor=c("alpha","beta[1:9]","sigd2","rho"),
                  data=data.for.jags,n.chains=4,method="parallel",
                  inits=list(inits1,inits2,inits3,inits4),thin=nthin,silent.jags=F,
                  plots=TRUE,burnin=nburn,sample=2500,psrf.target=1.05)
  gelman=gelman.diag(jagout)
  maxpsrf=max(gelman$psrf[,1])
  print(paste("maxpsrf =",maxpsrf))
} 
