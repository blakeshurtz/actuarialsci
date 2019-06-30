###JAGS script- running the Correlated Chain-Ladder (CCL) model
#estimate parameters alpha, beta, mu, logloss
#then predict incurred losses in year/development lag triangle


modelString = "model {
mu[1]<-alpha[w[1]]+beta[d[1]]
logloss[1]~dnorm(mu[1],1/sig2[1])
for (i in 2:length(w)){
mu[i]<-alpha[w[i]]+beta[d[i]]+rho*(logloss[i-1]-mu[i-1])*wne1[i]
logloss[i]~dnorm(mu[i],1/sig2[i])
}
#
# set up sig2
#
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
logelr~dunif(-1.5,0.5)
#
for (i in 1:9){
beta[i]~dunif(-5,5)
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
nburn=10000

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

#parameter estimates
summary(jagout)
b=as.matrix(jagout$mcmc) #extract information from jags output to process in R

#parameter vectors 
alpha=b[,1:10]
beta=cbind(b[,11:19],rep(0,dim(b)[1]))
rho=b[,20]
sigd2=b[,20:29]

#create loss incurred triangle as matrix
set.seed(12345)
Premium=subset(rdata,rdata$d==1)$premium
at.wd10=matrix(0,dim(b)[1],10) #almost all 0's
ss.wd10=rep(0,10) #all 0s
ms.wd10=rep(0,10) #all 0s
mu.wd11mw=matrix(0,dim(b)[1],10) #all 0s
mu.wd10=matrix(0,dim(b)[1],10) #all 0s
at.wd10[,1]=rep(rloss[55],dim(b)[1]) #almost all 0's

#fill in the matrix
mu.wd10[,1]=alpha[,1]+beta[,10]
for (w in 2:10){
  mu.wd10[,w]=alpha[,w]+beta[,10]+rho*(log(at.wd10[,w-1])-mu.wd10[,w-1])
  for (i in 1:dim(b)[1]){
    at.wd10[i,w]=rlnorm(1,mu.wd10[i,w],sqrt(sigd2[i,10]))  
  }
}

#calculate mean precited cost for year/lag 
ms.wd10[1]=mean(at.wd10[,1])
for (w in 2:10){
  ms.wd10[w]=mean(at.wd10[,w])
  ss.wd10[w]=sd(at.wd10[,w])
}

#generate summary stats
Pred.CCL=rowSums(at.wd10)
ms.td10=mean(Pred.CCL)
ss.td10=sd(Pred.CCL)
CCL.Estimate=round(ms.wd10)
CCL.S.E.=round(ss.wd10)
CCL.CV=round(CCL.S.E./CCL.Estimate,4)
act=sum(subset(aloss,adata$d==10)[1:10])
pct.CCL=sum(Pred.CCL<=act)/length(Pred.CCL)*100
pct.CCL

#create prediction table as data frame
W=c(1:10,"Total")
CCL.Estimate=c(CCL.Estimate,round(ms.td10))
CCL.S.E.=c(CCL.S.E.,round(ss.td10))
CCL.CV=c(CCL.CV,round(ss.td10/ms.td10,4))
Premium=c(Premium,sum(Premium))
Outcome=subset(aloss,adata$d==10)
Outcome=c(Outcome,sum(Outcome))
Group=rep(grpcode,11)
CCL.Pct=c(rep(NA,10),pct.CCL)
risk=data.frame(W,Premium,CCL.Estimate,CCL.S.E.,CCL.CV,Outcome,CCL.Pct)
print(risk)

#histograms
#hist(rho,main="",xlim=c(-1,1),xlab=expression(rho))
#hist(Pred.CCL,main="Predictive Distribution of Outcomes",xlab="")
#crosscorr.plot(jagout$mcmc)

#optional data export
#setwd("results"); write.csv(risk,file="outfilestan5.csv",row.names=F); setwd("..")

