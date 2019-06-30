###optional diagnostics
crosscorr.plot(jagout$mcmc)

###parameter estimates
summary(jagout)
#extract information from jags output to process in R
b=as.matrix(jagout$mcmc)
alpha=b[,1:10]
beta=cbind(b[,11:19],rep(0,dim(b)[1]))
rho=b[,20]
sigd2=b[,20:29]

###simulate loss statistics by accident year reflecting total risk for the JAGS model
set.seed(12345)
Premium=subset(rdata,rdata$d==1)$premium
at.wd10=matrix(0,dim(b)[1],10) #almost all 0's
ss.wd10=rep(0,10) #all 0s
ms.wd10=rep(0,10) #all 0s
mu.wd11mw=matrix(0,dim(b)[1],10) #all 0s
mu.wd10=matrix(0,dim(b)[1],10) #all 0s
at.wd10[,1]=rep(rloss[55],dim(b)[1]) #almost all 0's

###now, fill in the matrix
#warnings why? NA's produced! Maybe because of log...
mu.wd10[,1]=alpha[,1]+beta[,10]
for (w in 2:10){
  mu.wd10[,w]=alpha[,w]+beta[,10]+rho*(log(at.wd10[,w-1])-mu.wd10[,w-1])
  for (i in 1:dim(b)[1]){
    at.wd10[i,w]=rlnorm(1,mu.wd10[i,w],sqrt(sigd2[i,10]))  
  }
}

#take means
ms.wd10[1]=mean(at.wd10[,1])
for (w in 2:10){
  ms.wd10[w]=mean(at.wd10[,w])
  ss.wd10[w]=sd(at.wd10[,w])
}

Pred.CCL=rowSums(at.wd10)
ms.td10=mean(Pred.CCL)
ss.td10=sd(Pred.CCL)
CCL.Estimate=round(ms.wd10)
CCL.S.E.=round(ss.wd10)
CCL.CV=round(CCL.S.E./CCL.Estimate,4)
act=sum(subset(aloss,adata$d==10)[1:10])
pct.CCL=sum(Pred.CCL<=act)/length(Pred.CCL)*100


###put CCL accident year statistics into a data frame
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
par(mfrow=c(1,1))
hist(rho,main="",xlim=c(-1,1),xlab=expression(rho))
hist(Pred.CCL,main="Predictive Distribution of Outcomes",xlab="")

###optional data export
setwd("results")
write.csv(risk,file="outfile.csv",row.names=F)
setwd("..")
