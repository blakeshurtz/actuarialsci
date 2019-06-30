###is there a better way to import these libraries? Packrat?
library(coda)
library(ChainLadder)
library(runjags)

#models are executed on a per-insurer basis. Define variables for industry and group. See Appendix A for groups included in model
setwd("data")
insurer.data="comauto_pos.csv"
grpcode="353"

#import data
a=read.csv(insurer.data)
setwd("..")

#extract Schedule P triangle data
ins.line.data=function(g.code){
  b=subset(a,a$GRCODE==g.code)
  name=b$GRNAME; grpcode=b$GRCODE
  w=b$AccidentYear; d=b$DevelopmentLag
  cum_incloss=b[,6]; cum_pdloss=b[,7]; bulk_loss=b[,8]
  dir_premium=b[,9]; ced_premium=b[,10]; net_premium=b[,11]
  single=b[,12]; posted_reserve97=b[,13]  
  # get incremental paid losses - assume data is sorted by ay and lag
  inc_pdloss=numeric(0)
  for (i in unique(w)){
    s=(w==i)
    pl=c(0,cum_pdloss[s])
    ndev=length(pl)-1  
    il=rep(0,ndev)
    for (j in 1:ndev){            
      il[j]=pl[j+1]-pl[j]
    }
    inc_pdloss=c(inc_pdloss,il)
  }
  data.out=data.frame(grpcode,w,d,net_premium,dir_premium,ced_premium,
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,posted_reserve97)
  return(data.out)
}

#extract triangle
cdata=ins.line.data(grpcode)
w=cdata$w-1987; d=cdata$d #index year and development lag at 1
o1=100*d+w; o=order(o1); w=w[o]; d=d[o] #order the data by development, then year
premium=cdata$net_premium[o] #premium
cpdloss=cdata$cum_pdloss[o] #cumulative paid loss
cpdloss=pmax(cpdloss,1) #parallel max order
wne1=ifelse(w==1,0,1) #first year, yes or no

#create new variable incremental loss- incurred losses net of reinsurance, this is what shows up in triangle
incloss=cdata$cum_incloss[o]-cdata$bulk_loss[o] 
incloss=pmax(incloss,1) 

#data frame containing both current known and future known losses 
adata=data.frame(grpcode,w,d,premium,cpdloss,incloss,wne1)
write.csv(adata, "353.csv")

#construct SECOND data frame for model, containing current known outcomes only 
rdata=subset(adata,(adata$w+adata$d)<12)
numw=length(unique(rdata$w))
rloss=rdata$incloss; aloss=adata$incloss

