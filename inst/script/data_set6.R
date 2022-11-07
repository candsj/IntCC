n1=20
n2=20
n3=20
n = n1+n2+n3
p=500
p1 =83
p2 = 83
p3 = 83

true.class=c(rep(1,n1),rep(2,n2),rep(3,n3))

data_set6=alist()

for (i in 1:100) {
  ################# normal distribution ####################
  c1=matrix(rnorm(n1*p1,mean=1.36,sd=3), ncol=p1,nrow=n1)
  c2=matrix(rnorm(n2*p2,mean=1.18, sd=3), ncol=p2,nrow=n2)
  c3=matrix(rnorm(n3*p3,mean=1, sd=3), ncol=p3,nrow=n3)
  
  normData = matrix(rnorm(n*p, mean=0, sd=3),nrow=n, ncol=p)
  normData[1:n1,1:p1]= c1
  normData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  normData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3
  
  
  
  ######## binomial data #########################
  c1=matrix(rbinom(n1*p1,size=1, prob=0.38), ncol=p1,nrow=n1)
  c2=matrix(rbinom(n2*p2,size=1, prob=0.34), ncol=p2,nrow=n2)
  c3=matrix(rbinom(n3*p3,size=1, prob=0.3), ncol=p3,nrow=n3)
  
  binomData = matrix(rbinom(n*p, size=1, prob=0.2),nrow=n, ncol=p)
  binomData[1:n1,1:p1]= c1
  binomData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  binomData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3
  
  
  
  
  ######### poisson data ############################
  c1=matrix(rpois(n1*p1,lambda=1.3), ncol=p1,nrow=n1)
  c2=matrix(rpois(n2*p2,lambda=1.15), ncol=p2,nrow=n2)
  c3=matrix(rpois(n3*p3,lambda=1), ncol=p3,nrow=n3)
  
  poisData = matrix(rpois(n*p, lambda=0.8),nrow=n, ncol=p)
  poisData[1:n1,1:p1]= c1
  poisData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  poisData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3
  
  
  
  
  ################## multinomial distribution ##############
  c1=matrix(sample(1:3,size=n1*p1, replace=TRUE, prob=c(0.5,0.3,0.2)), ncol=p1,nrow=n1)
  c2=matrix(sample(1:3,size=n2*p2, replace=TRUE, prob=c(0.3,0.4,0.3)), ncol=p2,nrow=n2)
  c3=matrix(sample(1:3,size=n3*p3, replace=TRUE, prob=c(0.2,0.3,0.5)), ncol=p3,nrow=n3)
  
  multData = matrix(sample(1:3,size=n*p, replace=TRUE, prob=c(0.33,0.33,0.33)),nrow=n, ncol=p)
  multData[1:n1,1:p1]= c1
  multData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  multData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3
  
  
  
  
  
  #biID = seq(1,p,2)
  #multData[,biID] = binomData[,biID]
  
  ### some variables only have one category; need to remove these variables
  ncat = function(x){length(unique(x))}
  len = apply(binomData,2,ncat)
  unicat = which(len==1)
  if (length(unicat)>0) {
    binomData = binomData[,-unicat]
  }
  
  data_for_klic <- list()
  data_for_klic[[1]]=normData
  data_for_klic[[2]]=binomData
  data_for_klic[[3]]=poisData
  data_for_klic[[4]]=multData
  
  data_set6[[i]]=data_for_klic
  
}

save(data_set6,file="/home/mercury02/huchuang/wcc/data/data_set6.RData")

