#The following code is a derivative work of the test code from the iClusterPlus package.
#This simulation data contains a list of 4 datasets which are normal, binomial, Poisson and multinomial.

n1=20
n2=20
n3=20
n = n1+n2+n3
p=30
p1 =5
p2 = 5
p3 = 5

true.class=c(rep(1,n1),rep(2,n2),rep(3,n3))
exampleData=alist()

  ################# normal distribution ####################
  set.seed(1)
  c1=matrix(rnorm(n1*p1,mean=2.2,sd=1.7), ncol=p1,nrow=n1)
  set.seed(2)
  c2=matrix(rnorm(n2*p2,mean=1.6, sd=1.7), ncol=p2,nrow=n2)
  set.seed(3)
  c3=matrix(rnorm(n3*p3,mean=1, sd=1.7), ncol=p3,nrow=n3)

  set.seed(4)
  normData = matrix(rnorm(n*p, mean=0, sd=1.7),nrow=n, ncol=p)
  normData[1:n1,1:p1]= c1
  normData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  normData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

  ######## binomial data #########################
  set.seed(1)
  c1=matrix(rbinom(n1*p1,size=1, prob=0.54), ncol=p1,nrow=n1)
  set.seed(2)
  c2=matrix(rbinom(n2*p2,size=1, prob=0.42), ncol=p2,nrow=n2)
  set.seed(3)
  c3=matrix(rbinom(n3*p3,size=1, prob=0.3), ncol=p3,nrow=n3)

  set.seed(4)
  binomData = matrix(rbinom(n*p, size=1, prob=0.1),nrow=n, ncol=p)
  binomData[1:n1,1:p1]= c1
  binomData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  binomData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

  ######### poisson data ############################
  set.seed(1)
  c1=matrix(rpois(n1*p1,lambda=1.9), ncol=p1,nrow=n1)
  set.seed(2)
  c2=matrix(rpois(n2*p2,lambda=1.5), ncol=p2,nrow=n2)
  set.seed(3)
  c3=matrix(rpois(n3*p3,lambda=1), ncol=p3,nrow=n3)

  set.seed(4)
  poisData = matrix(rpois(n*p, lambda=0.6),nrow=n, ncol=p)
  poisData[1:n1,1:p1]= c1
  poisData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
  poisData[(n1+n2+1):n,(p1+p2+1):(p1+p2+p3)]= c3

  ################## multinomial distribution ##############
  set.seed(1)
  c1=matrix(sample(1:3,size=n1*p1, replace=TRUE, prob=c(0.8,0.15,0.05)), ncol=p1,nrow=n1)
  set.seed(2)
  c2=matrix(sample(1:3,size=n2*p2, replace=TRUE, prob=c(0.3,0.4,0.3)), ncol=p2,nrow=n2)
  set.seed(3)
  c3=matrix(sample(1:3,size=n3*p3, replace=TRUE, prob=c(0.05,0.15,0.8)), ncol=p3,nrow=n3)

  set.seed(4)
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

  exampleData[[1]]=normData
  exampleData[[2]]=binomData
  exampleData[[3]]=poisData
  exampleData[[4]]=multData

save(exampleData,file="inst/extdata/exampleData.RData")
