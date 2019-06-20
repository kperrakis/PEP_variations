library(rmutil)
library(twopiece)

included.vars<-c(1,3,4)
b<-c(0.3, 0.5, 1.0)
s<-1.0
s_laplace<-1/sqrt(2)
s_1<-(1-0.5)
s_2<-(1+0.5)

# choose error distribution

errors<-'normal'
errors<-'a-normal'
errors<-'laplace'
errors<-'h-normal'

sim1<-list()
N<-c(15,30,50,100,500,1000)
p<-10
R<-100

set.seed(1)

for(k in 1:length(N))
{
sim1[[k]]<-list()
n<-N[k]
x<-matrix(rnorm(n*p), ncol=p, nrow=n)
y<-matrix(nrow=n, ncol=R)
for (i in 1:R)
{
if(errors=='normal') { y[,i] <-  x[,included.vars]%*%b  + rnorm(n, 0, s)}
if(errors=='laplace') { y[,i] <-  x[,included.vars]%*%b  + rlaplace(n, 0, s_laplace)}
if(errors=='a-normal') {y[,i] <-  x[,included.vars]%*%b  + rtp3(n, 0, s_1, s_2, rnorm)}
  if (errors == 'h-normal') {
    n.errors <- rnorm(n, 0, s)
    h.errors <- exp(x[, included.vars] %*% b) * n.errors
    constant <- sqrt(c(var(n.errors) / var(h.errors)))
    h.errors <- h.errors * constant
    y[, i] <-  x[, included.vars] %*% b  + h.errors
  }
}
sim1[[k]]$x<-x
sim1[[k]]$y<-y
}

names(sim1)<-paste("n=", N, sep="")
sim1$N<-N
sim1$true.effects<-included.vars
sim1$true.model.index <- sum(2^(p-included.vars))+1
sim1$true.beta<-b
sim1$true.s <- s