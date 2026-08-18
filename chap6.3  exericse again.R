n<-10^6
myq1<-function(x){
res<-exp(-abs(x)^3/3)
return(res)}


la<-seq(-3,3,0.1)
myq1(la)/dnorm(la)  

min(myq1(la)/dnorm(la))
1/min(myq1(la)/dnorm(la)) #so the const must be at least 36

y<-rnorm(n)
u<-runif(n)
m<-ifelse(u<=(myq1(y)/(dnorm(y)*36)),1,0)
yy<-y[m==1]
mean(yy^2)

## the above envelope is not efficient.

mycheck<-function(sigma=1,c1=1){
	myseq<-seq(-3,3,0.01)
	env1<-c1*dnorm(myseq,0,sigma)
	plot(myseq,myq1(myseq),type="l",ylim=c(0,max(c(myq1(myseq),env1))))
	lines(myseq,c1*dnorm(myseq,0,sigma),type="l",col=2)
	res<-c1*dnorm(myseq,0,sigma)/myq1(myseq)
	return(min(res))
}

aa<-mycheck(1.5,4)

n<-10^6
y<-rnorm(n)
u<-runif(n)
m<-ifelse(u<=(myq1(y)/(dnorm(y)*36)),1,0)
yy<-y[m==1]
mean(yy^2)
length(yy)/n   #acceptance rate

n<-10^7
y<-rnorm(n,0,1.5)
u<-runif(n)
m<-ifelse(u<=(myq1(y)/(dnorm(y,0,1.5)*4)),1,0)
yy<-y[m==1]
mean(yy^2)
length(yy)/n   #acceptance rate

aa<-mycheck(1,1)
aa<-mycheck(1,2)
aa<-mycheck(1,3)  #looks ok

n<-10^7
y<-rnorm(n,0,1)
u<-runif(n)
m<-ifelse(u<=(myq1(y)/(dnorm(y,0,1)*3)),1,0)
yy<-y[m==1]
mean(yy^2)
length(yy)/n   #acceptance rate


## true value?
c1<-integrate(myq1,-Inf,Inf)$value



myq2<-function(x){
res<-x^2*exp(-abs(x)^3/3)/c1
return(res)}

integrate(myq2,-Inf,Inf)

n<-10^6
x<-rnorm(n)
x<-sort(x)
wstar<-myq1(x)/dnorm(x)
w<-wstar/sum(wstar)
sum(x^2*w)


x2<-x[1:(n-1)]
aa<-sum(diff(x)*x2^2*myq1(x2))
bb<-sum(diff(x)*myq1(x2))
aa/bb

