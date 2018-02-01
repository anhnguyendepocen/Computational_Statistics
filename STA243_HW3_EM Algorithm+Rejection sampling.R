

###  STA 243 HW3
# problem 2
alpha=1
d=1
iter=0
while (abs(d)>10^-8){
  iter=iter+1
  alpha1=sqrt((22.36+2*alpha^2)/4)
  d=abs((alpha1-alpha)/alpha)
  alpha=alpha1
}
alpha
iter

# STA243 HW3 Q3

u=runif(5000)
f = function(x){-log(1-x*(1-exp(-2)))}
x.in <- f(u)
x.true <- function(x){exp(-x)/(1-exp(-2))}

plot(density(x.in),main = 'Plot of estimated density v.s. true density')
curve(x.true, from = -.5, to = 2.5,add = TRUE, lty = 2, col = 2)

# problem 4
qx=function(x){
  exp(-x)/(1+x^2)
}
g1=function(x){exp(-x)}
g2=function(x){2/pi/(1+x^2)}
alpha1=1
x1=vector()
ptm<-proc.time()
for (i in 1:10000){
  u=runif(1,0,1)
  x=runif(1,0,1)
  x=-log(x)
  while (u>qx(x)/alpha1/g1(x)){
    u=runif(1,0,1)
    x=runif(1,0,1)
    x=-log(x)
  }
  x1[i]=x
}
proc.time()-ptm
X1=sample(x1,5000)
plot(density(X1),xlim=c(0,5),main='The estimated density of sample X with g1(x)')

alpha2=pi/2
x2=vector()
ptm2<-proc.time()
for (i in 1:10000){
  u=runif(1,0,1)
  x=runif(1,0,1)
  x=tan(pi*x/2)
  while (u>qx(x)/alpha2/g2(x)){
    u=runif(1,0,1)
    x=runif(1,0,1)
    x=tan(pi*x/2)
  }
  x2[i]=x
}
proc.time()-ptm2
X2=sample(x2,5000)
plot(density(X2),xlim=c(0,5),main='The estimated density of sample X with g2(x)')

integrate(qx,0,Inf)
