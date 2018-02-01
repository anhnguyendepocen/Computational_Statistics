#### STA 243 HW4
##Problem2
fi=function(x){
  exp(-x^2/2)/sqrt(2*pi)
}
integrate(fi,1,2)

hw=function(x,v)
{
  if(x>=1&x<=2)
  {
    v*exp(-x^2/2+(x-1.5)^2/2/v^2)
  }
  else
  {
    0
  }
}

mu_hat=function(v,N)
{
  x=rnorm(N,1.5,v)
  hw_vector=vector()
  for(i in 1:length(x))
  {
    hw_vector[i]=hw(x[i],v)
  }
  output=list()
  output$mu=hw_vector
  output$mu_hat=mean(hw_vector)
  output
}
N=100000
muoutput1=mu_hat(0.1,N)
muoutput2=mu_hat(1,N)
muoutput3=mu_hat(10,N)
hist(muoutput1$mu,main='histogram of importance sampling with v=0.1',xlab = 'sampled value')
muoutput1$mu_hat
hist(muoutput2$mu,main='histogram of importance sampling with v=1',xlab = 'sampled value')
muoutput2$mu_hat
hist(muoutput3$mu,main='histogram of importance sampling with v=10',xlab = 'sampled value')
muoutput3$mu_hat



###Problem6
b = c(seq(0.1,1,0.1), seq(1:20))
a = b
mean1 = matrix(0, nrow=30, ncol=30)
mean2 = matrix(0, nrow=30, ncol=30)
f = function(theta1,theta2,x)
{
  x^(-3/2)*exp(-theta1*x-theta2/x+2*sqrt(theta1*theta2)+log(sqrt(2*theta2)))
}

N=1000
for(i in a)
{
  for(j in b)
  {
    y = rgamma(N, shape = i, rate = j)
    
    u = runif(N,0,1)
    z = numeric(N)
    
    theta1=1.5
    theta2=2
    
    z[1] = sqrt(theta2/theta1)
    for(k in 2:N)
    {
      r=1
      rtemp = (f(theta1,theta2,y[k]) * dgamma(z[k-1], shape = i, rate = j))/(f(theta1,theta2,z[k-1]) * dgamma(y[k], shape = i, rate = j))
      
      r=min(r,rtemp)  
      
      
      if(r > u[k])
      { z[k] = y[k] }
      else
      { z[k] = z[k-1] }
      
    }
    mean1[i,j] = mean(z)
    mean2[i,j] = mean(1/z)
  }
}


truemean1 = sqrt(theta2/theta1)
truemean2 = sqrt(theta1/theta2) + 1/(2*theta2)

error1 = abs(mean1-truemean1)/truemean1
error2 = abs(mean2-truemean2)/truemean2
which((error1+error2)==min(error1+error2), arr.ind = T)
min(error1+error2) 