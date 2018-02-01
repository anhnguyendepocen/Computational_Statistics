#### STA 243 HW1
# Problem 1
### (1c)
x=c(-13.87,-2.53,-2.44,-2.40,-1.75,-1.34,-1.05,-0.23,-0.07,0.27,1.77,2.76,3.29,3.47,3.71,3.80,4.24,4.53,43.21,56.75)
n=length(x)
theta=seq(-50,50,0.01)
ltheta<-matrix(0,nrow=length(theta),ncol=1)
for (i in 1:length(theta)){
  ltheta[i]=-n*log(pi)-sum(log(1+(theta[i]-x)^2))
} 
plot(theta,ltheta,type='l',xlab=expression(theta),ylab='Log likelihood')

### (1d)
mle_newton<-function(thetaguess,X,maxiter,tol,convergencefactor){
  mle=matrix(0,nrow=length(thetaguess),ncol=1)
  for (i in 1:length(thetaguess)){
    theta0=thetaguess[i];
    iter=1;err=1;
    n=length(X)
    F=-2*sum((theta0-X)/(1+(theta0-X)^2))
    DF=-2*sum((1-(theta0-X)^2)/(1+(theta0-X)^2)^2)
    for (iter in 1:maxiter+1){
      theta1=theta0-F/DF;
      err=abs(theta1-theta0)/(abs(theta0)+convergencefactor)
      if((err<tol)||(is.na(err))) break
      theta0=theta1
      F=-2*sum((theta0-X)/(1+(theta0-X)^2))
      DF=-2*sum((1-(theta0-X)^2)/(1+(theta0-X)^2)^2)  
    }
    if ((iter<=maxiter)&&((is.na(err))<1)){
      mle[i]=theta1;
    }else{
      mle[i]='not converge'
    }  
  }
  return(mle)
}

### (1e)
mle_Fisher<-function(thetaguess,X,maxiter,tol,convergencefactor){
  mle=matrix(0,nrow=length(thetaguess),ncol=1)
  for (i in 1:length(thetaguess)){
    theta0=thetaguess[i];
    iter=1;err=1;
    n=length(X);
    F=-2*sum((theta0-X)/(1+(theta0-X)^2));
    FI=n/2;
    for (iter in 1:maxiter+1){
      theta1=theta0+F/FI;
      err=abs(theta1-theta0)/(abs(theta0)+convergencefactor)
      if((err<tol)||(is.na(err))) break
      theta0=theta1
      F=-2*sum((theta0-X)/(1+(theta0-X)^2))
    }
    if ((iter<=maxiter)&&((is.na(err))<1)){
      mle[i]=theta1;
    }else{
      mle[i]='not converge'
    }  
  }
  return(mle)
}


theta<-c(-11,-1,0,1.4,4.1,4.8,7,8,38)
mle1<-mle_newton(theta,x,1000,10^-6,10^-30)
result_1c<-cbind(theta,mle1)
mle2<-mle_Fisher(theta,x,1000,10^-6,10^-30)
mle2_2<-mle_newton(mle2,x,1000,10^-6,10^-30)
result_1d<-cbind(theta,mle2)
result_1d2<-cbind(theta,mle2_2)
write.table(result_1d,'~/Downloads/1d.txt',row.names=FALSE)
write.table(result_1d2,'~/Downloads/1d.txt',row.names=FALSE)


####### Problem 2
x<-c(0.52,1.96,2.22,2.28,2.28,2.46,2.50,2.53,2.54,2.99,3.47,3.53,3.70,3.88,3.91,4.04,4.06,4.82,4.85,5.46)
n=length(x)

### (2a)
theta=seq(-pi,pi,0.001)
ltheta<-matrix(0,nrow=length(theta),ncol=1)
for (i in 1:length(theta)){
  ltheta[i]=-n*log(2*pi)+sum(log(1-cos(x-theta[i])))
} 
plot(theta,ltheta,type='l',xlab=expression(theta),ylab='Log likelihood')

### (2b)
thetamom<-asin(mean(x)-pi)

### (2c)
mle_newton2<-function(thetaguess,X,maxiter,tol,convergencefactor){
  mle=matrix(0,nrow=length(thetaguess),ncol=1)
  for (i in 1:length(thetaguess)){
    theta0=thetaguess[i];
    iter=1;err=1;
    n=length(X)
    F=sum(-sin(x-theta0)/(1-cos(x-theta0)))
    DF=sum(1/(cos(x-theta0)-1))
    for (iter in 1:maxiter+1){
      theta1=theta0-F/DF;
      err=abs(theta1-theta0)/(abs(theta0)+convergencefactor)
      if((err<tol)||(is.na(err))) break
      theta0=theta1
      F=sum(-sin(x-theta0)/(1-cos(x-theta0)))
      DF=sum(1/(cos(x-theta0)-1))
    }
    if ((iter<=maxiter)&&((is.na(err))<1)){
      mle[i]=theta1;
    }else{
      mle[i]='not converge'
    }  
  }
  return(mle)
}
mle_newton2(thetamom,x,1000,10^-6,10^-30)

### (2d)
theta0<-c(-2.7,2.7)
mle_newton2(theta0,x,1000,10^-6,10^-30)

### (2e)
theta4<-seq(-pi,pi,by=2*pi/200)
mle4<-mle_newton2(theta4,x,1000,10^-6,10^-30)
plot(theta4,mle4,type='l',xlab=expression(theta),ylab='MLE')


# for (i in 1:length(mle4)-1){
#   d=mle4[i+1]-mle4[i]
#   if (abs(d)<10^-3){
#     mle4[i+1]=mle4[i]
#   }
# }
split(theta4,mle4)
result2d<-cbind(theta4,mle4)
write.table(result2d,'~/Downloads/2d.txt',row.names=FALSE)

######## Problem 3
x<-c(0.02,0.06,0.11,0.22,0.56,1.10,0.02,0.06,0.11,0.22,0.56,1.10)
# x<-rep(c(0.02,0.06,0.11,0.22,0.56,1.1),rep(2,6))
y<-c(47,97,123,152,191,200,76,107,139,159,201,207)
## (a)
ys<-1/y
u<-1/x
beta<-coef(lsfit(u,ys))
theta1<-1/beta[1];theta1
theta2<-beta[2]*theta1;theta2


## (b)
theta0<-c(theta1,theta2)
theta0<-as.matrix(theta0)

min_newton<-function(thetaguess,x,y,maxiter,tol,convergencefactor){
  min_theta=matrix(1,ncol=ncol(thetaguess),nrow=nrow(thetaguess))
  for (i in 1:ncol(thetaguess)){
    theta0=thetaguess[,i];
    iter=1;err=1;
    n=length(x);
    hessian<-matrix(1,2,2)
    DF<-matrix(1,2,1)
    hessian[1,1]=sum(2*x^2/((x+theta0[2])^2));
    hessian[2,2]=sum((-4*theta0[1]*(x*y)/(x+theta0[2])^3)+6*theta0[1]^2*x^2/((x+theta0[2])^4))
    hessian[1,2]=sum((2*x*y/(x+theta0[2])^2)-4*theta0[1]*x^2/((x+theta0[2])^3))
    hessian[2,1]=hessian[1,2]
    DF[1]=sum((-2*x*y/(x+theta0[2]))+2*theta0[1]*x^2/((x+theta0[2])^2))
    DF[2]=sum((2*theta0[1]*x*y/(x+theta0[2])^2)-2*theta0[1]^2*x^2/(x+theta0[2])^3)
#    g0=sum(y-theta0[1]*x/(x+theta0[2]))^2    
    for (iter in 1:maxiter+1){
      theta1=theta0-solve(hessian)%*%DF;
      err=sqrt(t(theta1-theta0)%*%(theta1-theta0))/(sqrt(t(theta0)%*%(theta0))+convergencefactor)
      if((err<tol)||(is.na(err))) break
      theta0=theta1;

      hessian[1,1]=sum(2*x^2/((x+theta0[2])^2))
      hessian[2,2]=sum((-4*theta0[1]*(x*y)/(x+theta0[2])^3)+6*theta0[1]^2*x^2/((x+theta0[2])^4))
      hessian[1,2]=sum((2*x*y/(x+theta0[2])^2)-4*theta0[1]*x^2/((x+theta0[2])^3))
      hessian[2,1]=hessian[1,2]
      DF[1]=sum((-2*x*y/(x+theta0[2]))+2*theta0[1]*x^2/((x+theta0[2])^2))
      DF[2]=sum((2*theta0[1]*x*y/(x+theta0[2])^2)-2*(theta0[1]^2)*(x^2/(x+theta0[2])^3)) 
    }
    if ((iter<=maxiter)&&((is.na(err))<1)){
      min_theta[,i]=theta1;
    }else{
      min_theta[,i]='not converge'
    }  
  }
  return(min_theta)
}
min_newton(theta0,x,y,1000,10^-6,10^-30)

####### 3(c)
l1 = function(t1,t2){sum(2*(y-(t1*x)/(x+t2))*(-x/(x+t2))) }
l2 = function(t1,t2){ sum(2*(y-(t1*x)/(x+t2))*((t1*x)/(x+t2)^2))}
g = function(start){
  t1 = start[1]
  t2 = start[2] 
  sum((y-(t1*x)/(x+t2))^2)
}
steepestdescent = function(start, tol, convergencefactor, maxiter){ 
  t1 = start[1]
  t2 = start[2]
  alpha=0.001
  new=matrix(,2,1)
  new = start - alpha*c(l1(t1,t2),l2(t1,t2)) 
  while (g(new)>g(start)) {
    alpha = alpha*.5
    new = start - alpha*c(l1(t1,t2),l2(t1,t2)) 
  }
  n =1
  #  while((sqrt((t(new-start))%*%(new-start)))/(sqrt((t(start)%*%(start))+delta)) > epsilon & n < iteration){
  while(sum(abs(new-start)/(abs(start)+convergencefactor)) > tol & n < maxiter){
    start = new
    t1 = start[1]
    t2 = start[2]
    alpha = 0.0001
    new = start - alpha*c(l1(t1,t2),l2(t1,t2)) 
    while (g(new)>g(start)) {
      alpha = alpha*.5
      new = start - alpha*c(l1(t1,t2),l2(t1,t2))
    }
    n=n+1
  }
  c(new,n)
}
steepestdescent(theta0,10^-6,10^-30,10000000)
# steepestdescent<-function(thetaguess,stepsize,x,y,maxiter,tol,convergencefactor){
#   min_theta=matrix(1,ncol=ncol(thetaguess),nrow=nrow(thetaguess))
#   for (i in 1:ncol(thetaguess)){
#     theta0=thetaguess[,i];
#     iter=1;err=1;
#     DF<-matrix(0,2,1);
#     DF[1]=sum((-2*x*y/(x+theta0[2]))+2*theta0[1]*x^2/((x+theta0[2])^2))
#     DF[2]=sum((2*theta0[1]*x*y/(x+theta0[2])^2)-2*(theta0[1]^2)*(x^2/(x+theta0[2])^3))
#     theta1=theta0-stepsize*DF;
#     for (iter in 1:maxiter+1){   
# #      err=sqrt(t(theta1-theta0)%*%(theta1-theta0))/(sqrt(t(theta0)%*%(theta0))+convergencefactor)
#       
#       g0=sum((y-theta0[1]*x/(x+theta0[2]))^2) 
#       g1=sum((y-theta1[1]*x/(x+theta1[2]))^2)
#       
# #      error=abs(g0-g1)/(abs(g0)+convergencefactor)
# 
#       if (g0<g1){
#         stepsize=stepsize/2
#         
#       }
#       DF[1]=sum((-2*x*y/(x+theta0[2]))+2*theta0[1]*x^2/((x+theta0[2])^2))
#       DF[2]=sum((2*theta0[1]*x*y/(x+theta0[2])^2)-2*(theta0[1]^2)*(x^2/(x+theta0[2])^3))
#       theta1=theta0-stepsize*DF;
#       err=sum(abs((theta1-theta0))/(abs((theta0))+convergencefactor))
# #err=sqrt(t(theta1-theta0)%*%(theta1-theta0))/(sqrt(t(theta0)%*%(theta0))+convergencefactor)
#       if((err<tol)||(is.na(err))) {
#         break} else{
#           theta0=theta1;
#         }          
# 
#     }
#     if ((iter<=maxiter)&&((is.na(err))<1)){
#       min_theta[,i]=theta1;
#     }else{
#       min_theta[,i]='not converge'
#     }  
#   }
#   return(list(min_theta,stepsize,iter))
# }
# steepestdescent(theta0,0.1,x,y,100000000,10^-6,10^-30)

######### 3(d)
gauss_newton<-function(thetaguess,x,y,maxiter,tol,convergencefactor){
  min_theta=matrix(1,ncol=ncol(thetaguess),nrow=nrow(thetaguess))
  for (i in 1:ncol(thetaguess)){
    theta0=thetaguess[,i];
    iter=1;err=1;
    n=length(x)
    A<-matrix(,n,2)
    A[,1]=x/(x+theta0[2])
    A[,2]=-theta0[1]*x/((x+theta0[2])^2)
   
    Z=y-theta0[1]*x/(x+theta0[2])
    Z<-matrix(Z,n,1)

    #    g0=sum(y-theta0[1]*x/(x+theta0[2]))^2    
    for (iter in 1:maxiter+1){
      theta1=theta0+(solve(t(A)%*%A))%*%(t(A))%*%Z
      err=sqrt(t(theta1-theta0)%*%(theta1-theta0))/(sqrt(t(theta0)%*%(theta0))+convergencefactor)
      if((err<tol)||(is.na(err))) break
      theta0=theta1;
      A<-matrix(,n,2)
      A[,1]=x/(x+theta0[2])
      A[,2]=-theta0[1]*x/((x+theta0[2])^2)
      
      Z=y-theta0[1]*x/(x+theta0[2])
      Z<-matrix(Z,n,1)

    }
    if ((iter<=maxiter)&&((is.na(err))<1)){
      min_theta[,i]=theta1;
    }else{
      min_theta[,i]='not converge'
    }  
  }
  return(min_theta)
}
gauss_newton(theta0,x,y,1000,10^-6,10^-30)
