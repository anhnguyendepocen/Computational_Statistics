##### STA243 HW5
### STA 243 HM 5 Tongbi Tu & Yuan Chen
## Problem 1
#(a)
#Implement the cross-validation and generalized cross-validation methods for choosing the smoothing parameter .
#(b)
# use the AICc criterion

library(MASS)

n = 200
xi = seq(0.5/n, 199.5/n, by = 1/n)
Xi = runif(n, 0, 1)
e = rnorm(n,0,1)


phi = function(x){
  exp(-x^2/2)/sqrt(2*pi)
}

f = function(x){
  1.5*phi((x - 0.35)/0.15) - phi((x - 0.8)/0.04)
}

sigmanj = function(j){
  0.02 + 0.04*(j-1)^2
}


F_inv = function(Xi, j = 1){
  qbeta(p = Xi, shape1 = (j+4)/5, shape2 = (11-j)/5)
}


fj = function(x, j = 1){
  sqrt(x*(1-x)) * sin( (2 * pi * (1 + 2^( (9 - 4 * j)/5) ))/ (x + 2^( (9 - 4 * j)/5) ))
}

vj = function(x, j = 1){
  (0.15 * (1 + 0.4*(2*j-7)*(x-0.5)))^2
}

yn = function(x, j = 1, epsilon = e){
  f(x) + sigmanj(j) * epsilon
}

yd = function(Xji, j = 1, epsilon = e){
  f(Xji) + 0.1 * epsilon
}

ys = function(x, j = 1, epsilon = e){
  fj(x, j) + 0.2 * epsilon
}

yv = function(x, j = 1, epsilon = e){
  f(x) + sqrt(vj(x, j)) * epsilon
}



t = seq(0, 1, length.out = 32)[2:31]

X = t(sapply(xi, function(x){
  xk = ifelse(x > t, x - t, 0)
  c(1, x, x^2, x^3, xk^3)
}))

D = diag(c(rep(0, 4), rep(1, 30)))

H_lambda = function(lambda = 0){
  X %*% solve(t(X) %*% X + lambda * D) %*% t(X)
}

H = X %*% solve(t(X) %*% X) %*% t(X)


cv = function(y, lambda){
  yhat = H_lambda(lambda) %*% y
  h = diag(H_lambda(lambda))
  sum(((y-yhat)/(1-h))^2)
}

gcv = function(y, lambda){
  yhat = H_lambda(lambda) %*% y
  h = sum(diag(H_lambda(lambda)))/200
  sum(((y-yhat)/(1-h))^2)
}

aicc = function(y, lambda){
  yhat = H_lambda(lambda) %*% y
  tr = sum(diag(H_lambda(lambda)))
  log(sum((y-yhat)^2)) + 2*(tr + 1)/(200 - tr - 2)
}

risk = function(y, lambda){
  yhat = H %*% y
  rss = sum((y - yhat)^2)
  sum((y - H_lambda(lambda) %*% y)^2) + rss/(200 - sum(diag(H))) * (2* sum(diag(H_lambda(lambda))) - 200)
}



mse = function(f, lambda, y){
  sum((f(xi) - H_lambda(lambda) %*% y)^2)
}



algorithm = function(j = 1, x = xi, f_y = yn, f = f){
  e = rnorm(200,0,1)
  y = f_y(x, j, e)
  lam = sapply(c(cv, gcv, aicc, risk), function(method) optimize(function(l) method(y, l), c(0, 10))$minimum)
  mse_min = optimize(function(l) mse(f, l, y), c(0, 10))$objective
  sapply(lam, function(l) mse(f, l, y))/mse_min
}


r1 = lapply(1:6, function(j) sapply(1:200, function(i) algorithm(j, xi, yn, f)))
r2 = lapply(1:6, function(j){
  xi = sort(F_inv(Xi, j))
  X = t(sapply(xi, function(x){
    xk = ifelse(x > t, x - t, 0)
    c(1, x, x^2, x^3, xk^3)
  }))
  H = X %*% ginv(t(X) %*% X) %*% t(X)
  sapply(1:200, function(i) algorithm(j, xi, yd, f))
})
r3 = lapply(1:6, function(j) sapply(1:200, function(i) algorithm(j, xi, ys, function(x) fj(x, j))))
r4 = lapply(1:6, function(j) sapply(1:200, function(i) algorithm(j, xi, yv, f)))





par(mar=c(1,1,1,1))
par(mfrow = c(3, 4))
sapply(1:6, function(j){
  plot(xi, yn(xi, j), pch = 1, main = paste0("j = ", j),cex = 0.1,ylab = '',xlab = '')
  lines(xi, f(xi))
  boxplot(t(log(r1[[j]])), xaxt = 'n',ylab = 'r', xlab = 'criteria')
  axis(1, at = 1:4 , labels = c('CV','GCV','AICc','Cp'))
})
sapply(1:6, function(j){
  plot(F_inv(Xi, j), yd(F_inv(Xi, j), j), pch = 1, main = paste0("j = ", j),cex = 0.1,ylab = '',xlab = '')
  lines(sort(F_inv(Xi, j)), f(sort(F_inv(Xi, j))))
  boxplot(t(log(r2[[j]])), xaxt = 'n',ylab = 'r', xlab = 'criteria')
  axis(1, at = 1:4 , labels = c('CV','GCV','AICc','Cp'))
})
sapply(1:6, function(j){
  plot(xi, ys(xi, j), pch = 1, main = paste0("j = ", j),cex = 0.1,ylab = '',xlab = '')
  lines(xi, fj(xi, j))
  boxplot(t(log(r3[[j]])), xaxt = 'n',ylab = 'r', xlab = 'criteria')
  axis(1, at = 1:4 , labels = c('CV','GCV','AICc','Cp'))
})
sapply(1:6, function(j){
  plot(xi, yv(xi, j), pch = 1, main = paste0("j = ", j),cex = 0.1,ylab = '',xlab = '')
  lines(xi, f(xi))
  boxplot(t(log(r4[[j]])), xaxt = 'n',ylab = 'r', xlab = 'criteria')
  axis(1, at = 1:4 , labels = c('CV','GCV','AICc','Cp'))
})





###generate data
setwd("~/Desktop/STA_2017_spring/STA243")
phi=function(x) {1/sqrt(2*pi)*exp(-x^2/2)}
f_fun=function(x,i,j)
{
  
  if(i==3)
  {
    sqrt(x*(1-x))*sin(2*pi*(1+2^((9-4*j)/5))/(x+2^((9-4*j)/5)))
  }
  else
  {
    1.5*phi((x-0.35)/0.15)-phi((x-0.8)/0.04) 
  }
}


generatenoisedata = function(j,err) {
  sigma = 0.02 + 0.04*(j-1)^2
  x = (c(1:200)-0.5)/200
  y = f_fun(x,1,j) + sigma*err
  return(data.frame(x,y))
}

generatedensitydata = function(j,err) {
  sigma = 0.1
  X = sort(runif(200,0,1))
  x= qbeta(X,(j+4)/5, (11-j)/5)
  y = f_fun(x,2,j) +sigma*err
  return(data.frame(x,y))
}


generatespatialdata=function(j,err)
{
  x=((1:200)-0.5)/200
  sigma=0.2
  y=f_fun(x,3,j)+sigma*err
  data.frame(x,y)
}

generatevariancedata=function(j,err)
{
  x=((1:200)-0.5)/200
  y=f_fun(x,4,j)+abs(0.15*(1+0.4*(2*j-7)*(x-0.5)))*err
  data.frame(x,y)
}
generatealldata=function()
{
  err=rnorm(200,0,1)
  noisedata = lapply(1:6,function(j) generatenoisedata(j,err))
  densitydata=lapply(1:6,function(j) generatedensitydata(j,err))
  spatialdata=lapply(1:6,function(j) generatespatialdata(j,err))
  variancedata=lapply(1:6,function(j) generatevariancedata(j,err))
  output=list()
  output$noisedata=noisedata
  output$densitydata=densitydata
  output$spatialdata=spatialdata
  output$variancedata=variancedata
  output
}

t=(1:30)/31
D=diag(c(rep(0,4),rep(1,30)))

H=function(data,lambda)
{
  sig=function(x)
  {
    x[x<0]=0
    x
  }
  x=data$x
  x_mat=data.frame(int=1,x,x^2,x^3)
  for(i in 1:length(t))
  {
    x_mat=cbind(x_mat,sig(x-t[i])^3)
  }
  x_mat=as.matrix(x_mat)
  x_mat%*%solve(t(x_mat)%*%x_mat+lambda*D)%*%t(x_mat)
}


CV=function(lambda,data)
{
  H_mat=H(data,lambda)
  f_hat=H_mat%*%data$y
  1/200*(sum(((data$y-as.numeric(f_hat))/(1-diag(H_mat)))^2))
}




GCV = function(lambda, data) {
  H_mat=H(data,lambda)
  f_hat=as.numeric(H_mat%*%data$y)
  1/200*sum((data$y-f_hat)^2)/((1-1/200*sum(diag(H_mat)))^2)
}




AIC=function(lambda,data)
{
  H_mat=H(data,lambda)
  f_hat=H_mat%*%data$y
  tr=sum(diag(H_mat))
  log(sum((data$y-f_hat)^2))+2*(tr+1)/(198-tr)
}

RISK=function(lambda,data)
{
  H_mat=H(data,lambda)
  f_hat=H_mat%*%data$y
  MSE=sum((data$y-f_hat)^2)/(200-34)
  sum((data$y-f_hat)^2)+MSE*(2*sum(diag(H_mat))-200)
}

deviate=function(lambda,data,ftrue)
{
  H_mat=H(data,lambda)
  f_hat=H_mat%*%data$y
  sum((ftrue-f_hat)^2)
}

output=data.frame()


for(k in 1:200)
{
  dataset=generatealldata()
  for(i in 1:4)
  {
    for(j in 1:6)
    {
      data=dataset[[i]][[j]]
      ftrue=f_fun(data$x,i,j)
      mindevi=optimize(deviate,c(0,1),data,ftrue,maximum=FALSE)$objective
      
      
      cv_opt=optimize(CV,c(0,1),data,maximum=FALSE)$minimum
      cv_f_hat=H(data,cv_opt)%*%data$y
      cv_r=sum((ftrue-cv_f_hat)^2)/mindevi
      print(data.frame(k,i,j,method="cv",r=cv_r))
      output=rbind(output,data.frame(k,i,j,method="cv",r=cv_r))
      
      
      gcv_opt=optimize(GCV,c(0,1),data,maximum=FALSE)$minimum
      gcv_f_hat=H(data,gcv_opt)%*%data$y
      gcv_r=sum((ftrue-gcv_f_hat)^2)/mindevi
      print(data.frame(k,i,j,method="gcv",r=gcv_r))
      output=rbind(output,data.frame(k,i,j,method="gcv",r=gcv_r))
      
      aic_opt=optimize(AIC,c(0,1),data,maximum=FALSE)$minimum
      aic_f_hat=H(data,aic_opt)%*%data$y
      aic_r=sum((ftrue-aic_f_hat)^2)/mindevi
      print(data.frame(k,i,j,method="aic",r=aic_r))
      output=rbind(output,data.frame(k,i,j,method="aic",r=aic_r))
      
      risk_opt=optimize(RISK,c(0,1),data,maximum=FALSE)$minimum
      risk_f_hat=H(data,risk_opt)%*%data$y
      risk_r=sum((ftrue-risk_f_hat)^2)/mindevi
      print(data.frame(k,i,j,method="risk",r=risk_r))
      output=rbind(output,data.frame(k,i,j,method="risk",r=risk_r))
    }
  }
}

write.csv(output,'output.csv')
output=read.csv('output.csv')

names=c("noise data","density data","spatial data","variance data")

dataset=generatealldata()
par(mfrow=c(3,4))
for(i in 1:4)
{
  for(j in 1:6)
  {
    print(plot(dataset[[i]][[j]],main=paste('j=',j)))
    print(lines(x=dataset[[i]][[j]]$x,y=f_fun(dataset[[i]][[j]]$x,i,j),col='red'))
    print(boxplot(log(r)~method,data=output[output$r>=1&output$i==i&output$j==j,]))
  }
  print(title(names[i],outer=TRUE,line=-1))
}
