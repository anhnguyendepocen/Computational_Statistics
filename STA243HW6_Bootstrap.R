## Problme 1
###(c)
X=runif(50,0,3)

p_bootstrap=lapply(1:5000,function(i) runif(50,0,max(X)))
ptheta_hat_bootstrap=sapply(p_bootstrap,max)
var_p_bootstrap=var(ptheta_hat_bootstrap)
var_true=50*9/((51)^2*52)
var_p_bootstrap
var_true

###(d)
np_bootstrap=lapply(1:5000, function(i) sample(X,50,replace=TRUE))
nptheta_hat_bootstrap=sapply(np_bootstrap,max)
var_np_bootstrap=var(nptheta_hat_bootstrap)
var_np_bootstrap


###(e)

hist(ptheta_hat_bootstrap,main='histgram of parametric bootstrap')
plot(density(ptheta_hat_bootstrap),main='density of parametric bootstrap')
hist(nptheta_hat_bootstrap,main='nonparametric bootstrap')
plot(density(nptheta_hat_bootstrap),main='density of nonparametric bootstrap')
```

###(f)

true_density=function(x)
{
  50*(x^49)/(3^50)
}
U=runif(5000,0,3)
hist(3*U^(1/50),main='true histgram')
plot((0:3000)/1000,sapply((0:3000)/1000,true_density),type='l',main='true density')


## Problem 2
###(a)

truefunction<-function(x)
{ 
  t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81) 
  h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2) 
  temp <- 0 
  for(i in 1:11) 
  { 
    temp <- temp + h[i]/2 * (1 + sign(x - t[i])) 
  } 
  return(temp) 
}

n<-512 
x<-(0:(n-1))/n 
f<-truefunction(x) 
set.seed(0401)
y<-f+rnorm(f)/3 
plot(x,y)
lines(x,f)
dataset1=data.frame(x,y,f)

truefunction2=function(x)
{
  4*x-2+2*exp(-16*(4*x-2)^2)
}

f2=truefunction2(x)
mod_f=optimize(truefunction2,c(0,1),maximum=TRUE)$objective
y2=f2+rnorm(length(f2),sd=mod_f/5)
dataset2=data.frame(x,y=y2,f=f2)




fhat=function(l,y)
{
  l[1]=1
  cut=c(which(l==1),length(l)+1)
  yhat=vector()
  for(i in 1:(length(cut)-1))
  {
    yhat[cut[i]:(cut[i+1]-1)]=mean(y[cut[i]:(cut[i+1]-1)])
  }
  yhat
}

MDL=function(l,y)
{
  l[1]=1
  cut=c(which(l==1),length(l)+1)
  Bhat=length(cut)-1
  nhat=vector()
  yhat=vector()
  for(i in 1:Bhat)
  {
    nhat[i]=cut[i+1]-cut[i]
    yhat[cut[i]:(cut[i+1]-1)]=mean(y[cut[i]:(cut[i+1]-1)])
  }
  n=length(y)
  MDL=Bhat*log(n)+1/2*sum(log(nhat))+n/2*log(1/n*sum((y-yhat)^2))
  MDL
}

AIC=function(l)
{
  l[1]=1
  cut=c(which(l==1),length(l)+1)
  Bhat=length(cut)-1
  yhat=vector()
  for(i in 1:Bhat)
  {
    yhat[cut[i]:(cut[i+1]-1)]=mean(y[cut[i]:(cut[i+1]-1)])
  }
  AIC=n*log(1/n*sum((y-yhat)^2))+log(n)*2*Bhat
  AIC
}

geneticalgoritm=function(x,y,Pcross,Pc,S,criterion,largestsame)
{
  originalset=list()
  for(i in 1:S)
  {
    originalset[[i]]=sample(c(0,1),length(y),replace=TRUE)
  }
  
  
  best_value=vector()
  best_chromosome=rep(0,length(y))
  #record_best_chromosome=list()
  
  parentset=originalset
  
  bestsamecount=1
  k=0
  while(bestsamecount<=largestsame)
  {
    k=k+1
    criterion_value=sapply(parentset,function(u) criterion(u,y))
    best_value[k]=min(criterion_value)
    if(sum(parentset[[which.min(criterion_value)[1]]]!=best_chromosome)==0)
    {
      bestsamecount=bestsamecount+1
    }
    else
    {
      bestsamecount=1
    }
    best_chromosome=parentset[[which.min(criterion_value)[1]]]
    #record_best_chromosome[[k]]=parentset[[which.min(criterion_value)[1]]]
    r=rank(-criterion_value)
    p=r/sum(r)
    childset=list()
    child=vector()
    for(j in 1:S)
    {
      coin=sample(c(0,1),1,prob=c(Pcross,1-Pcross))
      if(coin==0)#crossover
      {
        parentindicate=sample(1:S,2,replace=FALSE,prob=p)
        parent1=parentset[[parentindicate[1]]]
        parent2=parentset[[parentindicate[2]]]
        for(i in 1:length(parent1))
        {
          child[i]=sample(c(parent1[i],parent2[i]),1,prob=c(0.5,0.5))
        }
        childset[[j]]=child
      }
      else#mutate
      {
        parentindicate=sample(1:S,1,prob=p)
        parent1=parentset[[parentindicate]]
        for(i in 1:length(parent1))
        {
          child[i]=sample(c(1-parent1[i],parent1[i]),1,prob=c(Pc,1-Pc))
        }
        childset[[j]]=child
      }
    }
    parentset=childset
  }
  
  yhat=fhat(best_chromosome,y)
  record=data.frame(generation=1:k,min_value=best_value)
  output=list()
  output$best_criterion_value_each_generation=record
  output$yhat=yhat
  output$best_chromosome=best_chromosome
  best_chromosome[1]=1
  cut=c(which(best_chromosome==1),length(best_chromosome)+1)
  output$piecenum=length(cut)-1
  output
}

library(ggplot2)
makeplot1=function(output,x,y,f)
{
  origindata=data.frame(x=x,y=y,type='point')
  originline=data.frame(x=x,y=f,type='trueline')
  result=data.frame(x=x,y=output$yhat,type='estimateline')
  infosheet=rbind(origindata,originline,result)
  ggplot()+geom_line(data=infosheet[infosheet$type!='point',],aes(x,y,color=type))+geom_point(data=infosheet[infosheet$type=='point',],aes(x,y))+theme(legend.text = element_text(size = 13))
}
makeplot2=function(output)
{
  ggplot(output$best_criterion_value_each_generation,aes(x=generation,y=min_value))+geom_line()+theme(legend.text = element_text(size = 13))
}

output1=geneticalgoritm(dataset1$x,dataset1$y,0.9,0.1,300,MDL,5)
output2=geneticalgoritm(dataset2$x,dataset2$y,0.9,0.1,300,MDL,5)

makeplot1(output1,dataset1$x,dataset1$y,dataset1$f)
makeplot2(output1)
makeplot1(output2,dataset2$x,dataset2$y,dataset2$f)
makeplot2(output2)


###(b)


#dataset1_boot
dataset1$fhat=output1$yhat
dataset1$residual=dataset1$y-dataset1$fhat

boot_residual1=lapply(1:100,function(i)
  sample(dataset1$residual,length(dataset1$residual),replace=TRUE))

dataset1_boot=data.frame()
for(i in 1:100)
{
  dataset1_boot=rbind(dataset1_boot,data.frame(x=x,y=dataset1$fhat+boot_residual1[[i]],i=i))
}

write.csv(dataset1_boot,'dataset1_boot')

dataset=read.csv('dataset1_boot')

output=data.frame()
for(i in 1:10)
{
  output=rbind(output,data.frame(x=dataset$x[dataset$i==i],yhat=geneticalgoritm(dataset$x[dataset$i==i],dataset$y[dataset$i==i],0.9,0.1,300,MDL,5)$yhat,i))
}

write.csv(output,'output1')

setwd('D:/study/course/2016 spring/243/hw6/dataset1_boot_output/')
boot_result1=data.frame()
for(i in 1:100)
{
  boot_result1=rbind(boot_result1,read.csv(paste('output',i,sep='')))
}

boot_result1=boot_result1[,-1]
boot_result1$i=as.factor(boot_result1$i)
ggplot(boot_result1,aes(x,yhat,group=i))+geom_line()

yup=vector()
ydown=vector()
for(i in 1:512)
{
  ylist=sort(boot_result1$yhat[seq(i,i+99*512,512)])
  yup[i]=(ylist[2]+ylist[3])/2
  ydown[i]=(ylist[97]+ylist[98])/2
}
boot_band1=data.frame(x=boot_result1$x[boot_result1$i==1],y=dataset1$y,true=dataset1$f,yup,ydown)
ggplot(data=boot_band1)+geom_point(aes(x,y))+geom_line(aes(x,y=true),col='blue')+geom_line(aes(x,y=yup),col='red')+geom_line(aes(x,y=ydown),col='red')+geom_line(data=dataset1,aes(x,y=fhat),color='yellow')+theme(legend.text = element_text(size = 20))




#dataset1 pair
boot_pair=lapply(1:100,function(i) {pairindicate=sort(sample(1:512,replace=TRUE))
                                    data.frame(x=dataset1$x[pairindicate],y=dataset1$y[pairindicate] )})

dataset1_pair=data.frame()
for(i in 1:100)
{
  dataset1_pair=rbind(dataset1_pair,data.frame(boot_pair[[i]],i))
}

write.csv(dataset1_pair,'dataset1_pair')




setwd('D:/study/course/2016 spring/243/hw6/dataset1_pair_output/')

pair_result1=data.frame()
for(i in 1:10)
{
  pair_result1=rbind(pair_result1,read.csv(paste('output',i,sep='')))
}
pair_result1=pair_result1[,-1]
pair_result1$i=as.factor(pair_result1$i)
library(ggplot2)
ggplot(pair_result1,aes(x,y=yhat,group=as.factor(i)))+geom_line()

pair_band1=
  pair_result1 %>% 
  group_by(x) %>% 
  arrange(yhat) %>% 
  mutate(up=yhat[0.025*length(yhat)],down=yhat[0.975*length(yhat)]) %>% 
  filter(i==2)
pair_band1$y=dataset1$y
pair_band1$true=dataset1$f
ggplot(data=pair_band1)+geom_point(data=dataset1,aes(x,y))+geom_line(data=dataset1,aes(x,y=f),col='blue')+geom_line(aes(x,y=up),col='red')+geom_line(aes(x,y=down),col='red')+geom_line(data=dataset1,aes(x,y=fhat),col='yellow')


ggplot(pair_result1,aes(x,yhat,group=i))+geom_line()





#dataset2_boot
setwd('D:/study/course/2016 spring/243/hw6/')
dataset2$fhat=output2$yhat
dataset2$residual=dataset2$y-dataset2$fhat

boot_residual2=lapply(1:100,function(i)
  sample(dataset2$residual,length(dataset2$residual),replace=TRUE))

dataset2_boot=data.frame()
for(i in 1:100)
{
  dataset2_boot=rbind(dataset2_boot,data.frame(x=x,y=dataset2$fhat+boot_residual2[[i]],i=i))
}

write.csv(dataset2_boot,'dataset2_boot')



setwd('D:/study/course/2016 spring/243/hw6/dataset2_boot_output/')
boot_result2=data.frame()
for(i in 1:10)
{
  boot_result2=rbind(boot_result2,read.csv(paste('output',i,sep='')))
}
boot_result2=boot_result2[,-1]
boot_result2$i=as.factor(boot_result2$i)


yup=vector()
ydown=vector()
for(i in 1:512)
{
  ylist=sort(boot_result2$yhat[seq(i,i+99*512,512)])
  yup[i]=(ylist[2]+ylist[3])/2
  ydown[i]=(ylist[97]+ylist[98])/2
}
boot_band2=data.frame(x=boot_result2$x[boot_result2$i==1],y=dataset2$y,true=dataset2$f,yup,ydown)
ggplot(data=boot_band2)+geom_point(aes(x,y))+geom_line(aes(x,y=true),col='blue')+geom_line(aes(x,y=yup),col='red')+geom_line(aes(x,y=ydown),col='red')+geom_line(data=dataset2,aes(x,fhat),col='yellow')



ggplot(boot_result2,aes(x,yhat,group=i))+geom_line()


#dataset2 pair
setwd('D:/study/course/2016 spring/243/hw6/')
boot_pair2=lapply(1:100,function(i) {pairindicate=sort(sample(1:512,replace=TRUE))
                                     data.frame(x=dataset2$x[pairindicate],y=dataset2$y[pairindicate] )})

dataset2_pair=data.frame()
for(i in 1:100)
{
  dataset2_pair=rbind(dataset2_pair,data.frame(boot_pair2[[i]],i))
}


write.csv(dataset2_pair,'dataset2_pair')

setwd('D:/study/course/2016 spring/243/hw6/dataset2_pair_output/')

pair_result2=data.frame()
for(i in 1:10)
{
  pair_result2=rbind(pair_result2,read.csv(paste('output',i,sep='')))
}
pair_result2=pair_result2[,-1]
pair_result2$i=as.factor(pair_result2$i)
ggplot(pair_result2,aes(x,yhat,group=i))+geom_line()

pair_band2=
  pair_result2 %>% 
  group_by(x) %>% 
  arrange(yhat) %>% 
  mutate(up=yhat[0.025*length(yhat)],down=yhat[0.975*length(yhat)]) %>% 
  filter(i==2)
pair_band2$y=dataset2$y
pair_band2$true=dataset2$f
ggplot()+geom_point(data=dataset2,aes(x,y))+geom_line(data=dataset2,aes(x,y=f),col='blue')+geom_line(data=pair_band2,aes(x,y=yup),col='red')+geom_line(data=pair_band2,aes(x,y=ydown),col='red')+geom_line(data=dataset2,aes(x,y=fhat),col='yellow')

### Problem 3
#1
#(a)
#True value
(qt(0.75,3) - qt(0.25,3))/1.34 #1.14163
X = rt(25,3)
jack_t = function(X){
  n = length(X)
  theta = rep(NA,times = n) for (i in 1:n){
    sample = X[-i]
    theta[i] = (quantile(sample,0.75)-quantile(sample,0.25))/1.34 }
  theta_mean = mean(theta)
  se = sqrt((n-1)/n*sum((theta-theta_mean)^2))
  CI = c(theta_mean-1.645*se,theta_mean+1.645*se) return(CI)
}
#(b)
boot_t = function(X,B){
  n = length(X)
  theta = rep(NA,times = B) for (i in 1:B){
    s = sample(X,size = n,replace = TRUE)
    theta[i] = (quantile(s,0.75)-quantile(s,0.25))/1.34 }
  theta_mean = mean(theta)
  se = sd(theta)
  CI = c(theta_mean-1.645*se,theta_mean+1.645*se) return(CI)
}
#(c)
boot_q = function(X,B){
  n = length(X)
  theta = rep(NA,times = B) for (i in 1:B){
    s = sample(X,size = n,replace = TRUE)
    theta[i] = (quantile(s,0.75)-quantile(s,0.25))/1.34 }
  CI = c(quantile(theta,0.05),quantile(theta,0.95))
  return(CI) }
#to compare
boot_com = function(B,numint){
  CI1 = matrix(NA,nrow = numint,ncol = 2) CI2 = matrix(NA,nrow = numint,ncol = 2) CI3 = matrix(NA,nrow = numint,ncol = 2)
  for (i in 1:numint){ X = rt(25,3) CI1[i,]=jack_t(X) CI2[i,]=boot_t(X,B) CI3[i,]=boot_q(X,B)
  }
  CI =list(CI1,CI2,CI3)
  return(CI) }
l=boot_com(2000,50)