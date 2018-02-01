### STA243 HW2
############# Problem 2 Genetic Algorithm
truefunction=function(x)
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



MDL=function(chromo){
  chromo[1]=1 # the startin point needs to be obtained
  breakpoint=c(which(chromo==1),length(chromo)+1)
  Bhat=length(breakpoint)-1
  fhat=vector()
  nhat=vector()
  for (i in 1:Bhat){
    
    fhat[breakpoint[i]:(breakpoint[i+1]-1)]=mean(y[breakpoint[i]:(breakpoint[i+1]-1)]) 
    nhat[i]=breakpoint[i+1]-breakpoint[i]
  }
  n=length(y)
  MDL=Bhat*log(n)+1/2*sum(log(nhat))+n/2*log(1/n*sum((y-fhat)^2))
  MDL
}

AIC=function(chromo){
  chromo[1]=1 # the startin point needs to be obtained
  breakpoint=c(which(chromo==1),length(chromo)+1)
  Bhat=length(breakpoint)-1
  fhat=vector()
  nhat=vector()
  for (i in 1:Bhat){
    fhat[breakpoint[i]:(max((breakpoint[i+1]-1),1))]=mean(y[breakpoint[i]:(breakpoint[i+1]-1)]) 
    nhat[i]=max(breakpoint[i+1]-breakpoint[i],1)
  }
  n=length(y)
  AIC=2*Bhat*log(n)+n*log(1/n*sum((y-fhat)^2))
  AIC
}

GA=function(x,y,samplesize,pcross,pmutation,objfun,largestsame){
  ### initial population of chromosome
  parent=list()
  for (i in 1:samplesize){
    parent[[i]]=sample(c(0,1),length(y),replace=T)
  }
  bestsame=1
  gen=0
  bestchromo=vector()
  recordbestchromo=list() # record best chromo at each generation
  bestchromo=sample(c(0,0),length(y),replace=T)
  bestobj=vector() # store the best function result of each generation
  while (bestsame<largestsame){
    gen=gen+1
    
    
    obj=sapply(parent,objfun) # evaluate the object function for all the chromosomes
    
    bestobj[gen]=min(obj) # best result at each generation
    
    min_ind=which.min(obj)[1] # index of the best chromosome
    
    if (sum(parent[[min_ind]]!=bestchromo)==0){
      bestsame=bestsame+1
    }else{
      bestsame=1
    }
    bestchromo=parent[[min_ind]]  #update best_chromo
    
    
    recordbestchromo[[gen]]=bestchromo
    
    
    rankofchromo=rank(-obj)
    prob_chromo=rankofchromo/sum(rankofchromo)
    
    ##### generate offspring/child
    decision=sample(c(0,1),1,prob=c(pcross,1-pcross))
    offspring=list()
    child=vector()
    
    for (k in 1:samplesize){
      if (decision==0){
        ### perform cross
        # choose parent chromosomes for the child
        parent_chosen=sample(1:samplesize,2,replace=F,prob=prob_chromo)
        parent1=parent[[parent_chosen[1]]]
        parent2=parent[[parent_chosen[2]]]
        for (i in 1:length(y)){
          child[i]=sample(c(parent1[i],parent2[i]),1,prob=c(0.5,0.5))
        }
        child[1]=0
        offspring[[k]]=child
      }else{
        ### perform mutation
        # choose parent chromosome for the child
        parent_chosen=sample(1:samplesize,1,prob=prob_chromo)
        parent1=parent[[parent_chosen]]
        for (i in 1:length(y)){
          child[i]=sample(c(parent1[i],1-parent1[i]),1,prob=c(1-pmutation,pmutation))
        }
        child[1]=0  # take the first gene as 0
        
        offspring[[k]]=child
      }    
    }
    
    xxx=sample(1:samplesize,1)
    
    offspring[[xxx]]=parent[[min_ind]] # keep the best in the previous generation
    ##offspring[[1]]=parent[[min_ind]] # keep the best in the previous generation
    ##offspring[[samplesize]]=parent[[min_ind]] # keep the best in the previous generation
    
    parent=offspring # update parent choromosomes for the next generation
    
    
  }
  output=list()
  output$best_objective=bestobj
  output$best_chromosome=recordbestchromo
  output
}

############################################
### the following parameters can be changedd to test the sensitivity
samplesize=300
pcross=0.9
pmutation=0.05
objfun=MDL
largestsame=5
############################################
# ptm<-proc.time()
# #result=GA(x,y,samplesize,pcross,pmutation,objfun,largestsame)
# result1=GA(x,y,300,0.9,0.05,MDL,20) # this one no elist procedure
# proc.time()-ptm

ptm<-proc.time()
#result=GA(x,y,samplesize,pcross,pmutation,objfun,largestsame)
result4=GA(x,y,300,0.9,0.05,MDL,20)
proc.time()-ptm




fhat_result=function(chromo){
  chromo[1]=1 # the startin point needs to be obtained
  breakpoint=c(which(chromo==1),length(chromo)+1)
  Bhat=length(breakpoint)-1
  fhat=vector()
  nhat=vector()
  for (i in 1:Bhat){
    fhat[breakpoint[i]:(max((breakpoint[i+1]-1),1))]=mean(y[breakpoint[i]:(breakpoint[i+1]-1)]) 
    nhat[i]=max(breakpoint[i+1]-breakpoint[i],1)
  }
  #  result=list("fhat"=fhat,"Bhat"=Bhat,"nhat"=nhat)
  result=list()
  result$fhat=fhat
  result$Bhat=Bhat
  result$nhat=nhat
  result
}

MDL20=result4$best_chromosome[[length(result4$best_objective)]]
AIC20=result3$best_chromosome[[length(result3$best_objective)]]
fhat_result1=fhat_result(MDL20)
MDL20yhat=fhat_result1$fhat
MDL20bhat=fhat_result1$Bhat
MDL20nhat=fhat_result1$nhat
fhat_result2=fhat_result(AIC20)
AIC20yhat=fhat_result2$fhat
AIC20bhat=fhat_result2$Bhat
AIC20nhat=fhat_result2$nhat
plot(1:length(result4$best_objective),result4$best_objective,type='l',xlab='generation',ylab='min(MDL)')
#par(new=T)
plot(1:length(result3$best_objective),result3$best_objective,type='l',col='red',xlab='generation',ylab='min(AIC)')

#### fit for AIC
plot(x,y)
lines(x,AIC20yhat,col='red',lty=1,lwd = 2)

### fit for MDL
plot(x,y)
lines(x,f)
par(new=T)
lines(x,MDL20yhat,col='blue',lty=1,lwd=2)









############# Problem 1 Simulated Annealing Algorithm

distance=matrix(c(0,1,2,4,9,8,3,2,1,5,7,1,2,9,3,
                  1,0,5,3,7,2,5,1,3,4,6,6,6,1,9,
                  2,5,0,6,1,4,7,7,1,6,5,9,1,3,4,
                  4,3,6,0,5,2,1,6,5,4,2,1,2,1,3,
                  9,7,1,5,0,9,1,1,2,1,3,6,8,2,5,
                  8,2,4,2,9,0,3,5,4,7,8,3,1,2,5,
                  3,5,7,1,1,3,0,2,6,1,7,9,5,1,4,
                  2,1,7,6,1,5,2,0,9,4,2,1,1,7,8,
                  1,3,1,5,2,4,6,9,0,3,3,5,1,6,4,
                  5,4,6,4,1,7,1,4,3,0,9,1,8,5,2,
                  7,6,5,2,3,8,7,2,3,9,0,2,1,8,1,
                  1,6,9,1,6,3,9,1,5,1,2,0,5,4,3,
                  2,6,1,2,8,1,5,1,1,8,1,5,0,9,6,
                  9,1,3,1,2,2,1,7,6,5,8,4,9,0,7,
                  3,9,4,3,5,5,4,8,4,2,1,3,6,7,0),byrow=T,nrow=15)
rownames(distance)= c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O')
colnames(distance)= c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O')

### travelpath calculates the length of the travel path.
travelpath=function(theta){
  travel=0
  for (i in 1:(length(theta)-1)){
    travel=travel+distance[theta[i],theta[i+1]]
  }
  travel
}
alpha=function(tao,probab){
  alpha=probab*tao
  alpha
}
betam=function(m){
  betam=m
}


simulated_annealing=function(tao,p,maxj){
#   tao=400 # initial temperature
#   p=0.999
  m=100
  theta=sample(1:ncol(distance),replace=F) # initial path length
  
  j=0
  finalpath=vector()
  
  while (j<=maxj){
    j=j+1
    finalpath[j]=travelpath(theta)
    
    tao=alpha(tao,p)
    mj=betam(m)
    
    for (i in 1:mj){
      
      theta_candi=sample(1:ncol(distance),replace=F) # candidate path
      
      delta=travelpath(theta_candi)-travelpath(theta)
      if (delta<=0){
        theta=theta_candi
      }else{
        choice=sample(c(0,1),1,prob=c(exp(-delta/tao),1-exp(-delta/tao)))
        if (choice==0){
          theta=theta_candi
        }else{
          theta=theta
        }
      }
    }
        
  }
  result_SA=list() # output results
  result_SA$shortest_path=colnames(distance)[theta]
  result_SA$shortest_length=travelpath(theta)
  result_SA$every_final_path=finalpath     
  result_SA
}

ptm<-proc.time()
tao=400 # initial temperature
p=0.999

result_SA1=simulated_annealing(tao,p,30000)
proc.time()-ptm

