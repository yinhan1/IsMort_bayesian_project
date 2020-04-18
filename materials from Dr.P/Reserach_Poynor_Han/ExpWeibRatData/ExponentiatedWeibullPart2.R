##########################Inference for Functionals #################

#######################################################  Exponentiated Weibull Hazard Rate Function

exp.hr=function(x,alpha,theta,sigma){
  
  HR = exp.pdf(x,alpha,theta,sigma)/(1-exp.dis(x,alpha,theta,sigma))
  return(HR)}

############ Thinning the Posterior Draws #############

# Final.res=matrix(0,nrow=2001,ncol=3)
# Final.res[1,1]=Final_res[1,1] ### change mh6.ad
# for(i in 1:2000){
#   j= (20*(i)) #### change 15 if needed
#   Final.res[(i+1),]=Final_res[j,]
# }

Final.res=matrix(0,nrow=201,ncol=4)
Final.res[1,1]=Final_res[1,1] ### change mh6.ad
for(i in 1:200){
  j= (10*(i)) #### change 15 if needed
  Final.res[(i+1),]=Final_res[j,]
}
par(mfrow = c(4, 2))
plot(Final.res[-1,1])
acf(Final.res[-1,1])

plot(Final.res[-1,2])
acf(Final.res[-1,2])

plot(Final.res[-1,3])
acf(Final.res[-1,3])

plot(Final.res[-1,4])
acf(Final.res[-1,4])


##########################  Survival Function  ###
# grid=seq(0.01,max(y),.01)
# q=matrix(0,nrow=(length(grid)),ncol=3)
# sur=numeric(2000)
# 
# for(j in 1:length(grid)){
#   for(i in 1:2000){
#     sur[i]=exp.dis(grid[j],Final.res[(1+i),1],Final.res[(1+i),2],Final.res[(1+i),3])
#     #sur[i]=exp.dis(grid[j],Final.ad[(1+i),1],Final.ad[(1+i),2],Final.ad[(1+i),3])
#     
#     }
#   q[(j),]=quantile(sur,probs=c(0.025,0.5,0.975),names=F)
#   print(j)
# }

y = res
grid=seq(0.01,max(y),.01)
q=matrix(0,nrow=(length(grid)),ncol=3)
sur=numeric(200)

for(j in 1:length(grid)){
  for(i in 1:200){
    sur[i]=exp.dis(grid[j],
                   alpha = exp(Final.res[(1+i),1]),
                   theta = exp(Final.res[(1+i),2]),
                   sigma = exp(Final.res[(1+i),3] + Final.res[(1+i),4]))
    #sur[i]=exp.dis(grid[j],Final.ad[(1+i),1],Final.ad[(1+i),2],Final.ad[(1+i),3])
    
  }
  q[(j),]=quantile(sur, probs=c(0.025,0.5,0.975), names=F)
  print(j)
}

#grid.ad = grid
#sur.ad = q

grid.res = grid
sur.res = q


#plot((365*grid.ad),(1-sur.ad[,3]),type="l",lty="dashed",col=1,xlab="Time (in years)",ylab="Survival Function")
#lines((365*grid.ad),(1-sur.ad[,2]),type="l",col=3)
#lines((365*grid.ad),(1-sur.ad[,1]),type="l",lty="dashed",col=1)
#lines((365*grid.ad),(1-sur.ad[,3]),type="l",lty="dashed",col=1)


plot((365*grid.res),(1-sur.res[,3]),type="l",lty="dashed",col=1,xlab="Time (in years)",ylab="Survival Function")
lines((365*grid.res),(1-sur.res[,2]),type="l",col=3)
lines((365*grid.res),(1-sur.res[,1]),type="l",lty="dashed",col=1)
# lines((365*grid.res),(1-sur.res[,3]),type="l",lty="dashed",col=1)


legend(0,1,c("Ad libitum:","2.5 and 97.5 Quantiles","50 Quantile"),lty=c("blank","dashed","solid"),col=c(1,4,1),)


###################################################  Density  #### pdf.r ### pdf.a
grid=seq(0.01,max(y),0.01) 
q=matrix(0,nrow=(length(grid)),ncol=3)
den=numeric(200)

for(j in 1:length(grid)){
  for(i in 1:200){
    #den[i]=exp.pdf(grid[j],Final.ad[(1+i),1],Final.ad[(1+i),2],Final.ad[(1+i),3])
    den[i]=exp.pdf(grid[j],
                   alpha = exp(Final.res[(1+i),1]),
                   theta = exp(Final.res[(1+i),2]),
                   sigma = exp(Final.res[(1+i),3] + Final.res[(1+i),4]))
    }
  q[(j),]=quantile(den,probs=c(0.025,0.5,.975),names=F)
  
  print(j)
}


#pdf.ad=q
#grid.ad=grid
pdf.res=q
grid.res=grid

#plot(density(ad*365,from=0),lwd=2,col=3,xlab="Time (in days)",main="",ylab="Density")
#hist(y*365,freq=FALSE,xlab = "Time (in days)",ylab="Density")
#lines(grid.ad*365,pdf.ad[,2]/365,type="l",col=3)
#lines(grid.ad*365,pdf.ad[,1]/365,type="l",lty="dashed",col=1)
#lines(grid.ad*365,pdf.ad[,3]/365,type="l",lty="dashed",col=1)
#legend(0,0.004,c("Ad libitum:","1st and 99th Quantiles","50th Quantile","Data Density"),lty=c("blank","dashed","solid","solid"),lwd=c(1,1,1,2),col=c(1,1,3,3))

hist(y*365,freq=FALSE,xlab = "Time (in days)",ylab="Density")
lines(grid.res*365,pdf.res[,2]/365,type="l",col=3)
lines(grid.res*365,pdf.res[,1]/365,type="l",lty="dashed",col=1)
lines(grid.res*365,pdf.res[,3]/365,type="l",lty="dashed",col=1)




###################################################  Hazard Rate  ### hr.r ### hr.a
#grid=seq(0.01,max(ad),.1)
q=matrix(0,nrow=(length(grid)),ncol=3)
haz=numeric(200)

for(j in 1:length(grid)){
  for(i in 1:200){
    haz[i]=exp.hr(grid[j],
                  alpha = exp(Final.res[(1+i),1]),
                  theta = exp(Final.res[(1+i),2]),
                  sigma = exp(Final.res[(1+i),3] + Final.res[(1+i),4]))
  }
  q[(j),]=quantile(haz,probs=c(0.025,0.5,0.975),names=F)
  print(j)
}


hr.ad=q

plot(grid,hr.ad[,3],type="l",lty="dashed",xlim=c(0,4),col=2,xlab="Time (in years)",ylab="HR")
lines(grid,hr.ad[,2],type="l",col=3)
lines(grid,hr.ad[,1],type="l",lty="dashed",col=3)
lines(grid,hr.ad[,3],type="l",lty="dashed",col=3)

for(j in 1:length(grid)){
  for(i in 1:200){
    haz[i]=exp.hr(grid[j],
                  alpha = exp(Final.res[(1+i),1]),
                  theta = exp(Final.res[(1+i),2]),
                  sigma = exp(Final.res[(1+i),3] + Final.res[(1+i),4]))
  }
  q[(j),]=quantile(haz,probs=c(0.025,0.5,0.975),names=F)
  print(j)
}
hr.res=q

plot(grid.res,hr.res[,3],type="l",lty="dashed",xlim=c(0,4),col=2,xlab="Time (in years)",ylab="HR")
lines(grid.res,hr.res[,2],type="l",col=3)
lines(grid.res,hr.res[,1],type="l",lty="dashed",col=3)
lines(grid.res,hr.res[,3],type="l",lty="dashed",col=3)


##################################################  MRL   ## mrl.r ### mrl.a


m.expweibull=function(x,alpha,theta,sigma){
  M=0.01
  max=max(x)
  grid=seq(0.01,max,M)
  T=length(grid)
  m=numeric(T)
  mrl=matrix(0,nrow=length(alpha),ncol=T)
  
  for(j in 1:length(alpha)){
    exp.sur=function(t){1-exp.dis(t,alpha[j],theta[j],sigma[j])}
    mu=integrate(exp.sur,0,Inf,stop.on.error=FALSE)$value
    
    for(i in 1:T){
      m[i]=(exp.sur(grid[i]))^(-1)*(mu)  - ( (exp.sur(grid[i]))^(-1)*(integrate(exp.sur,0,grid[i])$value))
    }
    
    mrl[j,]=m
    print(j)
    
  }
  mrl=apply(mrl,2,quantile,probs=c(.025,.5,.975),names=FALSE)
  return(mrl)
}

# mrl.ad = m.expweibull(grid.ad,Final.ad[-1,1],Final.ad[-1,2], Final.ad[-1,3])
# 
# 
# 
# plot(grid.ad,mrl.ad[2,],type="l",lty="dashed",xlim=c(0,4),ylim=c(-0.1,3.1),xlab="Time (in years)",ylab="Mean Residual Life Function",col=4)
# lines(grid.ad,mrl.ad[2,],type="l",col=3)
# lines(grid.ad,mrl.ad[1,],type="l",lty="dashed",col=3)
# lines(grid.ad,mrl.ad[3,],type="l",lty="dashed",col=3)


mrl.res = m.expweibull(grid,
                       alpha = exp(Final.res[(1+i),1]),
                       theta = exp(Final.res[(1+i),2]),
                       sigma = exp(Final.res[(1+i),3] + Final.res[(1+i),4]))



plot(grid.res,mrl.res[2,],type="l",lty="dashed",xlim=c(0,4),ylim=c(-0.1,3.1),xlab="Time (in years)",ylab="Mean Residual Life Function",col=4)
lines(grid.res,mrl.res[2,],type="l",col=3)
lines(grid.res,mrl.res[1,],type="l",lty="dashed",col=3)
lines(grid.res,mrl.res[3,],type="l",lty="dashed",col=3)


################ Shape of MRL ###########

#For adlibitum group alpha > 1 always which means mrl is DCR or UBT (or Weibull for theta = 1)

UBT = DCR =0
prod = numeric()
for(i in 1:200){
  prod[i]=Final.ad[(i+1),1]*Final.ad[(i+1),2]
  if(prod[i]<1 && Final.ad[(i+1),2] <1 ){
   
    UBT=UBT+1
  }else{DCR=DCR+1}
  
}
