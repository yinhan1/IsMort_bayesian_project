############# Data #############
res=c(105,193,211,236, 302,363,389,390,391,403,530,604, 605, 630, 716, 718, 727, 731, 749, 769, 770, 789, 804, 810, 811, 833, 868, 871, 875, 893, 897, 901, 906, 907, 919, 923, 931, 940, 957, 958, 961, 962, 974, 979, 982, 1001, 1008, 1010, 1011, 1012, 1014, 1017, 1032, 1039, 1045, 1046, 1047, 1057, 1063, 1070, 1073, 1076, 1085, 1090, 1094, 1099, 1107, 1119, 1120, 1128, 1129, 1131, 1133, 1136, 1138, 1144, 1149, 1160, 1166, 1170, 1173, 1181, 1183, 1188, 1190,1203, 1206, 1209, 1218, 1220, 1221, 1228, 1230, 1231, 1233, 1239, 1244, 1258, 1268, 1294, 1316, 1327, 1328, 1369, 1393, 1435)

ad=c(89, 104, 387, 465, 479, 494, 496, 514, 532, 533, 536, 545, 547, 548, 582, 606, 609, 619, 620, 621, 630, 635, 639, 648, 652, 653, 654, 660, 665, 667, 668, 670, 675, 677, 678, 678, 681, 684, 688, 694, 695, 697, 698, 702,704, 710, 711, 712, 715, 716, 717, 720, 721,730,731, 732, 733, 735, 736, 738, 739, 741, 743, 746, 749, 751, 753, 764, 765, 768, 770, 773, 777, 779, 780, 788, 791, 794, 796, 799, 801, 806, 807, 815, 836, 838, 850, 859, 894, 963)

#scale by 365 days

res = res/365
ad = ad/365

###################################################### Inverse CDF method for sampling from the Exponentiated Weibull


r.expweibull = function(n,alpha, theta, sigma) { 
  #alpha, theata = shape; sigma = scale
  Y = numeric(n)
  Y = runif(n)
  t = sigma*(-log(1 - (Y^(1/theta))))^(1/alpha)
  return(t)}


#######################################################  Exponentiated Weibull Density Function

exp.pdf=function(x,alpha,theta,sigma){
  
  d=(alpha*theta/sigma)*((1-exp(-((x/sigma)^alpha)))^(theta-1))*exp(-((x/sigma)^alpha))*(((x/sigma)^(alpha-1)))
  return(d)}


#######################################################  Exponentiated Weibull Distribution Function

exp.dis=function(x,alpha,theta,sigma){
  
  D = ( 1-exp(-((x/sigma)^alpha)))^theta
  return(D)}


############################################################# Metropolis Hastings

likelihood=function(x,alpha,theta,sigma){
  n=length(x)
  plike=1
  if(theta < 1 && ((exp(-(x[1]/sigma)^alpha) ==1))){
    plike=0
  }
  else{
    for(i in 1:n){
      
      plike=sum(log(exp.pdf(x,alpha,theta,sigma)))}
    #plike=sum(log((theta*alpha/sigma)*((1-exp(-(x/sigma)^alpha) )^(theta-1))*(exp(-(x/sigma)^alpha))* ((x/sigma)^(alpha-1) )))}
    #plike=sum(log((theta*alpha/sigma)) + (theta-1)*log(1-exp(-(x/sigma)^alpha)) +  (-(x/sigma)^alpha)  + (alpha-1)*log(x/sigma) )}
  }
  return(plike)}

prior=function(alpha,theta,sigma){
  dexp(alpha,rate=(1/2))*dexp(theta,rate=(1/5))*dexp(sigma,rate=(1/2)) #res
  #dexp(alpha,rate=(1/4))*dexp(theta,rate=(1/1))*dexp(sigma,rate=(1/2)) #ad
}


library(MASS)
m=diag(1,nrow=3)

N=40000


mh=function(x,N){
  #m=diag(1,nrow=3)
  alpha.log=theta.log=sigma.log=q1=as.numeric(N+1)
  pm=matrix(0,nrow=N+1,ncol=3)
  
  #alpha.log[1]=log(8.7)
  #theta.log[1]=log(0.38)
  #sigma.log[1]=log(3.4)
  
  #alpha.log[1]=log(2) #res
  #theta.log[1]=log(5)
  #sigma.log[1]=log(2)
  
  
  #alpha.log[1]=log(4) #ad
  #theta.log[1]=log(1)
  #sigma.log[1]=log(2)
  alpha.log[1]=log(13.8) #ad
  theta.log[1]=log(0.34)
  sigma.log[1]=log(2.2)
  a1=0
  #m = matrix( c(0.04, 0,     0, 
  #              0,    0.07,  0,
  #              0,    0,      0.001),
  #            ncol = 3, nrow=3, byrow = TRUE)
  m = matrix( c(0.028, -0.039,     0.004, 
                -0.039,    0.067,  -0.006,
                0.004,    -0.006,      0.001),
              ncol = 3, nrow=3, byrow = TRUE)
  
  b=2
  
  s=m*b
  for(i in 2:(N+1)){
    p.star =  mvrnorm(1,c(alpha.log[i-1],theta.log[i-1],sigma.log[i-1]),Sigma=s)
    
    
    if(likelihood(x,exp(p.star[1]),exp(p.star[2]),exp(p.star[3])) == 0){
      q1[i]=0
      }else{
      q1[i]=exp(
        (likelihood(x,exp(p.star[1]),exp(p.star[2]),exp(p.star[3])) + log(prior(exp(p.star[1]),exp(p.star[2]),exp(p.star[3]))) ) - (likelihood(x,exp(alpha.log[i-1]),exp(theta.log[i-1]),exp(sigma.log[i-1])) + log(prior(exp(alpha.log[i-1]),exp(theta.log[i-1]),exp(sigma.log[i-1])))  )
      )*exp(p.star[1])*exp(p.star[2])*exp(p.star[3])/(exp(alpha.log[i-1])*exp(theta.log[i-1])*exp(sigma.log[i-1]) )
      }
    
    if(runif(1) < q1[i]){
      alpha.log[i]=p.star[1]
      theta.log[i]=p.star[2]
      sigma.log[i]=p.star[3]
      a1 = a1 + 1
      }else{
      alpha.log[i]=alpha.log[i-1] 
      theta.log[i]=theta.log[i-1]
      sigma.log[i]=sigma.log[i-1]}
    
    pm[i,]=c(exp(alpha.log[i]),exp(theta.log[i]),exp(sigma.log[i]))
    # Prints every 1000th iteration
    if(i %% 1000==0){
      cat(paste0("iteration: ", i, "\n"))}
  }
  pm[1,1]=c(a1/N)
  par(mfrow=c(3,2))
  plot(pm[2:N+1,1],pch=".")
  acf(pm[2:N+1,1])
  plot(pm[2:N+1,2],pch=".")
  acf(pm[2:N+1,2])
  plot(pm[2:N+1,3],pch=".")
  acf(pm[2:N+1,3])
  return(pm)}


test <- mh(res, N)
