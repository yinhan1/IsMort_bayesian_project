#### data ----------------------- ####
res = c(
  105,
  193,
  211,
  236,
  302,
  363,
  389,
  390,
  391,
  403,
  530,
  604,
  605,
  630,
  716,
  718,
  727,
  731,
  749,
  769,
  770,
  789,
  804,
  810,
  811,
  833,
  868,
  871,
  875,
  893,
  897,
  901,
  906,
  907,
  919,
  923,
  931,
  940,
  957,
  958,
  961,
  962,
  974,
  979,
  982,
  1001,
  1008,
  1010,
  1011,
  1012,
  1014,
  1017,
  1032,
  1039,
  1045,
  1046,
  1047,
  1057,
  1063,
  1070,
  1073,
  1076,
  1085,
  1090,
  1094,
  1099,
  1107,
  1119,
  1120,
  1128,
  1129,
  1131,
  1133,
  1136,
  1138,
  1144,
  1149,
  1160,
  1166,
  1170,
  1173,
  1181,
  1183,
  1188,
  1190,
  1203,
  1206,
  1209,
  1218,
  1220,
  1221,
  1228,
  1230,
  1231,
  1233,
  1239,
  1244,
  1258,
  1268,
  1294,
  1316,
  1327,
  1328,
  1369,
  1393,
  1435
)

ad = c(
  89,
  104,
  387,
  465,
  479,
  494,
  496,
  514,
  532,
  533,
  536,
  545,
  547,
  548,
  582,
  606,
  609,
  619,
  620,
  621,
  630,
  635,
  639,
  648,
  652,
  653,
  654,
  660,
  665,
  667,
  668,
  670,
  675,
  677,
  678,
  678,
  681,
  684,
  688,
  694,
  695,
  697,
  698,
  702,
  704,
  710,
  711,
  712,
  715,
  716,
  717,
  720,
  721,
  730,
  731,
  732,
  733,
  735,
  736,
  738,
  739,
  741,
  743,
  746,
  749,
  751,
  753,
  764,
  765,
  768,
  770,
  773,
  777,
  779,
  780,
  788,
  791,
  794,
  796,
  799,
  801,
  806,
  807,
  815,
  836,
  838,
  850,
  859,
  894,
  963
)

res = res / 365
ad = ad / 365
t = c(ad, res)
x = c(rep(0, 90), rep(1, 106))


#### functions ------------------- ####

#### 100,000 TRIALS!

###################################################### Inverse CDF method for sampling from the Exponentiated Weibull


r.expweibull = function(n, alpha, theta, sigma) {
  #alpha, theata = shape; sigma = scale
  Y = numeric(n)
  Y = runif(n)
  t = sigma * (-log(1 - (Y ^ (1 / theta)))) ^ (1 / alpha)
  return(t)
}


#######################################################  Exponentiated Weibull Density Function

exp.pdf = function(x, alpha, theta, sigma) {
  d = (alpha * theta / sigma) * ((1 - exp(-((
    x / sigma
  ) ^ alpha))) ^ (theta - 1)) * exp(-((x / sigma) ^ alpha)) * (((x / sigma) ^
                                                                  (alpha - 1)))
  return(d)
}


#######################################################  Exponentiated Weibull Distribution Function

exp.dis = function(x, alpha, theta, sigma) {
  D = (1 - exp(-((x / sigma) ^ alpha))) ^ theta
  return(D)
}


############################################################# Metropolis Hastings

likelihood = function(x, alpha, theta, sigma) {
  n = length(x)
  plike = 1
  if (theta < 1 && ((exp(-(x[1] / sigma) ^ alpha) == 1))) {
    plike = 0
  }
  else{
    for (i in 1:n) {
      plike = sum(log(exp.pdf(x, alpha, theta, sigma)))
    }
  }
  return(plike)
}

prior = function(alpha, theta, beta0, beta1) {
  dexp(alpha, rate = (1 / 3)) * dexp(theta, rate = (1 / 1)) * dnorm(beta0, 1, 100) *
    dnorm(beta1, 0, 100)
}



N = 200

mh = function(t, x, N) {
  #m=diag(1,nrow=3)
  
  q1 = as.numeric(N + 1)
  pm = matrix(0, nrow = N + 1, ncol = 4)
  
  
  # alpha.log.cur = log(3)
  # theta.log.cur = 0 #log(.2)
  # beta0.cur =   1 #-0.44 #-.45
  # beta1.cur =  0 #1.22#1.16
  
  alpha.log.cur = log(10)
  theta.log.cur = log(0.3) #log(.2)
  beta0.cur =   0.8 #-0.44 #-.45
  beta1.cur =  0.4 #1.22#1.16
  ######################## alpha/theta #################
  a1 = 0
  m = matrix(
    c(0.0136,-0.0155,-0.0155,    0.0257),
    ncol = 2,
    nrow = 2,
    byrow = TRUE
  )
  
  # b = .5
  b = 0.05
  s = m * b
  ######################## beta0/beta1 #################
  
  a2 = 0
  m2 = matrix(
    c(.00024, 0,
      0,  .00059),
    ncol = 2,
    nrow = 2,
    byrow = TRUE
  )
  
  # b2 = .5
  b2 = 0.05
  s2 = m2 * b2
  
  for (i in 2:(N + 1))
  {
    #p.star =  mvrnorm(1,c(alpha.log.cur,theta.log.cur, beta0.cur, beta1.cur),Sigma=s2)
    
    
    #    if(likelihood(t,exp(p.star[1]),exp(p.star[2]),exp(p.star[3] + p.star[4]*x)) == 0)
    #    {
    #      q1=0}
    
    #    else{
    
    #      q1=exp(
    
    #        (likelihood(t,exp(p.star[1]),exp(p.star[2]),exp(p.star[3] + p.star[4]*x)) + log(prior(exp(p.star[1]),exp(p.star[2]),p.star[3], p.star[4]) ) ) - (likelihood(t,exp(alpha.log.cur),exp(theta.log.cur),exp(beta0.cur + beta1.cur*x)) +log(prior(exp(alpha.log.cur),exp(theta.log.cur),beta0.cur,beta1.cur))  )
    
    #      )*exp(p.star[1])*exp(p.star[2])/(exp(alpha.log.cur)*exp(theta.log.cur) )
    #    }
    #    if(runif(1) < q1) {
    #      alpha.log.cur=p.star[1]
    #      theta.log.cur=p.star[2]
    #      beta0.cur = p.star[3]
    #      beta1.cur = p.star[4]
    #      a1 = a1 + 1}
    
    ######################## MH for alpha/theta #################
    p.star =  mvrnorm(1, c(alpha.log.cur, theta.log.cur), Sigma = s)
    
    if (likelihood(t, exp(p.star[1]), exp(p.star[2]), exp(beta0.cur + beta1.cur *
                                                          x)) == 0) {
      q1 = 0
    } else{
      q1 = exp((likelihood(
        t,
        exp(p.star[1]),
        exp(p.star[2]),
        exp(beta0.cur + beta1.cur * x)
      ) + log(
        prior(exp(p.star[1]), exp(p.star[2]), beta0.cur, beta1.cur)
      )) - (likelihood(
        t,
        exp(alpha.log.cur),
        exp(theta.log.cur),
        exp(beta0.cur + beta1.cur * x)
      ) + log(
        prior(
          exp(alpha.log.cur),
          exp(theta.log.cur),
          beta0.cur,
          beta1.cur
        )
      ))) * exp(p.star[1]) * exp(p.star[2]) / (exp(alpha.log.cur) * exp(theta.log.cur))
    }
    
    if (runif(1) < q1) {
      alpha.log.cur = p.star[1]
      theta.log.cur = p.star[2]
      a1 = a1 + 1
    }
    
    
    ######################## MH for beta0/beta1 #################
    
    beta.star =  mvrnorm(1, c(beta0.cur, beta1.cur), Sigma = s2)
    
    
    if (likelihood(t,
                   exp(alpha.log.cur),
                   exp(theta.log.cur),
                   exp(beta.star[1] + beta.star[2] * x)) == 0) {
      q2 = 0
    } else{
      q2 = exp((likelihood(
        t,
        exp(alpha.log.cur),
        exp(theta.log.cur),
        exp(beta.star[1] + beta.star[2] * x)
      )) + log(prior(
        exp(alpha.log.cur),
        exp(theta.log.cur),
        beta0.cur,
        beta1.cur
      )) - (likelihood(
        t,
        exp(alpha.log.cur),
        exp(theta.log.cur),
        exp(beta0.cur + beta1.cur * x)
      ) + log(
        prior(
          exp(alpha.log.cur),
          exp(theta.log.cur),
          beta0.cur,
          beta1.cur
        )
      )))
     }
    if (runif(1) < q2)
    {
      beta0.cur = beta.star[1]
      beta1.cur = beta.star[2]
      a2 = a2 + 1
    }
    
    pm[i, ] = c(alpha.log.cur, theta.log.cur, beta0.cur, beta1.cur)
    
    if (i %% 100 == 0) {
      cat(paste0("iteration: ", i, "\n"))
    }
    
  }
  
  pm[1, 1] = c(a1 / N)
  pm[1, 2] = c(a2 / N)
  par(mfrow = c(4, 2))
  plot(exp(pm[2:N + 1, 1]), pch = ".")
  acf(pm[2:N + 1, 1])
  plot(exp(pm[2:N + 1, 2]), pch = ".")
  acf((pm[2:N + 1, 2]))
  plot(pm[2:N + 1, 3], pch = ".")
  acf(pm[2:N + 1, 3])
  plot(pm[2:N + 1, 4], pch = ".")
  acf(pm[2:N + 1, 4])
  return(pm)
}
library(MASS)
Final_res = mh(t, x, N=2000)

