
cdf = function(t, alpha, theta, sigma){
  (1 - exp(-(t/sigma)^alpha))^theta
}

pdf = function(t, alpha, theta, sigma){
  (alpha*theta/sigma) * 
    ((1 - exp(-(t/sigma)^alpha))^(theta-1)) * 
    (exp(-(t/sigma)^alpha)) * 
    ((t/sigma)^(alpha-1))
}

survival_function = function(t, alpha, theta, sigma){
  1 - cdf(t, alpha, theta, sigma)
}

hazard_rate = function(t, alpha, theta, sigma){
  pdf(t, alpha, theta, sigma) / (1 - cdf(t, alpha, theta, sigma))
}

library(tidyverse)
library(ggsci)
library(latex2exp)

data.frame(t = seq(0, 120)) %>% 
  mutate(s1 = survival_function(t, alpha=0.5, theta=2, sigma=5),
         s2 = survival_function(t, alpha=2, theta=2, sigma=50),
         s3 = survival_function(t, alpha=0.55, theta=4, sigma=5),
         s4 = survival_function(t, alpha=4, theta=0.15, sigma=150)) %>% 
  reshape2::melt(id.vars = "t") %>% 
  ggplot(aes(x = t, y = value, color = variable)) + 
  geom_line(size = 1.5) +
  scale_color_jco(name = TeX("--------- $\\alpha$ ------ $\\theta$ ------ $\\sigma$"),
                     labels = c("    0.5           2            5",
                                "      2            2           50",
                                "   0.55          4            5",
                                "      4          0.15       150")) +
  theme_minimal() +
  labs(x = "Time t", y = "Survival function S(t)") 




