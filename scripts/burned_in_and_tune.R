

#### ACO -------------------------------- ####

aco_draws <- read.csv('data/posterior_aco_trial_2020-04-11.csv')
aco_burned <- burn_n_thin_draws(aco_draws, jump = 200, burn_at = 3e3)

aco_burned %>% 
  rowid_to_column('iteration') %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

cov(log(aco_burned[,1:2]))

cov(aco_burned[,-c(1:2)])


apply(aco_burned, 2, median)
apply(aco_burned, 2, mean)

coda::effectiveSize(aco_burned[,1:2])
coda::effectiveSize(aco_burned[,1])


# ggplot for autocorrelation graphs ---------------------------------------
alpha_theta_draws <- co_burned %>% 
  as_tibble() %>% 
  select(alpha, theta) #%>% 
  # burn_n_thin_draws(jump = 260, burn_at = 1e3)

coda::effectiveSize(outcomes[,1:2])
coda::effectiveSize(outcomes[,1])

alpha <- 0.95
alpha.post.acf <- acf(alpha_theta_draws, plot=FALSE, lag.max = 1000)
alpha.acf_df <- data.frame(corr = alpha.post.acf$acf) %>% 
  mutate(lags = 1:n())
conf.lims <- c(-1,1) * qnorm((1 + alpha)/2)/sqrt(alpha.post.acf$n.used)

alpha.acf_df %>% 
  ggplot(aes(x=lags, y = corr.1)) + 
  labs(y="Autocorrelations", x="Lag", title= "Autocorrelation for alpha") +
  geom_segment(aes(xend=lags, yend=0)) +
  # geom_point() +
  # ylim(c(-0.2,0.2)) +
  # xlim(c(800, 1000)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  theme_minimal()









#### CO -------------------------------- ####

co_draws <- read.csv('data/posterior_co.csv')
co_burned <- burn_n_thin_draws(co_draws, jump = 260, burn_at = 1e3)


co_burned %>% 
  rowid_to_column('iteration') %>% 
  # mutate(iteration = ) %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, color = "brown3", size = 0.5, linetype = "dashed") +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

cov(log(co_burned[,1:2]))

cov(co_burned[,-c(1:2)])


apply(co_burned, 2, median)
apply(co_burned, 2, mean)
