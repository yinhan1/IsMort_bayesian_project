


co_draws <- read.csv('data/posterior_co_trial_2020-04-12.csv')
co_burned <- burn_n_thin_draws(co_draws, jump = 1, burn_at = 1e5)


co_burned %>% 
  rowid_to_column('iteration') %>% 
  # mutate(iteration = ) %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  # geom_hline(yintercept = 0, color = "brown3", size = 0.5, linetype = "dashed") +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

cov(log(co_burned[,1:2]))

cov(co_burned[,-c(1:2)])
