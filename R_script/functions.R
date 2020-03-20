
# Last Edited: 03/25/2020


#### remove day 1 and convert data tpyes ---------------------------- ####

clean_data <- function(data, cols_to_numeric){
  data %>% 
    mutate_all(function(x) as.factor(x)) %>% 
    mutate_at(cols_to_numeric, function(x) as.numeric(as.character(x))) %>% 
    filter(`Day after spray` != 1 & Treatment != "Control") %>% 
    droplevels()
}


#### convert death counts to death status -------------------------- ####

convert_to_death_status <- function(data){
  data[rep(1:nrow(data), data$Deaths),] %>% 
    mutate(Death_status = 1) %>% 
    select(-Deaths)
}


#### get total time on test and quantile --------------------------- ####

get_ttt_quantile <- function(ttt, length_out){
  ttt %>% 
    quantile(probs = seq(0, 1, length.out = length_out)) %>% 
    scales::rescale(to = c(0,1))
}


#### plot total time on test and quantile -------------------------- ####

plot_ttt_quantile <- function(data){
  data %>% 
  ggplot(aes(x = u, y = ttt, color = tag)) +
    geom_line(size = 0.9) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.7, linetype = "dotted", size = 0.7) +
    scale_color_jco() +
    facet_grid(`Age at spray`~Sex) +
    labs(x = "Quantile", y = "Scaled total time on test", color = "") +
    theme_bw()
}












