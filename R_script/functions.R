
# First Edited: 03/19/2020


#### remove day 1 and convert data tpyes ---------------------------- ####

clean_data <- function(data, cols_to_numeric){
  data %>% 
    mutate_all(function(x) as.factor(x)) %>% 
    mutate_at(cols_to_numeric, function(x) as.numeric(as.character(x))) %>% 
    mutate(Sex = recode(Sex, "F" = "Female", "M" = "Male")) %>% 
    filter(`Day after spray` != 1 & Treatment != "Control") %>% 
    droplevels() %>% 
    setnames(old = c("Age at spray", "Day after spray"),
             new = c("Age","Day"))
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
    facet_grid(Age ~ Sex) +
    labs(x = "Quantile", y = "Scaled total time on test", color = "") +
    theme_bw()
}

#### build survreg model ------------------------------------------- ####

get_survreg_model <- function(data, dist, new_formula){
  new_formula <- paste("Surv(Day, Death_status) ~", new_formula) %>% as.formula()
  survreg(new_formula, dist = dist, data = droplevels(data))
}

#### get AIC score for a model ---------------------- ####

get_aic_score <- function(data, dist, new_formula){
  model <- get_survreg_model(data, dist, new_formula)
  AIC(model)
}

#### best subset algorithm ----------------------------------------- ####

get_best_model <- function(data, dist){
  
  # list out all possible subsets excluding intercept only
  variable_list <- c("Treatment","Age","Sex","Treatment:Age","Treatment:Sex","Age:Sex")
  all_combn <- 
    do.call(expand.grid,replicate(length(variable_list),c(TRUE,FALSE),simplify=FALSE)) %>% 
    filter(rowSums(.) != 0)
  
  # create formulas
  formula_list <- 
    apply(all_combn, 1, function(x) variable_list[x]) %>% 
    lapply(function(x) paste0(x, collapse = "+")) %>% unlist()
  
  # fit model and extract AIC score
  tb <- 
    data.frame(all_combn, row.names = NULL) %>% 
    set_colnames(variable_list) %>% 
    mutate_all(as.integer) %>% 
    mutate(AIC = sapply(formula_list, function(x) get_aic_score(data,dist,x)))
  
  # kable 
  tb[tb$AIC == min(tb$AIC),] %>% kable() %>% kable_styling("striped", full_width = F)
}








