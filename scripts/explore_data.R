library(tidyverse)

IS_df <- readxl::read_excel('data/IS_data.xlsx')
glimpse(IS_df)

skimr::skim(IS_df)

names(IS_df) <- names(IS_df) %>% stringr::str_to_lower() %>% 
  stringr::str_replace_all(' ', '_')

names(IS_df)
IS_df %>% 
  select_if(is.numeric) %>% 
  gather(cols, value) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram() + 
  facet_grid(.~cols, scales = 'free_x') +
  theme_classic()



IS_df %>% 
  ggplot() +
  geom_point(aes(x = day_after_spray, y = deaths))
