## This script is to take values_deriv2_2.csv and make a nice .html table with a
## legend and title to put in the repository.

## author: Julien Beaulieu


library(tidyverse)
library(gt)

## load data

val <- read.csv("time_series_outputs/values_deriv2_2.csv") %>% 
  select(-X)


### make the table

table_html <- val %>%
  gt() %>%
  tab_header(
    title = "Second derivative values for intertidal variables",
    subtitle = md("The regime shift 1 is the one in 1990, the 
                  regime shift 2 is the one in 2000, and the regime shift 3 is 
                  the one in 2010. Upper and Lower 95% CI are the upper and lower 
                  95% credible intervals of the second derivative when it is the
                  most different from zero during each regime shift. The upper_zi 
                  95% CI and lower_zi 95% CI are the same but for the zero inflated 
                  part of the model. Variables are the different species and site 
                  combinations. BC stands for the site Babb’s Cove, NE stands for 
                  NE Appledore, SE stands for SE Appledore, SW stands for SW 
                  Appledore, and SE stands for SE Appledore. Global stands for 
                  the whole island scale.")
  ) %>%
  fmt_number(
    columns = c(upper, lower, upper_zi, lower_zi),
    decimals = 4
  ) %>%
  cols_label(
    var = "Variable",
    regime_shift = "Regime Shift",
    upper = "Upper Bound",
    lower = "Lower Bound",
    upper_zi = "Upper Bound (ZI)",
    lower_zi = "Lower Bound (ZI)"
  )


gtsave(table_html, "Sec_Derivative_Values.html")


