#'@rdname setmytheme
#'@title Customize graphical settings
#'@description This function customizes the theme used by ggplot2. 
#' Loads the packages ggplot2, ggthemes, latex2exp
#'@export
setmytheme <- function(){
  library(ggplot2)
  library(ggthemes)
  library(latex2exp)
  theme_set(theme_tufte())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                               colour = "gray"), 
               panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                               colour = "gray"),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               legend.position = "bottom")
}
