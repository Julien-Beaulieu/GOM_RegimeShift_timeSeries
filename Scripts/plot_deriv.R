# This scrpt is to make a function to plot derivative trends. The input of this 
# function is a data frame given by the Deriv function


plot_deriv <- function(data,ylab,col = "darkblue"){
  
  p <- ggplot() +
    geom_line(data = data, aes(y = y_deriv_estim, x = year), color = col, size = 1)+
    geom_ribbon(data = data, aes(y = y_deriv_estim, x = year, ymin = lower_CI, ymax = upper_CI),
                fill = col, color = col, alpha = 0.15)+
    geom_rect(aes(xmin = c(1987,1999.4,2009.5), xmax = c(1992.7,2002,2011.7), 
                  ymin=-Inf, ymax=Inf), 
              fill = "red", alpha = 0.1)+
    ylab(ylab)+
    xlab("Time")+
    geom_hline(yintercept = 0) +
    theme_bw()+
    theme(panel.background = element_blank(), 
          #panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank(),
          axis.text=element_text(size=15),
          axis.title=element_text(size=17)
    )
}