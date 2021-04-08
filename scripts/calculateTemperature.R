library(tidyverse)
#Energy quantiles for positive energy:

#0%           5%          10%          25%          50%          75% 
 # 4.720390e-09 1.614352e-04 3.863012e-04 1.583512e-03 5.874370e-03 2.175870e-02 
#90%          95%         100% 
#7.180116e-02 1.400400e-01 4.113520e+00 

#Energy quantiles for negative energy
#quantile(e$V5, c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0))
#0%            5%           10%           25%           50% 
#-4.615070e+00 -9.190113e-02 -4.324130e-02 -1.367753e-02 -4.186400e-03 
#75%           90%           95%          100% 
#  -1.320780e-03 -3.820660e-04 -1.496433e-04 -2.445490e-08

#Calculated from dpygrad_mode326_0 (307_0 -like), with etournay1_strings8b, with a sample from the first 15 million trials or se  
f <- function(x,t){ifelse(x<0, exp(x/t),exp(-x/t) )}
f <- function(x,t){exp(-x/t) }

temps <- c(0.000001,0.00001,0.00005)
x <- seq(0, 0.14, 0.14/500)

res <- data.frame()
for(t in temps){
y <- f(x, t)
aux <- data.frame(energy_diff = x, temp = t, value = y)
res <- rbind(res, aux)
}

gg<-ggplot(res, aes(x=energy_diff, y=value))+
  geom_line(aes(col=as.factor(temp))) +
  ylim(0, 1) +
  xlim(0,0.001) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
ggsave("prob_acceptance_with_energy_diff.pdf", gg)
   

