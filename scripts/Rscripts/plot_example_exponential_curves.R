

library(tidyverse)
a <- seq(0, 1, by = 0.01)
es <- c(0.1,0.5,1,2,5,10)
c <- expand.grid(a, es)
names(c)<-c("x", "exp")
c$multipl <- 1 - exp(-c$exp*c$x)
c$power <- (1 - exp(-c$x^c$exp))/(1 - exp(-1))

d <- c %>% gather("formula", "result", power, multipl)
d$fact_exp <- as.factor(as.character(d$exp))
d$expression <- ifelse(d$formula == "power", "(1 - exp(-x^alpha))/(1 - exp(-1))", " 1 - exp(-alpha*x)")

ggplot(d, aes(x=x, y=result, fill=fact_exp, col=fact_exp)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(.~expression)+
  ggtitle("Gradient shape with different formulas and exponents (we used formula in the right side)")+
  ylab("Tension") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 
    

