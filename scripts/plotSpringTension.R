
#Plot spring tension
library(tidyverse)
setwd("/home/carmoma/vertex/Vertex-model/etournay1_moveX2")
f <- list.files() %>% subset(grepl(".sprtab",.))


#sprtension_tensionMode2_exp0.5_orderMode1
d <- read_tsv(f[1])
plot(d$x_static, d$tension, xlab = "Spring position in X", 
     ylab = "Spring tension constant", col = as.factor(d$compartment))
plot(d$y_static, d$tension, xlab = "Spring position in Y", 
     ylab = "Spring tension constant", col = as.factor(d$compartment))
