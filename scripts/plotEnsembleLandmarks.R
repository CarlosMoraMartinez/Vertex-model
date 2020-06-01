
library(tidyverse)
setwd("/home/carmoma/vertex/Vertex-model/rens2/wing2E_landmarks")
d <- read.table("landmarks.csv", sep=",", header=T, stringsAsFactors = FALSE)
d_long <- gather(d, landmark, value, x1:y12, factor_key=FALSE)
dx = d_long %>% filter(grepl("x", landmark)) %>% mutate(landmark = gsub("x", "", landmark))
dy = d_long %>% filter(grepl("y", landmark)) %>% mutate(landmark = gsub("y", "", landmark))
d_long <- merge(dx, dy, by = c("name", "landmark"))
plot(d_long$value.x, d_long$value.y)
d_means <- aggregate(cbind(value.x, value.y) ~ landmark, d_long, function(x) c(mean = mean(x), sd = sd(x))) #can't be used in geom_point
d_long <- d_long %>% group_by(landmark) %>% mutate(mean_x = mean(value.x), mean_y = mean(value.y), 
                                                   sd_x = sd(value.x), sd_y = sd(value.y))
plot(means$mean_x, means$mean_y)
gg<- ggplot(d_long, aes()) + 
  xlab("x") + ylab("y") +
  xlim(180, 500) + ylim(45, 180) +
  coord_fixed(ratio = 1, expand = TRUE, clip = "on") +
  geom_point(x=d_long$value.x, y=d_long$value.y, fill="black") +
  geom_point(x=d_long$mean_x, y=d_long$mean_y, color="red", size=2) + 
  geom_errorbarh(aes(xmin=d_long$mean_x - d_long$sd_x,
                     xmax=d_long$mean_x + d_long$sd_x,
                     y=d_long$mean_y), colour="red", stat = "identity", size=1) +
  geom_errorbar(aes(ymin=d_long$mean_y - d_long$sd_y,
                    ymax=d_long$mean_y + d_long$sd_y,
                    x=d_long$mean_x), colour="blue", stat = "identity", size=2, width=2)
gg
