
library(tidyverse)
#Load simulation data
setwd("/home/carmoma/vertex/Vertex-model/rens2/wing2E_landmarks")
d <- read.table("landmarks.csv", sep=",", header=T, stringsAsFactors = FALSE)
d_long <- gather(d, landmark, value, x1:y12, factor_key=FALSE)
dx = d_long %>% filter(grepl("x", landmark)) %>% mutate(landmark = gsub("x", "", landmark))
dy = d_long %>% filter(grepl("y", landmark)) %>% mutate(landmark = gsub("y", "", landmark))
d_long <- merge(dx, dy, by = c("name", "landmark"))
plot(d_long$value.x, d_long$value.y)
d_means <- aggregate(cbind(value.x, value.y) ~ landmark, d_long, function(x) c(mean = mean(x), sd = sd(x))) #can't be used in geom_point
d_long <- d_long %>% group_by(landmark) %>% mutate(mean_x = mean(value.x), 
                                                   mean_y = mean(value.y), 
                                                   sd_x = sd(value.x), 
                                                   sd_y = sd(value.y))

#Plot without procrustes
(gg<- ggplot(d_long, aes()) + 
  xlab("x") + ylab("y") +
  xlim(180, 500) + ylim(45, 180) +
  coord_fixed(ratio = 1, expand = TRUE, clip = "on") +
  geom_point(x=d_long$value.x, y=d_long$value.y, fill="black", size=0.1) +
  geom_point(x=d_long$mean_x, y=d_long$mean_y, color="red", size=2) + 
  geom_errorbarh(aes(xmin=d_long$mean_x - d_long$sd_x,
                     xmax=d_long$mean_x + d_long$sd_x,
                     y=d_long$mean_y), colour="red", stat = "identity", size=1) +
  geom_errorbar(aes(ymin=d_long$mean_y - d_long$sd_y,
                    ymax=d_long$mean_y + d_long$sd_y,
                    x=d_long$mean_x), colour="red", stat = "identity", size=1))

##With procrustes
library(vegan)

doProcrustes<-function(target, d){
  ff<-function(x){matrix(unlist(x), ncol=2, byrow = T) }
  rotate <- d %>% select(matches("([xy][1234]|[xy]12)$", perl=T)) # %>% as.list %>% pmap( ff)
  rotate2 <- lapply(1:nrow(rotate), function(x)ff(rotate[x,]))
  rotate3 <- lapply(rotate2, function(x){
    procrustes(target, x)$Yrot %>% as.data.frame %>% set_names("x", "y") #%>% add_column(landmark=c("1", "2", "3", "4", "12"))
  }) %>% bind_rows %>% add_column(name=rep(d$name, each=5), landmark=rep(c("1", "2", "3", "4", "12"), nrow(d))) %>% 
    group_by(landmark) %>% mutate(mean_x = mean(x), 
                                  mean_y = mean(y), 
                                  sd_x = sd(x), 
                                  sd_y = sd(y))
  return(rotate3)
}
## Load real data
setwd("/home/carmoma/vertex/base_population/")
pop <- list.files() %>% subset(grepl(pattern=".csv", .)) %>% 
  map(function(x)read.table(x,sep=",", header=T, stringsAsFactors=F ) %>% mutate(group=x)) %>%
  bind_rows %>% rename(name = File)
#Mean for each landmark in real Drosophila base population
means <- pop %>% select(matches("^(x|y)", perl=T)) %>%  summarise_all(list(mean=mean, sd=sd))
target <- matrix(means[,grep("([xy][1234]|[xy]12)_mean", colnames(means))] %>% unlist, ncol=2, byrow=T) 
# Load simulation data
setwd("/home/carmoma/vertex/Vertex-model/rens2/wing2E_landmarks")
d <- read.table("landmarks.csv", sep=",", header=T, stringsAsFactors = FALSE)

sims_rotated <- doProcrustes(target, d)
pop_rotated <- doProcrustes(target, pop)

(gg<- ggplot(sims_rotated, aes(x, y)) + 
    geom_point(size=0.5, alpha=0.2) +
    xlab("x") + ylab("y") +
    coord_fixed(ratio = 1, expand = TRUE, clip = "on") +
    geom_point(aes(mean_x, mean_y), color="red", size=2) + 
    geom_errorbarh(aes(xmin=mean_x - sd_x,
                       xmax=mean_x + sd_x,
                       y=mean_y), colour="red", stat = "identity") +
    geom_errorbar(aes(ymin=mean_y - sd_y,
                      ymax=mean_y + sd_y,
                      x=mean_x), colour="red", stat = "identity") +
  theme(element_blank()) + theme_bw()
        )

(gg<- ggplot(pop_rotated , aes(x, y)) + 
    geom_point(size=0.5, alpha=0.2) +
    xlab("x") + ylab("y") +
    coord_fixed(ratio = 1, expand = TRUE, clip = "on") +
    geom_point(aes(mean_x, mean_y), color="red", size=2) + 
    geom_errorbarh(aes(xmin=mean_x - sd_x,
                       xmax=mean_x + sd_x,
                       y=mean_y), colour="red", stat = "identity") +
    geom_errorbar(aes(ymin=mean_y - sd_y,
                      ymax=mean_y + sd_y,
                      x=mean_x), colour="red", stat = "identity") +
    theme(element_blank()) + theme_bw()
)
