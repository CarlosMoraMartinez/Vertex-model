
#Plot spring tension
library(tidyverse)
library(wesanderson)

HINGE_X_LIMITS <- c(373.602, 544.648) #min blade and max hinge coordinates in etournay1_strings8b initial condition


getSpringData <- function(this_f){
  spr <- read_tsv(this_f[grep(".sprtab", this_f)])
  spr$time <- i
  return(spr)
}

d <- "/home/carmoma/vertex/Vertex-model/dpygrad_mode273/etournay1_strings8b_all_final/"
setwd(d)
f <- list.files() %>% subset(grepl("dpygrad", .) & !grepl("png$|avi$", .))
time <- str_match(f, "_[0-9]+_") %>% gsub("_", "", .) %>% as.numeric
conditions <- paste("sim. ", as.character(0:5), sep="")

res<-data.frame();

#for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
for(i in unique(time)){
  this_f <- f[time==i] %>% subset(!is.na(.))
  edges <- getSpringData(this_f)
  edges2 <- edges %>% filter(type %in% EDGE_TYPES & ! is.na(tension))
  res <- rbind(res, edges2)
}  
data <- res %>% rename(condition = time) %>% 
  mutate(condition = conditions[condition + 1],
         compartment = ifelse(y_static - min(y_static) > 0.5*(max(y_static) - min(y_static)), "Anterior", "Posterior")
                                                    )


g1<- ggplot(data, aes(x=x_static, y = tension, fill=compartment, color=compartment))+
  geom_point(size=1) +
  geom_vline(xintercept = HINGE_X_LIMITS, linetype=2) +
  facet_wrap(condition~., nrow = 4)+
  #geom_text(c("hinge", "blade"), aes(x = c(100, 900), y = c(mean(data$tension), mean(data$tension))))+
  #theme_minimal() +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=9, face="bold") #,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values=wes_palette(n=2, name="Royal1")) 

g2 <- ggplot(data, aes(x=x_static, y = y_static, xend=x_movable, yend=y_movable, fill=compartment, color=compartment))+
  geom_segment() +
  geom_point() +
  facet_wrap(condition~., nrow = 4)+
  #theme_minimal() +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=9, face="bold") #,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values=wes_palette(n=2, name="Royal1")) 

pdf("../string_gradient.pdf")
print(g1)
print(g2)
dev.off()

