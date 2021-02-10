library(tidyverse)
library(wesanderson)

getEdgeData <- function(this_f){
  points <- read_tsv(this_f[grep(".ptab", this_f)])
  edges <- read_tsv(this_f[grep(".edges", this_f)])
  verts <- edges$vertices %>% strsplit(",")
  edges$v1 <- sapply(verts, FUN=function(x)x[1]) %>% as.numeric
  edges$v2 <- sapply(verts, FUN=function(x)x[2]) %>% as.numeric
  match1 <- match(edges$v1, points$ind)
  match2 <- match(edges$v2, points$ind)
  edges$x1 <- points$x[match1]
  edges$x2 <- points$x[match2]
  edges$y1 <- points$y[match1]
  edges$y2 <- points$y[match2]
  edges$normx <- edges$x1 - edges$x2
  edges$normy <- edges$y1 - edges$y2
  edges$origin <- 0
  edges$mean_x <- 0.5*(edges$x1 + edges$x2)
  edges$mean_y <- 0.5*(edges$y1 + edges$y2)
  edges$time <- i
  return(edges)
}


d <- "/home/carmoma/vertex/Vertex-model/dpygrad_mode247/etournay1_strings8b_all_final/"
setwd(d)
f <- list.files() %>% subset(grepl("dpygrad", .) & !grepl("png$|avi$", .))
time <- str_match(f, "_[0-9]+_") %>% gsub("_", "", .) %>% as.numeric

res<-data.frame();

#for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
for(i in unique(time)){
  this_f <- f[time==i] %>% subset(!is.na(.))
  edges <- getEdgeData(this_f)
  edges2 <- edges %>% filter(type %in% c(7) & ! is.na(tension))
  res <- rbind(res, edges2)
}  

conditions <- c("min = 0; max = 800; exp = 0.1",
                "min = 0; max = 800; exp = 0.75",
                "min = 0; max = 800; exp = 2.0",
                "min = 0; max = 800; exp = 5.0",
                "min = 1; max = 800; exp = 0.1",
                "min = 1; max = 800; exp = 0.75",
                "min = 1; max = 800; exp = 2.0",
                "min = 1; max = 800; exp = 5.0",
                "min = 10; max = 800; exp = 0.1",
                "min = 10; max = 800; exp = 0.75",
                "min = 10; max = 800; exp = 2.0",
                "min = 10; max = 800; exp = 5.0",
                "min = 100; max = 800; exp = 0.1",
                "min = 100; max = 800; exp = 0.75",
                "min = 100; max = 800; exp = 2.0",
                "min = 100; max = 800; exp = 5.0")
inds <- res %>% filter(tension==2000) %>% pull(ind) %>% min

data <- res %>% filter(tension !=2000) %>% mutate(compartment = ifelse(ind < inds, "Anterior", "Posterior"))
data <- data %>% rename(condition = time) %>% mutate(condition = conditions[condition + 1])

g1<- ggplot(data, aes(x=mean_x, y = tension, fill=compartment, color=compartment))+
  geom_line(size=1) +
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

pdf("CuticleTensionAll247.pdf", width = 10, height = 6.6)
print(g1)
dev.off()

