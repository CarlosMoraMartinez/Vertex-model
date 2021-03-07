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



dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode291/dpygrad_mode291/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode239/dpygrad_mode239_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode239/dpygrad_mode239_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode239/dpygrad_mode239_3/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode239/dpygrad_mode239_4/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode240/dpygrad_mode240_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode240/dpygrad_mode240_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode240/dpygrad_mode240_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode240/dpygrad_mode240_3/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode240/dpygrad_mode240_4/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode242/dpygrad_mode242_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode242/dpygrad_mode242_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode242/dpygrad_mode242_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode242/dpygrad_mode242_3/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode242/dpygrad_mode242_4/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode241/dpygrad_mode241_0/etournay1_strings8b"
)

for(d in dirs){
  setwd(d)
## READ CELL FILES
  f<-list.files() %>% subset(grepl("moved",.)) #%>%subset(grepl(".celltab", .)) 
  time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub("\\.[a-z]+", "", ., perl=T) %>% as.numeric

  res<-data.frame();

#for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
for(i in c(0, 1, 5, 10, 20, 40)){
this_f <- f[time==i] %>% subset(!is.na(.))
  edges <- getEdgeData(this_f)
  edges2 <- edges %>% filter(type %in% c(1) & ! is.na(tension))

#Full Wing
  g1<-ggplot(edges, aes(x=x1, y=y1, fill=tension, col=tension)) + 
   geom_segment(aes(xend=x2, yend=y2)) +
   coord_fixed(ratio = 1) +
   theme_minimal() +
   theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
    ) +
   scale_color_gradient(low="blue", high="red",limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish)
#Only caption
  g2<-ggplot(edges, aes(x=x1, y=y1, fill=tension, col=tension)) + 
   geom_segment(aes(xend=x2, yend=y2)) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
     plot.title = element_text(size=16, face="bold.italic"),
     axis.title.x = element_text( size=16, face="bold"),
     axis.title.y = element_text(size=16, face="bold"),
     strip.text = element_text(size=16, face="bold") ,
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank()
   ) +
    scale_color_gradient(low="blue", high="red",limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish) +
    xlim(200, 250) +
    ylim(200, 300)
  name1 <- paste("../../", str_extract(d, "[1-9]+_[1-9]"), "_time", as.character(i), "_fullWing.png", sep="", collapse="")
  name2 <- paste("../../", str_extract(d, "[1-9]+_[1-9]"), "_time", as.character(i), "_HingeCaption.png", sep="", collapse="")
  png(name1, width = 1500, height = 500)
  print(g1)
  dev.off()
  png(name2)
  print(g2)
  dev.off()
  }
}


ggplot(edges2, aes(x=origin, y=origin, fill=tension, col=tension)) + 
  geom_segment(aes(xend=normx, yend=normy)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 


ggplot(edges2, aes(x=normx, y=normy, fill=tension, col=tension, alpha=tension)) +
  geom_point() +
  scale_color_gradient(low="blue", high="red",limits = c(min(edges2$tension), 0.5),oob = scales::squish) +
  scale_fill_gradient(low="blue", high="red",limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish)
  





