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

#FILE_BASE_STRING <- "etournay1_strings9_moved_"
#FILE_BASE_STRING <- "21grid[0-9]_moved_"
FILE_BASE_STRING <- "cellcuticle1_moved_"
ZOOM_POINT <- c(300, 350, 320, 380)
CUTICLE_TYPE = 8 #7 for string cuticle, 8 for full cuticle
TYPES <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)

dirs_all <- c("/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_0/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_1/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_2/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_3/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_4/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_5/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_6/cellcuticle1",
              "/home/carmoma/vertex/Vertex-model/cellcuticle_11/cellcuticle_11_7/cellcuticle1"
          
)

dirs <- dirs_all[c(1,4,7,10)]
dirs <- dirs_all[c(1, 2)]
dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode337/dpygrad_mode337_0/etournay1_strings9",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_3/etournay1_strings8b"
)

TIME_POINTS =c(0, 1, 5, 10)
for(d in dirs){
  setwd(d)
## READ CELL FILES
  f<-list.files() %>% subset(grepl("moved",.) & !grepl(".avi$", .)) #%>%subset(grepl(".celltab", .)) 
  time <- gsub(FILE_BASE_STRING, "", f) %>% gsub("\\.[a-z]+", "", ., perl=T) %>% as.numeric

  res<-data.frame();

#for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
for(i in TIME_POINTS){
this_f <- f[time==i] %>% subset(!is.na(.))
  edges <- getEdgeData(this_f)
  edges2 <- edges %>% filter(type %in% TYPES & ! is.na(tension))

#Full Wing
  edges$log_tension <- ifelse(edges$tension < 0, 0, log10(edges$tension))
  g1<-edges %>% ggplot(aes(x=x1, y=y1, fill=tension, col=tension)) + 
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
   scale_color_gradient(low="blue", high="red",limits = c(min(edges$tension), max(edges$tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(edges$tension), max(edges$tension)),oob = scales::squish)

  g1_log<-edges %>% ggplot(aes(x=x1, y=y1, fill=log_tension, col=log_tension)) + 
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
    scale_color_gradient(low="blue", high="red",limits = c(min(edges$log_tension), max(edges$log_tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(edges$log_tension), max(edges$log_tension)),oob = scales::squish)
  
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
    scale_color_gradient(low="blue", high="red",limits = c(min(edges$tension), max(edges$tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(edges$tension), max(edges$tension)),oob = scales::squish) +
    xlim(ZOOM_POINT[1], ZOOM_POINT[2]) +
    ylim(ZOOM_POINT[3], ZOOM_POINT[4])
  
  #Cuticle
  e2 <- edges %>% filter(type==CUTICLE_TYPE)
  e3 <- e2 %>% gather("coord", "value", x1, x2, y1, y2) %>% mutate(coord=gsub("y", "y_", gsub("x", "x_", coord))) %>% separate(coord,c("coord", "coordnum"), sep="_") %>% spread(coord, value)
  
  (g3<-e2 %>% ggplot(aes(x=x1, y=y1, fill=tension, col=tension)) + 
    geom_segment(aes(xend=x2, yend=y2)) +
    geom_point(data=e3, aes(x=x, y=y)) +
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
    scale_color_gradient(low="blue", high="red",limits = c(min(e2$tension), max(e2$tension)),oob = scales::squish) +
    scale_fill_gradient(low="blue", high="red",limits = c(min(e2$tension), max(e2$tension)),oob = scales::squish)
  )
  name1 <- paste("../../", str_extract(d, "[1-9]+_[1-9]"), "_time", as.character(i), "_fullWing.png", sep="", collapse="")
  name1b <- paste("../../", str_extract(d, "[1-9]+_[1-9]"), "_time", as.character(i), "_fullWingLog.png", sep="", collapse="")
  name2 <- paste("../../", str_extract(d, "[1-9]+_[1-9]"), "_time", as.character(i), "_HingeCaption.png", sep="", collapse="")
  png(name1, width = 1500, height = 500)
  print(g1)
  dev.off()
  png(name1b, width = 1500, height = 500)
  print(g1_log)
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
  





