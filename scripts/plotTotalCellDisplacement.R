

library(tidyverse)
library(wesanderson)


FILE_BASE_STRING <- "etournay1_strings10_moved_"
DIF_TIMES <- c(0, 40)

getCellData <- function(this_f){
  cells <- read_tsv(this_f[grep(".celltab", this_f)])
  cells$time <- i
  return(cells)
} 


dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode337/dpygrad_mode337_0/etournay1_strings10",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode325/dpygrad_mode325_3/etournay1_strings8b"
)

#for(d in dirs){
d <- dirs[1]
  setwd(d)
  ## READ CELL FILES
  f<-list.files() %>% subset(grepl("moved",.) & !grepl(".avi$", .)) #%>%subset(grepl(".celltab", .)) 
  time <- gsub(FILE_BASE_STRING, "", f) %>% gsub("\\.[a-z]+", "", ., perl=T) %>% as.numeric
  
  res<-data.frame();
  
  #for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
  for(i in c(0, 40)){
    this_f <- f[time==i] %>% subset(!is.na(.))
    cells <- getCellData(this_f)
    res <- rbind(res, cells)
  }
  
  a <- res %>% filter(time==0)
  b <- res %>% filter(time==40)
  
  a$x2 <-b$centroid_x[match(a$ind, b$ind)]
  a$y2 <-b$centroid_y[match(a$ind, b$ind)]
  a$ind2 <- b$ind[match(a$ind, b$ind)]
a <- a %>% rename(x1 = centroid_x, y1 = centroid_y) %>% 
  mutate(distance = sqrt((x1 - x2)^2 + (y1 - y2)^2)) %>% 
  filter(x2>0 & y2 > 0)


(g1<-a %>% ggplot(aes(x=x1, y=y1, fill=distance, col=distance, alpha=0.6)) + 
  geom_segment(aes(xend=x2, yend=y2)) +
    geom_point(aes(x=x2, y=y2), fill="black", col="black", size=0.2) +
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
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")
  )

#,limits = c(min(edges2$tension), max(edges2$tension)),oob = scales::squish