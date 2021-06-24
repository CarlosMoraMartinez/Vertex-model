
library(tidyverse)
setwd("/home/carmoma/vertex/Vertex-model/etournay1_strings6/")
setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode167/dpygrad_mode167_0/etournay1_strings4")

interest_cells <- c(10989, 11460,11051,11052,10991)
interest_cells <- c(7998, 8397)
interest_cells <- c(2687, 2308, 536, 867, 1097, 538, 3077, 4136, 4945, 5351) #prox inf corner initial cond
interest_cells <- c(7530,7935,8727,9867,9117,7937, 5351, 7670, 7535, 4945, 2687) #prox top corner initial cond
interest_cells <- c(861, 975, 750, 639, 531, 426)
interest_cells <- c(9859, 9860, 9982)

f <- list.files() %>% subset(grepl("_0\\.",.))

f2 <- list.files() %>% subset(grepl("_22\\.",.))


points <-read.table(subset(f, grepl(".ptab", f)), sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(cells = strsplit(cells, ","),
         edges = strsplit(edges, ","),
         neighbour_vertices = strsplit(neighbour_vertices, ","))
cells <-read.table(subset(f, grepl(".celltab", f)), sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(edges = strsplit(edges, ","),
         vertices = strsplit(vertices, ","))
edges <-read.table(subset(f, grepl(".edges", f)), sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(cells = strsplit(cells, ","),
         vertices = strsplit(vertices, ","))
spr <-read.table(subset(f, grepl(".sprtab", f)), sep="\t", dec=".", header=T, stringsAsFactors = F)


points2 <-read.table(f2[5], sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(cells = strsplit(cells, ","),
         edges = strsplit(edges, ","),
         neighbour_vertices = strsplit(neighbour_vertices, ","))
cells2 <-read.table(f2[2], sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(edges = strsplit(edges, ","),
         vertices = strsplit(vertices, ","))
edges2 <-read.table(f2[3], sep="\t", dec=".", header=T, stringsAsFactors = F)%>% 
  mutate(cells = strsplit(cells, ","),
         vertices = strsplit(vertices, ","))
spr2 <-read.table(f2[7], sep="\t", dec=".", header=T, stringsAsFactors = F)



cint <- cells %>% filter(ind %in% interest_cells)

colors <- c("red", "blue", "green", "black", "green")
#vert2plot <- strsplit(cint$vertices, ",")
vert2plot <- cint$vertices
res <- data.frame()
for(i in 1:length(interest_cells)){
  vthis <- points[match(vert2plot[[i]], points$ind) ,] #%>% filter(ind %in% vert2plot[[i]])
  vthis$cell <- interest_cells[i]
  vthis$color <- colors[i]
  res <- rbind(res, vthis)
}

#plot border
border <- edges %>% filter(type %in% c(4)) 
border$x0 <- points$x[match(sapply(border$vertices, FUN=function(x){x[1]}) %>% as.integer(), points$ind)]
border$x1<- points$x[match(sapply(border$vertices, FUN=function(x){x[2]}) %>% as.integer(), points$ind)]
border$y0 <- points$y[match(sapply(border$vertices, FUN=function(x){x[1]}) %>% as.integer(), points$ind)]
border$y1 <- points$y[match(sapply(border$vertices, FUN=function(x){x[2]}) %>% as.integer(), points$ind)]


res$ind2<- as.character(res$ind)
p <- ggplot() +
  geom_polygon(data=res, aes(x = x, y = y, col="black", fill = color, group = cell, alpha=0.1)) + 
  geom_label(data=res,aes(x = x, y = y, label=ind, alpha = 0.1)) + 
  geom_segment(data=border, aes(x=x0, y=y0, xend=x1, yend=y1, col=as.character(type)))#+
  #xlim(689,691) + ylim(422.8,423.05)
p



edges[edges$ind %in% c(24925, 25323, 24528, 23729), ]