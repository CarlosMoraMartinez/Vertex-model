#plot movement of frontier

library(tidyverse)
library(wesanderson)

getEdgeData <- function(this_f){
  points <- read_tsv(this_f[grep(".ptab", this_f)])
  cells <- read_tsv(this_f[grep(".celltab", this_f)])
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
  
  edgecells <- edges$cells %>% strsplit(",")
  edges$cell1 <- sapply(edgecells, FUN=function(x)x[1]) %>% as.numeric
  edges$cell2 <- sapply(edgecells, FUN=function(x)x[2]) %>% as.numeric
  match1 <- match(edges$cell1, cells$ind)
  match2 <- match(edges$cell2, cells$ind)
  edges$celltype1 <- cells$type[match1]
  edges$celltype2 <- cells$type[match2]
  
  edges$origin <- 0
  edges$mean_x <- 0.5*(edges$x1 + edges$x2)
  edges$mean_y <- 0.5*(edges$y1 + edges$y2)
  edges$time <- i
  return(edges)
}


dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode313/dpygrad_mode313_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode239/dpygrad_mode239_1/etournay1_strings8b"
)
for(d in dirs){
  setwd(d)
  ## READ CELL FILES
  f<-list.files() %>% subset(grepl("moved",.)) #%>%subset(grepl(".celltab", .)) 
  time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub("\\.[a-z]+", "", ., perl=T) %>% as.numeric
  
  res<-data.frame();
  
  #for(i in min(time[!is.na(time)]):max(time[!is.na(time)])){
  for(i in c(0:60)){
    this_f <- f[time==i] %>% subset(!is.na(.))
    edges <- getEdgeData(this_f)
    res <- rbind(res, edges)
  }
}

EdgeType <- c("blade", "hinge", "vein_hinge", "vein_blade", "tissue_boundary", 
"spring", "vein", "stringedge")

res <- res %>% mutate(type = EdgeType[type+1])
r2 <- res %>% filter(type=="tissue_boundary" & (celltype1==1 | celltype2 ==1) )
r3 <- res %>% filter((celltype1==1 & celltype2 ==0) | (celltype1==0 & celltype2 ==1))


r4 <- r3 %>% filter(ind %% 100 == 0 ) #%in% sample(r3$ind[r3$time==0], 100)
(g2 <- ggplot(r3, aes(x=time, y=(mean_x), 
                       fill=mean_y, color=mean_y))+
  geom_line(aes(color = as.factor(ind)), size=0.5, stroke=0.5)+
  #xlim(0, 1100) +
  #facet_wrap(ind~.)+
  ggtitle("X position of hinge-blade frontier with time") +
  ylab("X position") +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) )
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")


g2 <- ggplot(r3, aes(x=time, y=(mean_x), 
                     fill=mean_y, color=mean_y))+
  geom_line(size=0.5, stroke=0.5)+
  #xlim(0, 1100) +
  facet_wrap(ind~.)+
  ggtitle("X position of hinge-blade frontier with time") +
  ylab("Y position") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) + 
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")


g3 <- ggplot(r3, aes(x=time, y=(mean_y), 
                     fill=mean_x, color=mean_x))+
  geom_point(size=0.5, stroke=0.5)+
  #xlim(0, 1100) +
  #facet_grid(coordinate~.)+
  ggtitle("Y position of hinge-blade frontier with time") +
  ylab("X position of border") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) + 
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")

g3 <- ggplot(r2, aes(x=time, y=(mean_y), 
                     fill=mean_x, color=mean_x))+
  geom_point(size=0.5, stroke=0.5)+
  #geom_line(aes(color=as.factor(ind)))+
  #xlim(0, 1100) +
  #facet_grid(coordinate~.)+
  ggtitle("Y position of hinge-blade frontier with time") +
  ylab("Y position of border") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) + 
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")


g4 <- ggplot(r3, aes(x=mean_x, y=mean_y, 
                     fill=time, color=time))+
  geom_point(size=0.5, stroke=0.5)+
  #xlim(0, 1100) +
  #facet_grid(coordinate~.)+
  ggtitle("X and Y position of hinge-blade frontier") +
  ylab("Y") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) + 
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")


g5 <- ggplot(r3, aes(x=mean_x, y=mean_y, 
                     fill=time, color=time))+
  geom_point(size=0.5, stroke=0.5)+
  #xlim(0, 1100) +
  #facet_grid(coordinate~.)+
  ggtitle("X and Y position of border") +
  ylab("Y") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) + 
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")


r2_reserva <- r2
r3_reserva <- r3
rmean <- r2_reserva %>% group_by(time) %>% summarise(mean_xpos=mean(mean_x),
                                             mean_ypos = mean(mean_y),
                                             min_xpos = min(mean_x),
                                             max_xpos = max(mean_x),
                                             min_ypos = min(mean_y),
                                             max_ypos = max(mean_y)) %>% as.data.frame()
for(t in rmean$time){
  r2 <- r2_reserva %>% filter(time==t)
  r3 <- r3_reserva %>% filter(time==t)
  anterior <- r2 %>% filter(mean_y > 0.5*(min(r2$mean_y)+max(r2$mean_y))  )
  posterior <- r2 %>% filter(mean_y < 0.5*(min(r2$mean_y)+max(r2$mean_y))  )
  
  rmean[rmean$time==t, "central_distal_ind"] <- r3$ind[which(abs(r3$mean_y - 0.5*(min(r3$mean_y)+max(r3$mean_y))) == min(abs(r3$mean_y - 0.5*(min(r3$mean_y)+max(r3$mean_y)))))]
  rmean[rmean$time==t, "central_proximal_ind"] <-  r2$ind[which(abs(r2$mean_y - 0.5*(min(r2$mean_y)+max(r2$mean_y))) == min(abs(r2$mean_y - 0.5*(min(r2$mean_y)+max(r2$mean_y)))))]
  rmean[rmean$time==t, "central_anterior_ind"] <- anterior$ind[which(abs(anterior$mean_x - 0.5*(min(anterior$mean_x)+max(anterior$mean_x))) == min(abs(anterior$mean_x - 0.5*(min(anterior$mean_x)+max(anterior$mean_x)))))]
  rmean[rmean$time==t, "central_posterior_ind"] <-  posterior$ind[which(abs(posterior$mean_x - 0.5*(min(posterior$mean_x)+max(posterior$mean_x))) == min(abs(posterior$mean_x - 0.5*(min(posterior$mean_x)+max(posterior$mean_x)))))]
  rmean[rmean$time==t, "distal_ypos"] <- r3$mean_y[r3$ind == rmean[rmean$time==t,]$central_distal_ind]
  rmean[rmean$time==t, "distal_xpos"] <- r3$mean_x[r3$ind == rmean[rmean$time==t,]$central_distal_ind]
  rmean[rmean$time==t, "proximal_ypos"] <- r2$mean_y[r2$ind == rmean[rmean$time==t,]$central_proximal_ind]
  rmean[rmean$time==t, "proximal_xpos"] <- r2$mean_x[r2$ind == rmean[rmean$time==t,]$central_proximal_ind]
  rmean[rmean$time==t, "anterior_ypos"] <- anterior$mean_y[anterior$ind == rmean[rmean$time==t,]$central_anterior_ind]
  rmean[rmean$time==t, "anterior_xpos"] <- anterior$mean_x[anterior$ind == rmean[rmean$time==t,]$central_anterior_ind]
  rmean[rmean$time==t, "posterior_ypos"] <- posterior$mean_y[posterior$ind == rmean[rmean$time==t,]$central_posterior_ind]
  rmean[rmean$time==t, "posterior_xpos"] <- posterior$mean_x[posterior$ind == rmean[rmean$time==t,]$central_posterior_ind]
  
  rmean[rmean$time==t, "AnteriorDistal_ypos"] <- anterior$mean_y[anterior$mean_x == max(anterior$mean_x)]
  rmean[rmean$time==t, "AnteriorDistal_xpos"] <- anterior$mean_x[anterior$mean_x == max(anterior$mean_x)]
  rmean[rmean$time==t, "AnteriorProximal_ypos"] <- anterior$mean_y[anterior$mean_x == min(anterior$mean_x)]
  rmean[rmean$time==t, "AnteriorProximal_xpos"] <- anterior$mean_x[anterior$mean_x == min(anterior$mean_x)]
  
  rmean[rmean$time==t, "PosteriorDistal_ypos"] <- posterior$mean_y[posterior$mean_x == max(posterior$mean_x)]
  rmean[rmean$time==t, "PosteriorDistal_xpos"] <- posterior$mean_x[posterior$mean_x == max(posterior$mean_x)]
  rmean[rmean$time==t, "PosteriorProximal_ypos"] <- posterior$mean_y[posterior$mean_x == min(posterior$mean_x)]
  rmean[rmean$time==t, "PosteriorProximal_xpos"] <- posterior$mean_x[posterior$mean_x == min(posterior$mean_x)]
                                           
}
rmean$central_length <- rmean$distal_xpos - rmean$proximal_xpos
rmean$central_height <- rmean$anterior_ypos - rmean$posterior_ypos
rmean$height <- rmean$central_height
rmean$central_HLratio <- rmean$height/rmean$central_length

rmean$anterior_length <- rmean$AnteriorDistal_xpos - rmean$AnteriorProximal_xpos
rmean$anterior_HLratio <- rmean$height/rmean$anterior_length

rmean$posterior_length <- rmean$PosteriorDistal_xpos - rmean$PosteriorProximal_xpos
rmean$posterior_HLratio <- rmean$height/rmean$posterior_length


plot(rmean$time, rmean$length)
plot(rmean$time, rmean$height)
plot(rmean$time, rmean$HLratio)

plot(rmean$time, rmean$anterior_length)
plot(rmean$time, rmean$anterior_HLratio)
plot(rmean$time, rmean$posterior_length)
plot(rmean$time, rmean$posterior_HLratio)

#Plot all ratios by compartment
allLengths <- rmean %>% gather("measure", "value", central_height, central_length, central_HLratio, 
                     anterior_length, anterior_HLratio, 
                     posterior_length, posterior_HLratio) %>% 
  separate(measure, c("region", "parameter"), sep="_")

ggplot(allLengths, aes(x=time, y=value, fill=region, color=region))+
  facet_wrap(parameter~., scales = "free", ncol=1) + 
  geom_point() + 
  geom_line() +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 

allpoints <- rmean %>% gather("point", "pos", distal_ypos, distal_xpos, proximal_xpos, proximal_ypos,
                              anterior_xpos, anterior_ypos, posterior_xpos, posterior_ypos,
                              AnteriorProximal_xpos, AnteriorProximal_ypos, AnteriorDistal_xpos, AnteriorDistal_ypos,
                              PosteriorDistal_xpos, PosteriorDistal_ypos, PosteriorProximal_xpos, PosteriorProximal_ypos) %>%  
                              #mean_xpos, mean_ypos, min_xpos, min_ypos, max_xpos, max_ypos) %>% 
  separate(point, c("point_type", "coord"), sep="_") %>% 
  spread(coord, pos)
  #mutate(coord = ifelse(grepl("xpos", point), "x", "y"),
  #       point = gsub("_[xy]pos", "", point, perl = T)) %>% 


ggplot(allpoints, aes(x=xpos, y=ypos, fill=point_type, col=point_type))+
  geom_point()+ #aes(fill=as.factor(time), col=as.factor(time))
  geom_line()+
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold")
    #strip.text = element_text(size=16, face="bold") 
  ) 
