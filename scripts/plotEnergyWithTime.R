library(tidyverse)
library(wesanderson)

setwd("/home/carmoma/vertex/Vertex-model/hexgrid2021_19/hexgrid2021_19_6/21grid5")

f <- list.files() %>% subset(grepl(".ptab", .))
x <- lapply(f, FUN=function(x){a<-read_tsv(x);a$time <- gsub("21grid5_moved_", "", x) %>% gsub("\\.ptab", "", .) %>% as.numeric;return(a)}) %>% bind_rows
x2 <- x%>%filter(spring==-999 & movable == 1)

##Get extra info
edges <- read_tsv("21grid5_moved_0.edges")
cells <- read_tsv("21grid5_moved_0.celltab")
points <- read_tsv("21grid5_moved_0.ptab")
pnei <- points %>% mutate(neighbour_vertices = gsub(",$", "", neighbour_vertices)) %>% 
  separate(neighbour_vertices, c("nei1", "nei2", "nei3"), ",") %>% 
  gather("nei_num", "nei_id", nei1, nei2, nei3)
  #Determine which edges are in the frontier
e2 <- edges %>% 
    mutate(vertices = gsub(",$", "", vertices)) %>% 
    separate(vertices, c("v1", "v2"), ",") %>%
    separate(cells, c("c1", "c2"), ",") %>%
    mutate(v1=as.numeric(v1), v2=as.numeric(v2))%>%
    gather("vert", "vert_id", v1, v2) %>% 
    select(-vert) %>% 
    mutate(celltype1 = cells$type[match(c1, cells$ind)],
           celltype2 = cells$type[match(c2, cells$ind)],
           newtype = ifelse( is.na(celltype1) | is.na(celltype2),"border" ,ifelse(celltype1 == celltype2, celltype1, "frontier")),
           xpos = points$x[match(vert_id, points$ind)],
           ypos = points$y[match(vert_id, points$ind)]
    )
e2$newtype[e2$newtype==0] <- "blade"
e2$newtype[e2$newtype==1] <- "hinge"
x2$type <- e2$newtype[match(x2$ind, e2$vert_id)]

frontinds <- e2 %>% filter(newtype=="frontier") %>% pull(vert_id) %>% unique
pnei <-pnei %>%  mutate(neighbour_frontier = ifelse(!(ind %in% frontinds) & nei_id %in% frontinds, T, F))
neighbors <- pnei %>% filter(neighbour_frontier) %>% pull(ind) %>% unique

x2$type[(x2$ind %in% neighbors) & x2$type == "hinge"] <- "hinge_frontier_neighbour"
x2$type[(x2$ind %in% neighbors) & x2$type == "blade"] <- "blade_frontier_neighbour"
x2 <- x2 %>% mutate(type = factor(type, levels = c("hinge", "hinge_frontier_neighbour", "frontier", "blade_frontier_neighbour", "blade", "border")))

x2 %>% filter(time == 0) %>%ggplot(aes(x=x, y=y, fill=type, color=type)) + geom_point()

x2 %>% ggplot(aes(x=x, y=y, fill=energy, color=energy))+
  geom_point(size=0.2)+
  facet_wrap(time~.) +
  scale_color_gradient(low="blue", high="red") +
scale_fill_gradient(low="blue", high="red") +
  ggtitle("Energy at each timestep") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

x2$ratio_accepted <- x2$moves_accepted/(x2$moves_accepted + x2$moves_rejected)
x2$ratio_accepted[x2$time==0]<-NA
x2 %>% ggplot(aes(x=x, y=y, fill=ratio_accepted, color=ratio_accepted))+
  geom_point(size=0.2)+
  facet_wrap(time~.) +
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red") +
  ggtitle("Proportion of moves accepted") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )


x2 %>% filter(time==0) %>%  ggplot(aes(x=x, y=y, fill=energy, color=energy))+
  geom_point()+ #size=0.2
  #facet_wrap(time~.) +
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

x2 %>% filter(time==1) %>%  ggplot(aes(x=x, y=y, fill=ratio_accepted, color=ratio_accepted))+
  geom_point()+ #size=0.2
  #facet_wrap(time~.) +
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red") +
  ggtitle("Proportion accepted moves at time = 1") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

x2 %>% filter(time %in% c(0,1,10,20,40)) %>% 
  mutate(time = factor(time, levels = c(0, 1, 10, 20, 40))) %>% 
  ggplot(aes(x=type, y = energy, fill=time)) + 
  #facet_grid(time~.)+
  #geom_violin() + 
  geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )


x2 %>% filter(time %in% c(0,1,10,20,40)) %>% 
  mutate(time = factor(time, levels = c(0, 1, 10, 20, 40))) %>% 
  ggplot(aes(x=type, y = ratio_accepted, fill=time)) + 
  #facet_grid(time~.)+
  #geom_violin() + 
  ylab("Proportion of accepted movements per vertex") +
  geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )


x2 %>% filter(time %in% c(0,1,10,20,40)) %>% 
  mutate(time = factor(time, levels = c(0, 1, 10, 20, 40))) %>% 
  ggplot(aes(x=x, y = ratio_accepted, fill=time, col=time)) + 
  geom_smooth() +
  geom_point() +
  #facet_grid(time~.)+
  #geom_violin() + 
  ylab("Proportion of accepted movements per vertex") +
  #geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

x2 %>% filter(time %in% c(0,1,10,20,40)) %>% 
  mutate(time = factor(time, levels = c(0, 1, 10, 20, 40))) %>% 
  ggplot(aes(x=x, y = energy, fill=time, col=time)) + 
  geom_smooth() +
  geom_point() +
  #facet_grid(time~.)+
  #geom_violin() + 
  ylab("Energy per vertex") +
  ggtitle("Energy along the X axis, colored by time") +
  #geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

x2 %>% 
  ggplot(aes(x=time, y = energy, fill=type, col=type)) + 
  geom_smooth() +
  geom_point(size = 0.2) +
  #facet_grid(time~.)+
  #geom_violin() + 
  ylab("Energy per vertex") +
  ggtitle("Energy with time, colored by type") +
  #geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

x2 %>% 
  ggplot(aes(x=time, y = ratio_accepted, fill=type, col=type)) + 
  geom_smooth() +
  geom_point(size = 0.2) +
  #facet_grid(time~.)+
  #geom_violin() + 
  ylab("Energy per vertex") +
  ggtitle("Proportion of accepted movements with time, colored by type") +
  #geom_boxplot(notch = T) +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    axis.text = element_text( size=10, face="bold")
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

x2 %>% mutate(time = factor(time)) %>% 
  ggplot(aes(x=time, y = energy)) + 
  #geom_violin() + 
  geom_boxplot() +
  theme_minimal()
x3 <- x2 %>% group_by(time) %>% summarise(energy = sum(energy))
plot(x3$time, x3$energy)
x3 %>% ggplot(aes(x=time, y = energy))+geom_point() +geom_line()+ theme_bw()

#########
setwd("/home/carmoma/vertex/Vertex-model/hexgrid2021_19/hexgrid2021_19_6/21grid5")

f <- list.files() %>% subset(grepl(".ptab", .))
x <- lapply(f, FUN=function(x){a<-read_tsv(x);a$time <- gsub("21grid5_moved_", "", x) %>% gsub("\\.ptab", "", .) %>% as.numeric;return(a)}) %>% bind_rows
x2 <- x%>%filter(spring==-999 & movable == 1)

x2 %>% ggplot(aes(x=x, y=y, fill=energy, color=energy))+
  geom_point(size=0.2)+
  facet_wrap(time~.) +
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )




#### Example of derivatives
setwd("/home/carmoma/vertex/Vertex-model/hexgrid2021_21/hexgrid2021_21_0")
files <- list.files() %>% subset(grepl("energy",.))
time <- files %>% gsub("energy_|\\.tsv", "", .) %>% as.numeric
##Assume that, since there are no T1s, divisions or T2s, this will be the same
edges <- read_tsv("21grid5/21grid5_moved_0.edges")
cells <- read_tsv("21grid5/21grid5_moved_0.celltab")
points <- read_tsv("21grid5/21grid5_moved_0.ptab")
pnei <- points %>% mutate(neighbour_vertices = gsub(",$", "", neighbour_vertices)) %>% 
  separate(neighbour_vertices, c("nei1", "nei2", "nei3"), ",") %>% 
  gather("nei_num", "nei_id", nei1, nei2, nei3)

for(i in length(files):1){
en <- read.table(files[i], header=F, stringsAsFactors = F)
basename <- paste("time", as.character(time[i]), sep="", collapse="")

#names(en)[5:6]<- c("y", "x") ##CAREFUL! FIXED IN VERTEX MODEL
nnames <- en[1, 1] %>% gsub(">", "", .) %>% strsplit(., ":" )
nnames <- c("name", nnames[[1]])
nnames <- c(nnames[1:11], "prob_accept", nnames[12])
names(en)<-nnames
en$favorable <- ifelse(en$prob_accept == 0.05, "unfavourable", "favourable")

#Determine which edges are in the frontier
e2 <- edges %>% 
  mutate(vertices = gsub(",$", "", vertices)) %>% 
  separate(vertices, c("v1", "v2"), ",") %>%
  separate(cells, c("c1", "c2"), ",") %>%
  mutate(v1=as.numeric(v1), v2=as.numeric(v2))%>%
  gather("vert", "vert_id", v1, v2) %>% 
   select(-vert) %>% 
  mutate(celltype1 = cells$type[match(c1, cells$ind)],
        celltype2 = cells$type[match(c2, cells$ind)],
        newtype = ifelse( is.na(celltype1) | is.na(celltype2),"border" ,ifelse(celltype1 == celltype2, celltype1, "frontier")),
        xpos = points$x[match(vert_id, points$ind)],
        ypos = points$y[match(vert_id, points$ind)]
       )
 #Plot edges by type just to check data make sense
#e2 %>% ggplot(aes(x=xpos, y=ypos, fill=newtype, col=newtype)) + 
#  geom_point() +
#  scale_x_continuous(breaks=(seq(0, 120, 5)))

en$type <- e2$type[match(en$id, e2$vert_id)]
en$newtype <- e2$newtype[match(en$id, e2$vert_id)]
en$newtype[en$newtype==0] <- "blade"
en$newtype[en$newtype==1] <- "hinge"
enred <- en %>% filter(x >57 & x < 67 & newtype != "border")

#Plot arrows
try({
ensum<- en %>% mutate(
  angle = atan2(ydif, xdif),
  x_component = -1*cos(angle)*edif,
  y_component = -1*sin(angle)*edif,
  new_angle = atan2(y_component, x_component)
) %>% group_by(id) %>% summarise(x=mean(x), y=mean(y), varx=var(x), vary=var(y),
                                 sum_vectorx = sum(x_component),
                                 sum_vectory = sum(y_component),
                                 mean_vectorx = mean(x_component),
                                 mean_vectory = mean(y_component),
                                 new_angle  =mean(new_angle),
                                 type = unique(newtype)) 
  ensum <- ensum %>% mutate(vectorend_x1 = x + 0.21*sum_vectorx,
         vectorend_y1 = y + 0.2*sum_vectory,
         vectorend_x2 = x + mean_vectorx,
         vectorend_y2 = y + mean_vectory)

g_vectors <- ensum %>% filter(type != "border") %>% 
  ggplot(aes(x=x, y=y, 
             xend=vectorend_x1, yend = vectorend_y1, fill=type, col=type))+
  geom_point()+
  geom_segment(arrow = arrow(length = unit(0.15,"cm"))) +
  ggtitle("Vertex position and direction calculated from E of move trials") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
fname <- paste(basename, "vectors.pdf", sep="_", collapse="_")
ggsave(fname, g_vectors)
fname <- paste(basename, "vectors.png", sep="_", collapse="_")
ggsave(fname, g_vectors)
cat(fname, "\n")
})
#x vs energy difference
try({
g_xdif1 <- en %>% filter(newtype != "border") %>% 
  ggplot(aes(x=x, y=edif, fill=xdif, col=xdif))+
  geom_point(size = 0.3) +
  xlab("Vertex position") +
  ylab("Increase in energy of move trials") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
g_xdif2 <-en %>% filter(newtype != "border") %>% 
  ggplot(aes(x=x, y=edif, fill=newtype, col=newtype))+
  geom_point(size = 0.3) +
  xlab("Vertex position") +
  ylab("Increase in energy of move trials") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "x-vs-Edif1.pdf", sep="_", collapse="_")
ggsave(fname, g_xdif1)
fname <- paste(basename, "x-vs-Edif2.pdf", sep="_", collapse="_")
ggsave(fname, g_xdif2)
fname <- paste(basename, "x-vs-Edif1.png", sep="_", collapse="_")
ggsave(fname, g_xdif1)
fname <- paste(basename, "x-vs-Edif2.png", sep="_", collapse="_")
ggsave(fname, g_xdif2)
cat(fname, "\n")
})
#vv <- en %>% filter(id==693)
v2select<- c(438, 396, 395,439, 440, 442, 441, 397, 398, 443, 399, 445)
vtypes <-c("hinge","hinge", "hinge", "hinge", "hinge", "frontier", "frontier", "frontier", "frontier", "blade", "blade", "blade")
try({
vv <- en %>% filter(id%in% v2select)
vv$type <- vtypes[match(vv$id, v2select)]
vv$name <- paste(as.character(vv$id), as.character(vv$type), sep=" ")
vv <- vv %>% mutate(type = factor(type, levels = c("hinge", "frontier", "blade")))
vv <- vv[order(vv$type) ,]
g_vsel1 <- vv %>% ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif))+
  geom_point() +
  facet_wrap(.~name)+
  scale_color_gradient2(low="blue", high="red", mid="yellow") +
  scale_fill_gradient2(low="blue", high="red", mid="yellow") +
  ggtitle("Energy landscape among selected points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleSelected.pdf", sep="_", collapse="_")
ggsave(fname, g_vsel1)
fname <- paste(basename, "EnergyCircleSelected.png", sep="_", collapse="_")
ggsave(fname, g_vsel1)
cat(fname, "\n")
})
#enred2 <- enred %>% filter(newtype != "border") %>% slice_sample(n=10000)
#minf <- enred %>% filter(newtype=="frontier") %>% pull(edif) %>% quantile(0.1)
#maxf <- enred %>% filter(newtype=="frontier") %>% pull(edif) %>% quantile(0.9)
try({
g_circleregion1 <- enred %>%ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  facet_wrap(.~newtype)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape among selected points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleCentralRegion.pdf", sep="_", collapse="_")
ggsave(fname, g_circleregion1)
fname <- paste(basename, "EnergyCircleCentralRegion.png", sep="_", collapse="_")
ggsave(fname, g_circleregion1)
cat(fname, "\n")
})
##Circles in frontier. First all, then filtering out the ones with high values
try({
enred_fron <- enred %>%filter(newtype == "frontier") 
})
try({
g_frontier1 <- enred_fron %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
g_frontier2 <- enred_fron %>% 
  filter(id != 45 & id != 861) %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleFrontier1.pdf", sep="_", collapse="_")
ggsave(fname, g_frontier1)
fname <- paste(basename, "EnergyCircleFrontier2.pdf", sep="_", collapse="_")
ggsave(fname, g_frontier2)
fname <- paste(basename, "EnergyCircleFrontier1.png", sep="_", collapse="_")
ggsave(fname, g_frontier1)
fname <- paste(basename, "EnergyCircleFrontier2.png", sep="_", collapse="_")
ggsave(fname, g_frontier2)
cat(fname, "\n")
})
##Circle of neighbours
try({
frontinds <- enred_fron %>% pull(id) %>% unique
pnei <-pnei %>%  mutate(neighbour_frontier = ifelse(!(ind %in% frontinds) & nei_id %in% frontinds, T, F))
neighbors <- pnei %>% filter(neighbour_frontier) %>% pull(ind) %>% unique
enred_nei <- enred %>% filter(id %in% neighbors)
})
try({
g_neighbors1 <- enred_nei %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in neighbours to frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

g_neighbors2 <- enred_nei %>% filter(! (id %in% c(228, 776, 188, 816))) %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in neighbours to frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

g_neighbors3 <- enred_nei %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  #facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in neighbours to frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

g_neighbors4 <- enred_nei %>% filter(! (id %in% c(228, 776, 188, 816))) %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.8))+
  geom_point() +
  #facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") + #,limits=c(minf, maxf) 
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in neighbours to frontier points") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleFrontNeighbors1.pdf", sep="_", collapse="_")
ggsave(fname, g_neighbors1)
fname <- paste(basename, "EnergyCircleFrontNeighbors2.pdf", sep="_", collapse="_")
ggsave(fname, g_neighbors2)
fname <- paste(basename, "EnergyCircleFrontNeighbors1.png", sep="_", collapse="_")
ggsave(fname, g_neighbors1)
fname <- paste(basename, "EnergyCircleFrontNeighbors2.png", sep="_", collapse="_")
ggsave(fname, g_neighbors2)
fname <- paste(basename, "EnergyCircleFrontNeighbors3.pdf", sep="_", collapse="_")
ggsave(fname, g_neighbors3)
fname <- paste(basename, "EnergyCircleFrontNeighbors4.pdf", sep="_", collapse="_")
ggsave(fname, g_neighbors4)
fname <- paste(basename, "EnergyCircleFrontNeighbors3.png", sep="_", collapse="_")
ggsave(fname, g_neighbors3)
fname <- paste(basename, "EnergyCircleFrontNeighbors4.png", sep="_", collapse="_")
ggsave(fname, g_neighbors4)
cat(fname, "\n")
cat(fname, "\n")
})
#Classify all vertices and pull derivatives
try({
g_x_xdif <- enred %>% ggplot(aes(x=x, y=xdif,fill=edif, col=edif, alpha=0.5)) + 
  geom_point(size=0.2) +
  scale_color_gradient2(low="blue", high="red", mid="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("X vs movement in X, colors according to Energy difference") +
  xlab("x position") + 
  ylab("Movement in x") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "XvsXdiff.pdf", sep="_", collapse="_")
ggsave(fname, g_x_xdif)
fname <- paste(basename, "XvsXdiff.png", sep="_", collapse="_")
ggsave(fname, g_x_xdif)
cat(fname, "\n")
})
try({
enfron <- en %>% filter(newtype=="frontier")
})
try({
g_frontierTogether1 <- enfron %>%ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.1))+
  geom_point() +
  #facet_wrap(.~newtype)+
  scale_color_gradient2(low="blue", high="red", mid="white", limits = c(-0.015, 0.015)) +
  scale_fill_gradient2(low="blue", high="red", mid="white", limits = c(-0.015, 0.015)) +
  ggtitle("Energy landscape among points in frontier") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

g_frontierTogether2 <- enfron %>%
  filter(id != 45 & id != 861) %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.1))+
  geom_point() +
  #facet_wrap(.~newtype)+
  scale_color_gradient2(low="blue", high="red", mid="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape among points in frontier (removed outliers)") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

fname <- paste(basename, "EnergyCircleFrontier3.pdf", sep="_", collapse="_")
ggsave(fname, g_frontierTogether1)
fname <- paste(basename, "EnergyCircleFrontier4.pdf", sep="_", collapse="_")
ggsave(fname, g_frontierTogether2)
fname <- paste(basename, "EnergyCircleFrontier3.png", sep="_", collapse="_")
ggsave(fname, g_frontierTogether1)
fname <- paste(basename, "EnergyCircleFrontier4.png", sep="_", collapse="_")
ggsave(fname, g_frontierTogether2)
cat(fname, "\n")
})
try({
en$bin <- as.integer(en$x/10) %>% as.factor()
en$biny <- as.integer(en$y/10) %>% as.factor()

g_Xbins <- en %>% filter(newtype!="border") %>% #slice_sample(n=10000) %>% 
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.1))+
  geom_point(size=0.2) +
  facet_wrap(.~bin)+
  scale_color_gradient2(low="blue", high="red", mid="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in 10 binned regions across the X axis") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

g_XYbins <- en %>% filter(newtype!="border") %>%
  ggplot(aes(x=xdif, y=ydif, fill=edif, col=edif, alpha=0.1))+
  geom_point(size=0.1) +
  facet_grid(biny~bin)+
  scale_color_gradient2(low="blue", high="red", mid="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape in 10x10 binned regions across XY axis") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleXbins.pdf", sep="_", collapse="_")
ggsave(fname, g_Xbins)
fname <- paste(basename, "EnergyCircleXYbins.pdf", sep="_", collapse="_")
ggsave(fname, g_XYbins)
fname <- paste(basename, "EnergyCircleXbins.png", sep="_", collapse="_")
ggsave(fname, g_Xbins)
fname <- paste(basename, "EnergyCircleXYbins.png", sep="_", collapse="_")
ggsave(fname, g_XYbins)
cat(fname, "\n")
})

try({
#energy difference for each point
en$name2 <- paste(en$id, en$newtype, sep=" ")

#option1 to localize points in space (no grid because bad result)
#en$extrax <- en$xdif + 0.2*as.integer(en$id/30)
#en$extray <- en$ydif + 0.2*as.integer(en$id)%%28
#option 2
en$extrax <- en$xdif + 6*(en$x - min(en$x))/(max(en$x) - min(en$x))
en$extray <- en$ydif + 6*(en$y - min(en$y))/(max(en$y) - min(en$y))

g_circleByPoint <- en %>% filter(newtype!="border") %>%
  ggplot(aes(x=extrax, y=extray, fill=edif, col=edif, alpha=0.2))+
  geom_point(size=0.1) +
  #facet_wrap(.~id)+
  scale_color_gradient2(low="blue", high="red", mid="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white") +
  ggtitle("Energy landscape for all points at scaled locations") +
  xlab("x movement") + 
  ylab("y movement") +
  theme_minimal() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
fname <- paste(basename, "EnergyCircleByPoint.pdf", sep="_", collapse="_")
ggsave(fname, g_circleByPoint)
fname <- paste(basename, "EnergyCircleByPoint.png", sep="_", collapse="_")
ggsave(fname, g_circleByPoint)
cat(fname, "\n")
})
cat("****************\n")
}





######### Analyze a particular timestep in more detail

setwd("/home/carmoma/vertex/Vertex-model/hexgrid2021_8/hexgrid2021_8_0/21grid5")

f <- list.files() %>% subset(grepl("_30\\.", .))
edges <- read_tsv(f[grepl("edges", f)])
cells <- read_tsv(f[grepl("celltab", f)])
points <- read_tsv(f[grepl("ptab", f)])
v2select<- c(438, 396, 395,439, 440, 442, 441, 397, 398, 443, 399, 445)
vtypes <-c("hinge","hinge", "hinge", "hinge", "hinge", "frontier", "frontier", "frontier", "frontier", "blade", "blade", "blade")

e2 <- edges %>% 
  mutate(vertices = gsub(",$", "", vertices)) %>% 
  separate(vertices, c("v1", "v2"), ",") %>%
  separate(cells, c("c1", "c2"), ",") %>%
  mutate(v1=as.numeric(v1), v2=as.numeric(v2))%>%
  gather("vert", "vert_id", v1, v2) %>% 
  select(-vert) %>% 
  mutate(celltype1 = cells$type[match(c1, cells$ind)],
         celltype2 = cells$type[match(c2, cells$ind)],
         newtype = ifelse( is.na(celltype1) | is.na(celltype2),"border" ,ifelse(celltype1 == celltype2, celltype1, "frontier")),
         xpos = points$x[match(vert_id, points$ind)],
         ypos = points$y[match(vert_id, points$ind)]
  )
e2$newtype[e2$newtype==0] <- "blade"
e2$newtype[e2$newtype==1] <- "hinge"



e2$energy <- e2$tension*e2$length
edges$energy <- edges$tension*edges$length


e441 <- points %>% filter(ind == 441) %>% pull(edges) %>% gsub(",$", "", .) %>% strsplit(",") %>% unlist
e440 <- points %>% filter(ind == 440) %>% pull(edges) %>% gsub(",$", "", .) %>% strsplit(",") %>% unlist
c441 <- points %>% filter(ind == 441) %>% pull(cells) %>% gsub(",$", "", .) %>% strsplit(",") %>% unlist
c440 <- points %>% filter(ind == 440) %>% pull(cells) %>% gsub(",$", "", .) %>% strsplit(",") %>% unlist

e3 <- e2 %>% filter(ind %in% c(e440, e441))
eaux <- edges %>% filter(ind %in% c(e440, e441))
caux <- cells %>% filter(ind %in% c(c440, c441))

term_coeffs <- c(0.1, 2, 1)
K <- caux$K[1]
C <- caux$perim_contract[1]
A0 <- caux$preferred_area[1]
#### After new simulations
setwd("~/vertex/Vertex-model/hexgrid2021_21/hexgrid2021_21_0")
f <- list.files() %>% subset(grepl("detail", .))
time <- c(0, 10, 20, 30, 50)
a <- data.frame()
for(i in 1:length(f)){
  aux <-read_tsv(f[i])
  aux$time2 <- time[i]
  a <- rbind(a, aux)
}
#a <- read_tsv("mv20_detail440-441.tsv")

a <- a %>% mutate(
  old_term2 = 2*(old_elen0*old_eten0 + old_elen1*old_eten1 + old_elen2*old_eten2),
  new_term2 = 2*(new_elen0*new_eten0 + new_elen1*new_eten1 + new_elen2*new_eten2),
  old_term1 = 0.1*(K*(old_cellarea0 - A0)^2 + K*(old_cellarea1 - A0) + K*(old_cellarea2 - A0)),
  new_term1 = 0.1*(K*(new_cellarea0 - A0)^2 + K*(new_cellarea1 - A0) + K*(new_cellarea2 - A0)),
  old_term3 = C*old_cellper0^2 + C*old_cellper1^2 + C*old_cellper2^2, 
  new_term3 = C*new_cellper0^2 + C*new_cellper1^2 + C*new_cellper2^2, 
  old_energy = old_term1 + old_term2 + old_term3,
  new_energy = new_term1 + new_term2 + new_term3,
  old_term1a = 0.1*K*(old_cellarea0 - A0)^2,
  old_term1b = 0.1*K*(old_cellarea1 - A0)^2,
  old_term1c = 0.1*K*(old_cellarea2 - A0)^2,
  old_term2a = 2*old_elen0*old_eten0,
  old_term2b = 2*old_elen1*old_eten1,
  old_term2c = 2*old_elen2*old_eten2,
  old_term3a = C*old_cellper0^2,
  old_term3b = C*old_cellper1^2,
  old_term3c = C*old_cellper2^2,
  term1_inc = old_term1 - new_term1,
  term2_inc = old_term2 - new_term2,
  term3_inc = old_term3 - new_term3,
  term1plus3_inc = term1_inc + term3_inc,
  energy_inc = old_energy - new_energy,
  favourable = factor(ifelse(energy_inc < 0, "accepted", "rejected")),
  energy_error = energy_inc - edif,
  time2 = factor(time2, levels = c(0, 10, 20, 30, 50)),
  edge630 = factor(ifelse(id == 440 & (new_elen0 > old_elen0), "gets longer", "gets shorter")),
  edge629 = factor(ifelse(id == 441 & (new_elen0 > old_elen0), "gets longer", "gets shorter")),
  edge631 = factor(ifelse(id == 441 & (new_elen1 > old_elen1), "gets longer", "gets shorter")),
  cell190 = factor(ifelse(new_cellarea0 < old_cellarea0, "gets smaller", "gets bigger")),
  cell191 = factor(ifelse(id == 441, ifelse(new_cellarea1 < old_cellarea1, "gets smaller", "gets bigger"),NA))
)  

b <- a %>% filter(time2==20)

g<-a %>% ggplot(aes(x = term2_inc, y = term3_inc, fill=favourable, col = favourable))+
  facet_grid(id~time2)+
  geom_point()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 3 (perimeter term)") +
  ggtitle("Relationship between differences in term 2 and 3 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )+
  scale_fill_discrete(wes_palette("Zissou1", 2, type = "discrete")) +
  scale_color_discrete(wes_palette("Zissou1", 2, type = "discrete"))

g<-a %>% ggplot(aes(x = term2_inc, y = term3_inc, fill=cell191:favourable, col = cell191:favourable))+
  facet_grid(id~time2)+
  geom_point(alpha = 0.5, size=0.2)+
  geom_rug() +
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 3 (perimeter term)") +
  ggtitle("Relationship between differences in term 2 and 3 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
ggsave("cell191getsLonger.pdf", g, width = 12, height = 7)
g<-a %>% ggplot(aes(x = term2_inc, y = term3_inc, fill=cell190:favourable, col = cell190:favourable))+
  facet_grid(id~time2)+
  geom_point(alpha = 0.5, size=0.2)+
  geom_rug() +
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 3 (perimeter term)") +
  ggtitle("Relationship between differences in term 2 and 3 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
ggsave("cell190getsLonger.pdf", g, width = 12, height = 7)
g<-a %>% ggplot(aes(x = cell190:favourable, fill=cell190:favourable, col = cell190:favourable))+
  facet_grid(id~time2)+
  #geom_point(alpha = 0.5, size=0.2)+
  #geom_rug() +
  geom_bar()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 3 (perimeter term)") +
  ggtitle("Frequencies of movement trials - cell 190") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
ggsave("cell190getsLonger_barplot.pdf", g, width=12, height = 10)

g<-a %>% filter(id==441)%>% ggplot(aes(x = cell191:favourable, fill=cell191:favourable, col = cell191:favourable))+
  facet_grid(id~time2)+
  #geom_point(alpha = 0.5, size=0.2)+
  #geom_rug() +
  geom_bar()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 3 (perimeter term)") +
  ggtitle("Frequencies of movement trials - cell 191") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )
ggsave("cell191getsLonger_barplot.pdf", g, width=12, height = 6)

b %>% ggplot(aes(x = term2_inc, y = term1_inc, fill=xdif, col = xdif))+
  facet_grid(id~.)+
  geom_point()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 1 (area term)") +
  
  ggtitle("Relationship between differences in term 2 and 1 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    legend.title = element_text("Movement in x"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

b %>% ggplot(aes(x = term2_inc, y = term1plus3_inc, fill=xdif, col = xdif))+
  facet_grid(id~.)+
  geom_point()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 1 + 3 (area + perim)") +
  
  ggtitle("Relationship between differences in term 2 and 1 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    legend.title = element_text("Movement in x"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

b %>% ggplot(aes(x = new_cellper0, y = new_cellper2, fill=xdif, col = xdif))+
  facet_grid(id~.)+
  geom_point()+
  xlab("Increment in term 2 (tension term)") + 
  ylab("Increment in term 1 + 3 (area + perim)") +
  
  ggtitle("Relationship between differences in term 2 and 1 in movement trials") +
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    legend.title = element_text("Movement in x"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

asum <- a %>% group_by(time2, id) %>% summarise(term1 = mean(old_term1),
                                                term2 = mean(old_term2),
                                                term3 = mean(old_term3),
                                                term1a = mean(old_term1a),
                                                term1b = mean(old_term1b),
                                                term1c = mean(old_term1c),
                                                term2a = mean(old_term2a),
                                                term2b = mean(old_term2b),
                                                term2c = mean(old_term2c),
                                                term3a = mean(old_term3a),
                                                term3b = mean(old_term3b),
                                                term3c = mean(old_term3c),
                                                energy = mean(old_energy)) %>% 
  gather("term", "value", term1, term2, term3, term1a, term1b, term1c, term2a, term2b, term2c, term3a, term3b, term3c, energy) %>% 
  rename(time=time2) %>% 
  mutate(term_type = gsub("a|b|c", "", term, perl = T))


g<-asum %>% ggplot(aes(x=time, y=value, fill=term, col=term))+
  facet_grid(id~.) +
  geom_line() + 
  geom_point(aes(shape=term_type))+
  theme_bw() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    legend.title = element_text("Movement in x"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

ggsave("LineEnergyTermsOld.pdf", g)

c <- a %>% gather(term, value, matches("inc")) %>% 
  select(id, term, value, time2) %>% rename(time=time2) %>% 
  filter(term != "term1plus3_inc")
g <- c %>% ggplot(aes(x=term, y=value, fill=term))+
  facet_grid(id~time) +
  geom_violin(scale="width")+
  theme_bw() +
  #geom_jitter() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.text.x = element_text( size=10, face="bold", angle=45, vjust = 0.5),
    axis.title.y = element_text(size=16, face="bold"),
    legend.title = element_text("Movement in x"),
    strip.text = element_text(size=16, face="bold") ,
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank()
  )

ggsave("violinPlotEnergyTerms.pdf", g)

