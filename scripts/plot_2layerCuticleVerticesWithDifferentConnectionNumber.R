library(tidyverse)


f<-list.files()%>%subset(grepl("_0.", .)) #files previous to any movement
e <- read_tsv(f[3]) #.edges file
p <- read_tsv(f[6]) #.ptab file
 
#Vertices in cuticle that are present in less than 5 vertices
vertcount <- e %>%filter(type==7) %>% 
	mutate(vertices = gsub(",$", "", vertices)) %>% 
	separate(vertices, c("v1", "v2"), sep=",") %>% 
	gather("vpos", vertex, v1, v2)%>%
	group_by(vertex)%>%
	summarise(n = n()) %>% filter(n != 5)

#Cuticle vertices with a factor specifying number of connections
vcut <- p %>% filter(type == 3) %>% 
	mutate(num_edges = factor(ifelse(ind %in% vertcount$vertex, vertcount$n[match(ind, vertcount$vertex)], 5)), 
		size=ifelse(num_edges ==5, 0.5, 1))

ggplot(vcut, aes(x=x, y=y, col=num_edges, size=size))+geom_point()

