
library(tidyverse)
library(wesanderson)

nums <- 268:265
wing <- "/etournay1_strings8b"
sims <- 0:5
timesteps <- c(0, 1, 5, 10, 20, 40)
dirbase <- "/home/carmoma/vertex/Vertex-model/dpygrad_mode%NUMBER%/dpygrad_mode%NUMBER%_"


for(num in nums){
  dirsims <- dirbase %>% gsub("%NUMBER%", as.character(num), .)
  cat(dirsims, "******")
  names <- paste(as.character(num), as.character(sims), sep="_")
  dirs <- paste(dirsims, as.character(sims), wing, sep="")
  plotAll(dirs) #Defined below
}
# names <- c("Only temporal gradient", "temporal + PD gradient", "temporal + AP gradient", 
#            "temporal + PD + AP gradients", "temporal + PD Isaac's way")
# #names <- c("No area, no perim grad", 
#           "No area, perim X grad", 
#           "No area, perim Y grad", 
#            "No area, perim XY grad", 
#            "Area X grad, no perim grad", 
#            "Area X grad, perim X grad", 
#            "Area X grad, perim Y grad", 
#            "Area X grad, perim XY grad", 
#            "Area Y grad, no perim grad", 
#            "Area Y grad, perim X grad", 
#            "Area Y grad, perim Y grad", 
#            "Area Y grad, perim XY grad", 
#            "Area XY grad, no perim grad", 
#            "Area XY grad, perim X grad", 
#            "Area XY grad, perim Y grad", 
#            "Area XY grad, perim XY grad")


plotAll <- function(dirs){
for(d in dirs){
  
  setwd(d)
  f<-list.files() %>%subset(grepl(".celltab", .)) %>% subset(grepl("moved",.))
  fpoints<-list.files() %>%subset(grepl(".ptab", .)) %>% subset(grepl("moved",.))
  time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub(".celltab", "", .)
  ind2read <- which(time %in% as.character(timesteps))
  time <- time[ind2read]
  f <- f[ind2read]
  res<-data.frame();
  for(i in 1:length(f)){

    aux<-read.table(f[i], head=T, sep="\t", dec=".", stringsAsFactors = F) %>% tibble() 
    #To update centroid positions 
    #points <- read.table(fpoints[i], head=T, sep="\t", dec=".", stringsAsFactors = F) %>% tibble()
    #aux <- aux %>% mutate(vertices = sapply(vertices, function(x)as.integer(strsplit(x, ",")[[1]])),
    #          current_x = sapply(vertices, function(indicesthiscell){
    #            points %>% filter(ind %in% indicesthiscell) %>% select(x) %>% pull %>% mean
    #            }),
    #         current_y = sapply(vertices, function(indicesthiscell){
    #           points %>% filter(ind %in% indicesthiscell) %>% select(y) %>% pull %>% mean
    # }))
    

    
   aux$time <- time[i]
   res<-rbind(res, aux)
  };
  #res$time<-as.character(res$time)
  res<-res %>% mutate(type = recode(type+1, "blade", "hinge", "vein blade", "vein hinge"))
  res2 <- res %>%
    filter(time %in% timesteps & preferred_area > 0 & type %in% c("blade", "hinge") ) %>% 
    mutate(preferred_area_norm = preferred_area*(num_divisions+1))
  res3 <- res2 %>% gather("coordinate", "centroid_position", centroid_x, centroid_y)
  res3$time <- as.numeric(res3$time)
  #b <- res2 %>% filter(type == "hinge")
  g0a <- ggplot(res2, aes(x=centroid_x, y=centroid_y, fill=perim_contract, col=perim_contract)) + 
    geom_point() + 
    facet_wrap(factor(time, levels = timesteps)~.) +
    xlim(0, 1100) +
    ylim(0, 450) +
    ggtitle("perim_contract")
  g0b <- ggplot(res2, aes(x=centroid_x, y=centroid_y, fill=preferred_area_norm, col=preferred_area_norm)) + 
    geom_point() + 
    facet_wrap(factor(time, levels = timesteps)~.)+
    xlim(0, 1100) +
    ylim(0, 450) +
    ggtitle("preferred_area")
  g0c <- ggplot(res2, aes(x=centroid_x, y=centroid_y, fill=base_eq_perimcontr, col=base_eq_perimcontr)) + 
    geom_point() + 
    facet_wrap(factor(time, levels = timesteps)~.)+
    xlim(0, 1100) +
    ylim(0, 450) +
    ggtitle("base_perim_contract")
  g0d <- ggplot(res2, aes(x=centroid_x, y=centroid_y, fill=base_eq_area, col=base_eq_area)) + 
    geom_point() + 
    facet_wrap(factor(time, levels = timesteps)~.)+
    xlim(0, 1100) +
    ylim(0, 450) +
    ggtitle("base_preferred_area")
  g1 <- ggplot(res3, aes(x=centroid_position, y=(preferred_area*(num_divisions+1)), 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    xlim(0, 1100) +
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("Equilibrium Area") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) 

  g2 <- ggplot(res3, aes(x=centroid_position, y=(perim_contract), 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    xlim(0, 1100) +
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("Perimeter contractility") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) + 
    scale_color_gradient(low="blue", high="red") +
    scale_fill_gradient(low="blue", high="red")
  g3 <- ggplot(res3, aes(x=centroid_position, y=K, 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5, alpha=0.5)+
    xlim(0, 1100) +
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("K") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) + 
    scale_color_gradient(low="green", high="purple") +
    scale_fill_gradient(low="green", high="purple")
  #scale_fill_grey()
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest"))
  
  res3$ratio_area <- res3$area / res3$preferred_area
  g4 <- ggplot(res3, aes(x=centroid_position, y=ratio_area, 
                 fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    geom_hline(yintercept=1, linetype = 2) +
    xlim(0, 1100) +
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    xlab("Centroid position") +
    ylab("Area / Eq. Area") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) 

  setwd("../../")
  pdf(paste(gsub(" ", "_", names[which(dirs == d)]), ".pdf", sep="", collapse = ""))
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  print(g0a)
  print(g0b)
  print(g0c)
  print(g0d)
  dev.off()
}##for directories
}

