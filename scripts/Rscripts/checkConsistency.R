
########################
#
##This script checks that a particular wing is consistent. 
# First it checks that each vertex only appears in 3 or less edges and cells,
#   and that no more than 3 vertices have it as neighbour. 
#   Checks that each edge does not appear in more than 2 vertices and cells. 
# Then checks reciprodicy:
#   Checks, for each vertex, that all its neighbours have it as neighbour, 
#   all its edges have it as vertex and all of its cells have it as a vertex.
#   For each edge, it checks that all of its vertices and cells have it as an edge
#   For each cell, it checks that all of its vertices and edges have it as a cell
# Prints the output as a list. 
#
# To run it:
#   open a terminal in the directory in which the wing is found
#   run, for example:
#   Rscript ../../../scripts/checkConsistency.R etournay1_strings3_moved_0
#
############################

library(tidyverse)

wingname <- "etournay1_strings3_not_consistent_15"
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) wingname <- args[1]
  
cat("\n\nWINGNAME: ", wingname, "\n\n***\n\n")
f <- paste(wingname, c(".ptab", ".celltab", ".edges", ".sprtab"), sep="")



points <-read.table(f[1], sep="\t", dec=".", header=T, stringsAsFactors = F) %>% 
  mutate(cells = strsplit(cells, ","),
         edges = strsplit(edges, ","),
         neighbour_vertices = strsplit(neighbour_vertices, ","))
cells <-read.table(f[2], sep="\t", dec=".", header=T, stringsAsFactors = F) %>% 
  mutate(edges = strsplit(edges, ","),
         vertices = strsplit(vertices, ","))
edges <-read.table(f[3], sep="\t", dec=".", header=T, stringsAsFactors = F) %>% 
  mutate(cells = strsplit(cells, ","),
         vertices = strsplit(vertices, ","))

spr <-read.table(f[4], sep="\t", dec=".", header=T, stringsAsFactors = F)


#1) Count apparition

verts_in_edges <- edges$vertices %>% unlist %>% table %>% subset(. > 3)
verts_in_cells <- cells$vertices %>% unlist %>% table %>% subset(. > 3)
verts_in_neighbours <- points$neighbour_vertices %>% unlist %>% table %>% subset(. > 3)

edges_in_cells <- cells$edges %>% unlist %>% table %>% subset(. > 2)
edges_in_verts <- points$edges %>% unlist %>% table %>% subset(. > 2)

countlist <- list("verts_in_edges" = verts_in_edges,
                  "verts_in_cells" = verts_in_cells,
                  "verts_in_neighbours" = verts_in_neighbours,
                  "edges_in_cells" = edges_in_cells,
                  "edges_in_verts" = edges_in_verts)

cat("Elements with more than maximum apparitions in given list:\n")
print(countlist)
#Check reciprocal neighbours
cat("\n******\n")


checkReciprocity <- function(tab1, var1, tab2, var2){
  res = c()
  for(i in 1:nrow(tab1)){
    inds <- tab1[i, var1] %>% unlist %>% as.integer() %>% subset(. != -999)
    reciprocal <- c()
    for(j in inds){
      neig2 <- tab2[tab2$ind==j, var2] %>% unlist %>% as.integer
      reciprocal <- c(reciprocal, ifelse(tab1[i, "ind"] %in% neig2, "", as.character(j) ))
    }

  res = c(res, paste(reciprocal, sep=",", collapse=","))
  }
  res <- gsub("^,|,$|,,", "", res)
  names(res) <- tab1$ind 
  res <- res[res != ""]
  return(res)
}

#Run all in parallel
library("future")
plan(multiprocess)

#points %>% filter(sapply(edges, FUN=function(x){22436 %in% x}))

rec1 %<-% checkReciprocity(points, "cells", cells, "vertices")
rec2 %<-% checkReciprocity(points, "edges", edges, "vertices")
rec3 %<-% checkReciprocity(points, "neighbour_vertices", points, "neighbour_vertices")
rec4 %<-% checkReciprocity(edges, "cells", cells, "edges")

#Now these would be repeated (just to make sure)
rec5 %<-% checkReciprocity(edges, "vertices", points, "edges")
rec6 %<-% checkReciprocity(cells, "vertices", points, "cells")
rec7 %<-% checkReciprocity(cells, "edges", edges, "cells")

reclist <- lapply(ls(pattern = "^rec[1-7]"), get)
names(reclist) <- c("Cells in Vertices and Vertices in Cells",
                    "Edges in Vertices and Vertices in Edges",
                    "Neighbours in Vertices",
                    "Cells in Edges and Edges in Cells",
                    "Vertices in Edges and Edges in Vertices (bis)",
                     "Vertices in Cells and Cells in Vertices",
                    "Edges in Cells and Cells in Edges")

cat("Non-reciprocal apparitions:\n")
print(reclist)
cat("\n\n********\n\n")


## Check consistency in cells

res <- c()
for(i in 1:nrow(cells)){
  thisv <- points[points$ind %in% as.integer(unlist(cells$vertices[i])), ]
  thise <-edges[edges$ind %in% as.integer(unlist(cells$edges[i])), ]
  num_ok <- (nrow(thisv) == cells$num_vertices[i]) &&  (nrow(thise) == cells$num_vertices[i]) 
  edge_verts <- thise$vertices %>% unlist %>% unique %>% as.integer %>% subset(. != -999)
  v_in_e <- all(thisv$ind %in% edge_verts) 
  v_in_e2 <- all(edge_verts %in% thisv$ind)
  vert_edges <- thisv$edges %>% unlist %>% unique %>% as.integer %>% subset(. != -999)
  e_in_v <- all(thise$ind %in% vert_edges) 
  this_res <- paste(
    ifelse(num_ok, "", as.character(c(cells$num_vertices[i],
                                      cells$vertices[i] %>% unlist %>% as.character %>% paste(sep=",", collapse=","), 
                                      cells$edges[i] %>% unlist %>% as.character %>% paste(sep=",", collapse=",")
                                      )) %>% paste(sep="|", collapse="|")),
    ifelse(v_in_e, "", ";Some verts in cell not in edges;"),
    ifelse(v_in_e2, "", "Some verts in cell edges not in cell;"),
    ifelse(e_in_v, "", "Some edges in cell not in vertices"),
    sep="", collapse=""
  )
  res <- c(res, this_res)
}

names(res) <- cells$ind
res<-res[res != ""]

cat("Cells with non-coincident edges and vertices:\n")
print(res)

