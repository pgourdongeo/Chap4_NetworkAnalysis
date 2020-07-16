############################### EUCICOP bipartite Networks STEP 1 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe EUCICOP (affiliation) biparti. 
#               Résumé numérique des dimensions des graphes (différents filtrages)
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/EUCICOP_Inter/")

#Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(gridExtra)
library(Matrix)


#####Data 

## Participation (partnership): table City-Project. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
ParticipationEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/Participations_All_Eucicop_Europe.RDS")

## Information on project (for nodes)
ProjectEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/ProjectsEucicop_all_noduplicated.RDS")

### ==== CREATE NETWORK DATA ====

## Create matrix
# edge list
edgelist <- ParticipationEUCICOP %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, ID_PROJECT) %>% 
  group_by(geonameId,ID_PROJECT)%>% 
  summarise(weight = n()) # some weight > 1 because several members of a project are (representing the same city) in the same locality 
#(eg. Municipality of Barcelona + Metropolitan area of Barcelona)


edgelistnw <- edgelist %>% mutate(weight = 1) # edgelist non weighted (transform weight > 1)


# Remove project with only one locality
Morethan1project <- edgelistnw %>% 
  group_by(ID_PROJECT) %>% 
  summarise(N= n()) %>% 
  filter(N> 1) %>% 
  select(ID_PROJECT) %>% 
  deframe()

1- length(Morethan1project)/nrow(ProjectEUCICOP) # 20% of projects with only one locality

## Filter the edgelist

edgelistnw <- edgelistnw %>% filter(ID_PROJECT %in% Morethan1project)

## Convert edgelist into a matrix

# with possible multiple affiliations of a city in the same organisation (weighted)
#em = edges matrix
em <- edgelist %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="ID_PROJECT") %>% 
  as.matrix()

# Remove duplicate, only unique localities in each project (unweighted)
emnw <- edgelistnw  %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="ID_PROJECT") %>% 
  as.matrix()


## Filter to keep only cities present in two project

# with possible multiple affiliations of a city in the same organisation (weighted)
em2 <- em[,colSums(em) > 1]

em2f <- round((ncol(em)- ncol(em2))/ncol(em)*100, 0) # remove 42% of cities. Remain 5346 cities in more than one project

# Remove duplicate, only unique localities in each project (unweighted)
em2nw <- emnw[,colSums(emnw) > 1]

em2nwf <- round((ncol(emnw)- ncol(em2nw))/ncol(emnw)*100, 0) # remove 44% of cities. Remain 5025 cities in more than one project



### Choice to work only on emnw and em2nw, that are binary matrices
### ==== COMPUTE BASIC GLOBAL INDEXES ====


graphlevelindex_severalcomp <- function(matrix){
  require(igraph)
  require(Matrix)
  
  #density
  density <- nnzero(matrix)/ (ncol(matrix)*nrow(matrix))
  #size
  nbliens <- nnzero(matrix)
  #order type 1
  ncities <- ncol(matrix)
  #order type 2
  nasso <- nrow(matrix)
  #Diameter
  g <- graph.incidence(matrix, directed = FALSE)
  diam <- diameter(g,directed = FALSE)
  # comp
  
  comp <- components(g)
  
  nbcomp <- comp$no
  
  sizecomp <- 100*comp$csize[[1]]/vcount(g)
  #degrees
  meanDegCities <- mean(colSums(matrix != 0))
  medDegCities <- median(colSums(matrix != 0))
  meanDegAsso <- mean(rowSums(matrix != 0))
  medDegAsso <- median(rowSums(matrix != 0))
  
  
  result <- data.frame(Densité = density, 
                       Diamètre = diam,
                       Taille = nbliens,
                       NbComp = nbcomp ,
                       SizeComp = sizecomp,
                       NbVilles = ncities,
                       NbAssos = nasso,
                       MeanDegreeCity =  meanDegCities,
                       MedDegreeCity =  medDegCities,
                       MeanDegreeAsso = meanDegAsso,
                       MedDegreeAsso = medDegAsso)
  return(result)
}



### apply on matrices

#full network with doublons in project
dimem <- graphlevelindex_severalcomp(em)

# doublons in project but cities with degree 1 removed
dimem2 <- graphlevelindex_severalcomp(em2)

#full network without doublons in project
dimemnw <- graphlevelindex_severalcomp(emnw)

# without doublons in project and cities with degree 1 removed
dimem2nw <- graphlevelindex_severalcomp(em2nw)


## Results in a df
dfnetworklevel <- rbind(dimemnw, dimem2nw) %>% mutate_all(~round(., digits = 5))

#¨Preparation dataframe for export in french



NameFilter <- c("Filtrage des doublons\ndans les villes membres\nde chaque projet",
                "Filtrage des doublons\ndans les villes membres\nde chaque projet\n&\nFiltrage des villes\nprésentes dans un seul\nprojet") 

TypeGraph <- c("Binaire", "Binaire")
TypeGraph

# Pct cities filtered

Citiesfilter <- c(0,em2nwf)

#Insert description of filtering
dfnetworklevel <- dfnetworklevel %>% mutate(TypeGraph = TypeGraph) %>%
  mutate(Filtre = NameFilter) %>% 
  mutate(Citiesfilter = Citiesfilter)%>%
  mutate(NomMatrix = c("emnw", "em2nw"))%>%
  select(Filtre,NomMatrix,TypeGraph,Densité, Diamètre,NbComp,SizeComp, Taille,NbVilles, Citiesfilter, everything())

Varnames <- c("Filtrage du graphe",
              "Nom de la matrice",
              "Type de matrice",
              "Densité",
              "Diamètre",
              "Nb Composantes Connexes",
              "% de Sommets de la 1ère composante",
              "Taille\n(nb de liens)",
              "Nb de Villes",
              "Pct villes filtrées",
              "Nb de projets",
              "Degré moyen\ndes villes",
              "Degré médian\ndes villes",
              "Degré moyen\ndes projets",
              "Degré médian\ndes projets")

colnames(dfnetworklevel) <- Varnames




## export df as pdf
dfexport <- t(dfnetworklevel)


pdf(file= "OUT/dfBipartiEucicop_index.pdf", height = 7, width =6.5 )
grid.table(dfexport)
dev.off()

### ==== KEEP ONLY THE FIRST COMPONENT ====

# emnw
emnw.g <- graph.adjacency(emnw, directed = FALSE)
Comp <- components(emnw.g)
dg <- decompose.graph(emnw.g)
Comp1_emnw <- dg[[1]]
Comp <- components(Comp1_emnw)$no

emnw_1stcomp <- as_edgelist(Comp1_emnw) %>% as.data.frame() %>% mutate(weight = 1)

emnw_1stcomp <- emnw_1stcomp  %>% 
  pivot_wider(names_from = V2, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="V1") %>% 
  as.matrix()


# em2nw
em2nw.g <- graph.adjacency(em2nw, directed = FALSE)
Comp <- components(em2nw.g)
dg <- decompose.graph(em2nw.g)
Comp1_em2nw <- dg[[1]]
Comp <- components(Comp1_em2nw)$no

em2nw_1stcomp <- as_edgelist(Comp1_em2nw) %>% as.data.frame() %>% mutate(weight = 1)

em2nw_1stcomp <- em2nw_1stcomp  %>% 
  pivot_wider(names_from = V2, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="V1") %>% 
  as.matrix()

### ==== COMPUTE BASIC GLOBAL INDEXES ON FIRST COMPONENT MATRIX ====

graphlevelindex <- function(matrix){
  require(igraph)
  require(Matrix)
  
  #density
  density <- nnzero(matrix)/ (ncol(matrix)*nrow(matrix))
  #size
  nbliens <- nnzero(matrix)
  #order type 1
  ncities <- ncol(matrix)
  #order type 2
  nasso <- nrow(matrix)
  #Diameter
  diam <- diameter(graph.incidence(matrix, directed = FALSE),directed = FALSE)
  #degrees
  meanDegCities <- mean(colSums(matrix != 0))
  medDegCities <- median(colSums(matrix != 0))
  meanDegAsso <- mean(rowSums(matrix != 0))
  medDegAsso <- median(rowSums(matrix != 0))
  
  
  result <- data.frame(Densité = density, 
                       Diamètre = diam,
                       Taille = nbliens,
                       NbVilles = ncities,
                       NbAssos = nasso,
                       MeanDegreeCity =  meanDegCities,
                       MedDegreeCity =  medDegCities,
                       MeanDegreeAsso = meanDegAsso,
                       MedDegreeAsso = medDegAsso)
  return(result)
}

# Summary on first component matrices
dimemnw_1stcomp <- graphlevelindex(emnw_1stcomp)

dimem2nw_1stcomp <- graphlevelindex(em2nw_1stcomp)

dfnetworklevel <- rbind(dimemnw_1stcomp, dimem2nw_1stcomp) %>% mutate_all(~round(., digits = 5))

#¨Preparation dataframe for export in french



NameFilter <- c("Plus grande\ncomposante connexe\nde emnw\n(98% des sommets)",
                "Plus grande\ncomposante connexe\nde em2nw\n(99% des sommets)") 

TypeGraph <- c("Binaire", "Binaire")
TypeGraph

# Pct cities filtered
em2nwf <- round((ncol(emnw_1stcomp)- ncol(em2nw_1stcomp))/ncol(emnw_1stcomp)*100, 0) 
Citiesfilter <- c(0,em2nwf)

#Insert description of filtering
dfnetworklevel <- dfnetworklevel %>% mutate(TypeGraph = TypeGraph) %>%
  mutate(Filtre = NameFilter) %>% 
  mutate(Citiesfilter = Citiesfilter)%>%
  mutate(NomMatrix = c("emnw_1stcomp", "em2nw_1stcomp"))%>%
  select(Filtre,NomMatrix,TypeGraph,Densité, Diamètre, Taille,NbVilles, Citiesfilter, everything())

Varnames <- c("Filtrage du graphe",
              "Nom de la matrice",
              "Type de matrice",
              "Densité",
              "Diamètre",
              "Taille\n(nb de liens)",
              "Nb de Villes",
              "Pct villes filtrées",
              "Nb de projets",
              "Degré moyen\ndes villes",
              "Degré médian\ndes villes",
              "Degré moyen\ndes projets",
              "Degré médian\ndes projets")

colnames(dfnetworklevel) <- Varnames




## export df as pdf
dfexport <- t(dfnetworklevel)


pdf(file= "OUT/dfBipartiEucicop_1stcomp_index.pdf", height = 7, width =6.5 )
grid.table(dfexport)
dev.off()

#### SAVE RDS the em2nw matrix for STEP2

saveRDS(em2nw, "DataProd/Eucicop_em2nw.rds")
saveRDS(emnw, "DataProd/Eucicop_emnw.rds")

saveRDS(em2nw_1stcomp, "DataProd/Eucicop_em2nw_1stcomp.rds")
saveRDS(emnw_1stcomp, "DataProd/Eucicop_emnw_1stcomp.rds")
