############################### URBACT bipartite Networks STEP 1 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe URBACT (affiliation) biparti. 
#               Résumé numérique des dimensions des graphes (différents filtrages)
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.5.URBACT/")

#Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(gridExtra)
library(Matrix)


#####Data 

## Membership : table City-Asso. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
MembershipUrbact <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/URBACT_Membership_GNidCorr_complete.csv", stringsAsFactors = F)

## Information on associations (for nodes)
ProjectURBACT <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/UrbactNetworks_complete.csv", 
                        stringsAsFactors = F)

### ==== CREATE NETWORK DATA ====
MembershipUrbact$label <- paste(MembershipUrbact$geonameId, MembershipUrbact$asciiName, sep = "_")
## Create matrix
# edge list
edgelist <- MembershipUrbact  %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network) %>% 
  group_by(geonameId,Code_Network)%>% 
  summarise(weight = n()) # some weight > 1 because several members of an association are (representing the same city) in the same locality 
#(eg. Municipality of Barcelona + Metropolitan area of Barcelona)

edgelistnw <- edgelist %>% mutate(weight = 1) # edgelist non weighted (transform weight > 1)

## Convert edgelist into a matrix

# with possible multiple affiliations of a city in the same organisation (weighted)
#em = edges matrix
em <- edgelist %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()

# Remove duplicate, only unique localities in each association (unweighted)
emnw <- edgelistnw  %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()


## Filter to keep only cities present in two asso

# with possible multiple affiliations of a city in the same organisation (weighted)
em2 <- em[,colSums(em) > 1]

em2f <- round((ncol(em)- ncol(em2))/ncol(em)*100, 0) # remove 58% of cities. Remain 165 cities in more than one association

# Remove duplicate, only unique localities in each association (unweighted)
em2nw <- emnw[,colSums(emnw) > 1]

em2nwf <- round((ncol(emnw)- ncol(em2nw))/ncol(emnw)*100, 0) # remove 58% of cities. Remain 165 cities in more than one association




### Choice to work only on emnw and em2nw, that are binary matrices
### ==== COMPUTE BASIC GLOBAL INDEXES ====


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
  g <- graph.incidence(matrix, directed = FALSE)
  diam <- diameter(g,directed = FALSE)
  
  # comp
  
  comp <- components(g)
  
  nbcomp <- comp$no
  #degrees
  meanDegCities <- mean(colSums(matrix != 0))
  medDegCities <- median(colSums(matrix != 0))
  meanDegAsso <- mean(rowSums(matrix != 0))
  medDegAsso <- median(rowSums(matrix != 0))
  
  
  result <- data.frame(Densité = density, 
                       Diamètre = diam,
                       Taille = nbliens,
                       NbComp = nbcomp,
                       NbVilles = ncities,
                       NbAssos = nasso,
                       MeanDegreeCity =  meanDegCities,
                       MedDegreeCity =  medDegCities,
                       MeanDegreeAsso = meanDegAsso,
                       MedDegreeAsso = medDegAsso)
  return(result)
}



### apply on matrices

#full network with doublons in asso
dimem <- graphlevelindex(em)

# doublons in asso but cities with degree 1 removed
dimem2 <- graphlevelindex(em2)

#full network without doublons in asso
dimemnw <- graphlevelindex(emnw)

# wihtout doublons in asso and cities with degree 1 removed
dimem2nw <- graphlevelindex(em2nw)


## Results in a df
dfnetworklevel <- rbind(dimemnw, dimem2nw) %>% mutate_all(~round(., digits = 3))

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
  select(Filtre,NomMatrix,TypeGraph,Densité, Diamètre,NbComp, Taille,NbVilles, Citiesfilter, everything())

Varnames <- c("Filtrage du graphe",
              "Nom de la matrice",
              "Type de matrice",
              "Densité",
              "Diamètre",
              "Nb composantes",
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


pdf(file= "OUT/dfBipartiURBACT_index.pdf", height = 7, width =6.5 )
grid.table(dfexport)
dev.off()

#### SAVE RDS the em2nw matrix for STEP2

saveRDS(em2nw, "DataProd/URBACT_em2nw.rds")
saveRDS(emnw, "DataProd/URBACT_emnw.rds")



### ==== COMPUTE BASIC GLOBAL INDEXES ON TEMPORAL Matrice ====

# Vectors of networks
Urbact1 <- ProjectURBACT %>% filter(Phase == "Urbact I") %>% select(Code_Network) %>% deframe()
Urbact2 <- ProjectURBACT %>% filter(Phase == "Urbact II") %>% select(Code_Network) %>% deframe()
Urbact3 <- ProjectURBACT %>% filter(Phase == "Urbact III") %>% select(Code_Network) %>% deframe()

# Edges matrix filter

#Filter out networks on the matrix that keep city of degree 1
#Urbact 1
emnw1 <- emnw[rownames(emnw) %in% Urbact1, ]
emnw1 <- emnw1[,colSums(emnw1) > 0]

# Urbact 2
emnw2 <- emnw[rownames(emnw) %in% Urbact2, ]
emnw2 <- emnw2[,colSums(emnw2) > 0]

# Urbact 3

emnw3 <- emnw[rownames(emnw) %in% Urbact3, ]
emnw3 <- emnw3[,colSums(emnw3) > 0]

# Compute global indexes for the 3 matrices

#Urbact 1
dimemnw1 <- graphlevelindex(emnw1)

# Urbact 2

dimemnw2 <- graphlevelindex(emnw2)

# Urbact 3

dimemnw3 <- graphlevelindex(emnw3)

## Results in a df
dfnetworklevel <- rbind(dimemnw1, dimemnw2,dimemnw3) %>% mutate_all(~round(., digits = 3))

#¨Preparation dataframe for export in french



NameFilter <- c("Urbact I",
                "Urbact II",
                "Urbact III") 

TypeGraph <- c("Binaire", "Binaire","Binaire")
TypeGraph



#Insert description of filtering
dfnetworklevel <- dfnetworklevel %>% mutate(TypeGraph = TypeGraph) %>%
  mutate(Filtre = NameFilter) %>% 
  mutate(NomMatrix = c("emnw1", "emnw2", "emnw3"))%>%
  mutate(Période = c("2000-2006", "2007-2013","2014-2020"))%>%
  select(Filtre,NomMatrix,Période,TypeGraph,Densité, Diamètre,NbComp, Taille,NbVilles, everything())

Varnames <- c("Filtrage du graphe",
              "Nom de la matrice",
              "Période",
              "Type de matrice",
              "Densité",
              "Diamètre",
              "Nb composantes",
              "Taille\n(nb de liens)",
              "Nb de Villes",
              "Nb de projets",
              "Degré moyen\ndes villes",
              "Degré médian\ndes villes",
              "Degré moyen\ndes projets",
              "Degré médian\ndes projets")

colnames(dfnetworklevel) <- Varnames




## export df as pdf
dfexport <- t(dfnetworklevel)


pdf(file= "OUT/dfBipartiURBACT_index_phases.pdf", height = 7, width =6.5 )
grid.table(dfexport)
dev.off()


#### SAVE RDS the phase matrices for STEP2

saveRDS(emnw1, "DataProd/URBACT_emnw1.rds")
saveRDS(emnw2, "DataProd/URBACT_emnw2.rds")
saveRDS(emnw3, "DataProd/URBACT_emnw3.rds")

