############################### ETMUN bipartite Networks STEP 1 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe ETMUN (affiliation) biparti. 
#               Résumé numérique des dimensions des graphes (différents filtrages)
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")

#Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(gridExtra)
library(Matrix)


#####Data 

## Membership : table City-Asso. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_europe.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", 
                       stringsAsFactors = F)

### ==== CREATE NETWORK DATA ====
MembershipEtmun$label <- paste(MembershipEtmun$geonameId, MembershipEtmun$asciiName,MembershipEtmun$CountryCode, sep = "_")
## Create matrix
# edge list
edgelist <- MembershipEtmun %>% 
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

em2f <- round((ncol(em)- ncol(em2))/ncol(em)*100, 0) # remove 85% of cities. Remain 1976 cities in more than one association

# Remove duplicate, only unique localities in each association (unweighted)
em2nw <- emnw[,colSums(emnw) > 1]

em2nwf <- round((ncol(emnw)- ncol(em2nw))/ncol(emnw)*100, 0) # remove 88% of cities. Remain 1545 cities in more than one association

# Remove WWCAM ASSO with only one european city

assoremove <- "18332"

#Filter out the association on the matrix that keep city of degree 1
emnw <- emnw[!rownames(emnw) %in% assoremove, ]
emnw <- emnw[,colSums(emnw) > 0]

#Filter out the association on the matrix that keep city of more than 1 degree (only cities involved at least in 2 associations)
em2nw <-  em2nw[!rownames(em2nw) %in% assoremove, ]
em2nw <- em2nw[,colSums(em2nw) > 1]


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



NameFilter <- c("Filtrage des doublons\ndans les villes membres\nde chaque association",
                "Filtrage des doublons\ndans les villes membres\nde chaque association\n&\nFiltrage des villes\nprésentes dans une seule\nassociation") 

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
              "Nb d'associations",
              "Degré moyen\ndes villes",
              "Degré médian\ndes villes",
              "Degré moyen\ndes associations",
              "Degré médian\ndes associations")

colnames(dfnetworklevel) <- Varnames




## export df as pdf
dfexport <- t(dfnetworklevel)


pdf(file= "OUT/dfBipartiEtmun_index.pdf", height = 7, width =6.5 )
grid.table(dfexport)
dev.off()

#### SAVE RDS the em2nw matrix for STEP2

saveRDS(em2nw, "DataProd/Etmun_em2nw.rds")
saveRDS(emnw, "DataProd/Etmun_emnw.rds")
