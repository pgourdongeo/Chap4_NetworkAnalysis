############################### ETMUN bipartites Network  ############
#                               
#                          
# DESCRIPTION : Travail sur le graphe ETMUN (affiliation) biparti
#
#
# 
############################################################################## PG mai 2020


### ==== Load Package and Data ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")

#Packages
library(tidyverse)
library(tidylog)
library(skimr)
library(igraph)
library(gridExtra)
library(Matrix)
library(igraph)
library(clue)
library(sf)
library(mapview)

#####Data 
## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_GNidCorr.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", stringsAsFactors = F)


##### edge list, creating network

## Create matrix
# edge list
edgelist <- MembershipEtmun %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network) %>% 
  group_by(geonameId,Code_Network)%>% summarise(weight = n()) # some weight > 1 because several members of an association are (representing the same city) in the same locality 
#(eg. Municipality of Barcelona + Metropolitan area of Barcelona)

edgelistnw <- edgelist %>% mutate(weight = 1) # edgelist non weighted (transform weight > 1)
# Convert edgelist into a matrix
em <- edgelist %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()

emnw <- edgelistnw  %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()


# Filter to keep only cities present in two asso

em2 <- em[,colSums(em) > 1]
em2f <- round((ncol(em)- ncol(em2))/ncol(em)*100, 0) # remove 85% of cities. Remain 1976 cities in more than one association

em2nw <- emnw[,colSums(emnw) > 1]
em2nwf <- round((ncol(emnw)- ncol(em2nw))/ncol(emnw)*100, 0) # remove 88% of cities. Remain 1545 cities in more than one association

sum(em2nw)/ (ncol(em2nw)*nrow(em2nw))
dim(em2nw)

# fiter covenant and wwcam and ewt
assoremove <- c("04080", "18332","07675")

## Create matrix
# edge list
edgelist <- MembershipEtmun %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network) %>% filter(!Code_Network %in% assoremove)%>%
  group_by(geonameId,Code_Network)%>% summarise(weight = n()) # some weight > 1 because several members of an association are (representing the same city) in the same locality 
#(eg. Municipality of Barcelona + Metropolitan area of Barcelona)

edgelistnw <- edgelist %>% mutate(weight = 1) # edgelist non weighted (transform weight > 1)
# Convert edgelist into a matrix

emnwF <- edgelistnw  %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()


# Filter to keep only cities present in two asso



em3<- emnwF[,colSums(emnwF) > 1]

sum(em2nw)/ (ncol(em2nw)*nrow(em2nw))
dim(em2nw)

#### SIMPLE CODE TO FILTER THE MATRIX by organisations
dim(em3)
dim(em4)
assoremove <- c("04080", "18332","07675")
em4 <- em2nw[!rownames(em2nw) %in% assoremove, ]
em4 <- em4[,colSums(em4) > 1]
all.equal(em4,em3)
####### compute indexes 


graphlevelindex <- function(matrix){
  
  #density
  density <- nnzero(matrix)/ (ncol(matrix)*nrow(matrix))
  #size
  nbliens <- nnzero(matrix)
  #order type 1
   ncities <- ncol(matrix)
   #order type 2
  nasso <- nrow(matrix)
 
  #degrees
  meanDegCities <- mean(colSums(matrix != 0))
  medDegCities <- median(colSums(matrix != 0))
  meanDegAsso <- mean(rowSums(matrix != 0))
  medDegAsso <- median(rowSums(matrix != 0))
  
  # weighted degrees
  
  meanwDegCities <- mean(colSums(matrix))
  medwDegCities <- median(colSums(matrix))
  meanwDegAsso <- mean(rowSums(matrix))
  medwDegAsso <- median(rowSums(matrix))
  
  result <- data.frame(Densité = density, 
                       Taille = nbliens,
                       NbVilles = ncities,
                       NbAssos = nasso,
                       MeanDegreeCity =  meanDegCities,
                       MedDegreeCity =  medDegCities,
                       MeanDegreeAsso = meanDegAsso,
                       MedDegreeAsso = medDegAsso,
                       MeanwDegreeCity = meanwDegCities,
                       MeanwDegreeAsso = meanwDegAsso
                       )
  return(result)
}

max(colSums(em != 0))

### on different matrix
#full network with doublons in asso
dimem <- graphlevelindex(em)

# doublons in asso but cities with degree 1 removed
dimem2 <- graphlevelindex(em2)

#full network without doublons in asso
dimemnw <- graphlevelindex(emnw)

# wihtout doublons in asso and cities with degree 1 removed
dimem2nw <- graphlevelindex(em2nw)


## Results in a df
dfnetworklevel <- rbind(dimem, dimem2, dimemnw, dimem2nw) %>% mutate_all(~round(., digits = 3))

#¨Preparation dataframe for export in french

Varnames <- c("Densité",
              "Taille\n(nb de liens)",
              "Nb de Villes",
              "Nb d'associations",
              "Degré moyen\ndes villes",
              "Degré médian\ndes villes",
              "Degré moyen\ndes associations",
              "Degré médian\ndes associations",
              "Degré pondéré moyen\ndes villes",
              "Degré pondéré moyen\ndes associations"
              )

colnames(dfnetworklevel) <- Varnames

NameFilter <- c("Aucun filtrage", 
                "Filtrage des villes\nprésentes dans une seule\nassociation",
                "Filtrage des doublons\ndans les villes membres\nde chaque association",
                "Filtrage des doublons\ndans les villes membres\nde chaque association\n&\nFiltrage des villes\nprésentes dans une seule\nassociation") 

TypeGraph <- c("Valué", "Valué", "Binaire", "Binaire")
TypeGraph
#Pct cities filtered
Citiesfilter <- c(0,em2f,0,em2nwf)

dfnetworklevel <- dfnetworklevel %>% mutate('Type de graphe' = TypeGraph) %>%
                  mutate(Filtre = NameFilter) %>% 
                  select(Filtre,'Type de graphe', everything())%>%
                  rename("Filtrage du graphe"= Filtre)

### Complete with diameter and clustering coef
clusteringCoef <- function(matrix){
  require(tnet)
  AssoNet <- as.tnet(matrix, type="binary two-mode tnet") ## tnet take the type of nodes of the rows
  CityNet <- as.tnet(t(matrix),type="binary two-mode tnet" )
 
  #reinforcement
  reinforcement <- reinforcement_tm(AssoNet)
   
 
  #clustering Asso
  ClusteringAsso <- clustering_tm(AssoNet)
  
  #clustering Asso
  ClusteringCity<- clustering_tm(CityNet)
  
  result <- data.frame(reinforcement = reinforcement, 
                       ClusteringAsso = ClusteringAsso,
                       ClusteringCity = ClusteringCity)
  return(result)
}

degreeCity <- data.frame(ND = colSums(em2nw)/nrow(em2nw), geonameId = colnames(em2nw))
degreeASSO <- data.frame(ND = rowSums(em2nw)/ncol(em2nw), Code = rownames(em2nw))

degreeCity2 <- data.frame(ND = colSums(em3)/nrow(em3), geonameId = colnames(em3))
degreeASSO2 <- data.frame(ND = rowSums(em3)/ncol(em3), Code = rownames(em3))


normalizedegreeV1 <- data.frame(ND = colSums(m)/nrow(m))
normalizedegreeV2 <- data.frame(ND = rowSums(m)/ncol(m))


clustringem2nw <- clusteringCoef(em2nw)
## export df as pdf
dfexport <- t(dfnetworklevel)


pdf(file= "OUT/dfBipartiEtmun_index.pdf", height = 7, width =10 )
grid.table(dfexport)
dev.off()

                 
#"ISA","SA", "linkage density", "H2"
## Entity level 
AssoIndexes <- specieslevel(em2nw, index=c("degree","species strength" ,"normalised degree","species specificity", "nestedrank", "betweenness", "closeness"), 
                            level="both", nested.method = "binmatnest",nested.weighted= FALSE)

cities <- AssoIndexes$`higher level`
asso <- AssoIndexes$`lower level`
ND <- ND(em2nw, normalised=TRUE)
class(ND$lower)
NDAsso<-ND$lower

citieswithName <- cbind(cities, MembershipEtmun[, "asciiName"][match(rownames(cities), MembershipEtmun$geonameId)])

BC <-BC(em2nw, rescale=TRUE, cmode="undirected", weighted=FALSE)
CC <- CC(em2, cmode="undirected", rescale=TRUE)
degreedistr(em2nw)
degreedistr(em3)
bipartite:::networklevel
#,  "species strength",  "nestedrank", "species specificity", "NS", "betweenness", "closeness", "diversity"
library(sna)
#citation
citation("bipartite")
#liste des jeux de données disponibles
data(package = "bipartite")
#liste des fonctions disponibles
library(help="bipartite")
#chargement d'un jeu de données
#classe et structure
data(Safariland)
class(Safariland) #matrix
edit(Safariland) #graphe bipartite valué
#mesures des caractéristiques d'ensemble
networklevel(Safariland)
#mesures des caractéristiques des sommmets
#nécessite le package sna
specieslevel(Safariland)
#graphique degrés - fréquence cumulée
#noter le message d'avertissement
#la droite de régression placée par défaut n'est
#pas pertinente vu le nombre de sommmets
degreedistr(Safariland)
#comparaison du graphe avec des graphes bipartites aléatoires
#fonction qui peut être très lente selon l'ordre du graphe
#tous les indicateurs ne peuvent être testés selon les caractéristiques
#du graphe de départ
#ici, sont testés le nombre de proie(s) par prédateur ("generality")
#le nombre de prédateur(s) par proie ("vulnerability")
#une mesure de spécialisation du graphe ("H2")
#la force des interactions ("ISA") et le niveau d'asymétrie ("SA")
#voir la documentation pour plus de détails sur les mesures disponibles
Nulltestem2nw <-null.t.test(em2nw, index=c("generality", "vulnerability", "connectance","links per species"), nrep=4, N=20)
NulltestCC <-null.t.test(em2nw, index=c("cluster coefficient"), nrep=4, N=20)

Nulltestem3 <-null.t.test(em3, index=c("generality", "vulnerability", "connectance","links per species"), nrep=4, N=20)
NulltestCC <-null.t.test(em3, index=c("cluster coefficient"), nrep=4, N=20)
#visualisation de la matrice d'adjacence
visweb(Safariland)
#avec des cercles proportionnels et de la couleur...
visweb(Safariland, circles=TRUE, boxes=TRUE, outerbox.col="green",
       labsize=1, circle.max=1.8, text="no", box.border="red")
#autre possibilité de présentation
visweb(Safariland,square="b",box.col="green",box.border="red")
#visualisation du graphe
plotweb(Safariland)
#en colorant liens et sommets et en orientant les étiquettes
plotweb(Safariland, text.rot=90, col.high="blue",
        col.low="green", col.interaction="orange")
#pour connaître toutes les options disponibles
?plotweb
#visualisation des interactions entre espèces
#les lignes symbolisent les pollinisateurs/prédateurs
plotPAC(Safariland)
#recherche des clusters à l'aide de la modularité
#liste des appartenances aux clusters
#et visualisation
cW <- computeModules(Safariland, deep = FALSE, deleteOriginalFiles = TRUE,
                     steps = 1000000)
listModuleInformation(cW)
plotModuleWeb(cW)
#étudier les conséquences de la disparition d'une espèce
#par défaut, l'espèce la moins nombreuse est éliminée en premier
second.extinct(Safariland, participant="low", method="r", nrep=50,
               details=TRUE)
#tranformation en graphe valué de coappartenance
#les objets créés sont de type "matrix"
transL <- as.one.mode(Safariland, project="lower") #origine - origine
transH <- as.one.mode(Safariland, project="higher")#destination - destination
#transformation en liste de liens valués pour une utilisation
#avec un autre logiciel (ex. Pajek)
#la fonction écrit par défaut deux fichiers sans extension sur le disque dur
#l'un avec les noms, l'autre aves les liens (origine, destination, poids)
web2edges(Safariland, out.files=c("edges", "names", "groups")[1:2])
#générer un graphe bipartite aléatoire
#Nx : nombre de sommets dans le groupe x
#dens : nombre de liens moyen pour chaque espèce (2 par défaut)
g <- genweb(N1 = 10, N2 = 30, dens=1)



library(tnet)
data(tnet)

assonumber <- AssoEtmun %>% select(Code)%>% mutate(id = rownames(.))
net <- edgelistnw %>% select(-weight)%>% left_join(assonumber, by = c("Code_Network"= "Code"))%>% select(-Code_Network)%>%
  select(id, geonameId)%>%as.data.frame()

net$geonameId <- as.integer(net$geonameId)
net$id <- as.integer(net$id)

net <- as.tnet(net, type="binary two-mode tnet")

# Calculate the reinforcement coefficient (Robins and Alexander, 2004)
reinforcement_tm(net)

# Calculate the global coefficient (Opsahl, 2012)
clustering_tm(net)

g2m<- rbind(
  c(1,1,3),
  c(2,1,2),
  c(1,2,2),
  c(3,2,1),
  c(3,3,2),
  c(4,1,2),
  c(4,3,1))
g2m <- g2m[,1:2]
g2 <- as.tnet(g2m, type="binary two-mode tnet")
distance_tm(g2)
clustering_tm(g2)
clustering_local_tm(g2)
reinforcement_tm(g2)


##########################"" Test with https://rpubs.com/pjmurphy/317838  methodo



library(ade4) # If you have not already done so


asso_jaccard <- dist.binary(emnw, method=1, upper=TRUE, diag = FALSE) # Method #1 is "Jaccard Index"
city_jaccard <- dist.binary(t(emnw), method=1, upper=TRUE, diag = FALSE) 

asso_jaccard <- as.matrix(asso_jaccard)   
diag(asso_jaccard)<-0
summary(nnzero(asso_jaccard))
# asso_jaccard          # Look at the matrix before you binarize
asso_jaccard <- ifelse(asso_jaccard>0.95, 1, 0)     # Binarize

# jaccard_asso     # Take a look at the matrix if you like.

jacc_asso <- graph_from_adjacency_matrix(asso_jaccard,    # Create an igraph network
                                          mode = "undirected")
plot(jacc_asso)

#### City
city_jaccard <- as.matrix(city_jaccard)   
diag(city_jaccard)<-0
summary(city_jaccard)
# women_jaccard          # Look at the matrix before you binarize
city_jaccard <- ifelse(city_jaccard>0.95, 1, 0)     # Binarize

# jaccard_women      # Take a look at the matrix if you like.

jacc_city<- graph_from_adjacency_matrix(city_jaccard,    # Create an igraph network
                                         mode = "undirected")
plot(jacc_city)









NMIs<-function(matrix){
 # a is the length of the first mode
  # b is the length of the second mode
  # NOTE: a should be less than b (if not, take the transpose)
  # den is the (asymptotic) density of the simulated network
  # p is the probability of a within-community tie
  
  A<-matrix
  B<-matrix(0,nr=nrow(matrix),nc=nrow(matrix))
  C<-matrix(0,nr=ncol(matrix),nc=ncol(matrix))
  A2<-cbind(B,A)
  A3<-cbind(t(A),C)
  X<-rbind(A2,A3)
  X2<-graph.adjacency(X,mode="undirected")
  k<-(is.connected(X2))*1

  C<-A%*%t(A)
  for(i in 1:dim(A)[1]){C[i,i]<-0}
  D<-t(A)%*%A
  for(i in 1:dim(A)[2]){D[i,i]<-0}
  if (k==1) y1<-walktrap.community(graph.adjacency(C,weighted=TRUE))$membership else y1<-0
  if (k==1) y1<-as.matrix(y1) else y1<-0
  if (k==1) Y1<-as.numeric(y1) else y1<-0
  if (k==1) y2<-walktrap.community(graph.adjacency(D,weighted=TRUE))$membership  else y2<-0
  if (k==1) y2<-as.matrix(y2) else y2<-0
  if (k==1) Y2<-as.numeric(y2) else y2<-0

 
  # x1 is the NMI for walktrap on the first mode
  # x2 is the NMI for walktrap on the second mode
  #
  # Y1 is the NMI for walktrap on the first mode (dual projection)
  # Y2 is the NMI for walktrap on the second mode (dual projection)
  NMIs<-list(Y1,Y2)
  NMIs<-t(as.matrix(NMIs))
  names<-c('DP_Walk_First','DP_Walk_Second')
  colnames(NMIs)<-names
  return(NMIs)
  }

clustering<-NMIs(em2nw)

ClusterAsso <- data.frame(Code = rownames(em2nw), clusterAsso = as.factor(clustering[[1]])) %>% 
  left_join(select(AssoEtmun,Acronym,Name, Code )) %>% 
  left_join(degreeASSO)

ClusterCities <- data.frame(geonameId = colnames(em2nw), ClusterCities = as.factor(clustering[[2]])) %>% 
  left_join(select(MembershipEtmun,asciiName, geonameId,lng_GN,lat_GN )) %>% 
  distinct()%>%
  left_join(degreeCity)

skim(ClusterCities$ClusterCities)
table(ClusterCities$ClusterCities)
saveRDS(ClusterCities, "clusterCities1.rds")
#citiesclusterSF <- st_as_sf(ClusterCities, coords = c("lng_GN", "lat_GN"), crs = 4326)


## filter
clustering2 <- NMIs(em3)

ClusterAsso2 <- data.frame(Code = rownames(em3), clusterAsso = as.factor(clustering2[[1]])) %>% 
  left_join(select(AssoEtmun,Acronym,Name, Code )) %>% 
  left_join(degreeASSO2)

table(ClusterAsso2$clusterAsso)

ClusterCities2 <- data.frame(geonameId = colnames(em3), ClusterCities = as.factor(clustering2[[2]])) %>% 
  left_join(select(MembershipEtmun,asciiName, geonameId,lng_GN,lat_GN )) %>% 
  distinct()%>%
  left_join(degreeCity2)

table(ClusterCities2$ClusterCities)

saveRDS(ClusterCities2, "clusterCities2.rds")
#citiesclusterSF2 <- st_as_sf(ClusterCities2, coords = c("lng_GN", "lat_GN"), crs = 4326)
