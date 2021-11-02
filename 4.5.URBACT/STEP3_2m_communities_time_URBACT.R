############################### URBACT bipartite Networks STEP 3 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe URBACT (affiliation) biparti. 
#               Détection de communautés / Mots clés / analyse temporelle
#
# 
############################################################################## PG juin 2020

### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.5.URBACT/")

### Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(GGally)
library(tidygraph)
library(ggraph)
library(clv)
library(class)
library(mcclust) 
library(scales)
library(sf)
library(patchwork)
options(scipen = 999)

### Data 

## Participation (partnership): table City-Project. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
ParticipationURBACT <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/URBACT_Membership_GNidCorr_complete.csv", stringsAsFactors = F)

## Information on project (for nodes)
ProjectURBACT <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/UrbactNetworks_complete.csv", 
                           stringsAsFactors = F)
## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")
DbCityInfos <- DbCity %>% st_drop_geometry()

# Infos on Nodes (STEP 2)

NodesInfos <- readRDS("DataProd/NodesInfos_URBACT.RDS")

## shapes

## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)

# Matrices from STEP1  

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/URBACT_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/URBACT_em2nw.rds")

## Temporal matrices

# Urbact 1
emnw1 <- readRDS("DataProd/URBACT_emnw1.rds")
# Urbact 2
emnw2 <- readRDS("DataProd/URBACT_emnw2.rds")
# Urbact 3
emnw3 <- readRDS("DataProd/URBACT_emnw3.rds")

### ==== IGRAPH TRANSFORMATION ====



emnw.g <- graph.incidence(emnw , directed = F)

em2nw.g <-  graph.incidence(em2nw, directed = F)





### ==== COMMUNITY DETECTIONS ON WHOLE GRAPH ====


clustering <- function(g){
  require(igraph)
  results <- list()
  
  #partition
  louvain <- cluster_louvain(g)
  infomap <- infomap.community(g, nb.trials = 20)
  greedy <- cluster_fast_greedy(g)
  walktrap <- walktrap.community(g, steps = 10)
  
  #membership
  V(g)$louvain <- louvain$membership
  V(g)$infomap <-  infomap$membership
  V(g)$fgreedy <-  greedy$membership
  V(g)$walktrap <-walktrap$membership
  
  memberships <- igraph::as_data_frame(g, what = "vertices")
  
  results[["memberships"]] <- memberships
  
  
  infoclustering <- data.frame(Algo = c("louvain", "infomap", "fgreedy", "walktrap"), 
                               modularity = c(max(louvain$modularity), infomap$modularity, modularity(greedy), modularity(walktrap) ),
                               nclust = c(length(unique(louvain$membership)), 
                                          length(unique(infomap$membership)), 
                                          length(unique(greedy$membership)), 
                                          length(unique(walktrap$membership))))
  
  results[["info_clustering"]] <- infoclustering 
  results[["igraph_object"]] <- g
  return(results)
}

emnw.clustering <- clustering(emnw.g)

em2nw.clustering <- clustering(em2nw.g)


emnw.clustering$info_clustering
em2nw.clustering$info_clustering

## Check modularity

Modularity <- bind_rows(emnw.clustering$info_clustering,
                        em2nw.clustering$info_clustering, .id = "Matrix")

Modularity <- Modularity %>% mutate(Matrix = case_when(Matrix == 1 ~ "emnw", Matrix == 2 ~ "em2nw"))



g1 <- ggplot(Modularity) + 
  geom_col(aes(x = Algo, y = modularity, fill = Algo)) + 
  facet_wrap(~Matrix) + 
  labs(x = "", y = "modularité", fill = "Algorithme", subtitle = "Modularité")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), text = element_text(size = 16))

g2 <- ggplot(Modularity) + 
  geom_col(aes(x = Algo, y = nclust, fill = Algo), show.legend = FALSE) + 
  facet_wrap(~Matrix)+ scale_y_log10(n.breaks =8)+
  labs(x = "", y = "nombre de communautés (log 10)", subtitle = "Nombre de communautés", 
       caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), text = element_text(size = 16))


gridMod <- g1 + g2                                                                                                                         
ggsave(gridMod, filename = "OUT/Modularity_CombinedApproach_4matrix.pdf", width = 11.7, height = 8.3, units = "in" )


# Save membership 


Memberem2nw <- em2nw.clustering$memberships 

Memberemnw <- emnw.clustering$memberships

### ==== PAIWISE COMPARISON OF CLUSTERING ====
### Dataframe Classification
dfclustering <- emnw.clustering$memberships
DfCluster <- dfclustering %>% select(-type)


### Compute simulations of each indicators for comparision between random clustering



### Jaccard  

Nclust <- ncol(DfCluster)
Cn <- colnames(DfCluster)[c(2:Nclust)]
Nrepet <- 50
ValJaccard <- list()

for(i in Cn ){
  RealCluster <- DfCluster[,i]
  for(j in 1:Nrepet){
    ShuffledCluster <- sample(RealCluster)
    ValJaccard[paste(i,j,sep="_")] <- clv.Jaccard(std.ext(RealCluster,ShuffledCluster))
  }
  
}
tabJaccardsSimul <- as.data.frame(unlist(ValJaccard, use.names = TRUE))
colnames(tabJaccardsSimul) <- "Jaccard"
tabJaccardsSimul<-tabJaccardsSimul %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J")) %>%
  spread(key = I, value = Jaccard)

ResultSimulJaccard <- tabJaccardsSimul %>% select(-J)%>% gather()%>%
  group_by(key)%>% summarise(mean = mean(value), sd = sd(value))%>% rename("Algo" = key)




### VI  




Nrepet <- 50
ValVI <- list()
for(i in Cn ){
  RealCluster <- DfCluster[,i]
  for(j in 1:Nrepet){
    ShuffledCluster <- sample(RealCluster)
    ValVI[paste(i,j,sep="_")] <- vi.dist(RealCluster,ShuffledCluster, parts = T, base = 3217)
  }
  
}
tabVISimul <- as.data.frame(unlist(ValVI, use.names = TRUE))
colnames(tabVISimul) <- "VI"
tabVISimul<-tabVISimul %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J")) %>%
  spread(key = I, value = VI)

ResultSimulVI <- tabVISimul %>% 
  select(-J)%>% gather()%>%
  group_by(key)%>% summarise(mean = mean(value), sd = sd(value))%>%
  rename("Algo" = key)

### Arandi


Nrepet <- 50
ValArandi <- list()
for(i in Cn ){
  RealCluster <- DfCluster[,i]
  for(j in 1:Nrepet){
    ShuffledCluster <- sample(RealCluster)
    ValArandi[paste(i,j,sep="_")] <- arandi(RealCluster,ShuffledCluster, adjust = T)
  }
  
}
tabArandiSimul <- as.data.frame(unlist(ValArandi, use.names = TRUE))
colnames(tabArandiSimul) <- "Arandi"
tabArandiSimul<-tabArandiSimul %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J")) %>%
  spread(key = I, value = Arandi)

ResultSimulArandi <- tabArandiSimul %>% 
  select(-J)%>% gather()%>%
  group_by(key)%>% summarise(mean = mean(value), sd = sd(value))%>%
  rename("Algo" = key)


#### COMPUTE THE 3 INDEXES FOR EACH PAIRWISE COMPARISONS



index <- list()

for (i in Cn) {
  Clust1 <- DfCluster[, i]
  
  for (j in Cn){
    
    Clust2 <- DfCluster[, j]
    
    
    
    ###Count pairs
    
    Pairs <- std.ext(Clust1, Clust2)
    
    ## Compute Index
    
    
    Jaccard <- clv.Jaccard(Pairs)
    
    index[paste(i,j,sep="_")]<- Jaccard
  }
  
}  


tabJaccards <- as.data.frame(unlist(index, use.names = TRUE))
colnames(tabJaccards) <- "Jaccard"
tabJaccards<-tabJaccards %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J"))

ResultSimulJaccard2 <-  ResultSimulJaccard%>%
  rename(I = Algo) %>%
  select(-sd)%>%
  mutate(J = "Z_Mean50RandomSimul")%>%
  rename(Jaccard = mean)%>% select(I,J,everything())

tabJaccards <- rbind(tabJaccards, ResultSimulJaccard2)


### Loop Adjusted Rand Index

index3 <- list()

for (i in Cn) {
  Clust1 <- DfCluster[, i]
  
  for (j in Cn){
    
    Clust2 <- DfCluster[, j]
    
    ## Compute Index
    
    
    Arandi <- arandi(Clust1,Clust2, adjust =  T)
    
    index3[paste(i,j,sep="_")]<- Arandi
  }
  
}  


tabArandi <- as.data.frame(unlist(index3, use.names = TRUE))
colnames(tabArandi) <- "Arandi"
tabArandi<-tabArandi %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J")) 
#spread(key = I, value = Arandi)


#Add simul
ResultSimulArandi2 <- ResultSimulArandi%>%
  rename(I = Algo) %>%
  select(-sd)%>%
  mutate(J = "Z_Mean50RandomSimul")%>%
  rename(Arandi = mean) %>% select(I,J,everything())

tabArandi <- rbind(tabArandi, ResultSimulArandi2)

### Loop Variation of Information

index2 <- list()

for (i in Cn) {
  Clust1 <- DfCluster[, i]
  
  for (j in Cn){
    
    Clust2 <- DfCluster[, j]
    
    ## Compute Index
    
    
    VI <- vi.dist(Clust1,Clust2, parts = T, base = 3217)
    
    index2[paste(i,j,sep="_")]<- VI
  }
  
}  


tabVI <- as.data.frame(unlist(index2, use.names = TRUE))
colnames(tabVI) <- "VI"
tabVI<-tabVI %>%
  rownames_to_column() %>%
  separate(rowname, sep="_", into = c("I", "J")) 
#spread(key = I, value = VI)

#Add simul
ResultSimulVI2 <- ResultSimulVI%>%
  rename(I = Algo) %>%
  select(-sd) %>%
  mutate(J = "Z_Mean50RandomSimul")%>%
  rename(VI = mean) %>% select(I,J,everything())

tabVI <- rbind(tabVI, ResultSimulVI2)



##### Plot heat map

ggplot(tabJaccards, aes(x = I, y= J, fill = Jaccard))+
  geom_tile()+
  scale_fill_gradient2(low ="grey" , high = "orange")+
  geom_text(aes(label = round(Jaccard, 3))) + 
  scale_y_discrete(labels=c(  "fgreedy", "infomap","louvain","walktrap" ,"Moyenne\n50 simulations\naléatoires") ) + 
  labs( subtitle = "Indice de Jaccard", caption = "P.Gourdon 2020",
        y= NULL, x = NULL)+ 
  theme(text = element_text(size = 16))

ggsave( filename = "OUT/Jaccard_Cluster_em2nw_biparti_EUCICOP.pdf", width = 11.7, height = 8.3, units = "in" )

ggplot(tabVI, aes(x = I, y= J, fill = VI))+
  geom_tile()+
  scale_fill_gradient2(low ="grey" , high = "orange", limit = c(0,1) )+
  scale_y_discrete(labels=c(  "fgreedy", "infomap","louvain","walktrap" ,"Moyenne\n50 simulations\naléatoires") ) + 
  geom_text(aes(label = round(VI, 3))) +
  labs(subtitle = "Variation de l'information (VI)", caption = "P.Gourdon 2020",
       y= NULL, x = NULL) +
  theme(text = element_text(size = 16))

ggsave(filename = "OUT/VI_Cluster_em2nw_biparti_eucicop.pdf", width = 11.7, height = 8.3, units = "in" )

ggplot(tabArandi, aes(x = I, y= J, fill = Arandi))+
  geom_tile()+
  scale_fill_gradient2(low ="grey" , high ="orange")+
  geom_text(aes(label = round(Arandi, 3))) + 
  labs(title = "Comparaison du partitionnement selon différents algorithmes", subtitle = "Ajusted Rand Index (Arandi)", caption = "P.Gourdon 2020",
       y= NULL, x = NULL)

meanJaccard <- tabJaccards %>% group_by(tabJaccards$I)%>% summarise_at(.,3,funs(mean(.)))                         
meanArandi <- tabArandi %>% group_by(tabArandi$I)%>% summarise_at(.,3,funs(mean(.)))
meanVI <- tabVI %>% group_by(tabVI$I)%>% summarise_at(.,3,funs(mean(.)))

rm(index, index2, index3, Pairs, ResultSimulArandi, ResultSimulArandi2, ResultSimulJaccard, ResultSimulJaccard2, ResultSimulVI, ResultSimulVI2,
   tabArandiSimul,tabJaccardsSimul,tabVISimul, ValVI,ValJaccard, ValArandi)


### ==== EXPLORATORY MAPPING of THE CLUSTERING OUTCOMES ====


clusterG <- emnw.clustering$igraph_object


# add geometry 

NodesInfos <- NodesInfos %>%   left_join(select(ParticipationURBACT, name = geonameId, lng_GN, lat_GN)) %>% distinct()
# Prepare Tidygraph object
tdClusterG <- as_tbl_graph(clusterG)
namefilter <- tdClusterG  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(name) %>% deframe()

#joint info nodes
tdClusterG <- tdClusterG %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")

## 1. Plot each community of the partitions
# #Louvain
# ggraph(tdClusterG, layout = "fr")+
#   geom_edge_link(alpha = 0.1)+
#   geom_node_point(aes(color = as.factor(V(tdClusterG)$louvain), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
#                   show.legend = FALSE) + 
#   geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
#   facet_nodes(~V(tdClusterG)$louvain)
# 
# 
# 
# #Fast Greedy
ggraph(tdClusterG, layout = "fr")+
  geom_edge_link(alpha = 0.1)+
  geom_node_point(aes(color = as.factor(V(tdClusterG)$fgreedy), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
                  show.legend = FALSE) +
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  facet_nodes(~V(tdClusterG)$fgreedy)



#### 2. PLot cities community with geographical coordinates

## Create an SF object
geoCluster <- dfclustering %>% left_join(NodesInfos)

geoCluster <- geoCluster %>% filter(type == TRUE) %>% select(-Acro, -Date)

geoCluster <- geoCluster %>% filter(!is.na(lng_GN))

geoClusterSF <- st_as_sf(geoCluster, coords = c("lng_GN","lat_GN"), crs = 4326) %>% st_transform(3035)
class(geoClusterSF)
geoClusterSF <- geoClusterSF %>% mutate_at(.vars=c(3:6), ~as.character(.))


### ==== CHOOSE AND MAP 1 PARTITION  ====
# Transfo population
# filter pop
geoClusterSF <- geoClusterSF %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 5000)

length(unique(geoClusterSF$fgreedy))

TopClusterNCities <- geoClusterSF %>% group_by(fgreedy)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopClusterNCities$N)

TopClusterNCities <- TopClusterNCities %>% filter(N>5) %>% select(fgreedy) %>% deframe()
## MAP community EM2NW
geoClusterTop <- geoClusterSF %>% filter(fgreedy %in% TopClusterNCities)

n<- length(unique(geoClusterTop$fgreedy))
# fgreedy with facet 

library(randomcoloR)

palette <- distinctColorPalette(n)


names(palette) <- sort(as.numeric(unique(geoClusterTop$fgreedy)))

palette
palette[16] <- "black"
Nfgreedy <- geoClusterTop %>% group_by(fgreedy) %>% summarise(Ncities = n()) %>% st_drop_geometry()

ggplot()+
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.2) +
  geom_sf(data =geoClusterTop, aes(fill = fgreedy), alpha= 0.8, color = "black", shape = 21, stroke = 0.2 ,show.legend = FALSE)+
  scale_color_manual(values =palette)+
  labs(y="", x ="", 
       caption = "Note : localités de plus de 5000 habitants seulement.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  # geom_sf_text(aes(label = label, check_overlap = TRUE), size = 2)+ 
  facet_wrap(~ fgreedy)+
  geom_text( data    = Nfgreedy,
             mapping = aes(x = Inf, y = Inf, label = paste("N = ", Ncities, sep = "")),
             hjust   = 1.2,
             vjust   = 2, size = 3)+
  theme(panel.grid = element_blank(), 
        line = element_blank(), 
        rect = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

ggsave( filename = "OUT/fgreedy_URBACTCities_FacetMap_emnw.pdf", width = 11.7, height = 8.3, units = "in" )

# Louvain on synthesis map
geoClusterSF2 <- geoClusterSF %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 20000)

TopClusterNCities2 <- geoClusterSF2 %>% group_by(fgreedy)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopClusterNCities2$N)

TopClusterNCities2 <- TopClusterNCities2 %>% filter(N>10) %>% select(fgreedy) %>% deframe()

palette2 <- palette[names(palette) %in% as.numeric(TopClusterNCities2)]
geoClusterTop2 <- geoClusterSF2 %>% filter(fgreedy %in% TopClusterNCities2)

myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(geoClusterTop2$PopAdmin11)
s[[1]]
bks <- c(s[[1]],s[[4]], 2000000, s[[6]] )
lbs <-  c("20 000",  "250 000",  "2M", "9 M" )
palcol <- palette2

citiesMembers2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "white", color = "#bfbfbf", size = 0.5) +
  geom_sf(data = geoClusterTop2 ,
          mapping = aes(size = PopAdmin11,  colour = fgreedy), show.legend = "point", alpha = 0.9) +
  scale_color_manual( values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.2, 15))+
  annotate("text", label = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = geoClusterTop2 %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Communautés (Fast Greedy)") +
  guides(colour=guide_legend(ncol=2),size=guide_legend(ncol=2) )+
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5))

pdf(file = "OUT/FgreedyCities_emnw_URBACT_20kHab.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembers2 
dev.off()


### ====DUAL PROJECTION ====

### FUNCTION 


DblDUALPROJECTION <- function(matrix){
  
  require(dplyr)
  require(igraph)
  # create a list to stock the results
  results <- list()
  
  # adjency matrix of the 2 types of nodes
  Orga <- matrix %*% t(matrix)
  Cities <- t(matrix) %*% matrix
  
  ## adjency matrix of group overlap (directed weighted)
  OrgaPct <- Orga/diag(Orga)
  diag(OrgaPct) <- 0
  
  CitiesPct <- Cities/diag(Cities)
  diag(CitiesPct) <- 0
  
  ## degree Bipartite 
  # Orga
  degreeOrga <- as.data.frame(diag(Orga))
  colnames(degreeOrga) <- "Degree"
  degreeOrga$name <- rownames(degreeOrga)
  
  # Cities
  degreeCities <- as.data.frame(diag(Cities))
  colnames(degreeCities) <- "Degree"
  degreeCities$name <- rownames(degreeCities)
  
  #Set diagonal to zero to avoid loop
  diag(Orga) <- 0
  diag(Cities) <- 0
  
  #Make simple weighted graph
  OrgaW <- graph_from_adjacency_matrix(Orga, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  CitiesW <- graph_from_adjacency_matrix(Cities, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  #######Normalized edge weight by theoretical maximum 
  
  ### get edge list to compute normalization of edge weight
  elistOrg <- as_data_frame(OrgaW, what = "edges")
  elistCities <- as_data_frame(CitiesW, what = "edges")
  
  
  
  ## compute weight normalised 
  
  # orga
  elistOrg <- merge(elistOrg, degreeOrga, all.x = TRUE, by.x = c("from"), by.y = c("name"))
  colnames(elistOrg)[4] <- "degreeFrom" 
  elistOrg <- merge(elistOrg, degreeOrga, all.x = TRUE, by.x = c("to"), by.y = c("name"))
  colnames(elistOrg)[5] <- "degreeTo" 
  
  elistOrg$nWeight <- elistOrg$weight /  apply(elistOrg[,4:5], 1, FUN=min)
  
  elistOrg <- elistOrg[, c(1,2,6)]
  
  # cities
  elistCities <- merge(elistCities, degreeCities, all.x = TRUE, by.x = c("from"), by.y = c("name"))
  colnames(elistCities)[4] <- "degreeFrom" 
  elistCities <- merge(elistCities, degreeCities, all.x = TRUE, by.x = c("to"), by.y = c("name"))
  colnames(elistCities)[5] <- "degreeTo" 
  
  elistCities$nWeight <- elistCities$weight /  apply(elistCities[,4:5], 1, FUN=min)
  
  elistCities <- elistCities[,c(1,2,6)]
  
  # create graph with normalized weight
  
  OrgaNW <- graph_from_data_frame(elistOrg, directed = FALSE)
  
  CitiesNW <- graph_from_data_frame(elistCities, directed = FALSE)
  
  ## Group overlap graph (directed weighted)
  
  OrgaWD <- graph_from_adjacency_matrix(OrgaPct, mode = "directed", weighted = TRUE, diag = FALSE)
  CitiesWD <- graph_from_adjacency_matrix(CitiesPct,  mode = "directed", weighted = TRUE, diag = FALSE)
  
  
  ### COMMUNITY DETECTION
  
  ## simple  Weighted graph
  louvainOrgaW <- cluster_louvain(OrgaW)
  walktrapOrgaW <- walktrap.community(OrgaW, weights = E(OrgaW)$weight)
  louvainCitiesW <- cluster_louvain(CitiesW)
  walktrapCitiesW <- walktrap.community(CitiesW, weights = E(CitiesW)$weight)
  
  ## Normalized weighted graph 
  
  louvainOrgaNW <- cluster_louvain(OrgaNW, weights = E(OrgaNW)$nWeight)
  walktrapOrgaNW <- walktrap.community(OrgaNW, weights = E(OrgaNW)$nWeight)
  louvainCitiesNW <- cluster_louvain(CitiesNW, weights = E(CitiesNW)$nWeight)
  walktrapCitiesNW <- walktrap.community(CitiesNW, weights = E(CitiesNW)$nWeight)
  
  # membership
  
  # membership orga
  V(OrgaW)$louvainW <- louvainOrgaW$membership
  V(OrgaW)$walktrapW <- walktrapOrgaW$membership
  
  membershipsOrgaW <- igraph::as_data_frame(OrgaW, what = "vertices")
  
  V(OrgaNW)$louvainNW <- louvainOrgaNW$membership
  V(OrgaNW)$walktrapNW <- walktrapOrgaNW$membership
  
  membershipsOrgaNW <- igraph::as_data_frame(OrgaNW, what = "vertices")
  
  membershipOrga <- membershipsOrgaW %>% left_join(membershipsOrgaNW)
  
  # membership cities
  V(CitiesW)$louvainW <- louvainCitiesW$membership
  V(CitiesW)$walktrapW <- walktrapCitiesW$membership
  
  membershipsCitiesW <- igraph::as_data_frame(CitiesW, what = "vertices")
  
  V(CitiesNW)$louvainNW <- louvainCitiesNW$membership
  V(CitiesNW)$walktrapNW <- walktrapCitiesNW$membership
  
  membershipsCitiesNW <- igraph::as_data_frame(CitiesNW, what = "vertices")
  
  membershipCities <- membershipsCitiesW %>% left_join(membershipsCitiesNW)
  
  ## Get summary of different clustering in one dataframe
  
  infoclustering <- data.frame( Proj = c(rep("Organisations",4), rep("Cities", 4)),
                                Graphe = c(rep("Weighted",2), rep("Normalized Weighted", 2), rep("Weighted",2), rep("Normalized Weighted", 2)),
                                Algo = c("louvain", "walktrap", "louvain", "walktrap", "louvain", "walktrap", "louvain", "walktrap"),
                                modularity = c(max(louvainOrgaW$modularity), modularity(walktrapOrgaW),
                                               max(louvainOrgaNW$modularity), modularity(walktrapOrgaNW),
                                               max(louvainCitiesW$modularity), modularity(walktrapCitiesW),
                                               max(louvainCitiesNW$modularity), modularity(walktrapCitiesNW)),
                                nclust = c(length(unique(louvainOrgaW$membership)), 
                                           length(unique(walktrapOrgaW$membership)),
                                           length(unique(louvainOrgaNW$membership)), 
                                           length(unique(walktrapOrgaNW$membership)),
                                           length(unique(louvainCitiesW$membership)), 
                                           length(unique(walktrapCitiesW$membership)),
                                           length(unique(louvainCitiesNW$membership)), 
                                           length(unique(walktrapCitiesNW$membership)))
  )
  
  results[["membership_Orga"]] <- membershipOrga 
  
  results[["membershipCities"]] <-   membershipCities
  
  results[["info_clustering"]] <- infoclustering 
  
  results[["igraph_objects"]] <- list(OrgaW, OrgaNW, OrgaWD, CitiesW, CitiesNW, CitiesWD)
  
  return(results)
  
  
  
}




em2nwDUAL <- DblDUALPROJECTION(em2nw)

em2nwDUAL$info_clustering

table(em2nwDUAL$membershipCities$louvainNW)
table(em2nwDUAL$membership_Orga$louvainNW)

saveRDS(em2nwDUAL, "DataProd/DualProj_em2nw_urbact.rds")

emnwDUAL <- DblDUALPROJECTION(emnw)

emnwDUAL$info_clustering

saveRDS(emnwDUAL, "DataProd/DualProj_emnw_urbact.rds")


##EXPLORE COMMUNITIES OF EM2NW DUAL PROJECTION

# Comparing classification for cities

Memberemnw <- readRDS("DataProd/Membership_emnw_EUCICOP.rds") 
Memberem2nwCities <- Memberem2nw %>% filter(type ==TRUE)
Memberem2nwCities <- Memberem2nwCities %>% left_join(NodesInfos)

MemberemnwCitiesDUAL <- emnwDUAL$membershipCities 

Memberem2nwCitiesALL <- Memberem2nwCities %>% left_join(Memberem2nwCitiesDUAL)
Memberem2nwCitiesALL <- Memberem2nwCitiesALL %>% left_join(NodesInfos)

## JACCARD

# between 2 dualproj Louvain
Pairs <- std.ext(Memberem2nwCitiesALL$louvainNW,Memberem2nwCitiesALL$louvainW)

clv.Jaccard(Pairs) # 0.42


# Between louvain NW dual proj and louvain on bipartite

Pairs <- std.ext(Memberem2nwCitiesALL$louvainNW,Memberem2nwCitiesALL$louvain)

clv.Jaccard(Pairs)# 0.38

# same computation but for cities with more than 5000 inh

foo <- Memberem2nwCitiesALL %>% filter(PopAdmin11>10000)
Pairs <- std.ext(foo$louvainNW,foo$louvain)
clv.Jaccard(Pairs)#0.41

### MAP and Esplore Cities Dual Proj Louvain NW (normalized Weighted)

table(Memberem2nwCitiesDUAL$louvainNW)



MemberemnwCitiesDUAL <- MemberemnwCitiesDUAL %>% left_join(NodesInfos) %>% filter(!is.na(lng_GN))

MemberemnwCitiesDUAL_sf <- st_as_sf(MemberemnwCitiesDUAL, coords = c("lng_GN", "lat_GN"), crs = 4326) %>% st_transform(3035)
class(MemberemnwCitiesDUAL_sf)
MemberemnwCitiesDUAL_sf  <- MemberemnwCitiesDUAL_sf  %>% mutate_at(.vars=c(2:5), ~as.character(.))

MemberemnwCitiesDUAL_sf  <- MemberemnwCitiesDUAL_sf  %>%
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

library(randomcoloR)

paletteDUAL <- distinctColorPalette(length(unique(MemberemnwCitiesDUAL_sf$louvainNW)))

names(paletteDUAL) <- sort(as.numeric(unique(MemberemnwCitiesDUAL_sf$louvainNW)))

## Map 



myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(MemberemnwCitiesDUAL_sf$PopAdmin11)
s[[1]]
bks <- c(5000,20000, 50000, 150000, 1000000, 9000000)
lbs <-  c("5 000" ,"20 000","50 000",  "150 000",  "1M", "9 M" )
palcol <- paletteDUAL

citiesMembersDUALprojNW <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = MemberemnwCitiesDUAL_sf ,
          mapping = aes(size = PopAdmin11,  colour = louvainNW), show.legend = "point", alpha = 0.9) +
  scale_color_manual( values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.05, 10))+
  annotate("text", label = "Note : Louvain sur le graphe ville-ville, valué et normalisé.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = MemberemnwCitiesDUAL_sf %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Communautés (Louvain)") +
  guides(colour=guide_legend(ncol=2),size=guide_legend(ncol=2) )+
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5))



pdf(file = "OUT/LouvainCitiesDUALnw_emnw_URBACT.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembersDUALprojNW
dev.off()

## Community with kernel density
MemberemnwCitiesDUAL_sf <- MemberemnwCitiesDUAL_sf %>% mutate(X = st_coordinates(.)[,1], Y = st_coordinates(MemberemnwCitiesDUAL_sf)[,2])
TopLouvainNCities <- MemberemnwCitiesDUAL_sf %>% group_by(louvainNW)%>% summarise(N= n()) %>%st_drop_geometry() 
TopLouvainNCities2 <- TopLouvainNCities %>% filter(N>3) %>% select(louvainNW) %>% deframe()

MemberemnwCitiesDUAL_sf_filtered <- MemberemnwCitiesDUAL_sf %>% filter(louvainNW %in% TopLouvainNCities2)


palcol <- paletteDUAL[names(paletteDUAL) %in% as.numeric(TopLouvainNCities2)]

citiesMembersDUALprojNW_kernel <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = MemberemnwCitiesDUAL_sf_filtered , show.legend = "point", alpha = 0) +
  geom_density_2d(data = MemberemnwCitiesDUAL_sf_filtered , aes(x = X, y = Y, color = louvainNW),size = 0.1) +
  scale_color_manual( values = palcol)+
  annotate("text", label = "Note : Louvain sur le graphe ville-ville, valué et normalisé. Les 4 communautés de moins de 3 villes ont été retirées.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = MemberemnwCitiesDUAL_sf_filtered %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Communautés (Louvain)\nDensité de Kernel du semis de villes") +
  guides(colour=guide_legend(ncol=2, override.aes = list(size=2)),
         size=guide_legend(ncol=2))+
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5))



pdf(file = "OUT/LouvainCitiesDUALnw_emnw_URBACT_Kernel.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembersDUALprojNW_kernel
dev.off()

## Facet representation

TopLouvainNCities <- MemberemnwCitiesDUAL_sf %>% group_by(louvainNW)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopLouvainNCities$N)


library(ggalt)
library(spats)
test <- ggplot()+
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data =MemberemnwCitiesDUAL_sf , aes(color = louvainNW), alpha= 0, show.legend = FALSE, size =0.1)+
  scale_color_manual( values = palcol)+
  labs(y="", x ="", 
       caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  # geom_sf_text(aes(label = label, check_overlap = TRUE), size = 2)+ 
  facet_wrap(~ louvainNW)+
  geom_text( data    = TopLouvainNCities,
             mapping = aes(x = Inf, y = Inf, label = paste("N = ", N, sep = "")),
             hjust   = 1.2,
             vjust   = 2, size = 3)


test+ 
  geom_density_2d(data = MemberemnwCitiesDUAL_sf, aes(x = X, y = Y, color = louvainNW)) 


################################### ===== Reconstructing the network for the 3 phase ====
## four configurations of network (3 phases and all phases)

# All period
emnwDUAL <- readRDS("DataProd/DualProj_emnw_urbact.rds")
gCities <-emnwDUAL$igraph_objects[[4]]
gCitiesTd <- as_tbl_graph(gCities)
#Phase I
emnwDUAL_I <- DblDUALPROJECTION(emnw1)
gCitiesI <- emnwDUAL_I$igraph_objects[[4]]
gCitiesITd <- as_tbl_graph(gCitiesI)
#Phase II
emnwDUAL_II <- DblDUALPROJECTION(emnw2)
gCitiesII <- emnwDUAL_II$igraph_objects[[4]]
gCitiesIITd <- as_tbl_graph(gCitiesII)
#Phase III
emnwDUAL_III <- DblDUALPROJECTION(emnw3)
gCitiesIII <- emnwDUAL_III$igraph_objects[[4]]
gCitiesIIITd <- as_tbl_graph(gCitiesIII)
Sweight= 1

Sweightall = 2
### Filter weight 1
# Label
NodesInfos <- NodesInfos %>% mutate(label = recode(label, "Mestre" = "Venice"))

gCitiesTd <- gCitiesTd %>% activate(edges)%>% filter(weight > Sweightall) %>% activate(nodes) %>%  filter(!node_is_isolated()) %>%
  left_join(NodesInfos) %>% mutate(louvainFilter = group_louvain(weights = weight), degreeFilterNW = centrality_degree(weights = NULL)) %>%
  activate(edges) %>% rename( "Poids\n(nb de projets en commun)" = weight)


gCitiesITd <- gCitiesITd %>% activate(edges)%>% filter(weight > Sweight) %>% activate(nodes) %>%  filter(!node_is_isolated()) %>%
  left_join(NodesInfos) %>% mutate(louvainFilter = group_louvain(weights = weight), degreeFilterNW = centrality_degree(weights = NULL)) %>%
  activate(edges) %>% rename( "Poids\n(nb de projets en commun)" = weight)


gCitiesIITd <- gCitiesIITd %>% activate(edges)%>% filter(weight > Sweight) %>% activate(nodes) %>%  filter(!node_is_isolated()) %>%
  left_join(NodesInfos) %>% mutate(louvainFilter = group_louvain(weights = weight), degreeFilterNW = centrality_degree(weights = NULL)) %>%
  activate(edges) %>% rename( "Poids\n(nb de projets en commun)" = weight)


gCitiesIIITd <- gCitiesIIITd %>% activate(edges)%>% filter(weight > Sweight) %>% activate(nodes) %>%  filter(!node_is_isolated()) %>%
  left_join(NodesInfos) %>% mutate(louvainFilter = group_louvain(weights = weight), degreeFilterNW = centrality_degree(weights = NULL)) %>%
  activate(edges) %>% rename( "Poids\n(nb de projets en commun)" = weight)


## Max degree I II and III

degMax <- max(V(gCitiesITd)$degreeFilterNW, V(gCitiesIITd)$degreeFilterNW, V(gCitiesIIITd)$degreeFilterNW)

degMin <- min(V(gCitiesITd)$degreeFilterNW, V(gCitiesIITd)$degreeFilterNW, V(gCitiesIIITd)$degreeFilterNW)

weightMax <- max(E(gCitiesITd)$`Poids\n(nb de projets en commun)`,
                 E(gCitiesIITd)$`Poids\n(nb de projets en commun)`, 
                 E(gCitiesIIITd)$`Poids\n(nb de projets en commun)`)

weightMin <- min(E(gCitiesITd)$`Poids\n(nb de projets en commun)`,
                 E(gCitiesIITd)$`Poids\n(nb de projets en commun)`, 
                 E(gCitiesIIITd)$`Poids\n(nb de projets en commun)`)
# change weight name (problem because impossible to change width lab in ggraph)

# FilterG <-  FilterG %>% activate(edges) %>% rename( "Poids\n(nb villes membres\nen commun)" = weight)

## Plotting networks


g1<-ggraph(gCitiesITd, layout = 'fr') + 
  geom_edge_hive( aes(width = `Poids\n(nb de projets en commun)`), alpha = 0.3) +
  geom_node_point(aes(size = degreeFilterNW, color = as.character(louvainFilter))) + 
  geom_node_text(aes(label = label),repel = TRUE, size = 3) +
  scale_size(range = c(1,8), limits = c(degMin, degMax))+
  scale_edge_width(range = c(0.5, 3), limits = c(weightMin, weightMax))+
  labs(title = "Phase I (2003 - 2006)",
       caption = "Algorithme de spatialisation : fr\nSous-graphe poids supérieur à 1",
       size =  "Degré non pondéré", color = "Communautés (louvain)")+
    theme_graph(base_family = 'Helvetica')  +
    theme(legend.position="bottom", legend.box = "vertical") + 
  annotate("label",
           label = paste("Densité : ", round(graph.density(gCitiesITd),3),
                         "\nNb Villes : ", vcount(gCitiesITd),
                         "\nNb liens : ", ecount(gCitiesITd),
                         "\nNb composantes : ", components(gCitiesITd)$no, 
                         "\nTransitivité : ", round(transitivity(gCitiesITd, type = "undirected", weights = NULL),3), 
                         "\nModularité : ", round(modularity(gCitiesITd, membership = V(gCitiesITd)$louvainFilter),2),
                         sep = ""), 
           x = Inf, 
           y = Inf,
           hjust = 1,
           vjust =1,
           alpha =0.4,
           size = 2.5)

g1 

palcol <- distinctColorPalette(length(unique(V(gCitiesIITd)$louvainFilter)))

g2<-ggraph(gCitiesIITd, layout = 'fr') + 
  geom_edge_hive( aes(width = `Poids\n(nb de projets en commun)`), alpha = 0.3) +
  geom_node_point(aes(size = degreeFilterNW, color = as.character(louvainFilter))) + 
  geom_node_text(aes(label = label),repel = TRUE, size = 3) +
  scale_size(range = c(1,8), limits = c(degMin, degMax))+
  scale_color_manual(values = palcol)+
  scale_edge_width(range = c(0.5, 3), limits = c(weightMin, weightMax))+
  labs(title = "Phase II (2007 - 2013)",
       caption = "Algorithme de spatialisation : fr\nSous-graphe poids supérieur à 1",
       size =  "Degré non pondéré", color = "Communautés (louvain)")+
  theme_graph(base_family = 'Helvetica')  +
  theme(legend.position="bottom", legend.box = "vertical") +
  annotate("label",
           label = paste("Densité : ", round(graph.density(gCitiesIITd),3),
                         "\nNb Villes : ", vcount(gCitiesIITd),
                         "\nNb liens : ", ecount(gCitiesIITd),
                         "\nNb composantes : ", components(gCitiesIITd)$no, 
                         "\nTransitivité : ", round(transitivity(gCitiesIITd, type = "undirected", weights = NULL),3), 
                         "\nModularité : ", round(modularity(gCitiesIITd, membership = V(gCitiesIITd)$louvainFilter),2),
                         sep = ""), 
           x = Inf, 
           y = Inf,
           hjust = 1,
           vjust =1,
           alpha =0.4,
           size = 2.5)
g2 




palcol <- distinctColorPalette(length(unique(V(gCitiesIIITd)$louvainFilter)))  
g3<-ggraph(gCitiesIIITd, layout = 'fr') + 
  geom_edge_hive( aes(width = `Poids\n(nb de projets en commun)`), alpha = 0.3) +
  geom_node_point(aes(size = degreeFilterNW, color = as.character(louvainFilter))) + 
  geom_node_text(aes(label = label),repel = TRUE, size = 3) +
  scale_size(range = c(1,8), limits = c(degMin, degMax))+
  scale_edge_width(range = c(0.5, 3), limits = c(weightMin, weightMax))+
  scale_color_manual(values = palcol)+
  labs(title = "Phase III (2014 - 2020)",
       caption = "Algorithme de spatialisation : fr\nSous-graphe poids supérieur à 1",
       size =  "Degré non pondéré", color = "Communautés (louvain)")+
  theme_graph(base_family = 'Helvetica')  +
  theme(legend.position="bottom", legend.box = "vertical") +
  annotate("label",
           label = paste("Densité : ", round(graph.density(gCitiesIIITd),3),
                         "\nNb Villes : ", vcount(gCitiesIIITd),
                         "\nNb liens : ", ecount(gCitiesIIITd),
                         "\nNb composantes : ", components(gCitiesIIITd)$no, 
                         "\nTransitivité : ", round(transitivity(gCitiesIIITd, type = "undirected", weights = NULL),3),
                         "\nModularité : ", round(modularity(gCitiesIIITd, membership = V(gCitiesIIITd)$louvainFilter),2),
                         sep = ""), 
           x = Inf, 
           y = Inf,
           hjust = 1,
           vjust =1,
           alpha =0.4,
           size = 2.5)
g3 

gall<-ggraph(gCitiesTd, layout = 'fr') + 
  geom_edge_hive( aes(width = `Poids\n(nb de projets en commun)`, alpha =  `Poids\n(nb de projets en commun)`))+
  geom_node_point(aes(size = degreeFilterNW, color = as.character(louvainFilter))) + 
  geom_node_text(aes(label = label),repel = TRUE, size = 3) +
  scale_size(range = c(1,8))+
  scale_edge_alpha(range = c(0.1, 0.5))+
  scale_edge_width(range = c(0.5, 3))+
  labs(title = "Ensemble des phase (2003-2019)",
       caption = "Algorithme de spatialisation : fr\nSous-graphe poids supérieur à 2",
       size =  "Degré non pondéré", color = "Communautés (louvain)")+
  theme_graph(base_family = 'Helvetica')  +
  theme(legend.position="bottom", legend.box = "horizontal")+
  guides(color=guide_legend(nrow=2), size=guide_legend(nrow=2))
gall 

## Grid with two plots

grid1 <- g1 + g2 + plot_layout(ncol = 1, heights = c(1.5, 2))

#save and export

ggsave(g1, filename = "OUT/URBACTLouvain_PhaseI_weight2.pdf", width = 8.3, height = 8.3, units = "in" )

ggsave(g2, filename = "OUT/URBACTLouvain_PhaseII_weight2.pdf", width = 8.3, height = 8.3, units = "in" )

ggsave(g3, filename = "OUT/URBACTLouvain_PhaseIII_weight2.pdf", width = 8.3, height = 8.3, units = "in" )

## basic properties
# ordre et taille
vcount(gCitiesITd)
ecount(gCitiesITd)

vcount(gCitiesIITd)
ecount(gCitiesIITd)

vcount(gCitiesIIITd)
ecount(gCitiesIIITd)

# density
graph.density(gCitiesITd)
gCitiesITd

# check density
(183*2)/ (65*64)
graph.density(gCitiesIITd)

graph.density(gCitiesIIITd)
modularity(gCitiesITd, membership = V(gCitiesITd)$louvainFilter)

# N components
components(gCitiesITd)$no
components(gCitiesIITd)$no
components(gCitiesIIITd)$no

# transitivity
transitivity(gCitiesITd, type = "undirected", weights = NULL)
transitivity(gCitiesIITd, type = "undirected", weights = NULL)
transitivity(gCitiesIIITd, type = "undirected", weights = NULL)


### Cities that are in the 3 subgraph

NodesI <- gCitiesITd %>% activate(nodes) %>% fortify.tbl_graph()

NodesII <- gCitiesIITd %>% activate(nodes) %>% fortify.tbl_graph()

NodesIII <- gCitiesIIITd %>% activate(nodes) %>% fortify.tbl_graph()


## Regional structure of communities
RegionLouvainI<-table(NodesI$louvainFilter,NodesI$subregion)
prop.table(RegionLouvainI,1)

RegionLouvainII <- table(NodesII$louvainFilter,NodesII$subregion)
prop.table(RegionLouvainII,1)


table(NodesIII$louvainFilter) #  check communities 1, 2, 3, 4, and 5
RegionLouvainIII <- table(NodesIII$louvainFilter,NodesIII$subregion)
prop.table(RegionLouvainIII,1)



## Admin status structure of communities
AdminLouvainI<-table(NodesI$louvainFilter,NodesI$adminLevel)
prop.table(AdminLouvainI,1)

AdminLouvainII <- table(NodesII$louvainFilter,NodesII$adminLevel)
prop.table(AdminLouvainII,1)

AdminLouvainIII <- table(NodesIII$louvainFilter,NodesIII$adminLevel)
prop.table(AdminLouvainIII,1)

### change in "backcloth" == cities that appear and disappear between each sub-graph

Nodes3sub <- bind_rows(NodesI,NodesII, NodesIII, .id = "Phase")

CitiesPresence <- Nodes3sub %>% group_by(label, Phase) %>% summarise(Nphase = n())%>% ungroup %>% 
  mutate(Phase = paste0("Phase_", Phase))%>%
  pivot_wider(names_from = Phase, values_from = Nphase) %>% 
  select(label, Phase_1, Phase_2, Phase_3) %>% 
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate(Totpresence = rowSums(across(where(is.numeric)))) %>%
  mutate(Transition1to2 = case_when(Phase_1 == 1 & Phase_2 == 1 ~ "Remain",
                                    Phase_1 == 0 & Phase_2 == 1 ~ "Appear",
                                    Phase_1 == 1 & Phase_2 == 0 ~ "Disappear",
                                    TRUE ~ "Still Missing"),
         Transition2to3 = case_when(Phase_2 == 1 & Phase_3 == 1 ~ "Remain",
                                    Phase_2 == 0 & Phase_3 == 1 ~ "Appear",
                                    Phase_2 == 1 & Phase_3 == 0 ~ "Disappear",
                                    TRUE ~ "Still Missing"))

CitiesPresence <- CitiesPresence %>% left_join(select(NodesInfos, label, Country,subregion, adminLevel, population, PopAdmin11))

CitiesPresence$population <- as.numeric(CitiesPresence$population)
CitiesPresence <- CitiesPresence %>% arrange(Phase_3, population)
# n cities with less than 100 k inh
CitiesPresence %>% filter(Phase_1 == 1 & PopAdmin11 < 100000) # phase II = 10 (pop)  and 8 (PopAdmin11)

CitiesPresence %>% filter(Phase_2 == 1 & PopAdmin11 < 100000) # phase II = 17(pop) and 18 (PopAdmin11)

smc<-CitiesPresence %>% filter(Phase_3 == 1 & population < 100000)  # phase III = 26 (pop) and 20 (PopAdmin11)
  
## Count change in backcloth between phases
table(CitiesPresence$Phase_3, CitiesPresence$subregion)
préstable(CitiesPresence$Phase_3, CitiesPresence$adminLevel)

PhaseItoII <- table(CitiesPresence$Transition1to2)

PhaseIItoIII <- table(CitiesPresence$Transition2to3)

table(CitiesPresence$Transition1to2, CitiesPresence$subregion)
table(CitiesPresence$Transition2to3, CitiesPresence$subregion)
CitiesPresenceAlltime <- CitiesPresence %>% filter(Totpresence == 3)
  
  
  
### Check assortativity
g <- as.igraph(gCitiesIIITd) 

assortativity_degree(g)


# strength
V(g)$strength <- strength(g, vids = V(g), loops = F )
degStrength <- assortativity(g , V(g)$strength)

## Pop 


popAssort  <- assortativity(g, V(g)$population, directed = FALSE)

# Pop discrete 
summary(NodesInfos$PopAdmin11)
NodesInfos <- NodesInfos %>% filter(!is.na(PopAdmin11))
labelClass <- c( "4.Small City", "3.Medium-sized Small City","2.Large City",  "1.Very Large City")



cesDonnees <- NodesInfos$PopAdmin11
min<-min(cesDonnees)
max<- max(cesDonnees)
valBreaks <- c(min, 50000,150000, 500000, max)

NodesInfos$ClassPop11 <- cut(cesDonnees,
                             breaks = valBreaks,
                             labels = labelClass,
                             include.lowest = TRUE,
                             right= FALSE)

NameVertex <- V(g)$name
NodesInfos2 <- NodesInfos %>% filter(name %in% NameVertex)
g <-  set_vertex_attr(g,"PopClass11", index = NodesInfos2$name, value = NodesInfos2$ClassPop11)



Matrixg <- get.adjacency(g) %>% as.matrix()
library(assortnet)
NameMatrix <- data.frame(name= rownames(Matrixg)) %>% left_join(select(NodesInfos, name, ClassPop11))
assortment.discrete(Matrixg ,NameMatrix$ClassPop11 ,weighted=TRUE)

## admin


adminAssort <- assortativity_nominal(g,as.numeric(as.factor(V(g)$adminLevel)))

Matrixg <- get.adjacency(g) %>% as.matrix()

NameMatrix <- data.frame(name= rownames(Matrixg)) %>% left_join(select(NodesInfos, name, adminLevel))
assortment.discrete(Matrixg ,NameMatrix$adminLevel ,weighted=FALSE)

## SubRegion 
subregionAssort <-  assortativity_nominal(g, as.numeric(as.factor(V(g)$subregion)))
Matrixg <- get.adjacency(g) %>% as.matrix()
NameMatrix <- data.frame(name= rownames(Matrixg)) %>% left_join(select(NodesInfos, name, subregion))
assortment.discrete(Matrixg ,NameMatrix$subregion ,weighted=FALSE)


#################################################################################################################################################"


### ==== TEMPORAL REPRESENTATIONS AND ANALYSIS ====


library(bmotif)

Motif_emwn <- mcount(t(emnw),six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)

# Dico top down number of nodes for each motif

V1V2df <- data.frame(motif = c(1:17), V1_V2 = c("11","12","21","31","22","22", "13", "41", "32", "32","32", "32","23", "23", "23", "23","14"))

Motif_emwn <- Motif_emwn %>% left_join(V1V2df)

ggplot(Motif_emwn) + 
  geom_col(aes(x = reorder(as.character(motif), -normalise_sum) , y = normalise_sum, fill = as.character(nodes)) )  + 
  labs(x = "Motifs", y = "Frequence", fill = "Nombre de noeuds des motifs")+ scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.3, by = 0.05))

ggplot(Motif_emwn %>% filter(!nodes == 2)) + 
  geom_col(aes(x = reorder(as.character(motif), -normalise_sizeclass) , y = normalise_sizeclass, fill = as.character(nodes)) )  + 
  labs(x = "Motifs", y = "Frequence", fill = "Nombre de noeuds des motifs")+ scale_x_discrete(guide = guide_axis(n.dodge = 2))+
facet_wrap(~nodes, scales = "free", ncol = 2, nrow = 2)+theme(legend.position =  c(1,0), legend.justification = c(1.2,-0.2))

ggplot(Motif_emwn %>% filter(!nodes == 2 & !nodes == 3 & V1_V2 == "23" |V1_V2 =="32")) + 
  geom_col(aes(x = reorder(as.character(motif), -normalise_levelsize) , y = normalise_levelsize) )  + 
  labs(x = "Motifs", y = "Frequence")+ scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  facet_wrap(~V1_V2, scales = "free")

# Position of cities

Position_emnw <- node_positions(M = t(emnw), six_node = FALSE, weights_method = "none", weights_combine = "none", level = "rows", normalisation = "sum")



TopPositionByCity_emnw <- Position_emnw %>% mutate(city= rownames(.))%>%  pivot_longer(-city,names_to = "Position", values_to = "Frequency")


TopPositionByCity_emnw <- TopPositionByCity_emnw %>% group_by(city) %>% top_n(1)
unique(TopPositionByCity_emnw$Position)

##Package for ecology so in row predator (polinisator) in my case cities, in columns plants (in my case projects)

Motif_emwn1 <- mcount(t(emnw1),six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)

Motif_emwn2 <- mcount(t(emnw2),six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)

Motif_emwn3 <- mcount(t(emnw3),six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)


MotifTemporal <- bind_rows(Motif_emwn1,Motif_emwn2,Motif_emwn3,.id= "Phase")

ggplot(MotifTemporal) + 
  geom_col(aes(x = reorder(as.character(motif), -normalise_sum) , y = normalise_sum, fill = as.character(nodes)) )  + 
  labs(x = "Motifs", y = "Frequence", fill = "Nombre de noeuds des motifs")+
  facet_wrap(~Phase, ncol = 2, nrow = 2,  scales= "free")+ theme(legend.position =  c(1,0), legend.justification = c(1.2,-0.2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggplot(MotifTemporal %>% filter(!nodes == 2)) + 
  geom_col(aes(x = reorder(as.character(motif), -normalise_sizeclass) , y = normalise_sizeclass, fill = as.character(Phase)) , position = "dodge")  + 
  labs(x = "Motifs", y = "Frequence", fill = "Phase du programme URBACT")+ scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  facet_wrap(~nodes, scales = "free", ncol = 2, nrow = 2)+theme(legend.position =  c(1,0), legend.justification = c(1,0))


## Evolution of the top position

# Urbact 1
Position_emnw1 <- node_positions(M = t(emnw1), six_node = FALSE, weights_method = "none", weights_combine = "none", level = "rows", normalisation = "sum")


TopPositionByCity_emnw1 <- Position_emnw1 %>% mutate(city= rownames(.))%>%  pivot_longer(-city,names_to = "Position", values_to = "Frequency")


TopPositionByCity_emnw1 <- TopPositionByCity_emnw1 %>% group_by(city) %>% top_n(1)
unique(TopPositionByCity_emnw1$Position)

#Urbact 2
Position_emnw2 <- node_positions(M = t(emnw2), six_node = FALSE, weights_method = "none", weights_combine = "none", level = "rows", normalisation = "sum")

TopPositionByCity_emnw2<- Position_emnw2 %>% mutate(city= rownames(.))%>%  pivot_longer(-city,names_to = "Position", values_to = "Frequency")


TopPositionByCity_emnw2<- TopPositionByCity_emnw2%>% group_by(city) %>% top_n(1)
unique(TopPositionByCity_emnw2$Position)

#Urbact 3

Position_emnw3 <- node_positions(M = t(emnw3), six_node = FALSE, weights_method = "none", weights_combine = "none", level = "rows", normalisation = "sum")


TopPositionByCity_emnw3<- Position_emnw3 %>% mutate(city= rownames(.))%>%  pivot_longer(-city,names_to = "Position", values_to = "Frequency")


TopPositionByCity_emnw3<- TopPositionByCity_emnw3%>% group_by(city) %>% top_n(1)
unique(TopPositionByCity_emnw3$Position)

# All phases
TopPositionCityTemporal <- bind_rows(TopPositionByCity_emnw1, TopPositionByCity_emnw2, TopPositionByCity_emnw3, .id= "Phase")

SequencePosition <- TopPositionCityTemporal %>% select(-Frequency) %>% group_by(city, Phase) %>%
  mutate(row = row_number())%>% ungroup %>%
  pivot_wider( names_from = Phase, values_from = Position, values_fill = list(Position = NA)) %>% select(-row)

colnames(SequencePosition)<- c("Name", "Phase1","Phase2", "Phase3")

SequencePosition <- SequencePosition %>% mutate_at(vars(2:4), ~replace_na(., "Aucune Participation")) 

SequencePositionCityAllPhase <- SequencePosition %>% filter_at(vars(2:4), all_vars(!.== "Aucune Participation"))
#" Sequence and optimal matching
library(TraMineR)


CityPos.seq <- seqdef(SequencePosition, var = 2:4)

CityPos.seq <- seqdef(SequencePositionCityAllPhase, var = 2:4)
seqstatl(CityPos.seq)
seqdim(CityPos.seq)



seqdplot(CityPos.seq, main  = "Distribution des états", with.legend = T) 
seqfplot(CityPos.seq, main = "Les 10 séquences les plus fréquentes", with.legend = T)




# Decription des classes avec la frequence en SPS

seqtab(CityPos.seq, tlim = 1:10, format="SPS")

# Transition rate

TransitionRate <- seqtrate(CityPos.seq)
TransitionRate



#####2.3.2 / OM CAH



########  CAH



#Definition des couts 
##(cout constant pour LevenshteinII et Trate pour Dynamic Hamming)

couts <- seqsubm(CityPos.seq,method="CONSTANT", cval=100)
#couts <- seqsubm(CityPos.seq,method="TRATE")

#matrice de distances_
##(changer Indel ou sm selon la distance utilisee : LevenshteinII ou Hamming)

seq.om <- seqdist(CityPos.seq, method="OM", indel=1, sm=couts)

# classification choix du nb de classes

seq.agnes <- agnes(as.dist(seq.om), method="ward", keep.diss=FALSE)
# seq.agnes2 <- agnes(as.dist(seq.om), method="single", keep.diss=FALSE)
# seq.agnes3 <- agnes(as.dist(seq.om), method="complete", keep.diss=FALSE)
# 


# 
par(mfrow = c(1, 1))
plot(as.dendrogram(seq.agnes), leaflab= "none")
plot(sort(seq.agnes$height, decreasing=TRUE)[1:20], 
     type="s", xlab="nb de classes", ylab="inertie")



nbcl <- 6
seq.part <- stats::cutree(seq.agnes, nbcl)

seq.part <- factor(seq.part,labels=paste("Classe",1:nbcl,sep='.'))

# library(NbClust)
# NbClust(data = CityPos.seq,diss = distom, distance = NULL, method = "ward.D2",min.nc = 2, max.nc = 15,index = "all")
### State distribution plot


seqplot(CityPos.seq, group=seq.part, type="d" ,
        border=NA, withlegend=T)



### index plot

#ordre <- cmdscale(as.dist(seq.om),k=1)
#ordre <- sort.list 
seqiplot(CityPos.seq, group=seq.part,sortv="from.start", title = "Index plot",
         tlim=0, space=0, border=NA, with.legend=T, yaxis=FALSE )

# distance moyenne des sequences au centre de la classe


meanDistClass <- round(aggregate(disscenter(as.dist(seq.om), group=seq.part), 
                                 list(seq.part),mean)[,-1],1)

#sequence frequency plot

#seqfplot(CityPos.seq[seq.part != "classe.2", ], group=seq.part[seq.part != "classe.2"],
#        title = "Sequence frequency plot", withlegend=T)
#
seqfplot(CityPos.seq, group=seq.part,  with.legend=T)

str(seq.part)
levels(seq.part)
# etat modal de chaque classe

seqmsplot(CityPos.seq, group=seq.part, with.legend=T, main ="Profil Modal")

# entropie par classe

seqHtplot(CityPos.seq, group=seq.part, with.legend=T, main = "Entropie intraclasse")



##### ==== Trajectories of role in the programme ==== #

# nb of projects by start year

ProjectURBACT %>% group_by(Start) %>% summarise(N= n())

# regroup 2007 and 2010 and 2016

library(questionr)
# irec(ParticipationURBACT)

## Recodage de ParticipationURBACT$Start en ParticipationURBACT$Start_rec
ParticipationURBACT$Start_rec <- as.character(ParticipationURBACT$Start)
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2015"] <- "2015-2016"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2008"] <- "2007-2008"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2009"] <- "2009-2010"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2010"] <- "2009-2010"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2007"] <- "2007-2008"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2016"] <- "2015-2016"
ParticipationURBACT$Start_rec[ParticipationURBACT$Start == "2000"] <- "2003-2006"


PartUrbactPhase <- ParticipationURBACT %>% group_by(City.Statut, Start_rec, geonameId) %>% summarise(N= n())

PartUrbactPhase <- PartUrbactPhase %>% group_by(Start_rec, geonameId) %>% mutate(Ndate = n())

PartUrbactPhase <- PartUrbactPhase %>% mutate(StatutPhase = ifelse(Ndate==2, "Lead Partner", City.Statut))

PartUrbactPhase <- PartUrbactPhase %>% select(-N,-Ndate, - City.Statut) %>% distinct() %>% pivot_wider(names_from = Start_rec, values_from = StatutPhase)

PartUrbactPhase <- PartUrbactPhase %>% ungroup %>%  replace(is.na(.), "Aucune Participation")


## Analyse de séquence

library(TraMineR)


CityStat.seq <- seqdef(PartUrbactPhase, var = 2:8)
colnames(CityStat.seq)
Date <- c("0306","0708","0910","13","1516","18", "19")
Date <- c("2003-\n2006", "2007-\n2008", "2009-\n2010" ,"2013",
          "2015-\n2016" ,"2018", "2019" )

colnames(CityStat.seq) <- Date

seqstatl(CityStat.seq)
seqdim(CityStat.seq)


pdf(file = "OUT/dplot_URBACT_seq.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
seqdplot(CityStat.seq, with.legend = TRUE, axes = TRUE)
dev.off()





pdf(file = "OUT/fplot_URBACT_seq.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
seqfplot(CityStat.seq,  with.legend = T, border = NA)
dev.off()

# Decription des classes avec la frequence en SPS

seqtab(CityStat.seq, tlim = 1:10, format="SPS")

# Transition rate

TransitionRate <- seqtrate(CityStat.seq)
TransitionRate

library(kableExtra)
kable(TransitionRate, booktabs = T, "latex") %>%
  # column_spec(column = c(3:11), width = "2cm") %>% 
  kable_styling(position = "center", latex_options = "striped")


#####2.3.2 / OM CAH

### filtrage des séquences ne comprenant qu'une seule participation.
ind2remove <- PartUrbactPhase %>% 
  pivot_longer(-geonameId, names_to = "period", values_to = "state") %>% 
  group_by(geonameId, state) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "state", values_from = "n") %>% 
  arrange(geonameId) %>% 
  filter(`Aucune Participation` == 6 & Partner == 1) %>% 
  select(geonameId) %>% 
  deframe

PartUrbactPhase2 <- PartUrbactPhase %>% filter(!geonameId %in% ind2remove)

CityStat.seq2 <- seqdef(PartUrbactPhase2, var = 2:8)
colnames(CityStat.seq2)
Date <- c("0306","0708","0910","13","1516","18", "19")
Date <- c("2003-\n2006", "2007-\n2008", "2009-\n2010" ,"2013",
          "2015-\n2016" ,"2018", "2019" )

colnames(CityStat.seq2) <- Date

seqstatl(CityStat.seq2)
seqdim(CityStat.seq2)


pdf(file = "OUT/dplot_URBACT_seqfilter.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
seqdplot(CityStat.seq2, with.legend = TRUE, axes = TRUE)
dev.off()





pdf(file = "OUT/fplot_URBACT_seqfilter.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
seqfplot(CityStat.seq2,  with.legend = T, border = NA)
dev.off()

# Decription des classes avec la frequence en SPS

seqtab(CityStat.seq2, tlim = 1:10, format="SPS")

# Transition rate

TransitionRate <- seqtrate(CityStat.seq2)
TransitionRate

library(kableExtra)
kable(TransitionRate, booktabs = T, "latex") %>%
  # column_spec(column = c(3:11), width = "2cm") %>% 
  kable_styling(position = "center", latex_options = "striped")
########  CAH



#Definition des couts 
##(cout constant pour LevenshteinII et Trate pour Dynamic Hamming)

# couts <- seqsubm(CityStat.seq2,method="CONSTANT", cval=50)
couts <- seqsubm(CityStat.seq,method="TRATE")

#matrice de distances_
##(changer Indel ou sm selon la distance utilisee : LevenshteinII ou Hamming)

seq.om <- seqdist(CityStat.seq2, method="OM", indel=1, sm=couts)

# classification choix du nb de classes

seq.agnes <- agnes(as.dist(seq.om), method="ward", keep.diss=FALSE)
# seq.agnes2 <- agnes(as.dist(seq.om), method="single", keep.diss=FALSE)
# seq.agnes3 <- agnes(as.dist(seq.om), method="complete", keep.diss=FALSE)
# 


# 
par(mfrow = c(1, 1))
plot(as.dendrogram(seq.agnes), leaflab= "none")
plot(sort(seq.agnes$height, decreasing=TRUE)[1:20], 
     type="s", xlab="nb de classes", ylab="inertie")



nbcl <- 7
seq.part <- stats::cutree(seq.agnes, nbcl)

seq.part <- factor(seq.part,labels=paste("Classe",1:nbcl,sep='.'))

# library(NbClust)
# NbClust(data = CityStat.seq,diss = distom, distance = NULL, method = "ward.D2",min.nc = 2, max.nc = 15,index = "all")
### State distribution plot



pdf(file = "OUT/dplot_URBACT_seqCAH_filter.pdf", width = 8.3, height = 8.3, pagecentre = FALSE)
seqplot(CityStat.seq2, group=seq.part, type = "d")
dev.off()

### index plot

#ordre <- cmdscale(as.dist(seq.om),k=1)
#ordre <- sort.list 
seqiplot(CityStat.seq2, group=seq.part,sortv="from.start", title = "Index plot",
         tlim=0, space=0, border=NA, with.legend=T, yaxis=FALSE )

# distance moyenne des sequences au centre de la classe


meanDistClass <- round(aggregate(disscenter(as.dist(seq.om), group=seq.part), 
                                 list(seq.part),mean)[,-1],1)

#sequence frequency plot

#seqfplot(CityStat.seq[seq.part != "classe.2", ], group=seq.part[seq.part != "classe.2"],
#        title = "Sequence frequency plot", withlegend=T)
#


pdf(file = "OUT/fplot_URBACT_seqCAH_filter.pdf", width = 8.3, height = 8.3, pagecentre = FALSE)
seqfplot(CityStat.seq2, group=seq.part,  with.legend=T, border = NA)
dev.off()

str(seq.part)
levels(seq.part)
# etat modal de chaque classe

 seqmsplot(CityStat.seq2, group=seq.part, with.legend=T, main ="Profil Modal")

# entropie par classe

seqHtplot(CityStat.seq2, group=seq.part, with.legend=T, main = "Entropie intraclasse")

## get the typo
PartUrbactPhase2$OMclustersCAH <- seq.part
PartUrbactPhase2$OMclustersCAH <- as.integer(PartUrbactPhase2$OMclustersCAH)


OMclusterLabel <- data.frame(OMclustersCAH = unique(PartUrbactPhase2$OMclustersCAH), 
                             OMclustLabel = c("1. Lead Partner en début de séquence,\nparticipations irrégulières", 
                                              "2. Etats de Lead Partner majoritaires\net consécutifs sur l'ensemble de la séquence",
                                              "3. Etat de Lead Partner en fin de séquence,\nparticipations irrégulières",
                                              "4. Participations consécutives et majoritaires",
                                              "5. Participations consécutives\nou répétées en fin de séquence",
                                              "6. Deux ou trois participation\nen milieu de séquence",
                                              "7. Participations en début de séquence,\npuis participations irrégulières"))

PartUrbactPhase2 <- PartUrbactPhase2 %>% left_join(OMclusterLabel)
### MAP and Esplore Cities Dual Proj Louvain NW (normalized Weighted)




MemberemnwCities <- Memberemnw %>% filter(type ==TRUE)
MemberemnwCities <- MemberemnwCities %>% left_join(NodesInfos)

MemberemnwCities <- MemberemnwCities %>%filter(!is.na(lng_GN))

MemberemnwCities_sf <- st_as_sf(MemberemnwCities, coords = c("lng_GN", "lat_GN"), crs = 4326) %>% st_transform(3035)

MemberemnwCities_sf  <- MemberemnwCities_sf  %>%
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

## add OM typo
MemberemnwCities_sf <- MemberemnwCities_sf %>% left_join(select(PartUrbactPhase2,name = geonameId, OMclustersCAH, OMclustLabel))

MemberemnwCities_sf_OM <- MemberemnwCities_sf %>% filter(!is.na(OMclustersCAH)) 



library(randomcoloR)
library(RColorBrewer)

paletteDUAL <- distinctColorPalette(7)

paletteDUAL <- brewer.pal(n = 7, name = "Set2")
names(paletteDUAL) <- unique(MemberemnwCities_sf_OM$OMclustLabel)

## Map 



myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(MemberemnwCities_sf_OM$PopAdmin11)
s[[1]]
bks <- c(2500, 50000, 150000, 300000,1000000, 9000000)
lbs <-  c("2 500" ,"50 000","150 000",  "300 000",  "1M", "9 M" )
palcol <- paletteDUAL

citiesMembersOM <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = MemberemnwCities_sf_OM ,
          mapping = aes(size = PopAdmin11,  colour = OMclustLabel), show.legend = "point", alpha = 0.9) +
  scale_color_manual( values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.05, 10))+
  annotate("text", label = "Distance utilisée pour la classification : Levenshtein I.\nSources : EUCICOP 2019 ; BD URBACT / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = MemberemnwCities_sf_OM %>% filter(PopAdmin11 > 300000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + 
  labs(color = "Types de séquences de participation\ndans le programme URBACT\nentre 2000 et 2019\n(Optimal Matching)") +
  guides(colour=guide_legend(ncol=1),size=guide_legend(ncol=2) )+
  theme_void() +
  theme(legend.position =  c(0.20, 0.60), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5))



pdf(file = "OUT/URBACTcities_OM_CAH7.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembersOM
dev.off()

table(MemberemnwCities_sf_OM$OMclustersCAH, MemberemnwCities_sf_OM$adminLevel)
