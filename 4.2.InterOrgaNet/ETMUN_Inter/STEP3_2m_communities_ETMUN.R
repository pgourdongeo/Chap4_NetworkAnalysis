############################### ETMUN bipartite Networks STEP 3 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe ETMUN (affiliation) biparti. 
#               Détection de communautés 
#
# 
############################################################################## PG juin 2020

### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")

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

## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_europe.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", 
                       stringsAsFactors = F)

## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")
DbCityInfos <- DbCity %>% st_drop_geometry()

# Infos on Nodes (STEP 2)

NodesInfos <- readRDS("DataProd/NodesInfos_ETMUN.RDS")

## shapes

## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)


# Get long_lat 
library(skimr)
NodesInfos <- NodesInfos %>% left_join(select(MembershipEtmun, name = geonameId, lng_GN,lat_GN))%>% distinct()%>%
  left_join(select(DbCityInfos,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11), by = c("name" = "geonameId"))

## Matrices from STEP1  (REMIND WWCAM is removed a priori because only one european city connected to the whole graph)

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/Etmun_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/Etmun_em2nw.rds")

## Matrices from STEP 2 (removed Covenant of Mayors, CLimate Alliance, and EWT)

# matrix with all edges but unweighted (unique member localitities in each association)
emnwF <- readRDS("DataProd/emnwF.RDS")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nwF <- readRDS("DataProd/em2nwF.RDS")


nrow(emnw)
ncol(emnw)
### ==== IGRAPH TRANSFORMATION ====



emnw.g <- graph.incidence(emnw, directed = F)

em2nw.g <-  graph.incidence(em2nw, directed = F)

emnwF.g <- graph.incidence(emnwF, directed = F)

em2nwF.g <- graph.incidence(em2nwF, directed = F)





### ==== COMMUNITY DETECTIONS ON WHOLE GRAPH ====


clustering <- function(g){
  require(igraph)
  results <- list()
  
  #partition
  louvain <- cluster_louvain(g)
  infomap <- infomap.community(g, nb.trials = 15)
  greedy <- cluster_fast_greedy(g)
  walktrap <- walktrap.community(g, steps = 6)
  
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

emnwF.clustering <-  clustering(emnwF.g)

em2nwF.clustering <- clustering(em2nwF.g)

emnw.clustering$info_clustering
em2nw.clustering$info_clustering
emnwF.clustering$info_clustering
em2nwF.clustering$info_clustering

## Check modularity

Modularity <- bind_rows(emnw.clustering$info_clustering,
                        em2nw.clustering$info_clustering,
                        emnwF.clustering$info_clustering,
                        em2nwF.clustering$info_clustering, .id = "Matrix")

Modularity <- Modularity %>% mutate(Matrix = case_when(Matrix == 1 ~ "emnw", Matrix == 2 ~ "em2nw" ,Matrix == 3 ~ "emnwF", Matrix == 4 ~ "em2nwF"))



g1 <- ggplot(Modularity) + 
  geom_col(aes(x = Algo, y = modularity, fill = Algo)) + 
  facet_wrap(~Matrix) + 
  labs(x = "", y = "modularité", fill = "Algorithme", subtitle = "Modularité")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), text = element_text(size = 16))

g2 <- ggplot(Modularity) + 
  geom_col(aes(x = Algo, y = nclust, fill = Algo), show.legend = FALSE) + 
  facet_wrap(~Matrix)+ 
  labs(x = "", y = "nombre de communautés", subtitle = "Nombre de communautés", caption = "Source : ETMUN 2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), text = element_text(size = 16))
 

gridMod <- g1 + g2                                                                                                                         
ggsave(gridMod, filename = "OUT/Modularity_CombinedApproach_4matrix.pdf", width = 11.7, height = 8.3, units = "in" )

### ==== PAIWISE COMPARISON OF CLUSTERING ====
### Dataframe Classification
dfclustering <- em2nwF.clustering$memberships
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

ggsave( filename = "OUT/Jaccard_Cluster_em2nw_bipati.pdf", width = 11.7, height = 8.3, units = "in" )

ggplot(tabVI, aes(x = I, y= J, fill = VI))+
  geom_tile()+
  scale_fill_gradient2(low ="grey" , high = "orange", limit = c(0,1) )+
  scale_y_discrete(labels=c(  "fgreedy", "infomap","louvain","walktrap" ,"Moyenne\n50 simulations\naléatoires") ) + 
  geom_text(aes(label = round(VI, 3))) +
  labs(subtitle = "Variation de l'information (VI)", caption = "P.Gourdon 2020",
                                            y= NULL, x = NULL) +
  theme(text = element_text(size = 16))

ggsave(filename = "OUT/VI_Cluster_em2nw_bipati.pdf", width = 11.7, height = 8.3, units = "in" )

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


clusterG <- em2nwF.clustering$igraph_object



# Prepare Tidygraph object
tdClusterG <- as_tbl_graph(clusterG)
namefilter <- tdClusterG  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(name) %>% deframe()

#joint info nodes
tdClusterG <- tdClusterG %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")

## 1. Plot each community of the partitions
#Louvain
ggraph(tdClusterG, layout = "fr")+
  geom_edge_link(alpha = 0.1)+
  geom_node_point(aes(color = as.factor(V(tdClusterG)$louvain), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
                  show.legend = FALSE) + 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  facet_nodes(~V(tdClusterG)$louvain)

#WalkTrap
ggraph(tdClusterG, layout = "fr")+
  geom_edge_link(alpha = 0.1)+
  geom_node_point(aes(color = as.factor(V(tdClusterG)$walktrap), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
                  show.legend = FALSE) + 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  facet_nodes(~V(tdClusterG)$walktrap)

#Fast Greedy
ggraph(tdClusterG, layout = "fr")+
  geom_edge_link(alpha = 0.1)+
  geom_node_point(aes(color = as.factor(V(tdClusterG)$fgreedy), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
                  show.legend = FALSE) + 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  facet_nodes(~V(tdClusterG)$fgreedy)

#Infomap
ggraph(tdClusterG, layout = "fr")+
  geom_edge_link(alpha = 0.1)+
  geom_node_point(aes(color = as.factor(V(tdClusterG)$infomap), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
                  show.legend = FALSE) + 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  facet_nodes(~V(tdClusterG)$infomap)


#### 2. PLot cities community with geographical coordinates

## Create an SF object
geoCluster <- dfclustering %>% left_join(NodesInfos)

geoCluster <- geoCluster %>% filter(type == TRUE) %>% select(-Acro, -Date)


geoClusterSF <- st_as_sf(geoCluster, coords = c("lng_GN", "lat_GN"), crs = 4326)

geoClusterSF <- geoClusterSF %>% mutate_at(.vars=c(3:6), ~as.character(.))

## mApping

# Louvain
ggplot(geoClusterSF)+
  geom_sf(aes(color = louvain))+ 
  # geom_sf_text(aes(label = label, check_overlap = TRUE), size = 2)+ 
  facet_wrap(~ louvain)

# Infomap
ggplot(geoClusterSF)+
  geom_sf(aes(color = infomap))+ 
  # geom_sf_text(aes(label = label, check_overlap = TRUE), size = 2)+ 
  facet_wrap(~ infomap)

# Walktrap
ggplot(geoClusterSF)+
  geom_sf(aes(color = walktrap))+ 
  # geom_sf_text(aes(label = label), size = 2)+ 
  facet_wrap(~ walktrap)

# Fast Greedy
ggplot(geoClusterSF)+
  geom_sf(aes(color = fgreedy))+ 
  # geom_sf_text(aes(label = label), size = 2)+ 
  facet_wrap(~ fgreedy)



### ==== CREATE A SYNTHESIS OF CLUSTERING WITH KMODES ====

library(klaR)

k <- round(median(Modularity$nclust),0)
k = 8
KmodesDf <- dfclustering %>% mutate_at(.vars=c(3:6), ~as.character(.))
KmodesResults <- kmodes(KmodesDf[,3:6], k)

KmodesDf$Kmodes <- KmodesResults$cluster
KmodesResults$withindiff
KmodesResults$size




## Check modularity of the result
clusterG  <-  set_vertex_attr(clusterG,"Kmodes", index = KmodesDf$name, value = KmodesDf$Kmodes)

modularity(clusterG, membership = V(clusterG)$Kmodes)


# Prepare Tidygraph object
tdClusterG <- as_tbl_graph(clusterG)
namefilter <- tdClusterG  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(name) %>% deframe()

#joint info nodes
tdClusterG <- tdClusterG %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")


#Kmodes
ggraph(tdClusterG, layout = "lgl")+
  geom_edge_link(alpha = 0.01)+
  geom_node_point(aes(color = as.character(Kmodes), shape = type, size = 1-type),
                  show.legend = FALSE) + 
  # facet_nodes(~ as.character(Kmodes))+
geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)


### ==== CHOOSE AND MAP 2 PARTITIONS  ====

#### MAP CITIES MEMBERSHIP #####
Memberem2nw <- em2nw.clustering$memberships 

Memberem2nwF <- em2nwF.clustering$memberships


geoCluster2 <- Memberem2nw %>% left_join(NodesInfos)
geoCluster2F <- Memberem2nwF %>% left_join(NodesInfos)

geoCluster2 <- geoCluster2 %>% filter(type == TRUE) %>% dplyr::select(-Acro, -Date)
geoCluster2F <- geoCluster2F %>% filter(type == TRUE) %>% dplyr::select(-Acro, -Date)


geoClusterSF2 <- st_as_sf(geoCluster2, coords = c("lng_GN", "lat_GN"), crs = 4326)
geoClusterSF2 <- geoClusterSF2 %>% mutate_at(.vars=c(3:6), ~as.character(.)) %>% st_transform(crs = 3035)

geoClusterSF2F <- st_as_sf(geoCluster2F, coords = c("lng_GN", "lat_GN"), crs = 4326)
geoClusterSF2F <- geoClusterSF2F %>% mutate_at(.vars=c(3:6), ~as.character(.)) %>% st_transform(crs = 3035)

# saveRDS(Memberem2nw, "DataProd/Membership_em2nw_ETMUN.rds")
# saveRDS(Memberem2nwF, "DataProd/Membership_em2nwF_ETMUN.rds")
# Transfo population

geoClusterSF2 <- geoClusterSF2 %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 5000)

geoClusterSF2F <- geoClusterSF2F %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 5000)


## MAP community EM2NW

# with louvain
myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(geoClusterSF2$PopAdmin11)
s[[1]]
bks <- c(s[[1]],s[[3]],s[[4]], 500000, 2000000, s[[6]] )
lbs <-  c("5 000", "50 000", "200 000", "500 000", "2M", "14 M" )
palcol <- brewer.pal(n = 11, name = "Paired")

citiesMembers2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = geoClusterSF2 ,
          mapping = aes(size = PopAdmin11,  colour = louvain), show.legend = "point", alpha = 0.7) +
  scale_color_manual(breaks = c(1:11), values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.5, 14))+
  annotate("text", label = "Source : ETMUN 2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = geoClusterSF2 %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Communautés (Louvain)") +
  theme_void() +
  theme(legend.position =  c(0.18, 0.48), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7.5))

pdf(file = "OUT/LouvainCities_em2nw_ETMUN.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembers2 
dev.off()

## MAP community EM2NWF

# with fast greedy 
myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(geoClusterSF2F$PopAdmin11)
s[[1]]
bks <- c(s[[1]],s[[2]],s[[3]], 500000, 2000000, s[[6]] )
lbs <-  c("5 000", "50 000", "100 000", "500 000", "2M", "14 M" )
palcol <- brewer.pal(n = 7, name = "Paired")

citiesMembers2F <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = geoClusterSF2F ,
          mapping = aes(size = PopAdmin11,  colour = fgreedy), show.legend = "point", alpha = 0.7) +
  scale_color_manual(breaks = c(1:7), values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.5, 14))+
  annotate("text", label = "Source : ETMUN 2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = geoClusterSF2F %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Communautés (Fast Greedy)") +
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7.5))

pdf(file = "OUT/FGreedyCities_em2nwF_ETMUN.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembers2F 
dev.off()


####  RELATION BETWEEN COMMUNITY MEMBERSHIP AND POPULATION / ADMIN STATUS#####
#### Perform ANOVA (communities / PopAdmin11)

## EM2NW
D2 <- geoClusterSF2 %>% st_drop_geometry() 
anovaTab <- AnovaTab(df = D2 , varx = c("louvain"), 
                     vary = c("PopAdmin11")) 
resume <- ComputeRegression(D2, vardep = "PopAdmin11", varindep = "louvain")
resume <- resume$TABCOEF ## R2 = 0.07

# export DF
library(gridExtra)
pdf(file= "OUT/LouvainEM2NW_PopAdmin11_anovaTab.pdf", height = 7, width =6.5 )
grid.table(anovaTab)
dev.off()

pdf(file= "OUT/LouvainEM2NW_PopAdmin11_anovaR2.pdf", height = 7, width =6.5 )
grid.table(resume)
dev.off()

## EM2NWF

D2F <- geoClusterSF2F %>% st_drop_geometry()

anovaTab <- AnovaTab(df = D2F , varx = c("fgreedy"), 
                     vary = c("PopAdmin11")) 
resume <- ComputeRegression(D2F, vardep = "PopAdmin11", varindep = "fgreedy")
resume <- resume$TABCOEF ## R2 = 0.10

# export DF
pdf(file= "OUT/FgreedyEM2NWF_PopAdmin11_anovaTab.pdf", height = 7, width =6.5 )
grid.table(anovaTab)
dev.off()

pdf(file= "OUT/FgreedyEM2NWF_PopAdmin11_anovaR2.pdf", height = 7, width =6.5 )
grid.table(resume)
dev.off()
rm(resume, anovaTab, D2F,D2)

#### Chi2

## EM2NW (louvain / admin level)

tableCluster <- table(geoClusterSF2$louvain, geoClusterSF2$adminLevel)

testChi2 <- chisq.test(tableCluster)
testChi2$statistic
testChi2$expected
testChi2$residuals
resChi2<-as.data.frame(testChi2$residuals)
maxFreq <- max(resChi2$Freq)
minFreq <- min(resChi2$Freq)
ggplot(resChi2, aes(x = Var1, y = Var2, fill= Freq ))+ geom_tile() + scale_fill_gradient2(low = "#67a9cf", high = "#ef8a62", mid = "#f7f7f7", 
                                                                                          midpoint = 0, limit = c(minFreq,maxFreq)) + 
  labs(title = "Résidus relation cmmunauté-statut administratif (em2nw)", 
       x = "Communautés (Louvain)", 
       y = "Niveau administratif des localités (geonames)",
       fill = "Résidus de Pearson")

ggsave(filename = "OUT/LouvainEM2NW_Adminlevel_ResChi2_HeatMap.pdf", width = 11.7, height = 8.3, units = "in" )


library(vcd)

cramer <-assocstats(tableCluster)
cramer  # 0.26




## EM2NWF (Fgreedy / admin level)

tableCluster <- table(geoClusterSF2F$fgreedy, geoClusterSF2F$adminLevel)

testChi2 <- chisq.test(tableCluster)
testChi2$statistic
testChi2$expected
testChi2$residuals
resChi2<-as.data.frame(testChi2$residuals)
maxFreq <- max(resChi2$Freq)
minFreq <- min(resChi2$Freq)
ggplot(resChi2, aes(x = Var1, y = Var2, fill= Freq ))+ geom_tile() + scale_fill_gradient2(low = "#67a9cf", high = "#ef8a62", mid = "#f7f7f7", 
                                                                                          midpoint = 0, limit = c(minFreq,maxFreq)) + 
  labs(title = "Résidus relation communauté-statut administratif (em2nwF)", 
       x = "Communautés (Fast Greedy)", 
       y = "Niveau administratif des localités (geonames)",
       fill = "Résidus de Pearson")

ggsave(filename = "OUT/FgreedyEM2NWF_Adminlevel_ResChi2_HeatMap.pdf", width = 11.7, height = 8.3, units = "in" )


library(vcd)

cramer <-assocstats(tableCluster)
cramer  # 0.26

rm(cramer, testChi2, tableCluster)


#### MAP ASSOCIATION MEMBERSHIP #####


##FUNCTION

plotBipartOrga <- function(graph, NodesInfos, colName){
  require(dplyr)
  require(tidygraph)
  require(ggraph)
  require(RColorBrewer)
  plots <- list()
  
  ## Transfo in td graph 
  
  mypart <- enquo(colName)
  Tdg <- as_tbl_graph(graph) 
  Tdg <-  Tdg %>% activate(nodes) %>% mutate(!!mypart == as.character(!!mypart))
  
  namefilter <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(name) %>% deframe()
  
  #joint info nodes
  Tdg <- Tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")
  com <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(!!mypart) %>% deframe() %>% unique() 
  nColor <- length(com)
  colvect <- brewer.pal(n = nColor, name = "Paired")

  
  
  for(i in com ){
    tdtemp <- Tdg %>% activate(nodes) %>% filter(!!mypart == i)
    
    g <- ggraph(tdtemp, layout = "bipartite")+
      geom_edge_link(alpha = 0.1)+
      geom_node_point(aes(shape =type, size = 1-type), color = colvect[i],
                      show.legend = FALSE) + 
      geom_node_label(aes(label = Acro),repel = TRUE, size = 2.5)+
      labs(subtitle = paste("Communauté", i, sep = " "))+
      theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
    
    plots[[i]] <- g
    
  }
  
  return(plots) 
  
}

# grid plo

Gridplot <- function(myplots, n){
  require(cowplot)
  splitted_plots <- split(myplots, ceiling(seq_along(myplots)/n))
  
  lapply(splitted_plots, function(x) plot_grid(plotlist = x))
  
}


### LOUVAIN EM2NW

Assoem2nwG <- em2nw.clustering$igraph_object


plotEM2NW <- plotBipartOrga( Assoem2nwG, NodesInfos, colName = louvain)

G1 <- Gridplot(plotEM2NW, 4)
# select interesting  ones
G1select <- plotEM2NW[-c(1, 5,8)]
G1 <- Gridplot(G1select, 8)
G1
ggsave( filename = "OUT/EtmunAssoLouvain_em2nwf.pdf", width = 11.7, height = 8.3, units = "in" )
### FGREEDY EM2NWF
Assoem2nwFG <- em2nwF.clustering$igraph_object


plotEM2NWF <- plotBipartOrga( Assoem2nwFG, NodesInfos, colName = fgreedy)

G2 <- Gridplot(plotEM2NWF, 6)

G2[1]



ggsave( filename = "OUT/EtmunAssoFgreedy_em2nwf.pdf", width = 11.7, height = 8.3, units = "in" )

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


em2nwFDUAL <- DblDUALPROJECTION(em2nwF)

em2nwFDUAL$info_clustering

em2nwDUAL <- DblDUALPROJECTION(em2nw)

em2nwDUAL$info_clustering

em2nwFDUAL$igraph_objects

emnwDUAL <- DblDUALPROJECTION(emnw)

emnwFDUAL <- DblDUALPROJECTION(emnwF)

saveRDS(em2nwFDUAL, "DataProd/DualProj_em2nwF.rds")

saveRDS(em2nwDUAL, "DataProd/DualProj_em2nw.rds")
saveRDS(emnwFDUAL, "DataProd/DualProj_emnwF.rds")
saveRDS(emnwDUAL, "DataProd/DualProj_emnw.rds" )


