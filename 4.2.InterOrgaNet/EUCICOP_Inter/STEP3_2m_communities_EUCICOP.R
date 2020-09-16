############################### EUCICOP bipartite Networks STEP 3 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe EUCICOP (affiliation) biparti. 
#               Détection de communautés 
#
# 
############################################################################## PG juin 2020

### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/EUCICOP_Inter")

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
ParticipationEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/Participations_All_Eucicop_Europe.RDS")

## Information on project (for nodes)
ProjectEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/ProjectsEucicop_all_noduplicated.RDS")

## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")
DbCityInfos <- DbCity %>% st_drop_geometry()

# Infos on Nodes (STEP 2)

  NodesInfos <- readRDS("DataProd/NodesInfos_EUCICOP.RDS")

## shapes

## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)

# Matrices from STEP1  

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/Eucicop_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/Eucicop_em2nw.rds")

# First component of emnw
emnw_1stComp <- readRDS("DataProd/Eucicop_emnw_1stcomp.rds")


# First component of emnw
em2nw_1stComp <- readRDS("DataProd/Eucicop_em2nw_1stcomp.rds")


### ==== IGRAPH TRANSFORMATION ====



emnw.g <- graph.incidence(emnw_1stComp , directed = F)

em2nw.g <-  graph.incidence(em2nw_1stComp, directed = F)





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

# saveRDS(Memberem2nw, "DataProd/Membership_em2nw_EUCICOP.rds")
# saveRDS(Memberemnw, "DataProd/Membership_emnw_EUCICOP.rds")


### ==== PAIWISE COMPARISON OF CLUSTERING ====
### Dataframe Classification
dfclustering <- em2nw.clustering$memberships
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


clusterG <- em2nw.clustering$igraph_object

#just create  a variable name specific to project such as Acro in ETMUN
NodesInfos <- NodesInfos %>% mutate(Acro = ifelse(str_detect(name, "P"), label, NA))

# add geometry 

NodesInfos <- NodesInfos %>% left_join(select(DbCity, name = geonameId, geometry))
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
# ggraph(tdClusterG, layout = "fr")+
#   geom_edge_link(alpha = 0.1)+
#   geom_node_point(aes(color = as.factor(V(tdClusterG)$fgreedy), shape = V(tdClusterG)$type, size = 1-V(tdClusterG)$type),
#                   show.legend = FALSE) + 
#   geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
#   facet_nodes(~V(tdClusterG)$fgreedy)



#### 2. PLot cities community with geographical coordinates

## Create an SF object
geoCluster <- dfclustering %>% left_join(NodesInfos)

geoCluster <- geoCluster %>% filter(type == TRUE) %>% select(-Acro, -Date)


geoClusterSF <- st_as_sf(geoCluster)
class(geoClusterSF)
geoClusterSF <- geoClusterSF %>% mutate_at(.vars=c(3:6), ~as.character(.))


### ==== CHOOSE AND MAP 1 PARTITION  ====
# Transfo population
# filter pop
geoClusterSF <- geoClusterSF %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 5000)

length(unique(geoClusterSF$louvain))

TopLouvainNCities <- geoClusterSF %>% group_by(louvain)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopLouvainNCities$N)

TopLouvainNCities <- TopLouvainNCities %>% filter(N>5) %>% select(louvain) %>% deframe()
## MAP community EM2NW
geoLouvainTop <- geoClusterSF %>% filter(louvain %in% TopLouvainNCities)

n<- length(unique(geoLouvainTop$louvain))
# Louvain with facet 

library(randomcoloR)

palette <- distinctColorPalette(n)

names(palette) <- sort(as.numeric(unique(geoLouvainTop$louvain)))

palette
palette[16] <- "black"
Nlouvain <- geoLouvainTop %>% group_by(louvain) %>% summarise(Ncities = n()) %>% st_drop_geometry()

ggplot()+
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data =geoLouvainTop, aes(color = louvain), alpha= 0.8, show.legend = FALSE)+
  scale_color_manual(values =palette)+
  labs(y="", x ="", 
       caption = "Note : localités de plus de 5000 habitants seulement.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  # geom_sf_text(aes(label = label, check_overlap = TRUE), size = 2)+ 
  facet_wrap(~ louvain)+
  geom_text( data    = Nlouvain,
    mapping = aes(x = Inf, y = Inf, label = paste("N = ", Ncities, sep = "")),
    hjust   = 1.2,
    vjust   = 2, size = 3)

ggsave( filename = "OUT/Louvain_EUCICOPCities_FacetMap_em2nw.pdf", width = 11.7, height = 8.3, units = "in" )

# Louvain on synthesis map
geoClusterSF2 <- geoClusterSF %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 20000)

TopLouvainNCities2 <- geoClusterSF2 %>% group_by(louvain)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopLouvainNCities2$N)

TopLouvainNCities2 <- TopLouvainNCities2 %>% filter(N>10) %>% select(louvain) %>% deframe()

palette2 <- palette[names(palette) %in% as.numeric(TopLouvainNCities2)]
geoLouvainTop2 <- geoClusterSF2 %>% filter(louvain %in% TopLouvainNCities2)

myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(geoLouvainTop2$PopAdmin11)
s[[1]]
bks <- c(s[[1]],s[[4]], 2000000, s[[6]] )
lbs <-  c("20 000",  "150 000",  "2M", "14 M" )
palcol <- palette2

citiesMembers2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "white", color = "#bfbfbf", size = 0.5) +
  geom_sf(data = geoLouvainTop2 ,
          mapping = aes(size = PopAdmin11,  colour = louvain), show.legend = "point", alpha = 0.9) +
  scale_color_manual( values = palcol)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.2, 12))+
  annotate("text", label = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = geoLouvainTop2 %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
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

pdf(file = "OUT/LouvainCities_em2nw_EUCICOP_20kHab.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembers2 
dev.off()

### ====EXPLORATION OF PROJECTS MEMBERSHIPS ====


ProjectMembership <- dfclustering  %>% filter(type == FALSE)

ProjectMembership <- ProjectMembership %>% left_join(select(ProjectEUCICOP, Programme, Project.Name, Period, Thematic1,Thematic2,Thematic3,name = ID_PROJECT))


ProgrammeLouvain <- ProjectMembership %>% group_by(Programme,louvain)%>%
  summarise(N = n()) %>% 
  group_by(louvain) %>% 
  mutate(ppct = prop.table(N)*100, Ntot = sum(N))

Top5LouvainProject <-  ProgrammeLouvain %>% group_by(louvain)%>% top_n(5, N)

  
NProjectLouvain <- ProgrammeLouvain %>% select(louvain, Ntot) %>% distinct()


library(ggforce)

#Manual selection of Communitites 

ProjectToPlot <- c(17,1,11,3)
palette3 <- palette[names(palette) %in% ProjectToPlot]
Top5LouvainProgPlot <- Top5LouvainProject %>% filter(louvain %in% ProjectToPlot)
tempN <- NProjectLouvain %>% filter(louvain %in% ProjectToPlot)

library(questionr)
#irec(Top5LouvainProgPlot)

## Recodage de Top5LouvainProgPlot$Programme en Top5LouvainProgPlot$Programme_rec
Top5LouvainProgPlot$Programme_rec <- Top5LouvainProgPlot$Programme
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Baltic Sea Region"] <- "2000 - 2006 Baltic Sea Region"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Cadses "] <- "2000 - 2006 Cadses"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 EUREGIO - Euregio Rhein-Waal and Euregio Rhein-Mass-Nord (NL-DE)"] <- "2000 - 2006 EUREGIO -\n Rhein-Waal and\nRhein-Mass-Nord\n(NL-DE)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Finland - Estonia (FI-EE)"] <- "2000 - 2006 Finland -\nEstonia (FI-EE)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Islands (IT-FR)"] <- "2000 - 2006 Islands\n(IT-FR)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Slovenia - Hungary - Croatia (SI-HU-HR)"] <- "2000 - 2006 Slovenia -\nHungary - Croatia (SI-HU-HR)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2000 - 2006 Western Mediterranean "] <- "2000 - 2006 Western\nMediterranean"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 2 SEAS (FR-UK-BE-NL)"] <- "2007 - 2013 2 SEAS\n(FR-UK-BE-NL)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Central Baltic (FI-SE-EE-LA)"] <- "2007 - 2013 Central Baltic\n(FI-SE-EE-LA)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Central Europe"] <- "2007-2013 Central Europe"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 France (Channel) - England (FR-UK)"] <- "2007 - 2013 France (Channel) -\nEngland (FR-UK)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 France - Wallonia - Flanders (BE-FR)"] <- "2007 - 2013 France - Wallonia -\nFlanders (BE-FR)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Hungary - Croatia (HU-HR)"] <- "2007 - 2013 Hungary -\nCroatia (HU-HR)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Italy - France Maritime (IT-FR)"] <- "2007 - 2013 Italy -\nFrance Maritime (IT-FR)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Latvia-Lithuania (LV-LT)"] <- "2007 - 2013 Latvia-Lithuania (LV-LT)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Mediterranean Sea Basin ENPI CBC"] <- "2007 - 2013 Mediterranean\nSea Basin\nENPI CBC"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 North West Europe"] <- "2007 - 2013\nNorth West Europe"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 Programme MED"] <- "2007 - 2013 Programme MED"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 South Baltic (PL-SE-DK-LT-DE)"] <- "2007 - 2013 South Baltic\n(PL-SE-DK-LT-DE)"
Top5LouvainProgPlot$Programme_rec[Top5LouvainProgPlot$Programme == "2007 - 2013 South East Europe"] <- "2007 - 2013 South East Europe"

ggplot(Top5LouvainProgPlot) + 
  geom_bar(aes(x= reorder(Programme_rec, desc(ppct)), y = ppct, fill = as.character(louvain)), stat = "identity", show.legend = FALSE) + 
  facet_wrap( ~ louvain, scales ="free") + 
  scale_fill_manual( values = palette3)+
  geom_label(aes(label = N, x= Programme_rec, y = ppct), position = position_stack(0.5), color = "black")+
  geom_text( data    = tempN,
             mapping = aes(x = Inf, y = Inf, label = paste("NProj = ", Ntot, sep = "")),
             hjust   = 1.2,
             vjust   = 1.2, size = 2.5)+
   coord_flip() + 
  labs(y = "Pct de projets", x = "Programme de coopération", 
                      caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")
  


ggsave( filename = "OUT/Louvain_EUCICOPProject_Programme_em2nw.pdf", width = 8.3, height = 5.8, units = "in" )

### exploration of main programme for 4 communities with strong  regional patterns

ProjectToPlot <- c(8,5,28,2)
palette3 <- palette[names(palette) %in% ProjectToPlot]
Top5LouvainProgPlot <- Top5LouvainProject %>% filter(louvain %in% ProjectToPlot)
tempN <- NProjectLouvain %>% filter(louvain %in% ProjectToPlot)

# irec(Top5LouvainProgPlot)

## Recodage de Top5LouvainProgPlot$Programme en Top5LouvainProgPlot$Programme_rec
Top5LouvainProgPlot$Programme_rec <- fct_recode(Top5LouvainProgPlot$Programme,
               "2000 - 2006 France - Spain" = "2000 - 2006 France - Spain (FR-ES)",
               "2000 - 2006 Ireland - Wales" = "2000 - 2006 Ireland - Wales (IE-UK)",
               "2000 - 2006 Kvarken -\nMittskandia (FI-SE-NO)" = "2000 - 2006 Kvarken - Mittskandia (FI-SE-NO)",
               "2000 - 2006 Poland -\nUkraine - Belarus" = "2000 - 2006 Poland - Ukraine - Belarus (PL-UA-BY)",
               "2007 - 2013\nBlack Sea Basin\nENPI CBC" = "2007 - 2013 Black Sea Basin ENPI CBC",
               "2007 - 2013\nBulgaria -FYROM\n IPA CBC" = "2007 - 2013 Bulgaria - Former Yugoslav Republic of Macedonia IPA CBC (BG-FYROM)",
               "2007 - 2013 Greece - Bulgaria" = "2007 - 2013 Greece - Bulgaria (EL-BG)",
               "2007 - 2013 Hungary - Romania" = "2007 - 2013 Hungary - Romania (HU-RO)",
               "2007 - 2013 HU-SK-RO-UA\nENPI CBC" = "2007 - 2013 Hungary-Slovakia-Romania-Ukraine ENPI CBC",
               "2007 - 2013 Karelia\nENPI CBC" = "2007 - 2013 Karelia ENPI CBC",
               "2007 - 2013 PO-BY-UA\nENPI CBC" = "2007 - 2013 Poland-Belarus-Ukraine ENPI CBC",
               "2007 - 2013\nPoland - Slovak Republic" = "2007 - 2013 Poland - Slovak Republic (PL-SK)",
               "2007 - 2013 RO-UA-MD\nENPI CBC" = "2007 - 2013 Romania-Ukraine-Moldova ENPI CBC",
               "2007 - 2013\nSouth West Europe" = "2007 - 2013 South West Europe",
               "2007 - 2013 Spain -\nFrance - Andorra" = "2007 - 2013 Spain - France - Andorra (ES-FR-AD)",
               "2007 - 2013 Spain - Portugal" = "2007 - 2013 Spain - Portugal (ES-PT)",
               "2014 - 2020 Interreg\nIPA CBC\nBulgaria - FYROM" = "2014 - 2020 Interreg IPA CBC Bulgaria - Former Yugoslav Republic of Macedonia",
               "2014 - 2020 INTERREG VB\nNorthern Periphery\nand Arctic" = "2014 - 2020 INTERREG VB Northern Periphery and Arctic")
Top5LouvainProgPlot$Programme_rec <- as.character(Top5LouvainProgPlot$Programme_rec)
## Recodage de Top5LouvainProgPlot$Programme en Top5LouvainProgPlot$Programme_rec

Top5LouvainProgPlot$Programme_rec <- as.character(Top5LouvainProgPlot$Programme_rec)
ggplot(Top5LouvainProgPlot) + 
  geom_bar(aes(x= reorder(Programme_rec, desc(ppct)), y = ppct, fill = as.character(louvain)), stat = "identity", show.legend = FALSE) + 
  facet_wrap( ~ louvain, scales ="free") + 
  scale_fill_manual( values = palette3)+
  geom_label(aes(label = N, x= Programme_rec, y = ppct), position = position_stack(0.5), color = "black")+
  geom_text( data    = tempN,
             mapping = aes(x = Inf, y = Inf, label = paste("NProj = ", Ntot, sep = "")),
             hjust   = 1.05,
             vjust   = 1.2, size = 2.5)+
  coord_flip() + 
  labs(y = "Pct de projets", x = "Programme de coopération", 
       caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

ggsave( filename = "OUT/LouvainRegional_EUCICOPProject_Programme_em2nw.pdf", width = 8.3, height = 5.8, units = "in" )



### exploration of main thematic object for the latter 4 communities with strong  regional patterns

# create a thematic variable 
ProjectMembership <- ProjectMembership %>% mutate(FirstTheme = ifelse(is.na(Thematic1), Thematic2, Thematic1)) %>%
                      mutate(FirstTheme = ifelse(is.na(FirstTheme), Thematic3, FirstTheme))


library(skimr)
skim(ProjectMembership)
ThematicLouvain <- ProjectMembership %>% group_by(FirstTheme,louvain)%>%
  summarise(N = n()) %>% 
  group_by(louvain) %>% 
  mutate(ppct = prop.table(N)*100, Ntot = sum(N))

Top5LouvainTheme <-  ThematicLouvain %>% group_by(louvain)%>% top_n(5, N)


NThemeLouvain <-ThematicLouvain %>% select(louvain, Ntot) %>% distinct()


# select the 4 communities 

Top5LouvainThemPlot <- Top5LouvainTheme %>% filter(louvain %in% ProjectToPlot)
tempN <- NThemeLouvain %>% filter(louvain %in% ProjectToPlot)

# recode theme 

#irec(Top5LouvainThemPlot)
## Recodage de Top5LouvainThemPlot$FirstTheme en Top5LouvainThemPlot$FirstTheme_rec
Top5LouvainThemPlot$FirstTheme_rec <- fct_recode(Top5LouvainThemPlot$FirstTheme,
               "Agriculture and fisheries\nand forestry" = "Agriculture and fisheries and forestry",
               "Institutional cooperation and\ncooperation networks" = "Institutional cooperation and cooperation networks",
               "Sustainable management\nof natural resources" = "Sustainable management of natural resources")
Top5LouvainThemPlot$FirstTheme_rec <- as.character(Top5LouvainThemPlot$FirstTheme_rec)
# Plot

ggplot(Top5LouvainThemPlot) + 
  geom_bar(aes(x= reorder(FirstTheme_rec, desc(ppct)), y = ppct, fill = as.character(louvain)), stat = "identity", show.legend = FALSE) + 
  facet_wrap( ~ louvain, scales ="free") + 
  scale_fill_manual( values = palette3)+
  geom_label(aes(label = N, x= FirstTheme_rec, y = ppct), position = position_stack(0.5), color = "black")+
  geom_text( data    = tempN,
             mapping = aes(x = Inf, y = Inf, label = paste("NProj = ", Ntot, sep = "")),
             hjust   = 1.05,
             vjust   = 1.2, size = 2.5)+
  coord_flip() + 
  labs(y = "Pct de projets", x = "Principaux objectifs thématiques", 
       caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

ggsave( filename = "OUT/LouvainRegional_EUCICOPProject_Theme_em2nw.pdf", width = 8.3, height = 5.8, units = "in" )

### ==== EXPLORE EDGES INTER-COMMUNITIES ====

## get communities 
Memberem2nw <- readRDS("DataProd/Membership_em2nw_EUCICOP.rds") 

## set communities in the original graph

em2nw.g  <-  set_vertex_attr(em2nw.g  ,"louvain", index = Memberem2nw$name, value = Memberem2nw$louvain)

# Compute crossing
E(em2nw.g)$crossing <- igraph::crossing(communities = V(em2nw.g)$louvain, graph =  em2nw.g)


XCE = which(V(em2nw.g)[ends(em2nw.g, E(em2nw.g))[,1]]$louvain !=
              V(em2nw.g)[ends(em2nw.g, E(em2nw.g))[,2]]$louvain)


Sub_em2nw_intercom <- subgraph.edges(em2nw.g, E(em2nw.g)[XCE], delete.vertices = TRUE)
# compare sub graph

ecount(Sub_em2nw_intercom)/ecount(em2nw.g) # 17% of edges

vcount(Sub_em2nw_intercom)/vcount(em2nw.g) # 31% of edges

graph.density(Sub_em2nw_intercom) # 0.00058 



# Add degree and get nodes df

V(Sub_em2nw_intercom)$degree <- degree(Sub_em2nw_intercom, loops = F)

Inter_Louvain_em2nwNodes <- as_data_frame(Sub_em2nw_intercom, what = "vertices")

## Get nodes info
Inter_Louvain_em2nwNodes <- Inter_Louvain_em2nwNodes %>% left_join(NodesInfos)

CitiesInterLouvain <- Inter_Louvain_em2nwNodes %>% filter(type == TRUE)
ProjectInterLouvain <- Inter_Louvain_em2nwNodes %>% filter(type == FALSE)

ProjectInterLouvain <- ProjectInterLouvain %>% left_join(select(ProjectEUCICOP, Programme, name = ID_PROJECT))

SumProjectInterLouvain <- ProjectInterLouvain %>% group_by(Programme) %>% 
  summarise(DegreeSum = sum(degree))

SumProjectInterLouvain <- SumProjectInterLouvain %>% mutate( PctEdges = (DegreeSum/sum(DegreeSum))*100)

TopProgrammeInterLouvain <- SumProjectInterLouvain %>%  top_n(10, DegreeSum)
sum(TopProgrammeInterLouvain$PctEdges) # 42% of edges
#Plot top programme by degree of inter communities subgraph

ggplot(TopProgrammeInterLouvain ) + 
  geom_bar(aes(x= reorder(Programme, desc(PctEdges)), y = PctEdges), stat = "identity", show.legend = FALSE) + 
  geom_label(aes(label =DegreeSum, x= Programme, y = PctEdges), position = position_stack(0.5), color = "black")+
  geom_text( mapping = aes(x = Inf, y = Inf, label = paste("Nb Total de liens inter-communautés = ", sum(SumProjectInterLouvain$DegreeSum), sep = "")),
             hjust   = 1.05,
             vjust   = 1.2, size = 3)+
  coord_flip() + 
  labs(y = "% de liens inter-communautaires", x = "Programme de coopération", 
       caption = "Note : à partir du sous-graphe comprenant uniquement les liens entre communautés (louvain sur em2nw_1stcomp).\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

ggsave( filename = "OUT/InterLouvainEdges_EUCICOPProject_TopProgramme_em2nw.pdf", width = 8.3, height = 5.8, units = "in" )

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




em2nwDUAL <- DblDUALPROJECTION(em2nw_1stComp)

em2nwDUAL$info_clustering

table(em2nwDUAL$membershipCities$louvainNW)
table(em2nwDUAL$membership_Orga$louvainNW)

saveRDS(em2nwDUAL, "DataProd/DualProj_em2nw_eucicop.rds")

emnwDUAL <- DblDUALPROJECTION(emnw_1stComp)

emnwDUAL$info_clustering

saveRDS(emnwDUAL, "DataProd/DualProj_emnw_eucicop.rds")


##EXPLORE COMMUNITIES OF EM2NW DUAL PROJECTION

# Comparing classification for cities

Memberem2nw <- readRDS("DataProd/Membership_em2nw_EUCICOP.rds") 
Memberem2nwCities <- Memberem2nw %>% filter(type ==TRUE)
Memberem2nwCities <- Memberem2nwCities %>% left_join(NodesInfos)

Memberem2nwCitiesDUAL <- em2nwDUAL$membershipCities 

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



Memberem2nwCitiesDUAL <- Memberem2nwCitiesDUAL %>% left_join(NodesInfos)

Memberem2nwCitiesDUAL_sf <- st_as_sf(Memberem2nwCitiesDUAL)
class(Memberem2nwCitiesDUAL_sf)
Memberem2nwCitiesDUAL_sf  <- Memberem2nwCitiesDUAL_sf  %>% mutate_at(.vars=c(2:5), ~as.character(.))

Memberem2nwCitiesDUAL_sf  <- Memberem2nwCitiesDUAL_sf  %>%
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

library(randomcoloR)

paletteDUAL <- distinctColorPalette(length(unique(Memberem2nwCitiesDUAL_sf$louvainNW)))

names(paletteDUAL) <- sort(as.numeric(unique(Memberem2nwCitiesDUAL_sf$louvainNW)))

## Map 



myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(Memberem2nwCitiesDUAL_sf$PopAdmin11)
s[[1]]
bks <- c(5000,20000, 50000, 150000, 1000000, 14000000)
lbs <-  c("5 000" ,"20 000","50 000",  "150 000",  "1M", "14 M" )
palcol <- paletteDUAL

citiesMembersDUALprojNW <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = Memberem2nwCitiesDUAL_sf ,
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
  geom_sf_text(data = Memberem2nwCitiesDUAL_sf %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
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



pdf(file = "OUT/LouvainCitiesDUALnw_em2nw_EUCICOP.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembersDUALprojNW
dev.off()

## Community with kernel density
Memberem2nwCitiesDUAL_sf <- Memberem2nwCitiesDUAL_sf %>% mutate(X = st_coordinates(.)[,1], Y = st_coordinates(Memberem2nwCitiesDUAL_sf)[,2])
TopLouvainNCities <- Memberem2nwCitiesDUAL_sf %>% group_by(louvainNW)%>% summarise(N= n()) %>%st_drop_geometry() 
TopLouvainNCities2 <- TopLouvainNCities %>% filter(N>3) %>% select(louvainNW) %>% deframe()

Memberem2nwCitiesDUAL_sf_filtered <- Memberem2nwCitiesDUAL_sf %>% filter(louvainNW %in% TopLouvainNCities2)


palcol <- paletteDUAL[names(paletteDUAL) %in% as.numeric(TopLouvainNCities2)]

citiesMembersDUALprojNW_kernel <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = Memberem2nwCitiesDUAL_sf_filtered , show.legend = "point", alpha = 0) +
  geom_density_2d(data = Memberem2nwCitiesDUAL_sf_filtered , aes(x = X, y = Y, color = louvainNW),size = 0.1) +
  scale_color_manual( values = palcol)+
  annotate("text", label = "Note : Louvain sur le graphe ville-ville, valué et normalisé. Les 4 communautés de moins de 3 villes ont été retirées.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = Memberem2nwCitiesDUAL_sf_filtered %>% filter(PopAdmin11 > 500000), aes(label = label), size = 2.2, color = "#4d4d4d",
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



pdf(file = "OUT/LouvainCitiesDUALnw_em2nw_EUCICOP_Kernel.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesMembersDUALprojNW_kernel
dev.off()

## Facet representation

TopLouvainNCities <- Memberem2nwCitiesDUAL_sf %>% group_by(louvainNW)%>% summarise(N= n()) %>%st_drop_geometry() 
summary(TopLouvainNCities$N)


library(ggalt)
library(spats)
test <- ggplot()+
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data =Memberem2nwCitiesDUAL_sf , aes(color = louvainNW), alpha= 0, show.legend = FALSE, size =0.1)+
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
  geom_density_2d(data = Memberem2nwCitiesDUAL_sf, aes(x = X, y = Y, color = louvainNW)) 



