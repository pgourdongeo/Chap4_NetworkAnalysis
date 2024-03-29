################# ETMUN 1 modes projections Networks STEP 4 ###################
#                               
#                          
# DESCRIPTION : Travail sur les projections 1-mode du graphe biparti ETMUN. 
#               Complément des graphes et assortativité 
#
# 
############################################################################## PG juin 2020

### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_inter/")

### Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(tidygraph)
library(ggraph)
library(sf)
library(patchwork)
library(RColorBrewer)
options(scipen = 999)


### DATA
## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_europe.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", 
                       stringsAsFactors = F)


# 1 mod projection for em2nw and em2nwf

em2nwFDUAL <- readRDS("DataProd/DualProj_em2nwF.rds")

em2nwDUAL <- readRDS("DataProd/DualProj_em2nw.rds")

emnwDual <- readRDS("DataProd/DualProj_emnw.rds")

# Centrality measures (STEP 2)

Centralities <- readRDS("DataProd/Centralities_ETMUN_enww_em2nw.rds")

CentralitiesFiltered <- readRDS("DataProd/Centralities_ETMUN_enwwF_em2nwF.rds")

## Community detection on bipartite (STEP 3)

Memberem2nw <- readRDS("DataProd/Membership_em2nw_ETMUN.rds")

Memberem2nwF <- readRDS("DataProd/Membership_em2nwF_ETMUN.rds")

# Infos on Nodes (STEP 2)

NodesInfos <- readRDS("DataProd/NodesInfos_ETMUN.RDS")


## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")
DbCityInfos <- DbCity %>% st_drop_geometry() %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

## get GN info for Nodes Info
NodesInfos <- NodesInfos %>% left_join(select(MembershipEtmun, name = geonameId, lng_GN,lat_GN))%>% distinct()%>%
  left_join(select(DbCityInfos,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11), by = c("name" = "geonameId"))



### ==== 1-mode Cities  ====

# basic weight
em2nwCities.g <- em2nwDUAL$igraph_objects[[4]]

em2nwFCities.g <- em2nwFDUAL$igraph_objects[[4]]


## Degree and strength

### compute degrees on 1-mode network

V(em2nwCities.g)$degree <- degree(em2nwCities.g)
V(em2nwCities.g)$strength <- strength(em2nwCities.g, vids = V(em2nwCities.g), loops = F )
graph.density(em2nwCities.g) # 0.70
em2nwDUAL$info_clustering

V(em2nwFCities.g)$degree <- degree(em2nwFCities.g)
V(em2nwFCities.g)$strength <- strength(em2nwFCities.g, vids = V(em2nwFCities.g), loops = F )
graph.density(em2nwFCities.g) # 0.37
em2nwFDUAL$info_clustering


## Complement


em2nw_compl <- complementer(em2nwCities.g)

em2nwF_compl <- complementer(em2nwFCities.g)


### ==== Plotting Network  ====

## EM2NW
Tdg <- as_tbl_graph(em2nw_compl)

# Add complete node info
namefilter <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(name) %>% deframe()

Tdg <- Tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")


# Add Centralities of bipartite graph

Central2 <- Centralities[[1]] %>% select(geonameId | ends_with("2")) %>% rename(name = geonameId) 

Tdg <- Tdg %>% activate(nodes) %>% left_join(Central2  %>% filter(name %in% namefilter), by = "name")

summary(V(Tdg)$D2)

# Add community membership on bipartite 

Member2 <- Memberem2nw %>% filter(type == TRUE)

Tdg <- Tdg %>% activate(nodes) %>% left_join(Member2 , by = "name")

# Filter graph to get top degree (in bipartite)

TdgF <- Tdg %>% activate(nodes) %>% filter(D2 > 6) %>% filter(!node_is_isolated())

Nk <- Tdg %>% activate(nodes) %>% filter(D2 > 6) %>% fortify.tbl_graph() %>% nrow()

Nk/vcount(Tdg)

vcount(TdgF)/Nk
## Deal with color
ComVec <- TdgF  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(louvain) %>% deframe() %>% unique() %>%sort()
pal <-  brewer.pal(n = 11, name = "Paired")
pal <- pal[ComVec]
## Create a layout with geo coordinates

coords <-  TdgF %>% activate(nodes) %>% select(lng_GN, lat_GN) %>% fortify.tbl_graph() %>% rename(x = lng_GN, y = lat_GN)


manual_layout <- create_layout(graph = TdgF,
                               layout = "manual", node.positions = coords)


##Plot the complement 
g1 <- ggraph(manual_layout) + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point( aes(color = as.character(louvain), size = PopAdmin11) )+ 
  geom_node_text(aes(label = label),repel = TRUE, size = 2.5)+
  scale_color_manual(values = pal, breaks = ComVec)+
  scale_size(range = c(1,8))+
  labs(x = "long", y = "lat",
       color = "Communautés dans le biparti\n(louvain sur em2nw)", size = "Population Administrative 2011",
       caption = "Complément du graphe em2nw (villes-villes) :\nchaque lien représente l'absence de lien direct.\nSpatialisation : coordonnées géographiques\nSource : ETMUN 2019/ PG 2020")
 
g1bis <- ggraph(manual_layout) + 
  geom_edge_link(alpha = 0.2, start_cap = circle(4, 'mm'),
                 end_cap = circle(4, 'mm'))+
   geom_node_point( aes(color = adminLevel, size = PopAdmin11), alpha= 0.7 )+ 
  geom_node_text(aes(label = label),repel = TRUE, size =3)+
  scale_color_brewer(palette = "Set1")+
  scale_size(range = c(2,15))+
  labs(x = "long", y = "lat",
       color = "Niveau Administratif (Geonames)", size = "Population Administrative 2011",
       caption = "Complément du graphe em2nw (villes-villes) :\nchaque lien représente l'absence de lien direct.\nSpatialisation : coordonnées géographiques\nSource : ETMUN 2019/ PG 2020")

## EM2NW ###
Tdg <- as_tbl_graph(em2nwF_compl)

# Add complete node info
namefilter <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(name) %>% deframe()

Tdg <- Tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")


# Add Centralities of bipartite graph

Central2 <- CentralitiesFiltered[[1]] %>% select(geonameId | ends_with("2")) %>% rename(name = geonameId) 

Tdg <- Tdg %>% activate(nodes) %>% left_join(Central2  %>% filter(name %in% namefilter), by = "name")

summary(V(Tdg)$D2)

# Add community membership on bipartite 

Member2 <- Memberem2nwF %>% filter(type == TRUE)

Tdg <- Tdg %>% activate(nodes) %>% left_join(Member2 , by = "name")


# Filter graph to get  degree

TdgF <- Tdg %>% activate(nodes) %>% filter(D2 > 6) %>% filter(!node_is_isolated()) %>% mutate(degreeCompl  = centrality_degree())

Nk <- Tdg %>% activate(nodes) %>% filter(D2 > 6) %>% fortify.tbl_graph() %>% nrow()

Nk/vcount(Tdg)

vcount(TdgF)/Nk

## Deal with color
ComVec <- TdgF  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(fgreedy) %>% deframe() %>% unique() %>%sort()
pal <-  brewer.pal(n = 7, name = "Paired")
pal <- pal[ComVec]

## Create a layout with geo coordinates

coords <-  TdgF %>% activate(nodes) %>% select(lng_GN, lat_GN) %>% fortify.tbl_graph() %>% rename(x = lng_GN, y = lat_GN)


manual_layout <- create_layout(graph = TdgF,
                               layout = "manual", node.positions = coords)
##Plot the completment 
g2<-ggraph(manual_layout) + 
  geom_edge_link(alpha = 0.2, start_cap = circle(3, 'mm'),
                  end_cap = circle(3, 'mm'))+
  geom_node_point( aes(color = as.character(fgreedy), size = PopAdmin11), alpha= 0.7  )+ 
  geom_node_text(aes(label = label),repel = TRUE, size = 3)+
  scale_color_manual(values = pal, breaks = ComVec)+
  scale_size(range = c(2,12))+
  labs(x = "long", y = "lat",
       color = "Communautés dans le biparti\n(Fast Greedy sur em2nwF)", size = "Population Administrative 2011",
       caption = "Complément du graphe  em2nwF (villes-villes):\nchaque lien représente l'absence de lien direct.\nSpatialisation : coordonnées géographiques\nSources : ETMUN 2019/ PG 2020")


g2bis<-ggraph(TdgF, layout = "kk") + 
  geom_edge_hive(alpha = 0.2, start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'))+
  geom_node_point( aes(color = as.character(fgreedy), size = degreeCompl), alpha= 0.7  )+ 
  geom_node_text(aes(label = label),repel = TRUE, size = 3)+
  scale_color_manual(values = pal, breaks = ComVec)+
  scale_size(range = c(2,12))+
  labs(x = "long", y = "lat",
       color = "Communautés dans le biparti\n(Fast Greedy sur em2nwF)", size = "Degré \n(dans le complément du graphe)",
       caption = "Complément du graphe  em2nwF (villes-villes):\nchaque lien représente l'absence de lien direct.\nSpatialisation : kk\nSources : ETMUN 2019/ PG 2020")


ggsave(g1bis, filename = "OUT/Compl_Cities_em2nw_ETMUN_deg6.pdf",width = 8.3, height = 5.8)

ggsave(g2bis, filename = "OUT/Compl_Cities_em2nwF_ETMUN_deg6_KK.pdf",width = 8.3, height = 8.3)


### ==== 1-mode Associations ====

# basic weight
em2nwAsso.g <- em2nwDUAL$igraph_objects[[1]]

em2nwFAsso.g <- em2nwFDUAL$igraph_objects[[1]]


## Degree and strength

### compute degrees on 1-mode network

V(em2nwAsso.g)$degree <- degree(em2nwAsso.g)
V(em2nwAsso.g)$strength <- strength(em2nwAsso.g, vids = V(em2nwAsso.g), loops = F )
graph.density(em2nwAsso.g) # 0.74
em2nwDUAL$info_clustering

V(em2nwFAsso.g)$degree <- degree(em2nwFAsso.g)
V(em2nwFAsso.g)$strength <- strength(em2nwFAsso.g, vids = V(em2nwFAsso.g), loops = F )
graph.density(em2nwFAsso.g) # 0.72
em2nwFDUAL$info_clustering


## Complement


em2nw_compl <- complementer(em2nwAsso.g)

em2nwF_compl <- complementer(em2nwFAsso.g)


### ==== Plotting Network ASSO  ====

## EM2NW
Tdg <- as_tbl_graph(em2nw_compl)

# Add complete node info
namefilter <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(name) %>% deframe()

Tdg <- Tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")


# Add Centralities of bipartite graph

Central2 <- Centralities[[2]] %>% select(Code | ends_with("2")) %>% rename(name = Code) 

Tdg <- Tdg %>% activate(nodes) %>% left_join(Central2  %>% filter(name %in% namefilter), by = "name")

summary(V(Tdg)$D2)

# Add community membership on bipartite 

Member2 <- Memberem2nw %>% filter(type == FALSE)

Tdg <- Tdg %>% activate(nodes) %>% left_join(Member2 , by = "name")

# Filter graph to get top degree (in bipartite)

TdgF <- Tdg %>% activate(nodes) %>% filter(D2 > 25) %>% filter(!node_is_isolated())

## Deal with color
ComVec <- TdgF  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(louvain) %>% deframe() %>% unique() %>%sort()
pal <-  brewer.pal(n = 11, name = "Paired")
pal <- pal[ComVec]


##Plot the complement 
g3 <- ggraph(TdgF, layout = "lgl") + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point( aes(color = as.character(louvain), size = D2) )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  scale_color_manual(values = pal, breaks = ComVec)+
  scale_size(range = c(1,12))+
  labs(color = "Communautés dans le biparti\n(louvain sur em2nw)", size = "Degrés biparti (nb villes membres)",
       caption = "Complément du graphe em2nw (villes-villes) :\nchaque lien représente l'absence de lien direct.\nSpatialisation : coordonnées géographiques\nSource : ETMUN 2019/ PG 2020")+
  theme_graph()

g3bis <- ggraph(TdgF, layout = "lgl") + 
  geom_edge_link(alpha = 0.2, start_cap = circle(4, 'mm'),
                 end_cap = circle(4, 'mm'))+
  geom_node_point( aes(color = Country, size = D2), alpha= 0.7 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size =3)+
  scale_size(range = c(2,15))+
  labs(color = "Niveau Administratif (Geonames)", size =  "Degrés biparti (nb villes membres)",
       caption = "Complément du graphe em2nw (villes-villes) :\nchaque lien représente l'absence de lien direct.\nSpatialisation : coordonnées géographiques\nSource : ETMUN 2019/ PG 2020")+ 
  theme_graph()

## EM2NWF ###
Tdg <- as_tbl_graph(em2nwF_compl)

# Add complete node info
namefilter <- Tdg  %>% activate(nodes) %>% fortify.tbl_graph()%>% dplyr::select(name) %>% deframe()

Tdg <- Tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name")


# Add Centralities of bipartite graph

Central2 <- CentralitiesFiltered[[2]] %>% select(Code | ends_with("2")) %>% rename(name = Code) 

Tdg <- Tdg %>% activate(nodes) %>% left_join(Central2  %>% filter(name %in% namefilter), by = "name")

summary(V(Tdg)$D2)

# Add community membership on bipartite 

Member2 <- Memberem2nwF %>% filter(type == FALSE)

Tdg <- Tdg %>% activate(nodes) %>% left_join(Member2 , by = "name")

# Filter graph to get top degree (in bipartite)

TdgF <- Tdg %>% activate(nodes) %>% filter(D2 > 20) %>% filter(!node_is_isolated())

## Deal with color
ComVec <- TdgF  %>% activate(nodes) %>% fortify.tbl_graph()%>% select(fgreedy) %>% deframe() %>% unique() %>%sort()
pal <-  brewer.pal(n = 7, name = "Paired")
pal <- pal[ComVec]

##Plot the completment 
g4 <- ggraph(TdgF, layout = "kk" ) + 
  geom_edge_hive(alpha = 0.2, start_cap = circle(4, 'mm'),
                 end_cap = circle(4, 'mm'))+
  geom_node_point( aes(color = as.character(fgreedy), size = D2), alpha= 0.7  )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 3)+
  scale_color_manual(values = pal, breaks = ComVec)+
  scale_size(range = c(2,12))+
  labs(color = "Communautés dans le biparti\n(Fast Greedy sur em2nwF)", size =  "Degrés biparti (nb villes membres)",
       caption = "NB : Complément du graphe em2nwF (asso-asso)\nSpatialisation :kk\nSources : ETMUN 2019/ PG 2020")+
  theme_graph(base_family = 'Helvetica')+
  theme(legend.position="bottom", legend.box = "vertical")

ggsave(g4, filename = "OUT/Compl_Asso_em2nwF_ETMUN_deg20.pdf",width = 8.3, height = 8.3)

### ====ASSORTATIVITY  ====


# EM2NW

## degree

degAssort <- assortativity_degree(em2nwCities.g)
degAssort
# strength
degStrength <- assortativity(em2nwCities.g, V(em2nwCities.g)$strength)
degStrength
## Pop 
NameVertex <- get.vertex.attribute(em2nwCities.g)$name
NodesInfos2 <- NodesInfos %>% filter(name %in% NameVertex)
em2nwCities.g <-  set_vertex_attr(em2nwCities.g,"PopAdmin11", index = NodesInfos2$name, value = NodesInfos2$PopAdmin11)

# keep only cities with more than 5000 

g <- delete_vertices(em2nwCities.g, V(em2nwCities.g)$PopAdmin11 < 5000)

popAssort  <- assortativity(g, as.numeric(V(g)$PopAdmin11), directed = FALSE)
popAssort
## admin

em2nwCities.g <-  set_vertex_attr(em2nwCities.g,"adminLevel", index = NodesInfos2$name, value = NodesInfos2$adminLevel)

adminAssort <- assortativity_nominal(em2nwCities.g,as.numeric(as.factor(V(em2nwCities.g)$adminLevel)))

adminAssort
## Country 

em2nwCities.g <-  set_vertex_attr(em2nwCities.g,"countryCode", index = NodesInfos2$name, value = NodesInfos2$countryCode)

countryAssort <-  assortativity_nominal(em2nwCities.g,as.numeric(as.factor(V(em2nwCities.g)$countryCode)))

countryAssort
## SubRegion 

library(countrycode)

NodesInfos2$subregion <- countrycode(NodesInfos2$countryCode, "iso2c", "region")

em2nwCities.g <-  set_vertex_attr(em2nwCities.g,"subregion", index = NodesInfos2$name, value = NodesInfos2$subregion)

subregionAssort <-  assortativity_nominal(em2nwCities.g, as.numeric(as.factor(V(em2nwCities.g)$subregion)))
subregionAssort

## TYpo science Po
library(assortnet)
TypoScPo <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CountryPolitical/CountryInfo_PoliticalTypo.csv", stringsAsFactors = F)

NodesInfos2 <- NodesInfos2 %>% left_join(select(TypoScPo, countryCode = iso_a2, 11:15))

em2nwCities.g <-  set_vertex_attr(em2nwCities.g,"TypoDecentre", index = NodesInfos2$name, value = NodesInfos2$LocGovType_MunicipalAdmin)

MatrixEM2NW <- get.adjacency(em2nwCities.g) %>% as.matrix()

NameMatrix <- data.frame(name= rownames(MatrixEM2NW)) %>% left_join(select(NodesInfos2, name, LocGovType_MunicipalAdmin))
assortment.discrete(MatrixEM2NW ,NameMatrix$LocGovType_MunicipalAdmin,weighted=FALSE)

# EM2NWF

## degree

degAssort <- assortativity_degree(em2nwFCities.g)

# strength
degStrength <- assortativity(em2nwFCities.g, V(em2nwFCities.g)$strength)

## Pop 
NameVertex <- get.vertex.attribute(em2nwFCities.g)$name
NodesInfos2 <- NodesInfos %>% filter(name %in% NameVertex)
em2nwFCities.g <-  set_vertex_attr(em2nwFCities.g,"PopAdmin11", index = NodesInfos2$name, value = NodesInfos2$PopAdmin11)
library(assortnet)

# keep only cities with more than 5000 

g <- delete_vertices(em2nwFCities.g, V(em2nwFCities.g)$PopAdmin11 < 5000)

popAssort  <- assortativity(g, as.numeric(V(g)$PopAdmin11), directed = FALSE)



## admin

em2nwFCities.g <-  set_vertex_attr(em2nwFCities.g,"adminLevel", index = NodesInfos2$name, value = NodesInfos2$adminLevel)

adminAssort <- assortativity_nominal(em2nwFCities.g,as.numeric(as.factor(V(em2nwFCities.g)$adminLevel)))

## Country 

em2nwFCities.g <-  set_vertex_attr(em2nwFCities.g,"countryCode", index = NodesInfos2$name, value = NodesInfos2$countryCode)

countryAssort <-  assortativity_nominal(em2nwFCities.g,as.numeric(as.factor(V(em2nwFCities.g)$countryCode)))



# EMNW

emnwCities.g <- emnwDual$igraph_object[[4]]
rm(emnwDual)
## degree

degAssort <- assortativity_degree(emnwCities.g )

V(emnwCities.g)$degree <- degree(emnwCities.g, loops = F )
# strength
V(emnwCities.g)$strength <- strength(emnwCities.g, vids = V(emnwCities.g), loops = F )
degStrength <- assortativity(emnwCities.g , V(emnwCities.g )$strength)

## Pop 
NameVertex <- get.vertex.attribute(emnwCities.g )$name
NodesInfos2 <- NodesInfos %>% filter(name %in% NameVertex)
emnwCities.g  <-  set_vertex_attr(emnwCities.g ,"PopAdmin11", index = NodesInfos2$name, value = NodesInfos2$PopAdmin11)


# keep only cities with more than 5000 

g <- delete_vertices(emnwCities.g , V(emnwCities.g )$PopAdmin11 < 5000)

popAssort  <- assortativity(g, as.numeric(V(g)$PopAdmin11), directed = FALSE)

# Pop discrete 
summary(NodesInfos$PopAdmin11)
NodesInfos <- NodesInfos %>% filter(!is.na(PopAdmin11))
labelClass <- c( "7.Very Small City", "6.Small City", "5.Medium-sized Small City", "4.Medium-sized City","3.Large City",  "2.Very Large City", "1.Metropolis")



cesDonnees <- NodesInfos$PopAdmin11
min<-min(cesDonnees)
max<- max(cesDonnees)
valBreaks <- c(min,10000,25000,50000,100000, 500000, 1000000, max)

NodesInfos$ClassPop11 <- cut(cesDonnees,
                             breaks = valBreaks,
                             labels = labelClass,
                             include.lowest = TRUE,
                             right= FALSE)

NodesInfos2 <- NodesInfos %>% filter(name %in% NameVertex)
emnwCities.g <-  set_vertex_attr(emnwCities.g,"PopClass11", index = NodesInfos2$name, value = NodesInfos2$ClassPop11)

popclassAssort <- assortativity_nominal(emnwCities.g, as.numeric(as.factor(V(emnwCities.g)$PopClass11)))


g <- delete_vertices(emnwCities.g , V(emnwCities.g )$PopAdmin11 < 5000)

popclassAssort <- assortativity_nominal(g, as.numeric(as.factor(V(g)$PopClass11)))

## admin

emnwCities.g <-  set_vertex_attr(emnwCities.g,"adminLevel", index = NodesInfos2$name, value = NodesInfos2$adminLevel)

adminAssort <- assortativity_nominal(emnwCities.g,as.numeric(as.factor(V(emnwCities.g)$adminLevel)))

## Country 

emnwCities.g <-  set_vertex_attr(emnwCities.g,"countryCode", index = NodesInfos2$name, value = NodesInfos2$countryCode)

countryAssort <-  assortativity_nominal(emnwCities.g,as.numeric(as.factor(V(emnwCities.g)$countryCode)))

## SubRegion 

library(countrycode)

NodesInfos2$subregion <- countrycode(NodesInfos2$countryCode, "iso2c", "region")

emnwCities.g <-  set_vertex_attr(emnwCities.g,"subregion", index = NodesInfos2$name, value = NodesInfos2$subregion)

subregionAssort <-  assortativity_nominal(emnwCities.g, as.numeric(as.factor(V(emnwCities.g)$subregion)))

##### Degree distribution


CentralitiesEmnwCities <- igraph::as_data_frame(emnwCities.g, what = "vertices")

CentralitiesEmnwCities %>% keep(is.numeric) %>% select(degree, strength) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram() + 
  scale_y_log10() + 
   scale_x_log10()+
  labs(title = "Distributions des degrés des villes ETMUN",
       subtitle = "graphe villes-villes (à partir de emnw)", 
       y = "Nombre de villes (log10)", x = "Degrés (log10)",
       caption = "Sources : ETMUN 2019 / PG 2020")
