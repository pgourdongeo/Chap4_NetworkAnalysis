############################### ETMUN bipartite Networks STEP 2 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe ETMUN (affiliation) biparti. 
#               Indices locaux,distributions des degrés, ego network
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")

### Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(Matrix)
library(patchwork)
library(GGally)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
library(ggrepel)

options(scipen = 999)

### Data 

## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_europe.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", 
                       stringsAsFactors = F)

## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")


# Matrix from STEP1  (REMIND WWCAM is removed a priori because only one european city connected to the whole graph)

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/Etmun_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/Etmun_em2nw.rds")



### ==== DISTRIBUTION OF DEGREES and Centrality Measures ====

# Function centalities indexes

Centralitydf <- function(matrix1, matrix2){
  require(igraph)
  
  matrix1.g <- graph.incidence(matrix1, directed = FALSE)
  city <- data.frame(D = colSums(matrix1), 
                     ND = colSums(matrix1)/nrow(matrix1), 
                     geonameId = colnames(matrix1), 
                     stringsAsFactors = FALSE)
  asso <- data.frame(D = rowSums(matrix1), 
                     ND = rowSums(matrix1)/ncol(matrix1), 
                     Code = rownames(matrix1),
                     stringsAsFactors = FALSE)
  
  V(matrix1.g)$B <- betweenness(matrix1.g, normalized = FALSE, directed = FALSE)
  V(matrix1.g)$NB <- betweenness(graph = matrix1.g, normalized = TRUE, directed = FALSE)
  V(matrix1.g)$C <- closeness(graph = matrix1.g, normalized = FALSE)
  V(matrix1.g)$NC <- closeness(graph = matrix1.g, normalized = TRUE)
  V(matrix1.g)$N_Eigen_C <- eigen_centrality(graph = matrix1.g, scale = TRUE, directed = FALSE)$vector

  
  df1 <- as.data.frame(vertex_attr(matrix1.g), stringsAsFactors =FALSE)
  cityBC <- df1[df1$type== TRUE,]
  assoBC <- df1[df1$type== FALSE,]
  
  city <- city %>% left_join(cityBC, by= c("geonameId" = "name")) %>% select(-type)
  asso <- asso %>% left_join(assoBC, by= c("Code" = "name")) %>% select(-type)
  
  
  matrix2.g <- graph.incidence(matrix2, directed = FALSE)
  city2 <- data.frame(D2 = colSums(matrix2), 
                      ND2 = colSums(matrix2)/nrow(matrix2), 
                      geonameId = colnames(matrix2),
                      stringsAsFactors = FALSE)
  asso2 <- data.frame(D2 = rowSums(matrix2), 
                      ND2 = rowSums(matrix2)/ncol(matrix2), 
                      Code = rownames(matrix2), 
                      stringsAsFactors = FALSE)
  
  
  V(matrix2.g)$B2 <- betweenness(matrix2.g, normalized = FALSE, directed = FALSE)
  V(matrix2.g)$NB2 <- betweenness(graph = matrix2.g, normalized = TRUE, directed = FALSE)
  V(matrix2.g)$C2 <- closeness(graph = matrix2.g, normalized = FALSE)
  V(matrix2.g)$NC2 <- closeness(graph = matrix2.g, normalized = TRUE)
  V(matrix2.g)$N_Eigen_C2 <- eigen_centrality(graph = matrix2.g, scale = TRUE, directed = FALSE)$vector
  
  df2 <- as.data.frame(vertex_attr(matrix2.g), stringsAsFactors = FALSE)
  cityBC2 <- df2[df2$type== TRUE,]
  assoBC2 <- df2[df2$type== FALSE,]
  
  city2 <- city2 %>% left_join(cityBC2, by= c("geonameId" = "name")) %>% select(-type)
  asso2 <- asso2 %>% left_join(assoBC2, by= c("Code" = "name")) %>% select(-type)
  
  city <- merge(city, city2, by = "geonameId", all.x = T)
  asso <- merge(asso, asso2, by = "Code", all.x = T)

  
  result <- list(city, asso)
  return(result)
}

Centralities <- Centralitydf(emnw, em2nw)

saveRDS(Centralities, "DataProd/Centralities_ETMUN_enww_em2nw.rds")
### plot of degree distributions

# hist

Centralities[[1]] %>% keep(is.numeric) %>% select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram() + 
  scale_y_log10() + 
  # scale_x_log10()+
  labs(title = "Distributions des mesures de centralité des villes ETMUN",
       subtitle = "Sur les matrices emnw et em2nw", 
       y = "Nombre de villes (log10)", x = "",
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSources : ETMUN 2019 / PG 2020")

ggsave(filename = "OUT/CentralDistrib_2matrix_Cities_ETMUN.pdf",width = 8.3, height = 5.8 )

Centralities[[2]] %>% keep(is.numeric) %>% select(starts_with("N")) %>% 
  select(-N_Eigen_C , -N_Eigen_C2 ) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram()+ 
  # scale_y_log10() + 
  scale_x_log10(n.breaks = 6) +
  labs(title = "Distributions des mesures de centralité des associations ETMUN",
       subtitle = "Sur les matrices emnw et em2nw", 
       y = "Nombre d'associations", x = "Valeurs normalisées (échelle log 10)",
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation,  2 = em2nw\nSources : ETMUN 2019 / PG 2020")

ggsave(filename = "OUT/CentralDistrib_2matrix_Asso_ETMUN.pdf",width = 8.3, height = 5.8 )

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes ETMUN (emnw)")

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(ends_with("2") & starts_with("N" )) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes ETMUN (em2nw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (emnw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>%
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (em2nw)")


## cumulative normalized degree
ytext <- "Fréquence cumulée (log10)"
xtext <- "Degré normalisé (log10)"

g1 <-ggplot(Centralities[[1]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des villes ETMUN (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g2 <- ggplot(Centralities[[1]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des villes ETMUN (em2nw)", y = ytext,  x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g3 <- ggplot(Centralities[[2]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des associations ETMUN (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g4 <- ggplot(Centralities[[2]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des associations ETMUN (em2nw)", y = ytext, x = xtext,
       caption = "Sources : ETMUN 2019\nPG 2020") + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()


grid1 <- (g1 | g2) / (g3 | g4)
grid1
rm(g1,g2,g3,g4)


### ==== POWER LAW OF DEGREE DISTRIBUTION ====

# 
# fit1 <- fit_power_law(Centralities[[1]]$ND, implementation = c("plfit"))
# fit1
# ND2 <- Centralities[[1]] %>% filter(!is.na(ND2))%>% select(ND2)%>% deframe()
# fit2 <- fit_power_law(ND2, implementation = c("plfit"))
# fit2
# fit3 <- fit_power_law(Centralities[[2]]$ND, implementation = c("plfit"))
# fit3
# 
# ND2 <- Centralities[[2]] %>% filter(!is.na(ND2))%>% select(ND2)%>% deframe()
# fit4 <- fit_power_law(ND2, implementation = c("plfit"))
# fit4
# 
# g <- graph.incidence(emnw, directed = FALSE)
# d <- degree(g, normalized = FALSE) %>% as.data.frame()
# colnames(d) <- "degree"
# ggplot(d)+  geom_histogram(aes(degree))+ scale_x_log10()+ scale_y_log10()
# fit5 <- fit_power_law(d, implementation = c("plfit"))
# fit5
# 
# fit6 <- fit_power_law(append(Centralities[[2]]$ND, Centralities[[1]]$ND), implementation = c("plfit"))
# fit6
# 
# fit7 <-  fit_power_law(append(Centralities[[2]]$D, Centralities[[1]]$D), implementation = c("plfit"))
# fit7
# rm(ND2, fit1,fit2,fit3,fit4)
# 
# library(poweRlaw)
# 
# d2 <- append(Centralities[[2]]$ND, Centralities[[1]]$ND) %>% as.data.frame() 
# colnames(d2) <- "degBi"
# 
# m_pl =conpl$new(d2$degBi)         # discrete power law fitting
# est = estimate_pars(m_pl)# get xmin and alpha
# m_pl$setPars(est)
# bs_p = bootstrap_p(m_pl, no_of_sims = 10)  
# bs_p$p
# 
# 
# 
# m_pl = displ$new(d$degree)         # discrete power law fitting
# est = estimate_pars(m_pl)# get xmin and alpha
# m_pl$setPars(est)
# 
# # here we have the goodness-of-fit test p-value
# # as proposed by Clauset and al. (2009)
# bs_p = bootstrap_p(m_pl)  
# bs_p$p
# 
# m_ln = dislnorm$new(d$degree)
# m_ln$setXmin(m_pldeg$getXmin())
# est = estimate_pars(m_ln)
# m_ln$setPars(est)
# 
# m_ex <- disexp$new(d$degree)
# m_ex$setXmin(m_pldeg$getXmin())
# est = estimate_pars(m_ex)
# m_ex$setPars(est)
# 
# comp = compare_distributions(m_pldeg, m_ex)





#Set city name
namecity <- MembershipEtmun %>% 
  select(geonameId, asciiName) %>% 
  distinct() %>% 
  left_join(select(DbCity,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11))
Centralities[[1]] <- Centralities[[1]] %>% left_join(namecity)

# set names of association
Centralities[[2]] <- Centralities[[2]] %>% left_join(select(AssoEtmun, Code,  AssoName = Name))


### ==== TOP CENTRALITIES ====

### top N for each type of nodes
Top50Cities <- Centralities[[1]] %>% arrange_if(is.numeric, desc)%>% head(50)
Top10Asso <- Centralities[[2]] %>% arrange_if(is.numeric, desc)%>% head(10)

### Rank 
CityCentralities <- Centralities[[1]] %>% arrange(desc(D)) %>%  select( geonameId | asciiName | starts_with("N"))

RankCities <- CityCentralities %>% 
  mutate_if(is.numeric, 
            list( rank = ~ rank(desc(.), ties.method = 'min'))) 

# cities with a rank < 5 in any index
TopRank <- RankCities %>% filter_at(.vars = c(11:18), any_vars(.<5))

# tibble for ggplot rank
TopRankgg <- TopRank %>% select(geonameId | asciiName | ends_with("rank"))%>% 
  pivot_longer(-c(geonameId,asciiName), names_to = "indice", values_to = "rang")%>% 
  mutate(indice = str_remove(indice,pattern = "_rank")) %>% 
  mutate(indice = str_remove(indice,pattern = "N_")) %>%
  mutate(lograng = log10(rang))

#plot in coord polar
s <- summary(TopRankgg$lograng)
bks <- c(s[[1]], log10(2), s[[2]], s[[3]], s[[5]], log10(30), log10(50), s[[6]])
lbs <-  round(10^bks,0)

TopRankgg %>% 
  ggplot(aes(x = indice, y = lograng, group = asciiName, label = rang) ) + 
  geom_line() +
  geom_point()+
  geom_label_repel(aes( fill = indice),size = 3, label.padding = 0.2)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap(~asciiName) + 
  scale_y_continuous(trans = "reverse", 
                     breaks = bks, 
                    labels = lbs)+
  coord_polar() + 
  labs(x = "", y = "Rang (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSources : ETMUN 2019 / PG 2020") + 
  theme(axis.text.y=element_text(size=7),
        axis.text.x=element_blank(), legend.position="bottom")

ggsave(filename = "OUT/RankClock_TopAny4_ETMUN.pdf",width = 8.3, height = 7 )

rm(CityCentralities,TopRank,TopRankgg , s, bks, lbs)
#saveRDS(Top50Cities, "DataProd/Top50Cities_CentralitiesETMUN.rds")

#### ====ARTICULATION POINTS AND IGRAPH OBJECT ====

emnw.g <- graph.incidence(emnw, directed = F)

em2nw.g <-  graph.incidence(em2nw, directed = F)

#Info all nodes 
NameAsso <- AssoEtmun %>% select(name = Code,label = Name, Country=Country..secretariat., Acro =Acronym, Date )
nameCitiesjoint <- namecity %>% select(name= geonameId, label = asciiName, Country = countryCode, everything())
NodesInfos <- bind_rows(NameAsso,nameCitiesjoint) %>% filter(!duplicated(name)) ## Get the real info from the GN table (for country)
#saveRDS(NodesInfos, "DataProd/NodesInfos_ETMUN.RDS")
rm(nameCitiesjoint)


## Articulation points
#emnw
art <-articulation_points(emnw.g)
vecArt <- V(emnw.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)## 46 articulations points (only associations)

#em2nw
articulation_points(em2nw.g)# none

#components
components(em2nw.g)$no  # 1
components(emnw.g)$no  # 1

rm(art,vecArt)

#### ==== FILTER OUT SOME ASSOCIATIONS ====
## show outliers by degrees (2 associations)


#filter some asso outliers (visual detection) according to the ggpairs on Associations

Centralities[[2]]%>% filter(ND>0.10)# Covenant of mayors, Climate Alliance, WWCAM

Centralities[[2]] %>% keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% filter(ND<0.10 )%>%
  ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (emnw)")


# Covenant of mayors, Climate Alliance, EWT (subset of another asso)
assoremove <- c("04080",  "03241", "07675")
#assoremove <- c("Covenant of Mayors", "WWCAM","Climate Alliance")

#Filter out the association on the matrix that keep city of degree 1
emnwF <- emnw[!rownames(emnw) %in% assoremove, ]
emnwF <- emnwF[,colSums(emnwF) > 0]

#Filter out the association on the matrix that keep city of more than 1 degree (only cities involved at least in 2 associations)
em2nwF <-  em2nw[!rownames(em2nw) %in% assoremove, ]
em2nwF <- em2nwF[,colSums(em2nwF) > 1]


#saveRDS(emnwF, "DataProd/emnwF.RDS")
#saveRDS(em2nwF, "DataProd/em2nwF.RDS")

### Graph and Articulations points
#emnwF

emnwF.g <- graph.incidence(emnwF, directed = F)
components(emnwF.g)$no


art <-articulation_points(emnwF.g)
vecArt <- V(emnwF.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)## 49 articulations points (only associations)

#em2nwF
em2nwF.g <- graph.incidence(em2nwF, directed = F)
articulation.points(em2nwF.g)#none


### ==== FILTERED MATRIX DISTRIBUTION OF DEGREES and Centrality Measures ====

# Function from STEP1

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

# Summary of filtered matrices
dimemnwF <- graphlevelindex(emnwF)

dimem2nwF <- graphlevelindex(em2nwF)

#Centralities Filtered

CentralityFiltered <- Centralitydf(emnwF, em2nwF)

saveRDS(CentralityFiltered, "DataProd/Centralities_ETMUN_enwwF_em2nwF.rds")
## Histogram centralities

CentralityFiltered[[1]] %>% 
  keep(is.numeric) %>% 
  select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des degrés des villes ETMUN",
       subtitle = "Sur les matrices emnwF et em2nwF", 
       y = "Nombre de villes (log10)", x = "",
       caption =  "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSources : ETMUN 2019 / PG 2020")

ggsave(filename = "OUT/CentralDistrib_2matrixF_Cities_ETMUN.pdf",width = 8.3, height = 5.8 )

CentralityFiltered[[2]] %>% 
  keep(is.numeric) %>% 
  select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  # scale_y_log10() + 
  labs(title = "Distributions des degrés des Associations ETMUN",
       subtitle = "Sur les matrices emnwF et em2nwF", 
       y = "Nombre de villes (log10)", x = "",
       caption =  "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSources : ETMUN 2019 / PG 2020")

ggsave(filename = "OUT/CentralDistrib_2matrixF_Asso_ETMUN.pdf",width = 8.3, height = 5.8 )

#### ==== EGO NETWORK  =====

# Make ego network
# ego.emnw <- make_ego_graph(emnw.g, order = 2, nodes = V(emnw.g))
# ego.emnwF <- make_ego_graph(emnwF.g, order = 2, nodes = V(emnwF.g))
ego.em2nwF <- make_ego_graph(em2nwF.g, order = 2, nodes = V(em2nwF.g))

## Function to compute numeric summary on each ego network
egoDim <- function(egoglist, graph){
require(igraph)
  require(tnet)
  d <-NULL

for (i in 1:length(egoglist)){    #boucle qui a le nombre d'itérations correspondant à la longueur de la liste des ego
  idList <- i 
  name <- V(graph)$name[[i]]  #id
  type <- V(graph)$type[[i]]
  egi <- delete.vertices(egoglist[[i]], V(egoglist[[i]])[name == V(graph)$name[[i]]]) #remove ego
  nbAsso <- length(V(egi)[type == FALSE]) # Vcount type 1
  nbVilles <-length(V(egi)[type == TRUE])# Vcount type 2
  nbliens <- ecount(egi)                #edges count
  densite <- nbliens/(nbVilles*nbAsso)         #bipartite density
  comp <- components(egi)               # connex composante
  nbcomp <- comp$no                     #nb components
  sizecomp <- 100*comp$csize[[1]]/vcount(egi) #% vertices of the larger component
  
  elist <- get.edgelist(egi) %>% as.data.frame() %>% select(V2,V1)# get edge list to create a tnet object (inverse V1 and V2 to get a city-association list)
  tn <- as.tnet(elist, type="binary two-mode tnet") # tn object
  reinf <- reinforcement_tm(tn)
  # clust4cycleCities <- clustering_tm(tn) #clustering coeff of cities (take a huge amount of time)
  
  d <- rbind(d, data.frame(name,idList, type, nbAsso,nbVilles,
                           nbliens, densite, nbcomp, sizecomp, reinf)) 
}
return(d)
}

# result in df (keep only cities)

#ego.emnw.df <-  egoDim(ego.emnw , emnw.g)
#ego.emnwF.df <-  egoDim(ego.emnwF, emnwF.g)
ego.em2nwF.df <- egoDim(ego.em2nwF, em2nwF.g)

#
# egocity.emnw <- ego.emnw.df %>% filter(type == TRUE) %>% 
#   left_join(NodesInfos)
# 
# egocity.emnwF <- ego.emnwF.df %>% filter(type == TRUE) %>% 
#   left_join(NodesInfos)

egocity.em2nwF <- ego.em2nwF.df %>% filter(type == TRUE) %>% 
  left_join(NodesInfos)
saveRDS(egocity.em2nwF, "DataProd/egocity_em2nwF.RDS")
rm(ego.emnw, ego.emnw.df, ego.emnwF.df, ego.em2nwF.df)

### Exploration of ego network properties

egocity.em2nwF %>% keep(is.numeric) %>% select(-idList, -Date) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des propriétés des ego-network de chaque ville ETMUN",
       subtitle = "Sur la matrice em2nwF", 
       y = "Nombre de villes (log10)",
       caption = "ego-network d'ordre 2")

egocity.em2nwF %>% 
  keep(is.numeric) %>% select(-idList,  -Date) %>%
  select(-nbcomp, -sizecomp)%>%filter(PopAdmin11 < 7000000)%>%
  ggpairs(title = "Relations entre propriétés des ego-network de chaque ville ETMUN (em2nwF)")

#### ==== EGO NETWORK PROPERTIES CLUSTERING  =====


### PCA 
library(FactoMineR)
library(explor)

DfACP <- egocity.em2nwF %>% 
          select(c(1,4:10)) %>% ## keep properties of ego net
          mutate(reinf = replace(reinf, is.nan(reinf), 0))
          



rownames(DfACP) <- DfACP$name
DfACP <- DfACP %>% select(-name)
res.pca <- PCA(DfACP,
               scale.unit = TRUE,
               #ind.sup = 7,
               #quanti.sup = 7,
               # quali.sup = 5,
               graph = FALSE)



explor(res.pca)

CoordPCA<- as.data.frame(res.pca$ind$coord)

DfACP <- cbind(DfACP , CoordPCA)

DfACP <- DfACP %>% mutate(name = rownames(.)) 

## add PCA results

egocity.em2nwF <- egocity.em2nwF %>% left_join(select(DfACP, name, Dim.1, Dim.2, Dim.3))


## CAH on dim ACP 
library(cluster)
library(ggdendro)
library(reshape2)

suppressWarnings(source("../../multi_bivariateFunctions_exploratR/Functions_anova_cah.R"))



myVar <- c("Dim.1","Dim.2", "Dim.3")
cah <- ComputeClassif(df = egocity.em2nwF,
                      varquanti = myVar, method = "ward", stand = FALSE)


dendro <- PlotDendro(classifobj = cah)
dendro
inert <- PlotHeight(classifobj = cah)
inert
myProfiles <- PlotProfile(classifobj = cah, nbclus = 5)
myProfiles$PROFILE
table(myProfiles$CLUSID)
egocity.em2nwF <- egocity.em2nwF %>% mutate(CAH = myProfiles$CLUSID)



### Explore clustering outcomes in relation to population and administrative status of cities



### Test anova Classe and Pop

# correct population variable
egocity.em2nw <- egocity.em2nw %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

# Perform AnonVa

anovaTab <- AnovaTab(df = egocity.em2nw , varx = c("CAH"), 
                     vary = c("PopAdmin11")) 
resume <- ComputeRegression(egocity.em2nw, vardep = "PopAdmin11", varindep = "CAH")
resume$TABCOEF ## R2 = 0.0
rm(resume, anovaTab)
# ---- CHI2 ----

tableCAH <- table(egocity.em2nw$CAH, egocity.em2nw$adminLevel)



testChi2 <- chisq.test(tableCAH)
testChi2$statistic
testChi2$expected
testChi2$residuals
resChi2<-as.data.frame(testChi2$residuals)
maxFreq <- max(resChi2$Freq)
minFreq <- min(resChi2$Freq)
ggplot(resChi2, aes(x = Var1, y = Var2, fill= Freq ))+ geom_tile() + scale_fill_gradient2(low = "#67a9cf", high = "#ef8a62", mid = "#f7f7f7", 
                                                                                          midpoint = 0, limit = c(minFreq,maxFreq)) + 
   labs(title = "Résidus de Pearson du Test Chi 2", 
        x = "Classification des ego-networks", 
        y = "Niveau administratif des localités (geonames)",
        fill = "Résidus de Pearson")
  
ggsave(filename = "OUT/HeatMap_PearsonChi2_CAHegoNetandAdmin_ETMUN.pdf", device = "pdf",width = 8.3, height = 5.8 , units = "in")
class(resChi2)
testChi2$stdres

library(vcd)

cramer <-assocstats(tableCAH)
cramer  # 0.275

#### ==== MAPPING SOME EGO NETWORKS  =====


# create a sample of egonetworks

sampleEgo <- egocity.em2nw  %>% group_by(CAH) %>% sample_n(1, replace = T)
ego.em2nw[[1]]
# function to map all ego-networks

EgoMapping <- function(df, egolist, NodesInfos){
  require(ggraph)
  require(tidygraph)
  require(RColorBrewer)
  
  listplot <- list()
  
  for (i in 1:nrow(df)){
    # get the ego graph from the ego list
    
    local({    # add local because on saving plots in a list otherwise ggplot evaluate the global environnement and there is a mismatch between data and plot specifications
      i <- i;
      g <- egolist[[df[i,]$idList]] ;
      
      # get info of ego from df
      strInfo <- df[i,] ;
      
      # transfor into tidygraph
      tdg <-  as_tbl_graph(g) ;
      
      #get the names of all nodes
      namefilter <- tdg %>% activate(nodes) %>% fortify.tbl_graph()%>% select(name) %>% deframe();
      
      #joint info nodes
      tdg <- tdg %>% activate(nodes) %>% left_join(NodesInfos  %>% filter(name %in% namefilter), by = "name");
      
      # prepare type of nodes for aes
      tdg <- tdg %>% activate(nodes)%>%
        mutate(type2 = case_when(type == FALSE ~ "Asso", type == TRUE ~ "City"))%>% 
        mutate(type2 = ifelse(name == df[i,]$name, "ego", type2))%>%
        mutate(size = case_when(type2 == "Asso" ~ 4,type2 == "City" ~ 2, type2 == "ego" ~ 3 ));
      
      # plot
      
     
      plot <- ggraph(tdg, layout = "kk")+
        geom_edge_link(alpha = 0.2)+
        geom_node_point(aes(color = V(tdg)$type2, shape = V(tdg)$type, size = V(tdg)$size),
                        show.legend = FALSE, alpha = 0.8) + 
        scale_color_manual( values =c("lightskyblue", "mediumvioletred", "orange"), breaks = c("City", "ego", "Asso"), labels = c("Villes", "Ego", "Projets"))+
        scale_size(range = c(1.5,6))+
        labs(title = paste("Ego Network de ", strInfo$label,", ", strInfo$Country, sep = ""),
             subtitle = paste(strInfo$CAH), color = "Type de sommets",
             caption = "Algorithme de spatialisation : kk")+
        annotate("label",
                 label = paste("Densité : ", round(strInfo$densite,3),
                               "\nNb Villes : ", strInfo$nbVilles,
                               "\nNb composantes : ", strInfo$nbcomp, 
                               "\nReinforcement : ", round(strInfo$reinf,4), 
                               sep = ""), 
                 x = min(layout_with_kk(tdg)[,1]), 
                 y = min(layout_with_kk(tdg)[,1])+2.5,
                 hjust = 0, 
                 alpha =0.4,
                 size = 2.5)+
        theme_graph(base_family = 'Helvetica');
      
      
      Node <- strInfo$name;
      
      
      listplot[[Node]] <<- plot
      #print(plot)
    }) 
  }
  
  return(listplot)
  
}

#Apply function egoMapping
egoplot.em2nw <- EgoMapping(df = sampleEgo,egolist = ego.em2nw,NodesInfos =NodesInfos)

egoplot.em2nw[[5]]




# Grid plot function

Gridplot <- function(myplots, n){
  require(cowplot)
  splitted_plots <- split(myplots, ceiling(seq_along(myplots)/n))
  
  lapply(splitted_plots, function(x) plot_grid(plotlist = x))

}


#Try
gEgo <- Gridplot(egoplot.em2nw, 4)

gEgo[[1]]
ggsave( filename = "OUT/ETMUNEGONetwork_cluster5.pdf", width = 11.7, height = 8.3, units = "in" )

gEgo[[2]]

ggsave( filename = "OUT/ETMUNEGONetwork_cluster5_UniqueClasse5.pdf", width = 8.3, height = 5.8 , units = "in" )
