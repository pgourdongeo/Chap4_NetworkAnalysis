############################### URBACT bipartite Networks STEP 2 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe URBACT (affiliation) biparti. 
#               Indices locaux, distributions des degrés, ego network
#
# 
############################################################################## PG Aout 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.5.URBACT/")

### Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(Matrix)
library(patchwork)
library(GGally)
library(tidygraph)
library(ggraph)
library(skimr)
library(ggrepel)

options(scipen = 999)

### Data 

## Participation (partnership): table City-Project. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
ParticipationURBACT <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/URBACT_Membership_GNidCorr_complete.csv", stringsAsFactors = F)

## Information on project (for nodes)
ProjectURBACT <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/URBACT/UrbactNetworks_complete.csv", 
                            stringsAsFactors = F)

## GN info for cities
DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")

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
bibi1 <-Centralities[[1]]
bibi2 <-Centralities[[2]]
saveRDS(Centralities, "DataProd/Centralities_URBACT_enww_em2nw.rds")
### plot of degree distributions

# hist
Centralities[[1]] %>% keep(is.numeric) %>% select(starts_with("N")&- ends_with("2")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram() + 
  scale_y_log10() + 
  # scale_x_log10()+
  labs(title = "Distributions des mesures de centralité des villes URBACT",
       subtitle = "Sur la matrice emnw", 
       y = "Nombre de villes (log10)", x = "",
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")


summary(Centralities[[1]]$ND)
summary(Centralities[[1]]$D)
ggsave(filename = "OUT/CentralDistrib_Cities_URBACT.pdf",width = 8.3, height = 5.8 )

Centralities[[2]] %>% keep(is.numeric) %>% select(starts_with("N")&- ends_with("2")) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram()+ 
  scale_y_log10() + 
  # scale_x_log10(n.breaks = 6) +
  labs(title = "Distributions des mesures de centralité des projets URBACT",
       subtitle = "Sur la matrice emnw", 
       y = "Nombre de projets (échelle log 10)", x = "Valeurs normalisées",
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

summary(Centralities[[2]]$ND)
summary(Centralities[[2]]$D)
 ggsave(filename = "OUT/CentralDistrib_Project_URBACT.pdf",width = 8.3, height = 5.8 )

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes URBACT (emnw)")

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes URBACT (em2nw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets URBACT (emnw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>%
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets URBACT (em2nw)")


#Set city name

namecity <- ParticipationURBACT %>% select(geonameId, asciiName) %>% distinct()%>% 
  left_join(select(DbCity,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11))
Centralities[[1]] <- Centralities[[1]] %>% left_join(namecity)

# set names of association
Centralities[[2]] <- Centralities[[2]] %>% left_join(select(ProjectURBACT, Code = Code_Network,  ProjectName = Name))

### ==== TOP CENTRALITIES ====

### top N for each type of nodes
Top50Cities <- Centralities[[1]] %>% arrange_if(is.numeric, desc)%>% head(50)
Top50Projects <- Centralities[[2]] %>% arrange_if(is.numeric, desc)%>% head(50)



saveRDS(Top50Cities, "DataProd/Top50Cities_Centralities_EUCICOP.rds")

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
bks <- c(s[[1]], log10(2), s[[2]], s[[3]], s[[5]] , s[[6]])
lbs <-  round(10^bks,0)

TopRankgg %>% 
  ggplot(aes(x = indice, y = lograng, group = asciiName, label = rang) ) + 
  geom_line() +
  geom_point()+
  geom_label_repel(aes( fill = indice),size = 3, label.padding = 0.2)+
  scale_fill_brewer(palette = "Paired")+
  facet_wrap(~asciiName) + 
  scale_y_continuous(trans = "reverse", 
                     breaks = bks, 
                     labels = lbs)+
  coord_polar() + 
  labs(x = "", y = "Rang (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020") + 
  theme(axis.text.y=element_text(size=7),
        axis.text.x=element_blank(), legend.position="bottom")

ggsave(filename = "OUT/RankClock_TopAny5_URBACT.pdf",width = 8.3, height = 7 )

#### ====ARTICULATION POINTS AND IGRAPH OBJECT ====

emnw.g <- graph.incidence(emnw, directed = F)

em2nw.g <-  graph.incidence(em2nw, directed = F)

#Info all nodes 
NameProjects <- ProjectURBACT %>% select(name = Code_Network,label = Name, Date = Start )
nameCitiesjoint <- namecity %>% select(name= geonameId, label = asciiName, Country = countryCode, everything())

NodesInfos <- bind_rows(NameProjects,nameCitiesjoint) 
 # saveRDS(NodesInfos, "DataProd/NodesInfos_URBACT.RDS")
rm(nameCitiesjoint)
## Articulation points
#emnw
art <-articulation_points(emnw.g)
vecArt <- V(emnw.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)

skim(art)## 119 articulations points (1 is a city : Kavala, GR)
#em2nw

art <-articulation_points(em2nw.g)
vecArt <- V(em2nw.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)
skim(art) ## 1 (1 is a city : Kavala, GR)


#### ==== EGO NETWORK  =====

# Make ego network
ego.emnw <- make_ego_graph(emnw.g, order = 2, nodes = V(emnw.g))
ego.em2nw <- make_ego_graph(em2nw.g, order = 2, nodes = V(em2nw.g))

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
    nbProjects <- length(V(egi)[type == FALSE]) # Vcount type 1
    nbVilles <-length(V(egi)[type == TRUE])# Vcount type 2
    nbliens <- ecount(egi)                #edges count
    densite <- nbliens/(nbVilles*nbProjects)         #bipartite density
    comp <- components(egi)               # connex composante
    nbcomp <- comp$no                     #nb components
    sizecomp <- 100*comp$csize[[1]]/vcount(egi) #% vertices of the larger component
    
    # with ego for reinforcement
    elist <- get.edgelist(egoglist[[i]]) %>% as.data.frame() %>% select(V2,V1)# get edge list to create a tnet object (inverse V1 and V2 to get a city-association list)
    tn <- as.tnet(elist, type="binary two-mode tnet") # tn object
    reinf <- reinforcement_tm(tn)
    # clust4cycleCities <- clustering_tm(tn) #clustering coeff of cities (take a huge amount of time)
    
    d <- rbind(d, data.frame(name,idList, type, nbProjects, nbVilles,
                             nbliens, densite, nbcomp, sizecomp, reinf)) 
  }
  return(d)
}

# result in df (keep only cities)

ego.emnw.df <-  egoDim(ego.emnw , emnw.g)

ego.em2nw.df <- egoDim(ego.em2nw, em2nw.g)

#
egocity.emnw <- ego.emnw.df %>% filter(type == TRUE) %>% 
  left_join(NodesInfos)
# 
# egocity.emnwF <- ego.emnwF.df %>% filter(type == TRUE) %>% 
#   left_join(NodesInfos)

egocity.em2nw <- ego.em2nw.df %>% filter(type == TRUE) %>% 
  left_join(NodesInfos)

saveRDS(egocity.emnw, "DataProd/Df_egocity_emnw.RDS")
# saveRDS(egocity.em2nw, "DataProd/Df_egocity_em2nw.RDS")
rm(ego.emnw, ego.emnw.df, ego.emnwF.df, ego.em2nwF.df)


egocity.emnw <- readRDS("DataProd/Df_egocity_emnw.RDS")
egocity.em2nw <-readRDS("DataProd/Df_egocity_em2nw.RDS")
### Exploration of ego network properties

egocity.emnw %>% keep(is.numeric) %>% select(-idList, -Date) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des propriétés des ego-network de chaque ville EUCICOP",
       subtitle = "Sur la matrice emnw", 
       y = "Nombre de villes (log10)",
       caption = "ego-network d'ordre 2")

egocity.emnw %>% 
  keep(is.numeric) %>% select(-idList, -Date) %>%
  select(-nbcomp, -sizecomp)%>%
  ggpairs(title = "Relations entre ropriétés des ego-network de chaque ville EUCICOP (emnw)")


egocity.em2nw %>% keep(is.numeric) %>%  select(-idList, -Date) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des propriétés des ego-network de chaque ville EUCICOP",
       subtitle = "Sur la matrice em2nw", 
       y = "Nombre de villes (log10)",
       caption = "ego-network d'ordre 2")

egocity.em2nw %>% 
  keep(is.numeric) %>%  select(-idList, -Date) %>%
  select(-nbcomp, -sizecomp)%>%
  ggpairs(title = "Relations entre propriétés des ego-networks de chaque ville EUCICOP (em2nw)")


#### ==== EGO NETWORK PROPERTIES CLUSTERING  =====


### PCA 
library(FactoMineR)
library(explor)

DfACP <- egocity.emnw %>% 
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

#egocity.em2nw <- egocity.em2nw %>% left_join(select(DfACP, name, Dim.1, Dim.2, Dim.3))

egocity.emnw <- egocity.emnw %>% left_join(select(DfACP, name, Dim.1, Dim.2, Dim.3))

## CAH on dim ACP 
library(cluster)
library(ggdendro)
library(reshape2)

# compute classification ----

suppressWarnings(source("../multi_bivariateFunctions_exploratR/Functions_anova_cah.R"))


myVar <- c("Dim.1","Dim.2", "Dim.3")
cah <- ComputeClassif(df = egocity.emnw,
                      varquanti = myVar, method = "ward", stand = FALSE)


dendro <- PlotDendro(classifobj = cah)
dendro
inert <- PlotHeight(classifobj = cah)
inert
myProfiles <- PlotProfile(classifobj = cah, nbclus = 4)
myProfiles$PROFILE
table(myProfiles$CLUSID)
#egocity.em2nw <- egocity.em2nw %>% mutate(CAH = myProfiles$CLUSID)
egocity.emnw <- egocity.emnw %>% mutate(CAH = myProfiles$CLUSID)

## plot profile of CAH with original variables
PlotCAHProfile <- egocity.emnw %>% select(c(4:10)| CAH) %>%  mutate(reinf = replace(reinf, is.nan(reinf), 0)) %>% pivot_longer(-CAH)
library(questionr)
irec(PlotCAHProfile)
## Recodage de PlotCAHProfile$name en PlotCAHProfile$name_rec
PlotCAHProfile$name_rec <- fct_recode(PlotCAHProfile$name,
                                      "Nb Projets" = "nbProjects",
                                      "Nb Villes" = "nbVilles",
                                      "Nb Liens" = "nbliens",
                                      "Densité" = "densite",
                                      "Nb Composantes\nConnexes" = "nbcomp",
                                      "Taille de la\nplus grande\ncomposante (%)" = "sizecomp",
                                      "Reinforcement" = "reinf")
PlotCAHProfile$name_rec <- as.character(PlotCAHProfile$name_rec)
ggplot(PlotCAHProfile)+ 
  geom_boxplot(aes(x =CAH, fill = CAH, y = value), axes = TRUE, outlier.colour="black", outlier.size=0.5, outlier.alpha = 0.6)+ 
  facet_wrap(~ name_rec, scales = "free_y")+ 
  labs(y = "", x = "", fill = "Typologies\ndes ego-networks\ndes villes URBACT\n(CAH)", 
       caption = "NB : matrice emnw\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 12))

ggsave( filename = "OUT/URBACTEGONetwork_CAH4Profiles_BoxPlot_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )
# saveRDS(cah , "DataProd/cah_ego_eucicop_em2nw.rds")
# saveRDS(cah , "DataProd/cah_ego_eucicop_emnw.rds")

## plot profile of CAH with centralities indexes

# Add centralities mesures

CentralCities_emnw <- Centralities[[1]] %>% 
  select(-ends_with("2"))%>% 
  as.data.frame()%>% select(1:8) %>% rename(name = geonameId)

egocity.emnw <- egocity.emnw %>% left_join(CentralCities_emnw)

PlotCAHProfile <- egocity.emnw %>% select(CAH, ND,NB,NC,N_Eigen_C) %>% pivot_longer(-CAH)



ggplot(PlotCAHProfile)+ 
  geom_boxplot(aes(x =CAH, fill = CAH, y = value), axes = TRUE, outlier.colour="black", outlier.size=0.5, outlier.alpha = 0.6)+ 
  facet_wrap(~ name, scales = "free_y")+ 
  labs(y = "", x = "", fill = "Typologies\ndes ego-networks\ndes villes URBACT\n(CAH)", 
       caption = "Note : ces variables de centralités\nne sont pas prises en compte\ndans la CAH.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 12))

ggsave( filename = "OUT/URBACTEGONetwork_CAH4Profiles_BoxPlotCentralities_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )
### Explore clustering outcomes in relation to population and administrative status of cities



### Test anova Classe and Pop

# correct population variable
egocity.emnw <- egocity.emnw %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>% filter(!is.na(PopAdmin11))%>%
  mutate(CAH = as.character(CAH))

# Perform AnoVa Pop/typo

anovaTab <- AnovaTab(df = egocity.emnw %>% filter(!label == "London") , varx = c("CAH"), 
                     vary = c("PopAdmin11")) 
resume <- ComputeRegression(egocity.emnw, vardep = "PopAdmin11", varindep = "CAH")
resume$TABCOEF ## R2 = 0.19
rm(resume, anovaTab)


require(RColorBrewer)
require(ggrepel)
anovaPlot <- AnovaPlot(df = egocity.emnw %>% filter(!label == "London") , varx = c("CAH"), 
                       vary = c("PopAdmin11"),
                       tx = "Classe Ego-Network",
                       ty = "Population Administrative 2011",
                       source = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

### save plot
pdf(file = "OUT/ANOVAboxplot_EgoNet_Pop11__URBACT.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
anovaPlot
dev.off()
# ---- CHI2 ----

tableCAH <- table(egocity.emnw$CAH, egocity.emnw$adminLevel)



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

ggsave(filename = "OUT/HeatMap_PearsonChi2_CAHegoNetandAdmin_URBACT_emnw.pdf", device = "pdf",width = 8.3, height = 5.8 , units = "in")
class(resChi2)
testChi2$stdres

library(vcd)

cramer <-assocstats(tableCAH)
cramer  # 0.26


### Test with lead partner status

LeadPartnerCount <- ParticipationURBACT %>% 
                    group_by(geonameId, City.Statut)%>% 
                    summarise(nbLead= n()) %>% 
                    filter(City.Statut == "Lead Partner") %>% 
                    select(-City.Statut) %>% 
                    rename(name = geonameId)

egocity.emnw <- egocity.emnw %>% left_join(LeadPartnerCount)

# Perform AnoVa Lead partners count/typo


anovaTab <- AnovaTab(df = egocity.emnw %>% filter(!is.na(nbLead)) , varx = c("CAH"), 
                     vary = c("nbLead")) 
resume <- ComputeRegression(egocity.emnw %>% filter(!is.na(nbLead)), vardep = "nbLead", varindep = "CAH")
resume$TABCOEF ## R2 = 0.22
rm(resume, anovaTab)


require(RColorBrewer)
require(ggrepel)
anovaPlot <- AnovaPlot(df =egocity.emnw %>% filter(!is.na(nbLead)) , varx = c("CAH"), 
                       vary = c("nbLead"),
                       tx = "Classe Ego-Network",
                       ty = "Nombre de participation en tant que Lead Partner",
                       source = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

### save plot
pdf(file = "OUT/ANOVAboxplot_EgoNet_nbLead__URBACT.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
anovaPlot
dev.off()


## Chi 2 lead partner / typo ego

summary(egocity.emnw %>% filter(!is.na(nbLead)) %>% select(nbLead))

# discretize
egocity.emnw <- egocity.emnw %>% mutate(LeadPartClass = case_when(is.na(nbLead)~ "Jamais", 
                                                                  nbLead == 1 ~ "Une fois",
                                                                  nbLead == 2 ~ "Deux fois",
                                                                  nbLead > 2  ~ "Multiple Lead Partner"))


tableCAH <- table(egocity.emnw$CAH, egocity.emnw$LeadPartClass)



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
       y = "Nombre de participations en tant que lead partner",
       fill = "Résidus de Pearson")

ggsave(filename = "OUT/HeatMap_PearsonChi2_CAHegoNetandNbLead_URBACT_emnw.pdf", device = "pdf",width = 8.3, height = 5.8 , units = "in")
class(resChi2)
testChi2$stdres

library(vcd)

cramer <-assocstats(tableCAH)
cramer  # 0.30


## Small cities of type 1 and 4

SmallClass1 <- egocity.emnw %>% filter(CAH == "CLASSE 1" & PopAdmin11 < 100000) %>% 
  select(Ville = label, Pays = Country, Région = subregion,  "Niveau Administratif" = adminLevel, "Population Administrative 2011" = PopAdmin11,
         "Nombre de projets" = nbProjects,"Nombre de villes partenaires" = nbVilles, "Densité egonetwork" = densite,"Reinforcement" = reinf, "Nombre de composantes" = nbcomp,
          "Participation en tant que lead partner" = LeadPartClass)


SmallClass4 <- egocity.emnw %>% filter(CAH == "CLASSE 4" & PopAdmin11 < 100000) %>% 
  select(Ville = label, Pays = Country, Région = subregion,  "Niveau Administratif" = adminLevel, "Population Administrative 2011" = PopAdmin11,
         "Nombre de projets" = nbProjects,"Nombre de villes partenaires" = nbVilles, "Densité egonetwork" = densite, "Reinforcement" = reinf, "Nombre de composantes" = nbcomp,
         "Participation en tant que lead partner" = LeadPartClass)



# Export for latex
library(kableExtra)

kable(SmallClass1, booktabs = T, "latex") %>% 
  column_spec(column = c(3:11), width = "2cm") %>% 
  kable_styling(position = "center", latex_options = "striped")





#### ==== MAPPING SOME EGO NETWORKS  =====

#just create  a variable name specific to project such as Acro in ETMUN

NodesInfos <- NodesInfos %>% mutate(Acro = ifelse(is.na(Date), NA, label))
# create a sample of egonetworks

sampleEgo <- egocity.emnw  %>% group_by(CAH) %>% sample_n(1, replace = T)

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
      
      # Layout
      l <- layout_with_kk(tdg)
      l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
      l<- as.data.frame(l) %>% rename(x= V1, y = V2)
      manual_layout <- create_layout(tdg, layout = "manual", node.positions = l)
      # plot
      
      
      plot <- ggraph(manual_layout)+
        geom_edge_link(alpha = 0.2)+
        geom_node_point(aes(color = V(tdg)$type2, shape = V(tdg)$type, size = V(tdg)$size),
                        show.legend = TRUE, alpha = 0.8) + 
        geom_node_label(aes(label = Acro), repel = TRUE, size = 2, alpha = 0.6)+
        scale_color_manual( values =c("lightskyblue", "mediumvioletred", "orange"), breaks = c("City", "ego", "Asso"), labels = c("Villes", "Ego", "Projets"))+
        scale_size(range = c(1.5,6))+
        labs(title = paste("Ego Network de ", strInfo$label,", ", strInfo$Country, sep = ""),
             subtitle = paste(strInfo$CAH), color = "Type de sommets",
             caption = "Algorithme de spatialisation : kk")+
        annotate("label",
                 label = paste("Densité : ", round(strInfo$densite,3),
                               "\nNb Villes : ", strInfo$nbVilles,
                               "\nNb Projets : ", strInfo$nbProjects,
                               "\nNb composantes : ", strInfo$nbcomp, 
                               "\nReinforcement : ", round(strInfo$reinf,4), 
                               sep = ""), 
                 x = 1, 
                 y = 1,
                 hjust = 1,
                 vjust = 1,
                 alpha =0.4,
                 size = 2.5)+
        guides(size=FALSE, shape = FALSE) +
        theme_graph(base_family = 'Helvetica');
      
      
      Node <- strInfo$name;
      
      
      listplot[[Node]] <<- plot
      #print(plot)
    }) 
  }
  
  return(listplot)
  
}

#Apply function egoMapping

egoplot.emnw <- EgoMapping(df = sampleEgo,egolist = ego.emnw,NodesInfos =NodesInfos)
# Grid plot function

Gridplot <- function(myplots, n){
  require(cowplot)
  splitted_plots <- split(myplots, ceiling(seq_along(myplots)/n))
  
  lapply(splitted_plots, function(x) plot_grid(plotlist = x))
  
}


#Try
gEgo <- Gridplot(egoplot.emnw, 4)

gEgo[[1]]
ggsave( filename = "OUT/URBACTEGONetwork_cluster5.pdf", width = 11.7, height = 8.3, units = "in" )



gEgo <- Gridplot(egoplot.emnw, 4)

gEgo[[1]]
ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5_emnw.pdf", width = 11.7, height = 8.3, units = "in" )

gEgo[[2]]

ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5_UniqueClasse5_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )


#### ==== Geo Mapping ego network typology ====
library(sf)
## sf 

sf_egocity.emnw <- egocity.emnw %>% dplyr::left_join(select(ParticipationURBACT, name = geonameId, lng_GN,lat_GN))%>% distinct()

sf_egocity.emnw <- st_as_sf(sf_egocity.emnw, coords = c("lng_GN", "lat_GN"), crs = 4326)%>% st_transform(crs = 3035)

# filter class 3 and london
sf_egocity.emnw <- sf_egocity.emnw %>% filter(!CAH == "CLASSE 3"& !label == "London") 

# Deal with  mestre == Venice
sf_egocity.emnw <- sf_egocity.emnw %>% mutate(label = recode(label, "Mestre" = "Venice"))
## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)

myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(sf_egocity.emnw$PopAdmin11)
s
bks <- c(5000, 50000, 200000,  500000,1000000, 9000000 )
lbs <-  c("5000",  "50 000",   "200 000", "500 000", "1 M", "9 M")



citiesEgoTypo <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.2) +
  geom_sf(data = sf_egocity.emnw ,
          mapping = aes(size = PopAdmin11, color = CAH), show.legend = "point", alpha = 0.8) +
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.5, 6))+
  scale_color_manual(values = c("#F8766D" ,"#7CAE00","#C77CFF"),labels = c("Classe 1 :  Très haute participation, faible densité,\nnombre élevé de villes partenaires", 
                                  "Classe 2 : Participation moyenne,\nensemble de villes partenaires très diversifié", 
                                  "Classe 4 :  Participation faible, faible densité,\ngroupes de travail récurrents"))+
  labs(x = "", y = "", caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020", color = "Typologie des ego-networks") +
  geom_sf_text(data = sf_egocity.emnw %>% filter(PopAdmin11 > 100000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0.2, vjust = -0.5, 
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + 
  guides(size=guide_legend(ncol=2) )+
  theme_void() +  
  facet_wrap(~CAH,  nrow = 2,ncol = 2)+
  theme(legend.position =  c(1,0), legend.justification = c(1.2,-0.2),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5),
        strip.text = element_blank(), plot.caption = element_text(size = 7)) 

pdf(file = "OUT/EgoNetTypo_Cities_emnw_URBACT.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesEgoTypo
dev.off()

### ==== TEMPORAL DISTRIBUTION OF DEGREES and Centrality Measures  ====

# Function centalities indexes

CentralitydfTemp <- function(matrix1, matrix2, matrix3){
  require(igraph)
  
  # Matrix 1
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
  
  # Matrix 2
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
  
  # Matrix 3
  matrix3.g <- graph.incidence(matrix3, directed = FALSE)
  city3 <- data.frame(D3 = colSums(matrix3), 
                      ND3 = colSums(matrix3)/nrow(matrix3), 
                      geonameId = colnames(matrix3),
                      stringsAsFactors = FALSE)
  asso3 <- data.frame(D3 = rowSums(matrix3), 
                      ND3 = rowSums(matrix3)/ncol(matrix3), 
                      Code = rownames(matrix3), 
                      stringsAsFactors = FALSE)
  
  
  V(matrix3.g)$B3 <- betweenness(matrix3.g, normalized = FALSE, directed = FALSE)
  V(matrix3.g)$NB3 <- betweenness(graph = matrix3.g, normalized = TRUE, directed = FALSE)
  V(matrix3.g)$C3 <- closeness(graph = matrix3.g, normalized = FALSE)
  V(matrix3.g)$NC3 <- closeness(graph = matrix3.g, normalized = TRUE)
  V(matrix3.g)$N_Eigen_C3 <- eigen_centrality(graph = matrix3.g, scale = TRUE, directed = FALSE)$vector
  
  df3 <- as.data.frame(vertex_attr(matrix3.g), stringsAsFactors = FALSE)
  cityBC3 <- df3[df3$type== TRUE,]
  assoBC3 <- df3[df3$type== FALSE,]
  
  city3 <- city3 %>% left_join(cityBC3, by= c("geonameId" = "name")) %>% select(-type)
  asso3 <- asso3 %>% left_join(assoBC3, by= c("Code" = "name")) %>% select(-type)
  
  city <- merge(city, city2, by = "geonameId", all = T)
  asso <- merge(asso, asso2, by = "Code", all = T)
  
  city <- merge(city, city3, by = "geonameId", all = T)
  asso <- merge(asso, asso3, by = "Code", all = T)
  
  result <- list(city, asso)
  return(result)
}





Centralities <- CentralitydfTemp(matrix1 = emnw1, matrix2 =  emnw2, matrix3 =  emnw3)
bibi <- Centralities[[1]]
saveRDS(Centralities, "DataProd/Centralities_URBACT_phases.rds")
### plot of degree distributions

# hist
Centralities[[1]] %>% keep(is.numeric) %>% select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 3) +
  geom_histogram() + 
  scale_y_log10() + 
  # scale_x_log10()+
  labs(title = "Distributions des mesures de centralité des villes URBACT",
       subtitle = "Selon les 3 phases du programme", 
       y = "Nombre de villes (log10)", x = "",
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 1,2,3 = phase\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")


summary(Centralities[[1]]$ND)
summary(Centralities[[1]]$D)
ggsave(filename = "OUT/CentralDistrib_3phases_Cities_URBACT.pdf",width = 8.3, height = 5.8 )

Centralities[[2]] %>% keep(is.numeric) %>% select(starts_with("N")) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 3) +
  geom_histogram()+ 
  # scale_x_log10(n.breaks = 6) +
  labs(title = "Distributions des mesures de centralité des projets URBACT",
       subtitle = "Selon les 3 phases du programme", 
       y = "Nombre de projets", x = "Valeurs normalisées",
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation,  1,2,3 = phase\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

summary(Centralities[[2]]$ND)
summary(Centralities[[2]]$D)
ggsave(filename = "OUT/CentralDistrib_3phases_Project_URBACT.pdf",width = 8.3, height = 5.8 )

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & -ends_with("3") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes URBACT (emnw1)")

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes URBACT (emnw2)")

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(ends_with("3") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes URBACT (emnw3)")

Centralities[[2]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & -ends_with("3") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets URBACT (emnw1)")

Centralities[[2]] %>% 
  keep(is.numeric) %>%
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets URBACT (emnw2)")

Centralities[[2]] %>% 
  keep(is.numeric) %>%
  select(ends_with("3") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets URBACT (emnw3)")

#Set city name

Centralities[[1]] <- Centralities[[1]] %>% left_join(namecity)

bibi = Centralities[[1]]
# set names of association
Centralities[[2]] <- Centralities[[2]] %>% left_join(select(ProjectURBACT, Code = Code_Network,  ProjectName = Name))

### ==== TOP CENTRALITIES ====

### Rank 
CityCentralities <- Centralities[[1]] %>% arrange(desc(D)) %>%  select( geonameId | asciiName | starts_with("N"))

RankCities <- CityCentralities %>% 
  mutate_if(is.numeric, 
            list( rank = ~ rank(desc(.), ties.method = 'min', na.last = "keep"))) 

# cities with a rank < 5 in any index
TopRank <- RankCities %>% filter_at(.vars = c(15:26), any_vars(.<5))

# tibble for ggplot rank
TopRankgg <- TopRank %>% select(geonameId | asciiName | ends_with("rank"))%>% 
  pivot_longer(-c(geonameId,asciiName), names_to = "indice", values_to = "rang")%>% 
  mutate(indice = str_remove(indice,pattern = "_rank")) %>% 
  mutate(indice = str_remove(indice,pattern = "N_")) %>%
  mutate(phase = case_when(str_detect(indice, pattern = "2") ~ 2,
                           str_detect(indice, pattern = "3") ~ 3,
                           TRUE ~ 1))%>%
  mutate(indice = str_remove_all(indice,"[23]")) %>%
  mutate(lograng = log10(rang))

#plot in coord polar
s <- summary(TopRankgg$lograng)
bks <- c(s[[1]], log10(2), s[[2]], s[[3]], s[[5]] , s[[6]])
lbs <- round(10^bks,0)


toprankggphase <-TopRankgg %>% 
  ggplot(aes(x = phase, y = lograng, label = rang) ) + 
  geom_line(aes(col = indice)) +
  geom_point(aes(col = indice))+
  facet_wrap(~asciiName) + 
  scale_y_continuous(trans = "reverse", 
                     breaks = bks, 
                     labels = lbs)+
  scale_x_continuous(breaks = c(1,2,3), labels = c("2000-2006", "2007-2013", "2014-2020"), guide = guide_axis(n.dodge = 3))+
  labs(x = "Phases du Programme URBACT", y = "Rang (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020") + 
  theme(axis.text.y=element_text(size=7), legend.position="bottom", axis.text.x=element_text(size=8))

pdf(filename = "OUT/Rank_TopAny4_URBACT_Phase.pdf",width =11.7, height = 8.3)
toprankggphase
dev.off()

#### ==== TEMPORAL ARTICULATION POINTS AND IGRAPH OBJECT ====

emnw1.g <- graph.incidence(emnw1, directed = F)

emnw2.g <-  graph.incidence(emnw2, directed = F)

emnw3.g <-  graph.incidence(emnw3, directed = F)


## Articulation points
#emnw1
art <-articulation_points(emnw1.g)
vecArt <- V(emnw1.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)

skim(art)## 34 articulations points (all projects)
#emnw2

art <-articulation_points(emnw2.g)
vecArt <- V(emnw2.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)
skim(art) ## 56 (4 are cities : Soedertaelje, SE ; Lille, FR ; Albacete, ES ; Kavala, GR)

bc <- biconnected_components(emnw2.g)
#emnw3

art <-articulation_points(emnw3.g)
vecArt <- V(emnw3.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)
skim(art) ## 72 (One is a city : 	Loule, PT)
