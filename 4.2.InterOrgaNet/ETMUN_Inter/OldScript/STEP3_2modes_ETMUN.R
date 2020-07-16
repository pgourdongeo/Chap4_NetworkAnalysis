############################### ETMUN bipartite Networks  ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe ETMUN (affiliation) biparti. 
#               Indices locaux,distributions des degrés, détection de communautés
#
# 
############################################################################## PG juin 2020
### BECOME STEP 2

### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")

#Packages
library(tidyverse)
library(tidylog)
library(igraph)
library(gridExtra)
library(Matrix)
library(patchwork)
library(gghighlight)
library(GGally)

options(scipen = 999)
#####Data 

## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_GNidCorr.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", 
                       stringsAsFactors = F)


# Matrix from STEP1

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/Etmun_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/Etmun_em2nw.rds")


### ==== DISTRIBUTION OF DEGREES and Centrality Measures ====

# Function degrees

Centralitydf <- function(matrix1, matrix2){
  require(igraph)
  
  matrix1.g <- graph.incidence(matrix1, directed = FALSE)
  city <- data.frame(D = colSums(matrix1), 
                     ND = colSums(matrix1)/nrow(matrix1), 
                     geonameId = colnames(matrix1))
  asso <- data.frame(D = rowSums(matrix1), 
                     ND = rowSums(matrix1)/ncol(matrix1), 
                     Code = rownames(matrix1))
  
  V(matrix1.g)$B <- betweenness(matrix1.g, normalized = FALSE, directed = FALSE)
  V(matrix1.g)$NB <- betweenness(graph = matrix1.g, normalized = TRUE, directed = FALSE)
  V(matrix1.g)$C <- closeness(graph = matrix1.g, normalized = FALSE)
  V(matrix1.g)$NC <- closeness(graph = matrix1.g, normalized = TRUE)
  
  df1 <- as.data.frame(vertex_attr(matrix1.g))
  cityBC <- df1[df1$type== TRUE,]
  assoBC <- df1[df1$type== FALSE,]
  
  city <- city %>% left_join(cityBC, by= c("geonameId" = "name")) %>% select(-type)
  asso <- asso %>% left_join(assoBC, by= c("Code" = "name")) %>% select(-type)
  
  
  matrix2.g <- graph.incidence(matrix2, directed = FALSE)
  city2 <- data.frame(D2 = colSums(matrix2), ND2 = colSums(matrix2)/nrow(matrix2), geonameId = colnames(matrix2))
  asso2 <- data.frame(D2 = rowSums(matrix2), ND2 = rowSums(matrix2)/ncol(matrix2), Code = rownames(matrix2))
  
  
  V(matrix2.g)$B2 <- betweenness(matrix2.g, normalized = FALSE, directed = FALSE)
  V(matrix2.g)$NB2 <- betweenness(graph = matrix2.g, normalized = TRUE, directed = FALSE)
  V(matrix2.g)$C2 <- closeness(graph = matrix2.g, normalized = FALSE)
  V(matrix2.g)$NC2 <- closeness(graph = matrix2.g, normalized = TRUE)
  
  df2 <- as.data.frame(vertex_attr(matrix2.g))
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

### plot of degree distributions

# hist

Centralities[[1]] %>% keep(is.numeric) %>% select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des mesures de centralité des villes ETMUN",
       subtitle = "Sur les matrices emnw et em2nw", 
       y = "Nombre de villes (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nw")

Centralities[[2]] %>% keep(is.numeric) %>% select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram()+ 
  scale_y_log10() + 
  scale_x_log10() +
  labs(title = "Distributions des mesures de centralité des associations ETMUN",
       subtitle = "Sur les matrices emnw et em2nw", 
       y = "Nombre d'associations (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nw")

Centralities[[1]] %>% keep(is.numeric) %>% select(-ends_with("2") & starts_with("N")) %>% ggpairs(title = "Relations entre les mesures de centralité des villes ETMUN (emnw)")

Centralities[[1]] %>% keep(is.numeric) %>% select(ends_with("2") & starts_with("N")) %>% ggpairs(title = "Relations entre les mesures de centralité des villes ETMUN (em2nw)")

Centralities[[2]] %>% keep(is.numeric) %>% select(-ends_with("2") & starts_with("N")) %>% ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (emnw)")

Centralities[[2]] %>% keep(is.numeric) %>% select(ends_with("2") & starts_with("N")) %>% ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (em2nw)")


#Set city name
namecity <- MembershipEtmun %>% select(geonameId, asciiName) %>% distinct()
Centralities[[1]] <- Centralities[[1]] %>% left_join(namecity)

# set names of association
Centralities[[2]] <- Centralities[[2]] %>% left_join(select(AssoEtmun, Code,  AssoName = Name))


## Top by centralities

Top50Cities <- Centralities[[1]] %>% arrange_if(is.numeric, desc)%>% head(50)
Top10Asso <- Centralities[[2]] %>% arrange_if(is.numeric, desc)%>% head(10)



## cumulative normalized degree
ytext <- "Cumulative Frequency (log10)"
xtext <- "Normalized Degree (log10)"

g1 <-ggplot(Centralities[[1]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Normalized Degrees' Distribution of Cities (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g2 <- ggplot(Centralities[[1]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Normalized Degrees' Distribution of Cities (em2nw)", y = ytext,  x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g3 <- ggplot(Centralities[[2]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Normalized Degrees' Distribution of Associations (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g4 <- ggplot(Centralities[[2]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Normalized Degrees' Distribution of Associations (em2nw)", y = ytext, x = xtext,
       caption = "Sources : ETMUN 2019\nPG 2020") + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()


grid1 <- (g1 | g2) / (g3 | g4)

##########
emnw.g <- graph.incidence(emnw, directed = F)

em2nw.g <-  graph.incidence(em2nw, directed = F)

articulation_points(emnw.g)
articulation_points(em2nw.g)

#### ==== FILTER OUT SOME ASSOCIATIONS ====
## show outliers by degrees (2 associations)


#filter some asso outliers (visual detection) according to the ggpairs on Associations

Centralities[[2]]%>% filter(ND>0.10 | NC<0.2)# Covenant of mayors, Climate Alliance, WWCAM

Centralities[[2]] %>% keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% filter(ND<0.10 & NC>0.2)%>%
  ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (em2nw)")


# Covenant of mayors, Climate Alliance, WWCAM
assoremove <- c("04080", "18332", "03241")
#assoremove <- c("Covenant of Mayors", "WWCAM","Climate Alliance")

#Filter out the association on the matrix that keep city of degree 1
emnwF <- emnw[!rownames(emnw) %in% assoremove, ]
emnwF <- emnwF[,colSums(emnwF) > 0]

#Filter out the association on the matrix that keep city of more than 1 degree (only cities involved at least in 2 associations)
em2nwF <-  em2nw[!rownames(em2nw) %in% assoremove, ]
em2nwF <- em2nwF[,colSums(em2nwF) > 1]

## Graph

emnwF.g <- graph.incidence(emnwF, directed = F)

em2nwF.g <- graph.incidence(em2nwF, directed = F)

articulation_points(emnwF.g)
articulation.points(em2nwF.g)#none
### ==== FILTERED MATRIX DISTRIBUTION OF DEGREES and Centrality Measures ====

#Centralities Filtered

CentralityFiltered <- Centralitydf(emnwF, em2nwF)

## Histogram centralities

CentralityFiltered[[1]] %>% keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des degrés des villes ETMUN",
       subtitle = "Sur les matrices emnwF et em2nwF", 
       y = "Nombre de villes (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nwF")


CentralityFiltered[[2]] %>% keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram() + 
  scale_y_log10() + 
  labs(title = "Distributions des degrés des Associations ETMUN",
       subtitle = "Sur les matrices emnwF et em2nwF", 
       y = "Nombre de villes (log10)", 
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nwF")



#### ==== COMMUNITY DETECTION ON FILTERED MATRIX =====


# Function

DualProjectionCommunity<-function(matrix){
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
  Result<-list(Y1,Y2)
  Result<-t(as.matrix(Result))
  names<-c('DP_Walk_First','DP_Walk_Second')
  colnames(Result)<-names
  return(Result)
}

clustering <- DualProjectionCommunity(emnwF)

clustering2 <- DualProjectionCommunity(em2nwF)

ClusterAsso <- data.frame(Code = rownames(emnwF), clusterAsso = as.factor(clustering[[1]])) %>% 
  left_join(select(AssoEtmun,Acronym,Name, Code ))

ClusterCities <- data.frame(geonameId = colnames(emnwF), ClusterCities = as.factor(clustering[[2]])) %>% 
  left_join(namecity)
 
clustering[[2]]

ClusterAsso2 <- data.frame(Code = rownames(em2nwF), clusterAsso = as.factor(clustering2[[1]])) %>% 
  left_join(select(AssoEtmun,Acronym,Name, Code )) 

ClusterCities2 <- data.frame(geonameId = colnames(em2nwF), ClusterCities = as.factor(clustering2[[2]])) %>% 
  left_join(namecity)

table(ClusterCities2$ClusterCities)


### try to plot the network 
# get a df with nodes info

emnwFCentralities <- CentralityFiltered[[1]] %>% select(- ends_with("2")) %>% rename("name"= geonameId)
emnwFCentralitiesAsso <- CentralityFiltered[[2]] %>% select(- ends_with("2")) %>% rename("name" = Code)

emnwFCentralities <- bind_rows(emnwFCentralities, emnwFCentralitiesAsso)

ClusterAssoBind <- ClusterAsso %>% mutate(clusterAsso = paste0("A", clusterAsso)) %>% rename("clusterDP" = clusterAsso,"name" = Code, "label" = Name )

ClusterCitiesBind <- ClusterCities %>% mutate(ClusterCities = paste0("C", ClusterCities)) %>% rename("clusterDP" = ClusterCities,"name" = geonameId, "label" = asciiName )

CLusterEmnwF <- bind_rows(ClusterCitiesBind, ClusterAssoBind)

emnwF_NodesInfo <- emnwFCentralities %>% left_join(CLusterEmnwF)

#set all variables of the df as vertex attributes....

# Set full name
emnwF.g  <-  set_vertex_attr(emnwF.g ,"Degree", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$D)## same as a joint

# Set Short Name
emnwF.g  <-  set_vertex_attr(emnwF.g ,"ND", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$ND)

# Set Year of creation
emnwF.g   <- set_vertex_attr(emnwF.g ,"B", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$B)

# Set country of the seat
emnwF.g   <- set_vertex_attr(emnwF.g ,"NB", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$NB)

emnwF.g  <-  set_vertex_attr(emnwF.g ,"C", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$C)

emnwF.g  <-  set_vertex_attr(emnwF.g ,"NC", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$NC)

emnwF.g  <-  set_vertex_attr(emnwF.g ,"ClustDP", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$clusterDP)

emnwF.g  <-  set_vertex_attr(emnwF.g ,"Label", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$label)

emnwF.g  <-  set_vertex_attr(emnwF.g ,"Acro", index = emnwF_NodesInfo$name, value = emnwF_NodesInfo$Acronym)

unique(emnwF_NodesInfo$clusterDP)

V(emnwF.g)$shape <- c("square", "circle")[V(emnwF.g)$type+1]
V(emnwF.g)$size <- c(5, 1)[V(emnwF.g)$type+1]
library(RColorBrewer)
palcol <- c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
                    "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
                    "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
                    "springgreen2", "yellowgreen", "palegreen4",
                    "wheat2", "tan", "tan2", "tan3", "brown",
                    "grey70")# "grey50", "grey30


 palcol[V(emnwF.g)$ClustDP]
 l1 <- layout_with_graphopt(emnwF.g, charge=0.02)
 l1 <- layout_with_kk(emnwF.g, kkconst = vcount(emnwF.g)*0.004)
 
 
 l <- norm_coords(l1, ymin=-1, ymax=1, xmin=-1, xmax=1)
 par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(emnwF.g, vertex.color=palcol[as.numeric(as.factor(vertex_attr(emnwF.g, "ClustDP")))], vertex.label = V(emnwF.g)$Acro, layout = l) 
