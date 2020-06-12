############################### ETMUN Inter Organisations Network  ############
#                               
#                          
# DESCRIPTION : Création des données relationnelles à partir de la table
#               des adhésions ETMUN. Projection 1mode du biparti
#               Graphes asso = nombres de villes membres en commun
#
#
# 
############################################################################## PG mai 2020


### ==== Load Package and Date ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")
#Packages
library(tidyverse)
library(tidylog)
library(skimr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(flows)
library(GGally)
library(patchwork)
library(gghighlight)

#Data 
## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_GNidCorr.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", stringsAsFactors = F)


### ==== Create table Asso-Asso by shared members ====

# simple graph


edgelist <- MembershipEtmun %>% filter(!is.na(geonameId)) %>% select( geonameId, Code_Network)
edgelist <- MembershipEtmun %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network) %>% 
  group_by(geonameId,Code_Network)%>% summarise(weight = n())
edgelistnw <- edgelist %>% mutate(weight = 1) 
#create network from edge list

emnw <- edgelistnw  %>% 
  pivot_wider(names_from = geonameId, values_from = weight, values_fill = list(weight = 0)) %>% 
  column_to_rownames(var="Code_Network") %>% 
  as.matrix()
em2nw <- emnw[,colSums(emnw) > 1]
Network2modes <- graph.data.frame(edgelist)
Network2modes <- graph.incidence(em2nw , directed = F)
#check

is_bipartite(Network2modes)#FALSE

# Bipartite network

V(Network2modes)$type <- bipartite_mapping(Network2modes)$type 

#check

is_bipartite(Network2modes)#TRUE

## bi partite graph

NetProj <- bipartite.projection(Network2modes)

# select the inter-asso projection 

InterCities <- NetProj$proj2
InterOrga <- NetProj$proj1

## Add some info on nodes (just the ones we want)

# Set full name
InterOrga  <-  set_vertex_attr(InterOrga,"Label", index = AssoEtmun$Code, value = AssoEtmun$Name)## same as a joint

# Set Short Name
InterOrga <-  set_vertex_attr(InterOrga,"Acro", index = MembershipEtmun$Code_Network, value = MembershipEtmun$Network_Name)

# Set Year of creation
InterOrga  <- set_vertex_attr(InterOrga,"Year", index = AssoEtmun$Code, value = AssoEtmun$Date)

# Set country of the seat
InterOrga  <- set_vertex_attr(InterOrga,"CountrySeat", index = AssoEtmun$Code, value = AssoEtmun$Country..secretariat.)

# Set AssoSize of asso 

SizeAsso <- edgelistnw %>% group_by(Code_Network)%>% summarise(nMembers = n())
InterOrga  <- set_vertex_attr(InterOrga,"AssoSize", index = SizeAsso$Code_Network, value = SizeAsso$nMembers)
# compute degree

V(InterOrga)$degree <- degree(InterOrga)
V(InterOrga)$wdegree <- strength(InterOrga, vids = V(InterOrga), loops = F )


# transform into tidy graph and df

TdInterOrga <- as_tbl_graph(InterOrga)
DfInterOrga <- fortify.tbl_graph(TdInterOrga)

# #compute relative wdegree 


CompRelative_wdegree <- function(df, size , wdegree){
  df <- as.data.frame(df)
  df <- df[order(df[size]),]
  df$cumsumSize <- unlist(cumsum(df[size]))
  
  V <- nrow(df)
  secn <- seq(1:V)
  df$max_thwdegree <- NA
  for (i in secn){
    tmpCumSize <- df[i-1, "cumsumSize", drop = FALSE]
    tmpCumSize <- ifelse(nrow(tmpCumSize) == 0, 0, tmpCumSize)
    tmpSize <- df[i,size] *(V-i)
    df[i, "max_thwdegree"] <- sum(tmpCumSize[[1]], tmpSize, na.rm = TRUE)
  }
  
  df$norm_wdegree <- unlist(df[wdegree]/df$max_thwdegree)
  
  Wdensity <- paste0("Density of weigted graph (proj 1-mode of a 2-modes graph) : " ,
                     round(sum(df[wdegree])/sum(df$max_thwdegree),4))
  print(Wdensity)

  
  df <- df %>% select(-cumsumSize)
 
  return(df)
  
}


DfInterOrga2  <- CompRelative_wdegree(df = DfInterOrga,
                            size = "AssoSize",
                            wdegree = "wdegree")#"Density of weigted graph (proj one-mode of a 2modes graph) : 0.1603"

str(DfInterOrga2)
#bivariate overlook of new variable
DfInterOrga2 %>% keep(is.numeric) %>% filter(wdegree < 2000 & degree >2 & AssoSize < 1500)%>% select(-Year, - max_thwdegree) %>% ggpairs()

# set normalized wdegree on the original graph
InterOrga <- set_vertex_attr(InterOrga,"norm_wdegree", index = DfInterOrga2$name, value = unlist(DfInterOrga2$norm_wdegree))

# Update tidy and df
TdInterOrga <- as_tbl_graph(InterOrga)
DfInterOrga <- fortify.tbl_graph(TdInterOrga) %>% mutate(norm_wdegree = unlist(norm_wdegree) )


## Normalized egdes the same way
#divide the weight by the smallest AssoSize of the two vertices
# get edges
edgesInterOrga <- TdInterOrga %>% activate(edges) %>% fortify.tbl_graph()

#add vertices info
DfInterOrga <- DfInterOrga %>% mutate(id.edge = as.numeric(rownames(.)))
#from
edgesInterOrga <- edgesInterOrga %>% 
        left_join(select(DfInterOrga, id.edge, Labelfrom = Label, AssoSizefrom = AssoSize), 
                  by = c("from" = "id.edge"))
#to
edgesInterOrga <- edgesInterOrga %>% 
  left_join(select(DfInterOrga, id.edge, Labelto = Label, AssoSizeto = AssoSize), 
            by = c("to" = "id.edge"))

#get minimum assoSize

edgesInterOrga <- edgesInterOrga %>% rowwise()%>% mutate(MinAsso = min(AssoSizefrom, AssoSizeto))

# Normalized weight
edgesInterOrga <- edgesInterOrga %>% mutate(norm_weight = weight/ MinAsso)

# Set normalized weight in orignal graph
InterOrga <- set_edge_attr(InterOrga,"norm_weight", value = edgesInterOrga$norm_weight)

# Update tidy and df
TdInterOrga <- as_tbl_graph(InterOrga) %>% mutate(norm_wdegree = unlist(norm_wdegree) )
DfInterOrga <- fortify.tbl_graph(TdInterOrga)
### ==== Graph representation ====

## ggraph
## plot the graph unweighted
g1<-ggraph(TdInterOrga, layout = 'kk') + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point(color = "orange", aes(size = degree)) + 
  scale_size(range = c(0.01,5))+
  labs(title = "A", size = "Degré\n(nb d'associations\navec au moins\nune ville membre\nen commun)")+
  theme_graph(base_family = 'Helvetica')
 g1 

 
 ## plot the graph weighted
 # filter it for representation
 FilterG <- TdInterOrga %>% 
   activate(edges) %>% 
   filter(weight > 20 ) %>% 
   activate(nodes) %>% 
   filter(!node_is_isolated())
 
 # change weight name (problem because impossible to change width lab in ggraph)
 
 FilterG <-  FilterG %>% activate(edges) %>% rename( "Poids\n(nb villes membres\nen commun)" = weight)
 #layout
 l <- layout_with_fr( FilterG, niter = 10000, weights = sqrt(E(FilterG)$`Poids\n(nb villes membres\nen commun)`)/1000)

  g2 <-ggraph(FilterG, layout = "kk") + 
   geom_edge_link(aes(width = `Poids\n(nb villes membres\nen commun)`), alpha = 0.2)+
   geom_node_point( color = "indianred3", aes(size = wdegree)) + 
   geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
   scale_edge_width(range = c(0.1, 5)) +# control size
   scale_size(range = c(1,8))+
   labs(title = "B" ,size = "Degré pondéré\n(sommme du poids\ndes liens pour\nchaque sommet)",
   caption = "Note : Seuls les liens supérieurs\nà 20 sont représentés\n\nSources : ETMUN 2019. PG, 2020")+
   theme_graph(base_family = 'Helvetica')
 g2 
 


## Grid with two plots
 
 grid1<-g1 + g2 + plot_layout(ncol = 1, heights = c(1.5, 2))
#save
ggsave(grid1, filename = "OUT/GraphInterOrgaETMUN.pdf", width = 8.3, height = 8.3, units = "in" )

## plot the graph normalized
FilterG <- TdInterOrga %>% activate(nodes)%>% filter(wdegree < 2000 & degree >2 & AssoSize < 1500)
skim(E(FilterG)$norm_weight)
# discretize edges norm weight
FilterG <- FilterG %>% 
  activate(edges) %>% 
  mutate(strength_normw = case_when(norm_weight <= 0.05 ~ "1.Very Weak < 5%",
                                    norm_weight > 0.05 & norm_weight <= 0.10 ~ "2.Weak < 15%",
                                    norm_weight >0.10 & norm_weight <= 0.25 ~ "3.Strong < 30%",
                                    norm_weight >0.25 & norm_weight <= 0.50 ~ "4.Very Strong < 50%",
                                    norm_weight > 0.50  ~ "5.Exceptionnal > 50%"))


l <- layout_with_fr( FilterG, niter = 100, weights = sqrt(E(FilterG)$norm_weight))
g3<-ggraph(FilterG, layout = l) + 
  geom_edge_link(aes(width = norm_weight),alpha = 0.1)+
  geom_node_point( aes(color = norm_wdegree), size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  scale_color_gradient(low = "lightgoldenrod1", high = "indianred3" )+
  scale_edge_width(range = c(0.1, 5)) +
  labs(title = "All Edges")+
  theme_graph(base_family = 'Helvetica')
  
g3 

# facet the graph by class of normalized weight
l <- layout_with_fr( FilterG, niter = 1000, weights = sqrt(E(FilterG)$norm_weight))
g4<-ggraph(FilterG, layout = l) + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point( aes(color = norm_wdegree), size = 3 )+ 
  scale_color_gradient(low = "lightgoldenrod1", high = "indianred3" )+
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(size = "Degré Normalisé")+
  theme_graph(base_family = 'Helvetica')+
  facet_edges(~ strength_normw)
g4

# Create a facet grid manually (to recompute layout and remove isolated nodes each time)
l = "graphopt"
#1


fg1 <- g3

#2

f2 <- FilterG %>% 
  activate(edges)%>% 
  filter(norm_weight > 0.15)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg2 <- ggraph(f2, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Edges > 15%")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg2

#3

f3 <- FilterG %>% 
  activate(edges)%>% 
  filter(norm_weight > 0.30)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg3 <- ggraph(f3, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Edges > 30%")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg3

#4

f4 <- FilterG %>% 
  activate(edges)%>% 
  filter(norm_weight > 0.50)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg4 <- ggraph(f4, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Egdes > 50%")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg4

# 5
f5 <- FilterG %>% 
  activate(edges)%>% 
  filter(norm_weight > 0.70)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg5 <- ggraph(f5, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Edges > 70%",  
       caption = "Sources : ETMUN 2019. PG, 2020")+
  theme_graph(base_family = 'Helvetica',foreground = "grey60", border = TRUE)
fg5

# grid

grid2 <- fg3 / fg4 /fg5

#save
ggsave(grid2, filename = "OUT/EtmunInterOrga_Norm_grid.pdf", width = 8.3, height = 11.7, units = "in" )

### ==== Main numerical outputs ====
#Create df for classic ploting

Fgraph2 <- TdInterOrga %>% activate(nodes) %>% filter(AssoSize < 500) %>% filter(degree>2) # remove Covenant of Mayors (size= 10 000) and Climate Alliance (size =1742) + WWCAM (degree =1)

Fgraph2 <- fortify.tbl_graph(Fgraph2) 

## order and size
InterOrga
#Order = 59  , size 1230

## density

graph.density(InterOrga)# 0.72

## degree distribution

summary(V(InterOrga)$degree)

hist(V(InterOrga)$degree,breaks=10, xlab = "Degré des associations ETMUN", ylab = "Nombre d'associations", main = NA)

## Univariate overlook

Fgraph2 %>% keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

## bivariate overlook

Fgraph2 %>% keep(is.numeric) %>% select(-Year) %>% ggpairs()

#log
skim(Fgraph2)
Fgraph2Log10  <- Fgraph2  %>% 
  keep(is.numeric) %>%
  select(-Year) %>% 
  mutate_all( list(log10 = ~log10(.)))

ggpairs(Fgraph2Log10)
 
####### Relation between node variables ####


## wdegree with degree and asso size

# degree
p1  <- ggplot(DfInterOrga, aes(y= wdegree, x = degree, color = Acro)) + 
  geom_point()+
  labs(title = "A1", y = "Degré pondéré", x = "Degré")+
  gghighlight(degree <2 | wdegree > 2000 ,  label_key = Acro)  +
  theme_light()

# asso size
p2 <- ggplot(DfInterOrga %>% filter(AssoSize < 2000), aes(y= wdegree, x = AssoSize, color = Acro)) + 
  geom_point()+
  labs(title = "B1", y = "Degré pondéré", x = "Taille de l'association (nb membres)", 
       caption = "Note : Covenant of Mayors a été retiré\n(Taille = 10 000) " )+
  gghighlight(AssoSize> 1500 ,  label_key = Acro)  +
  theme_light()

### Regressions

  #1. regression weighted degree / degree
reg <- lm(log10(wdegree) ~ degree , data = Fgraph2)
summary(reg)


## Equation de la droite de regression :
eq = paste0("y = ", round(10^reg$coefficients[2],3), "^x * ", round(10^reg$coefficients[1],3)) #log10(y) = ax + b -> y = 10^(ax+b)= 10^(ax) * 10^b= (10^a)^x * 10^b



## Residuals
### add residuals and standart residuals to df
Fgraph2 <- Fgraph2 %>% 
  mutate(rezStand = residuals(reg, type = "pearson")) %>% 
  ungroup()

sdRez <- sd(Fgraph2$rezStand)

## Outliers
is_outlier <- function(x) {
  
  return(x < -1.5 * sdRez | x > 1.5 * sdRez)
  
}

### Ajout d'une variable Outlier au DF
Fgraph2 <- Fgraph2 %>%
  mutate(outlier_rezStand = ifelse(is_outlier(rezStand), 
                                   rezStand, 
                                   as.numeric(NA)))

## ggplot


regplot1 <- ggplot(Fgraph2, aes(x = degree, y = wdegree)) +
  geom_point () +
  labs(x = "Taille de l'association (nb membres)", 
       y = "Degré pondéré") + 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="#E69F00",
              linetype = "dashed", size = 1.5) +
  geom_label_repel(data = Fgraph2 %>% filter(!is.na(outlier_rezStand)), 
                   aes(label = Acro),
                   na.rm = TRUE, nudge_y = 0.05, color = "black", size = 2.5) +
  annotate(geom = "label", x = 25, y = 1000, label= paste0(eq, "\nR2 = ", 
                                                           round(summary(reg)$r.squared, 2)), hjust = 0, fill = "#E69F00", size = 3) +
  labs(title = "A2", y = "Degré pondéré (log10)", x = "Degré", 
       caption = "Note : Pour A2 et B2,\nles 3 associations outliers\nont été retirées")+
  scale_y_continuous(trans= 'log10')+
  theme_light()

pdf(file = "OUT/reg_wdgree_degree_ETMUN.pdf",width = 8.3, height = 5.8, pagecentre =FALSE)
regplot1
dev.off()

  #2. regression weighted degree / Size Asso

reg <- lm(wdegree ~ AssoSize , data = Fgraph2)
summary(reg)


## Equation de la droite de regression :
eq = paste0("y = ", round(reg$coefficients[2],2), " * x + ", round(reg$coefficients[1],2))

## Residuals
### add residuals and standart residuals to df
Fgraph2 <- Fgraph2 %>% 
  mutate(rezStand = residuals(reg, type = "pearson")) %>% 
  ungroup()

sdRez <- sd(Fgraph2$rezStand)


### Ajout d'une variable Outlier au DF
Fgraph2 <- Fgraph2 %>%
  mutate(outlier_rezStand = ifelse(is_outlier(rezStand), 
                                   rezStand, 
                                   as.numeric(NA)))

## ggplot


regplot2 <- ggplot(Fgraph2, aes(x = AssoSize, y = wdegree)) +
  geom_point () +
  theme_light()+
  labs(x = "Taille de l'association (nb membres)", 
       y = "Degré pondéré") + 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="#E69F00",
              linetype = "dashed", size = 1.5) +
  geom_label_repel(data = Fgraph2 %>% filter(!is.na(outlier_rezStand)), 
                   aes(label = Acro),
                   na.rm = TRUE, nudge_y = 0.05, color = "black", size = 2.5) +
  annotate(geom = "label", x = 310, y = 200, label= paste0(eq, "\nR2 = ", 
                                                           round(summary(reg)$r.squared, 2)), hjust = 0, fill = "#E69F00", size = 3) +
  labs(title= "B2", caption = "Sources : ETMUN 2019 / PG 2020", size = 2.5) 

pdf(file = "OUT/reg_wdgree_size_ETMUN.pdf",width = 8.3, height = 5.8, pagecentre =FALSE)
regplot2
dev.off()

## Graphiques regressions outliers filtrés

p1bis <- regplot1

p2bis <- regplot2

# Planche grahiques relation wdegree avec size et degree


grid3 <- (p1 | p2 ) / (p1bis| p2bis)

ggsave(grid3, filename = "OUT/interOrgaETMUN_wDbyD_andSize.pdf", width = 8.3, height = 8.3, units = "in" )



 

####====Compute same variables for a sub graph (with 3 asso filtered)====

subgraph1 <- TdInterOrga %>% activate(nodes) %>% filter(AssoSize < 500) %>% filter(degree>2) 

V(subgraph1)$degree <- degree(subgraph1)

V(subgraph1)$wdegree <- strength(subgraph1, vids = V(subgraph1), loops = F )

subgraph1 <- subgraph1 %>% 
              activate(nodes) %>%
                mutate(relative_wdegree = wdegree/(AssoSize*(graph_order()-1)))%>%
                mutate(louvain = group_louvain())
                
modularity(subgraph1, membership = V(subgraph1)$louvain)
graph.density(subgraph1)


#Plot the subgraph
filterG2 <- subgraph1 %>% 
           activate(edges) %>% filter(weight > 20) %>% 
          activate(nodes) %>% filter(!node_is_isolated())

filterG2 <- filterG2 %>% mutate(louvain = group_louvain( ))
modularity(filterG2, membership = V(filterG2)$louvain)


l <- layout_with_fr(filterG2, niter = 10000, weights = sqrt(E(filterG2)$weight)/1000)
l <- layout_with_kk(filterG2,weights = sqrt(E(filterG2)$weight)/10)

gsub2<-ggraph(filterG2, layout = l) + 
  geom_edge_link(aes(width = weight), alpha = 0.2)+
  geom_node_point( aes(size = relative_wdegree, color = as.factor(louvain))) + 
  geom_node_text(aes(label = Acro), repel = T, size = 2.5)+
  scale_edge_width(range = c(0.5, 8)) +# control size
  scale_size(range = c(1,8))+
  labs(size = "Degré\n(nb d'associations\navec au moins\nune ville membre\nen commun)",
       caption = "Note : Seuls les liens supérieurs\nà 20 sont représentés\n\nSources : ETMUN 2019. PG, 2020")+
  theme_graph(base_family = 'Helvetica')
gsub2 
####==== Dominant flows =====##




