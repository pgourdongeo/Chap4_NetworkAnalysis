############################### ETMUN Inter-Organisations Network  ############
#                               
#                          
# DESCRIPTION :  Projection 1-mode du biparti
#               Graphes asso = nombres de villes membres en commun.
#               Représentations graphiques + calcul des indicateurs relatifs
#
# 
############################################################################## PG mai-juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/ETMUN_Inter")
#Packages
library(tidyverse)
library(skimr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
#library(flows)
library(GGally)
library(patchwork)
library(gghighlight)
library(tidylog)

#Data 

## Membership : table City-Asso. (for edges)
MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_GNidCorr.RDS")

## Information on associations (for nodes)
AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", stringsAsFactors = F)

# Matrix from STEP2

em2nw <- readRDS("DataProd/Etmun_em2nw.rds")


### ==== CREATE NETWORK DATA ====

Network2modes <- graph.incidence(em2nw , directed = FALSE)

#check

is_bipartite(Network2modes)#TRUE

### 1-mode projections

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

# Set AssoSize of asso (equivalent of Asso degrees in the non-filtered and unweighted matrix = emnw in STEP1 )
##
edgelist <- MembershipEtmun %>% filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network)

edgelist <- MembershipEtmun %>% 
  filter(!is.na(geonameId)) %>% 
  select( geonameId, Code_Network) %>% 
  group_by(geonameId,Code_Network) %>% 
  summarise(weight = n())

edgelistnw <- edgelist %>% mutate(weight = 1) 

SizeAsso <- edgelistnw %>% group_by(Code_Network) %>% 
  summarise(nMembers = n())

InterOrga  <- set_vertex_attr(InterOrga,"AssoSize", index = SizeAsso$Code_Network, value = SizeAsso$nMembers)

### compute degrees

V(InterOrga)$degree <- degree(InterOrga)
V(InterOrga)$wdegree <- strength(InterOrga, vids = V(InterOrga), loops = F )

# transform into tidy graph and df

TdInterOrga <- as_tbl_graph(InterOrga)
DfInterOrga <- fortify.tbl_graph(TdInterOrga)

### ==== COMPUTE RELATIVE WEIGHTED DEGREES ====

# NOTE : the idea is to make the weighted variables (weighted degrees for nodes and weights for edges) 
# taking the theoretical maximum in case the original bipartite graph is complete. 
# This theoretical maximum thus depends on the size (number of member cities) of each association. 
# Example: The maximum weight of a edges between two associations A and B of sizes 10 and 20 
# is 10 (the 10 member cities of association A are also members of association B).

# Function

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

# Apply the function to the interorga DF

DfInterOrga2  <- CompRelative_wdegree(df = DfInterOrga,
                                      size = "AssoSize",
                                      wdegree = "wdegree")#"Density of weigted graph (proj one-mode of a 2modes graph) : 0.1603"

#bivariate overlook of new variable
DfInterOrga2 %>% 
  keep(is.numeric) %>% 
  filter(wdegree < 2000 & degree >2 & AssoSize < 1500) %>% 
  select(-Year, - max_thwdegree) %>% 
  ggpairs()

# set relative wdegree on the original graph
InterOrga <- set_vertex_attr(InterOrga,"relative_wdegree", index = DfInterOrga2$name, value = unlist(DfInterOrga2$norm_wdegree))

# Update tidy and df
TdInterOrga <- as_tbl_graph(InterOrga)

DfInterOrga <- fortify.tbl_graph(TdInterOrga) %>% 
  mutate(relative_wdegree = unlist(relative_wdegree ))


### ==== COMPUTE RELATIVE EDGES WEIGHT ==== 

#divide the weight by the smallest AssoSize of the two vertices
# get edges
edgesInterOrga <- TdInterOrga %>% activate(edges) %>% 
  fortify.tbl_graph()

#add vertices info
DfInterOrga <- DfInterOrga %>% 
  mutate(id.edge = as.numeric(rownames(.)))
#from
edgesInterOrga <- edgesInterOrga %>% 
  left_join(select(DfInterOrga, id.edge, Labelfrom = Label, AssoSizefrom = AssoSize), 
            by = c("from" = "id.edge"))
#to
edgesInterOrga <- edgesInterOrga %>% 
  left_join(select(DfInterOrga, id.edge, Labelto = Label, AssoSizeto = AssoSize), 
            by = c("to" = "id.edge"))

#get minimum assoSize

edgesInterOrga <- edgesInterOrga %>% 
  rowwise() %>%
  mutate(MinAsso = min(AssoSizefrom, AssoSizeto))

# Normalized weight
edgesInterOrga <- edgesInterOrga %>% 
  mutate(norm_weight = weight/ MinAsso)

# Set normalized weight in orignal graph
InterOrga <- set_edge_attr(InterOrga,"relative_weight", value = edgesInterOrga$norm_weight)

# Update tidy and df
TdInterOrga <- as_tbl_graph(InterOrga) %>% 
  mutate(relative_wdegree = unlist(relative_wdegree))

DfInterOrga <- fortify.tbl_graph(TdInterOrga)


### ==== FIG 4.3 GRAPH MAPPING ====

####### BASIC GRAPHS UNWEIGHTED AND WEIGHTED
## ggraph
## plot the graph unweighted
g1<-ggraph(TdInterOrga, layout = 'kk') + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point(color = "orange", aes(size = degree)) + 
  scale_size(range = c(0.01,5))+
  labs(title = "A", size = "Degré\n(nb d'associations\navec au moins\nune ville membre\nen commun)",
       caption = "Algorithme de spatialisation : kk")+
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
       caption = "Algorithme de spatialisation : FR\nNote : Seuls les liens supérieurs\nà 20 sont représentés\n\nSources : ETMUN 2019. PG, 2020")+
  theme_graph(base_family = 'Helvetica')
g2 



## Grid with two plots

grid1 <- g1 + g2 + plot_layout(ncol = 1, heights = c(1.5, 2))

#save and export
ggsave(grid1, filename = "OUT/GraphInterOrgaETMUN.pdf", width = 8.3, height = 8.3, units = "in" )

rm(grid1, g1,g2, FilterG)
### ==== EXPLORING NODE DEGREES and CHOICE TO FILTER OUT SOME ASSOCIATION ====
### ==== Main numerical outputs ====

## order and size
InterOrga
#Order = 59  , size 1230

## density (1-mode)

graph.density(InterOrga)# 0.72

## degree distribution

summary(V(InterOrga)$degree)

hist(V(InterOrga)$degree,breaks=10, xlab = "Degré des associations ETMUN", ylab = "Nombre d'associations", main = NA)

## Univariate overlook

DfInterOrga %>% keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

## bivariate overlook

DfInterOrga %>% keep(is.numeric) %>% select(-Year) %>% ggpairs()


####### Relation between node variables ####


## wdegree with degree and asso size

# degree
p1  <- ggplot(DfInterOrga, aes(y= wdegree, x = degree, color = Acro)) + 
  geom_point()+
  labs(title = "A1", y = "Degré pondéré", x = "Degré")+
  gghighlight(degree <2 | wdegree > 2000 ,  label_key = Acro)  +
  theme_light()
p1
# asso size
p2 <- ggplot(DfInterOrga %>% filter(AssoSize < 2000), aes(y= wdegree, x = AssoSize, color = Acro)) + 
  geom_point()+
  labs(title = "B1", y = "Degré pondéré", x = "Taille de l'association (nb membres)", 
       caption = "Note : Covenant of Mayors a été retiré\n(Taille = 10 000) " )+
  gghighlight(AssoSize> 1500 ,  label_key = Acro)  +
  theme_light()
p2

### ==== FIG 4.4. FILTER OUT 3 ASSOCIATION AND EXPLORE BIVARIATE RELATIONS====

TDInterOrgaFiltered <- TdInterOrga %>% activate(nodes) %>% filter(AssoSize < 500) %>% filter(degree>2) 
# remove Covenant of Mayors (size= 10 000) and Climate Alliance (size =1742) + WWCAM (degree =1)

DfInterOrgaFiltered <- fortify.tbl_graph(TDInterOrgaFiltered) 



### Regressions on filtered Data relation wdegree with assoSize and degree

#1. regression weighted degree / degree
reg <- lm(log10(wdegree) ~ degree , data = DfInterOrgaFiltered)
summary(reg)


## equation :
eq = paste0("y = ", round(10^reg$coefficients[2],3), "^x * ", round(10^reg$coefficients[1],3)) #log10(y) = ax + b -> y = 10^(ax+b)= 10^(ax) * 10^b= (10^a)^x * 10^b



## Residuals
### add residuals and standart residuals to df
DfInterOrgaFiltered <- DfInterOrgaFiltered %>% 
  mutate(rezStand = residuals(reg, type = "pearson")) %>% 
  ungroup()

sdRez <- sd(DfInterOrgaFiltered$rezStand)

## Outliers
is_outlier <- function(x) {
  
  return(x < -1.5 * sdRez | x > 1.5 * sdRez)
  
}

### adding an Outlier variable in the DF
DfInterOrgaFiltered <- DfInterOrgaFiltered %>%
  mutate(outlier_rezStand = ifelse(is_outlier(rezStand), 
                                   rezStand, 
                                   as.numeric(NA)))

## ggplot


regplot1 <- ggplot(DfInterOrgaFiltered, aes(x = degree, y = wdegree)) +
  geom_point () +
  labs(x = "Taille de l'association (nb membres)", 
       y = "Degré pondéré") + 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="#E69F00",
              linetype = "dashed", size = 1.5) +
  geom_label_repel(data = DfInterOrgaFiltered %>% filter(!is.na(outlier_rezStand)), 
                   aes(label = Acro),
                   na.rm = TRUE, nudge_y = 0.05, color = "black", size = 2.5) +
  annotate(geom = "label", x = 25, y = 1000, label= paste0(eq, "\nR2 = ", 
                                                           round(summary(reg)$r.squared, 2)), hjust = 0, fill = "#E69F00", size = 3) +
  labs(title = "A2", y = "Degré pondéré (log10)", x = "Degré", 
       caption = "Note : Pour A2 et B2,\nles 3 associations outliers\nont été retirées")+
  scale_y_continuous(trans= 'log10')+
  theme_light()

regplot1

#export
pdf(file = "OUT/reg_wdgree_degree_ETMUN.pdf",width = 8.3, height = 5.8, pagecentre =FALSE)
regplot1
dev.off()

#2. regression weighted degree / Size Asso

reg <- lm(wdegree ~ AssoSize , data = DfInterOrgaFiltered)
summary(reg)


## equation :
eq = paste0("y = ", round(reg$coefficients[2],2), " * x + ", round(reg$coefficients[1],2))

## Residuals
### add residuals and standardized residuals to df
DfInterOrgaFiltered <- DfInterOrgaFiltered %>% 
  mutate(rezStand = residuals(reg, type = "pearson")) %>% 
  ungroup()

sdRez <- sd(DfInterOrgaFiltered$rezStand)


### adding an Outlier variable in the DF
DfInterOrgaFiltered <- DfInterOrgaFiltered %>%
  mutate(outlier_rezStand = ifelse(is_outlier(rezStand), 
                                   rezStand, 
                                   as.numeric(NA)))

## ggplot


regplot2 <- ggplot(DfInterOrgaFiltered, aes(x = AssoSize, y = wdegree)) +
  geom_point () +
  theme_light()+
  labs(x = "Taille de l'association (nb membres)", 
       y = "Degré pondéré") + 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="#E69F00",
              linetype = "dashed", size = 1.5) +
  geom_label_repel(data = DfInterOrgaFiltered %>% filter(!is.na(outlier_rezStand)), 
                   aes(label = Acro),
                   na.rm = TRUE, nudge_y = 0.05, color = "black", size = 2.5) +
  annotate(geom = "label", x = 310, y = 200, label= paste0(eq, "\nR2 = ", 
                                                           round(summary(reg)$r.squared, 2)), hjust = 0, fill = "#E69F00", size = 3) +
  labs(title= "B2", caption = "Sources : ETMUN 2019 / PG 2020", size = 2.5) 

regplot2

#export
pdf(file = "OUT/reg_wdgree_size_ETMUN.pdf",width = 8.3, height = 5.8, pagecentre =FALSE)
regplot2
dev.off()

## Plot regression on filtered data (3 associations removed)

p1bis <- regplot1

p2bis <- regplot2

# Making a grid of plot 


grid2 <- (p1 | p2 ) / (p1bis| p2bis)

ggsave(grid2, filename = "OUT/interOrgaETMUN_wDbyD_andSize.pdf", width = 8.3, height = 8.3, units = "in" )

rm(p1,p2,p1bis,p2bis, regplot1,regplot2)
### ==== GRAPH REPRESENTATION USING RELATIVE VALUE ====

skim(E(TDInterOrgaFiltered)$relative_weight)

# discretize edges relative weight
TDInterOrgaFiltered <- TDInterOrgaFiltered %>% 
  activate(edges) %>% 
  mutate(relative_weight_class = case_when(relative_weight <= 0.05 ~ "1.Très Faibles, moins de 5%",
                                           relative_weight > 0.05 & relative_weight <= 0.10 ~ "2.Faibles, entre 5 et 10%",
                                           relative_weight >0.10 & relative_weight <= 0.25 ~ "3.Forts, entre 10 et 25%",
                                           relative_weight >0.25 & relative_weight <= 0.50 ~ "4.Très Forts, entre 25 et 50%",
                                           relative_weight > 0.50  ~ "5.Exceptionels, plus de 50%"))


#l <- layout_with_fr( TDInterOrgaFiltered, niter = 100, weights = sqrt(E(TDInterOrgaFiltered)$relative_weight))
l <- layout_with_lgl(TDInterOrgaFiltered)

#Change name relative weight in french before plotting

g3df <-  TDInterOrgaFiltered  %>% 
  activate(edges) %>% rename( "Poids relatif (%)" = relative_weight)

#Plot
g3<- ggraph(g3df, layout = l) + 
  geom_edge_link(aes(width = `Poids relatif (%)`, alpha = `Poids relatif (%)`))+
  geom_node_point( aes(color = relative_wdegree), size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  scale_color_continuous(low = "lightgoldenrod1", high = "firebrick3" )+
  scale_edge_width(range = c(0.1, 5)) +
  scale_edge_alpha(range = c(0.05,0.2))+
  labs(title = "Ensemble des liens",color = "Degré pondéré relatif (%)", 
  caption = "Poids relatif % (liens) = part de villes en commun observée\npar rapport au maximun théorique.\n \nAlgorithme de spatialisation : lgl")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)

g3 

# facet the graph by class of normalized weight


g4<-ggraph(TDInterOrgaFiltered, layout = l) + 
  geom_edge_link(alpha = 0.2)+
  geom_node_point( aes(color = relative_wdegree), size = 3 )+ 
  scale_color_gradient(low = "lightgoldenrod1", high = "indianred3" )+
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(color = "Degré pondéré relatif (%)")+
  theme_graph(base_family = 'Helvetica')+
  facet_edges(~ relative_weight_class)
g4

# Create a facet grid manually (to recompute layout and remove isolated nodes each time)
l = "graphopt"

skim(E(TDInterOrgaFiltered)$relative_weight)
#1 all edges


fg1 <- g3

#2

f2 <- TDInterOrgaFiltered %>% 
  activate(edges)%>% 
  filter(relative_weight > 0.10)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg2 <- ggraph(f2, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Liens > 10% (médiane)")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg2

#3

f3 <- TDInterOrgaFiltered %>% 
  activate(edges)%>% 
  filter(relative_weight > 0.25)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

fg3 <- ggraph(f3, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Liens > 25% (3ème quartile)")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg3

#4

f4 <- TDInterOrgaFiltered %>% 
  activate(edges)%>% 
  filter(relative_weight > 0.50)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

l <- layout_with_graphopt(f4, charge = 0.05,mass =10)

fg4 <- ggraph(f4, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Liens > 50%")+
  theme_graph(base_family = 'Helvetica', foreground = "grey60", border = TRUE)
fg4

# 5
f5 <- TDInterOrgaFiltered %>% 
  activate(edges)%>% 
  filter(relative_weight > 0.70)%>%
  activate(nodes) %>% 
  filter(!node_is_isolated())

l <- layout_with_graphopt(f5, charge = 0.05,mass =10)

fg5 <- ggraph(f5, layout = l) + 
  geom_edge_link(alpha = 0.3)+
  geom_node_point(color = "indianred3", size = 3 )+ 
  geom_node_text(aes(label = Acro),repel = TRUE, size = 2.5)+
  labs(title = "Liens > 70%",  
       caption = "Algorithme de spatialisation : graphopt\n \nSources : ETMUN 2019. PG, 2020")+
  theme_graph(base_family = 'Helvetica',foreground = "grey60", border = TRUE)
fg5

# grid

grid3 <- fg1 + (fg4 / fg5) +  plot_layout(ncol = 2, width =  c(1.5, 1))

#save
ggsave(grid3, filename = "OUT/EtmunInterOrga_Norm_grid.pdf", width = 11.7, height = 8.3, units = "in" )



