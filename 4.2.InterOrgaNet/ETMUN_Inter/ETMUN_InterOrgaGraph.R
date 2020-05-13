############################### ETMUN Inter Organisations Network  ############
#                               
#                          
# DESCRIPTION : Création des données relationnelles à partir de la table
#               des adhésions ETMUN.
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
library(visNetwork)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(flows)
#Data 

MembershipEtmun <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/ETMUN_Membership_GNidCorr.RDS")

AssoEtmun <- read.csv2("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/ETMUN/BD_ETMUN_OrganizationsWithMembersCities.csv", stringsAsFactors = F)


### ==== Create table Asso-Asso by shared members ====

# simple graph


edgelist <- MembershipEtmun %>% filter(!is.na(geonameId)) %>% select( geonameId, Code_Network)

Network2modes <- graph.data.frame(edgelist)

V(Network2modes)$type <- V(Network2modes)$name %in% edgelist[,1]

# bi partite

NetProj <- bipartite.projection(Network2modes)

InterOrga <- NetProj$proj1

# Set full name
InterOrga  <-  set_vertex_attr(InterOrga,"Label", index = AssoEtmun$Code, value = AssoEtmun$Name)
# Set Short Name
InterOrga <-  set_vertex_attr(InterOrga,"Acro", index = MembershipEtmun$Code_Network, value = MembershipEtmun$Network_Name)
# Set Year of creation
InterOrga  <- set_vertex_attr(InterOrga,"Year", index = AssoEtmun$Code, value = AssoEtmun$Date)
# Set country of the seat
InterOrga  <- set_vertex_attr(InterOrga,"CountrySeat", index = AssoEtmun$Code, value = AssoEtmun$Country..secretariat.)

#degree

V(InterOrga)$degree <- degree(InterOrga)
V(InterOrga)$wdegree <- strength(InterOrga, vids = V(InterOrga), loops = F )




# Set size asso

SizeAsso <- edgelist %>% group_by(Code_Network)%>% summarise(nMembers = n())
InterOrga  <- set_vertex_attr(InterOrga,"SizeAsso", index = SizeAsso$Code_Network, value = SizeAsso$nMembers)



# transform into tidy graph
TdInterOrga <- as_tbl_graph(InterOrga)

### ==== Graph representation ====

## HeatMap

palf <- colorRampPalette(c("yellow", "red")) 

heatmap(netm, Rowv = NA, Colv = NA, col = palf(100),
        scale="none", margins=c(10,10) )
legend(x="bottomright", legend=c(min(netm), round(mean(netm),2), max(netm)), 
       fill=colorRampPalette(c("yellow", "red"))(3))




## ggraph
FilterG <- TdInterOrga %>% 
  activate(edges) %>% 
  filter(weight > 20 ) %>% 
  activate(nodes) %>% 
  filter(!node_is_isolated())
# change weight name
FilterG <- FilterG %>% activate(edges) %>% rename( "Poids\n(nb villes membres\nen commun)" = weight)
#plot
g1<-ggraph(FilterG, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(width = `Poids\n(nb villes membres\nen commun)`), alpha = 0.2)+
  geom_node_point(color = "orange", aes(size = degree)) + 
  geom_node_text(aes(label = Acro), repel = T, size = 2.5)+
  scale_edge_width(range = c(0.1, 10)) +# control size
  scale_size(range = c(1,8))+
  labs(size = "Degré\n(nb d'associations\navec au moins\nune ville membre\nen commun)",
       caption = "Note : Seuls les liens supérieurs\nà 20 sont représentés\n\nSources : ETMUN 2019. PG, 2020")+
  theme_graph(base_family = 'Helvetica')
 g1 

ggsave(g1, filename = "OUT/GraphInterOrgaETMUN.pdf", width = 8.3, height = 5.8, units = "in" )




### ==== Main numerical outputs ====

## order and size
InterOrga
#Order = 59  , size 1230

# density

graph.density(InterOrga)# 0.72

# degree distribution


hist(V(InterOrga)$degree,breaks=10, xlab = "Degré des associations ETMUN", ylab = "Nombre d'associations", main = NA)

#### Size Asso and Weigthed degree

Fgraph2 <- TdInterOrga %>% activate(nodes) %>% filter(SizeAsso < 500) # remove Covenant of Mayors (n= 10 000) and Climate Alliance (n =1742)

Fgraph2 <- fortify.tbl_graph(Fgraph2)
ggplot(Fgraph2, aes(x = SizeAsso, y = wdegree))+
  geom_point()

# regression

reg <- lm(wdegree ~ SizeAsso , data = Fgraph2)
summary(reg)


## Equation de la droite de regression :
eq = paste0("y = ", round(reg$coefficients[2],2), " * x + ", round(reg$coefficients[1],2))

## Residuals
### add residuals and standart residuals to df
Fgraph2 <- Fgraph2 %>% 
  mutate(rezStand = residuals(reg, type = "pearson")) %>% 
  ungroup()

sdRez <- sd(Fgraph2$rezStand)

## Outliers
is_outlier <- function(x) {
  
  return(x < -1 * sdRez | x > 1 * sdRez)
  
}

### Ajout d'une variable Outlier au DF
Fgraph2 <- Fgraph2 %>%
  mutate(outlier_rezStand = ifelse(is_outlier(rezStand), 
                                   rezStand, 
                                   as.numeric(NA)))

## ggplot


regplot <- ggplot(Fgraph2, aes(x = SizeAsso, y = wdegree)) +
  geom_point () +
  theme_light() +
  labs(x = "Taille de l'association (nb membres)", 
       y = "Degré pondéré\n(nb de villes membres en commun avec les autres associations)") + 
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="#E69F00",
              linetype = "dashed", size = 1.5) +
  geom_label_repel(data = Fgraph2 %>% filter(!is.na(outlier_rezStand)), 
                   aes(label = Acro),
                   na.rm = TRUE, nudge_y = 0.05, color = "black", size = 2.5) +
  annotate(geom = "label", x = 350, y = 200, label= paste0(eq, "\nR2 = ", 
                                                          round(summary(reg)$r.squared, 2)), hjust = 0, fill = "#E69F00", size = 3) +
  labs(caption = "Sources : ETMUN 2019 / PG 2020", size = 2.5) 

pdf(file = "OUT/reg_wdgree_size_ETMUN.pdf",width = 8.3, height = 5.8, pagecentre =FALSE)
regplot
dev.off()


### Dominant flows


