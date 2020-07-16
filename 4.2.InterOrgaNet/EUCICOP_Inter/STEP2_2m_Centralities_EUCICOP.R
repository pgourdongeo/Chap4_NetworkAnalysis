############################### EUCICOP bipartite Networks STEP 2 ###################
#                               
#                          
# DESCRIPTION : Travail sur le graphe EUCICOP (affiliation) biparti. 
#               Indices locaux,distributions des degrés, ego network
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/EUCICOP_Inter")

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

options(scipen = 999)

### Data 

## Participation (partnership): table City-Project. (for edges). EUROPE frame filtered (check STEP0 in Data folder)
ParticipationEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/Participations_All_Eucicop_Europe.RDS")

## Information on project (for nodes)
ProjectEUCICOP <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/EUCICOP/ProjectsEucicop_all_noduplicated.RDS")

## GN info for cities
DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")

# Matrices from STEP1  

# matrix with all edges but unweighted (unique member localitities in each association)
emnw <- readRDS("DataProd/Eucicop_emnw.rds")

# Matrix unweighted and with only member cities more than 1 degree (at least in 2 associations)
em2nw <- readRDS("DataProd/Eucicop_em2nw.rds")

# First component of emnw
emnw_1stComp <- readRDS("DataProd/Eucicop_emnw_1stcomp.rds")


# First component of emnw
em2nw_1stComp <- readRDS("DataProd/Eucicop_em2nw_1stcomp.rds")


### NB : WE WORK ON THE 1ST COMPONENT ONLY 

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





Centralities <- Centralitydf(emnw_1stComp, em2nw_1stComp)

saveRDS(Centralities, "DataProd/Centralities_EUCICOP_enww_em2nw_1stComp.rds")
### plot of degree distributions

# hist
Centralities[[1]] %>% keep(is.numeric) %>% select(starts_with("N")) %>%
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram() + 
  scale_y_log10() + 
  # scale_x_log10()+
  labs(title = "Distributions des mesures de centralité des villes EUCICOP",
       subtitle = "Premières composantes de emnw et em2nw", 
       y = "Nombre de villes (log10)", x = "",
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")


summary(Centralities[[1]]$ND)
summary(Centralities[[1]]$D)
# ggsave(filename = "OUT/CentralDistrib_2matrix_Cities_EUCICOP.pdf",width = 8.3, height = 5.8 )

Centralities[[2]] %>% keep(is.numeric) %>% select(starts_with("N")) %>% 
  gather() %>% 
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free", ncol = 2) +
  geom_histogram()+ 
  scale_y_log10() + 
  # scale_x_log10(n.breaks = 6) +
  labs(title = "Distributions des mesures de centralité des projets EUCICOP",
       subtitle = "Premières composantes de emnw et em2nw", 
       y = "Nombre de projets", x = "Valeurs normalisées (échelle log 10)",
       caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation,  2 = em2nw\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")

summary(Centralities[[2]]$ND)
summary(Centralities[[2]]$D)
# ggsave(filename = "OUT/CentralDistrib_2matrix_Project_EUCICOP.pdf",width = 8.3, height = 5.8 )

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes EUCICOP (emnw)")

Centralities[[1]] %>% 
  keep(is.numeric) %>% 
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des villes EUCICOP (em2nw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>% 
  select(-ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets EUCICOP (emnw)")

Centralities[[2]] %>% 
  keep(is.numeric) %>%
  select(ends_with("2") & starts_with("N")) %>% 
  ggpairs(title = "Relations entre les mesures de centralité des projets EUCICOP (em2nw)")

## cumulative normalized degree
ytext <- "Fréquence cumulée (log10)"
xtext <- "Degré normalisé (log10)"

g1 <-ggplot(Centralities[[1]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des villes EUCICOP (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g2 <- ggplot(Centralities[[1]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des villes EUCICOP (em2nw)", y = ytext,  x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g3 <- ggplot(Centralities[[2]], aes(ND)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des associations EUCICOP (emnw)", y = ytext, x = xtext) + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()

g4 <- ggplot(Centralities[[2]], aes(ND2)) + 
  stat_ecdf(geom = "line", pad = FALSE) + 
  labs(title = "Distribution des degrés des associations EUCICOP (em2nw)", y = ytext, x = xtext,
       caption = "Sources : ETMUN 2019\nPG 2020") + 
  scale_x_log10( n.breaks =10) + 
  scale_y_log10()


grid1 <- (g1 | g2) / (g3 | g4)
grid1

#Clean objects
rm(g1,g2,g3,g4)

### ==== POWER LAW OF DEGREE DISTRIBUTION ====


fit1 <- fit_power_law(Centralities[[1]]$ND, implementation = c("plfit"))
fit1
ND2 <- Centralities[[1]] %>% filter(!is.na(ND2))%>% select(ND2)%>% deframe()
fit2 <- fit_power_law(ND2, implementation = c("plfit"))
fit2
fit3 <- fit_power_law(Centralities[[2]]$ND, implementation = c("plfit"))
fit3

ND2 <- Centralities[[2]] %>% filter(!is.na(ND2))%>% select(ND2)%>% deframe()
fit4 <- fit_power_law(ND2, implementation = c("plfit"))
fit4

g <- graph.incidence(emnw_1stComp, directed = FALSE)
d <- degree(g, normalized = TRUE) %>% as.data.frame()
colnames(d) <- "degree"
ggplot(d)+  geom_histogram(aes(degree))+ scale_x_log10()+ scale_y_log10()
fit5 <- fit_power_law(d$degree, implementation = c("plfit"))
fit5

fit6 <- fit_power_law(append(Centralities[[2]]$ND, Centralities[[1]]$ND), implementation = c("plfit"))
fit6


fit7 <-  fit_power_law(append(Centralities[[2]]$D, Centralities[[1]]$D), implementation = c("plfit"))
fit7


d2 <- append(Centralities[[2]]$ND, Centralities[[1]]$ND) %>% as.data.frame()


colnames(d2) <- "degBi"
ggplot(d2)+  geom_histogram(aes(degBi))+ scale_x_log10(n.breaks = 10)+ scale_y_log10()

nunderminx <- d2 %>% filter(degBi < fit6$xmin) %>% nrow()
1- nunderminx/nrow(d2)# powerlaw estimated on only 2% of the vertices...

rm(ND2, fit1,fit2,fit3,fit4, fit5, fit7, d, d2)
rm(g)

# library(bipartite)
# 
# test <-degreedistr(emnw_1stComp)
# test$`lower level dd fits`
# library(poweRlaw)
# m_pl =conpl$new(d2$degBi)         # discrete power law fitting
# est = estimate_pars(m_pl)# get xmin and alpha
# m_pl$setPars(est)
# bs_p = bootstrap_p(m_pl, no_of_sims = 100)
# bs_p$p


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

namecity <- ParticipationEUCICOP %>% select(geonameId, asciiName) %>% distinct()%>% 
  left_join(select(DbCity,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11))
Centralities[[1]] <- Centralities[[1]] %>% left_join(namecity)

# set names of association
Centralities[[2]] <- Centralities[[2]] %>% left_join(select(ProjectEUCICOP, Code = ID_PROJECT,  ProjectName = Acronym))

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
TopRank <- RankCities %>% filter_at(.vars = c(11:18), any_vars(.<6))

# tibble for ggplot rank
TopRankgg <- TopRank %>% select(geonameId | asciiName | ends_with("rank"))%>% 
  pivot_longer(-c(geonameId,asciiName), names_to = "indice", values_to = "rang")%>% 
  mutate(indice = str_remove(indice,pattern = "_rank")) %>% 
  mutate(indice = str_remove(indice,pattern = "N_")) %>%
  mutate(lograng = log10(rang))%>% mutate(asciiName = recode(asciiName, "Mestre"= "Venice"))

#plot in coord polar
s <- summary(TopRankgg$lograng)
bks <- c(s[[1]], log10(2), s[[2]], s[[3]], s[[5]] , s[[6]])
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
       caption = "D = Degré, B = Betweeness, C = Closeness, Eigen_C = Eigen vector centrality, N = normalisation, 2 = em2nw\nSource : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020") + 
  theme(axis.text.y=element_text(size=7),
        axis.text.x=element_blank(), legend.position="bottom")

ggsave(filename = "OUT/RankClock_TopAny5_EUCICOP.pdf",width = 8.3, height = 7 )



##Look project for top
ProjectsTopCities <- ParticipationEUCICOP %>% filter(geonameId %in% TopRank$geonameId)

ProjectsTopCities <- ProjectsTopCities %>% group_by(Programme, asciiName)%>% summarise(n= n()) 
villes <- c("Ljubljana", 	
            "Bologna", 	
            "Vienna", "Riga")
ProjectsTopCities <- ProjectsTopCities %>% group_by(asciiName) %>% top_n(5) %>% filter(asciiName %in% villes )

library(questionr)

irec(ProjectsTopCities, Programme)
## Recodage de ProjectsTopCities$Programme en ProjectsTopCities$Programme_rec
ProjectsTopCities$Programme_rec <- fct_recode(ProjectsTopCities$Programme,
               "2000 - 2006\nBaltic Sea Region" = "2000 - 2006 Baltic Sea Region",
               "2000 - 2006\nCadses" = "2000 - 2006 Cadses ",
               "2000 - 2006\nWestern Mediterranean" = "2000 - 2006 Western Mediterranean ",
               "2007 - 2013\nAlpine Space" = "2007 - 2013 Alpine Space",
               "2007 - 2013\nAustria - Czech Republic\n(AT-CZ)" = "2007 - 2013 Austria - Czech Republic (AT-CZ)",
               "2007 - 2013\nBaltic Sea Region" = "2007 - 2013 Baltic Sea Region",
               "2007 - 2013\nCentral Baltic\n(FI-SE-EE-LA)" = "2007 - 2013 Central Baltic (FI-SE-EE-LA)",
               "2007 - 2013\nCentral Europe" = "2007 - 2013 Central Europe",
               "2007 - 2013\nEstonia - Latvia (EE-LV)" = "2007 - 2013 Estonia - Latvia (EE-LV)",
               "2007 - 2013\nInterreg IVC" = "2007 - 2013 Interreg IVC",
               "2007 - 2013\nItaly - Slovenia (IT-SI)" = "2007 - 2013 Italy - Slovenia (IT-SI)",
               "2007 - 2013\nProgramme MED" = "2007 - 2013 Programme MED",
               "2007 - 2013\nSlovak Republic - Austria\n(SK-AT)" = "2007 - 2013 Slovak Republic - Austria (SK-AT)",
               "2007 - 2013\nSouth East Europe" = "2007 - 2013 South East Europe",
               "2014 - 2020\nINTERREG V-A\nFinland - Estonia - Latvia -\nSweden (Central Baltic)" = "2014 - 2020 INTERREG V-A Finland - Estonia - Latvia - Sweden (Central Baltic)")
ProjectsTopCities$Programme_rec <- as.character(ProjectsTopCities$Programme_rec)


ggplot(ProjectsTopCities) + 
  geom_bar(aes(x= reorder(Programme_rec, desc(n)), y = n, fill = asciiName), stat = "identity") + 
  facet_wrap(~ asciiName, scales ="free" ) + 
  coord_flip() + labs(y = "Nombre de projets", fill = "Ville", x = "Programme de coopération", 
                      caption = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")


ggsave(filename = "OUT/TopCitiesEUCICOP_programme_projects.pdf",width = 8.3, height = 5.8 )
rm(CityCentralities,TopRank,TopRankgg , s, bks, lbs)

#### ====ARTICULATION POINTS AND IGRAPH OBJECT ====

emnw.g <- graph.incidence(emnw_1stComp, directed = F)

em2nw.g <-  graph.incidence(em2nw_1stComp, directed = F)

#Info all nodes 
NameProjects <- ProjectEUCICOP %>% select(name = ID_PROJECT,label = Acronym, Date = Period )
nameCitiesjoint <- namecity %>% select(name= geonameId, label = asciiName, Country = countryCode, everything())

NodesInfos <- bind_rows(NameProjects,nameCitiesjoint) 
# saveRDS(NodesInfos, "DataProd/NodesInfos_EUCICOP.RDS")
rm(nameCitiesjoint)
## Articulation points
#emnw
art <-articulation_points(emnw.g)
vecArt <- V(emnw.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)

skim(art)## 3256  articulations points (597 are Cities)
#em2nw

art <-articulation_points(em2nw.g)
vecArt <- V(em2nw.g)[art]$name
art <- NodesInfos %>% filter(name %in% vecArt)
skim(art) ## 744 articulations points (597 are cities)

##Plot cities articulations

citiesArt <- art %>% filter(!str_detect(name, "P"))
citiesArt <- citiesArt %>% mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))%>%
  filter(PopAdmin11> 0)
citiesArt <- st_as_sf(citiesArt)

##Map

## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)

myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(citiesArt$PopAdmin11)
s[[1]]
bks <- c(5000,10000, 50000, 100000,  500000, 3000000 )
lbs <-  c("5000",  "10 000",  "50 000" , "100 000", "500 000", "3 M")


citiesArt2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = citiesArt ,
          mapping = aes(size = PopAdmin11), color = "#3384b0", show.legend = "point", alpha = 0.8) +
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(0.5, 8))+
  annotate("text", label = "Sources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = citiesArt %>% filter(PopAdmin11 > 100000), aes(label = label), size = 2.2, color = "#4d4d4d",
               check_overlap = TRUE) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + 
  guides(size=guide_legend(ncol=2) )+
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7.5))

pdf(file = "OUT/ArtPoint_Cities_em2nw_EUCICOP.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesArt2
dev.off()

#components
components(em2nw.g)$no  # 1
components(emnw.g)$no  # 1

rm(art,vecArt)

#### ==== FILTER OUT SOME ASSOCIATIONS ====
# ## show outliers by degrees (2 associations)
# 
# 
# #filter some asso outliers (visual detection) according to the ggpairs on Associations
# 
# Centralities[[2]]%>% filter(ND>0.10)# Covenant of mayors, Climate Alliance, WWCAM
# 
# Centralities[[2]] %>% keep(is.numeric) %>% 
#   select(-ends_with("2") & starts_with("N")) %>% filter(ND<0.10 )%>%
#   ggpairs(title = "Relations entre les mesures de centralité des associations ETMUN (em2nw)")
# 
# 
# # Covenant of mayors, Climate Alliance, EWT (subset of another asso)
# assoremove <- c("04080",  "03241", "07675")
# #assoremove <- c("Covenant of Mayors", "WWCAM","Climate Alliance")
# 
# #Filter out the association on the matrix that keep city of degree 1
# emnwF <- emnw[!rownames(emnw) %in% assoremove, ]
# emnwF <- emnwF[,colSums(emnwF) > 0]
# 
# #Filter out the association on the matrix that keep city of more than 1 degree (only cities involved at least in 2 associations)
# em2nwF <-  em2nw[!rownames(em2nw) %in% assoremove, ]
# em2nwF <- em2nwF[,colSums(em2nwF) > 1]
# 
# 
# saveRDS(emnwF, "DataProd/emnwF.RDS")
# saveRDS(em2nwF, "DataProd/em2nwF.RDS")
# ### Graph and Articulations points
# #emnwF
# 
# emnwF.g <- graph.incidence(emnwF, directed = F)
# components(emnwF.g)$no
# 
# art <-articulation_points(emnwF.g)
# vecArt <- V(emnwF.g)[art]$name
# art <- NodesInfos %>% filter(name %in% vecArt)## 49 articulations points (only associations)
# 
# #em2nwF
# em2nwF.g <- graph.incidence(em2nwF, directed = F)
# articulation.points(em2nwF.g)#none
# 
# 
# ### ==== FILTERED MATRIX DISTRIBUTION OF DEGREES and Centrality Measures ====
# 
# # Function from STEP1
# 
# graphlevelindex <- function(matrix){
#   require(igraph)
#   require(Matrix)
#   
#   #density
#   density <- nnzero(matrix)/ (ncol(matrix)*nrow(matrix))
#   #size
#   nbliens <- nnzero(matrix)
#   #order type 1
#   ncities <- ncol(matrix)
#   #order type 2
#   nasso <- nrow(matrix)
#   #Diameter
#   diam <- diameter(graph.incidence(matrix, directed = FALSE),directed = FALSE)
#   #degrees
#   meanDegCities <- mean(colSums(matrix != 0))
#   medDegCities <- median(colSums(matrix != 0))
#   meanDegAsso <- mean(rowSums(matrix != 0))
#   medDegAsso <- median(rowSums(matrix != 0))
#   
#   
#   result <- data.frame(Densité = density, 
#                        Diamètre = diam,
#                        Taille = nbliens,
#                        NbVilles = ncities,
#                        NbAssos = nasso,
#                        MeanDegreeCity =  meanDegCities,
#                        MedDegreeCity =  medDegCities,
#                        MeanDegreeAsso = meanDegAsso,
#                        MedDegreeAsso = medDegAsso)
#   return(result)
# }
# 
# # Summary of filtered matrices
# dimemnwF <- graphlevelindex(emnwF)
# 
# dimem2nwF <- graphlevelindex(em2nwF)
# 
# #Centralities Filtered
# 
# CentralityFiltered <- Centralitydf(emnwF, em2nwF)
# 
# ## Histogram centralities
# 
# CentralityFiltered[[1]] %>% keep(is.numeric) %>% 
#   gather() %>% 
#   ggplot(aes(value)) + 
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() + 
#   scale_y_log10() + 
#   labs(title = "Distributions des degrés des villes ETMUN",
#        subtitle = "Sur les matrices emnwF et em2nwF", 
#        y = "Nombre de villes (log10)", 
#        caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nwF")
# 
# 
# CentralityFiltered[[2]] %>% keep(is.numeric) %>% 
#   gather() %>% 
#   ggplot(aes(value)) + 
#   facet_wrap(~ key, scales = "free") +
#   geom_histogram() + 
#   scale_y_log10() + 
#   labs(title = "Distributions des degrés des Associations ETMUN",
#        subtitle = "Sur les matrices emnwF et em2nwF", 
#        y = "Nombre de villes (log10)", 
#        caption = "D = Degré, B = Betweeness, C = Closeness, N = normalisation, 2 = em2nwF")
# 


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

 saveRDS(egocity.emnw, "DataProd/Df_egocity_emnw_1stcomp.RDS")
# saveRDS(egocity.em2nw, "DataProd/Df_egocity_em2nw_1stcomp.RDS")
rm(ego.emnw, ego.emnw.df, ego.emnwF.df, ego.em2nwF.df)


egocity.emnw <- readRDS("DataProd/Df_egocity_emnw_1stcomp.RDS")
egocity.em2nw <-readRDS("DataProd/Df_egocity_em2nw_1stcomp.RDS")
### Exploration of ego network properties

egocity.emnw %>% keep(is.numeric) %>% select(-idList) %>%
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
  keep(is.numeric) %>% select(-idList) %>%
  select(-nbcomp, -sizecomp)%>%
  ggpairs(title = "Relations entre ropriétés des ego-network de chaque ville EUCICOP (emnw)")


egocity.em2nw %>% keep(is.numeric) %>% select(-idList) %>%
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
  keep(is.numeric) %>% select(-idList) %>%
  select(-nbcomp, -sizecomp)%>%
  ggpairs(title = "Relations entre propriétés des ego-networks de chaque ville EUCICOP (em2nw)")

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

egocity.em2nw <- egocity.em2nw %>% left_join(select(DfACP, name, Dim.1, Dim.2, Dim.3))

egocity.emnw <- egocity.emnw %>% left_join(select(DfACP, name, Dim.1, Dim.2, Dim.3))

## CAH on dim ACP 
library(cluster)
library(ggdendro)
library(reshape2)

# compute classification ----

suppressWarnings(source("../../multi_bivariateFunctions_exploratR/Functions_anova_cah.R"))


myVar <- c("Dim.1","Dim.2", "Dim.3")
cah <- ComputeClassif(df = egocity.emnw,
                      varquanti = myVar, method = "ward", stand = FALSE)


dendro <- PlotDendro(classifobj = cah)
dendro
inert <- PlotHeight(classifobj = cah)
inert
myProfiles <- PlotProfile(classifobj = cah, nbclus = 5)
myProfiles$PROFILE
table(myProfiles$CLUSID)
egocity.em2nw <- egocity.em2nw %>% mutate(CAH = myProfiles$CLUSID)
egocity.emnw <- egocity.emnw %>% mutate(CAH = myProfiles$CLUSID)

## plot profile of CAH with original variables
PlotCAHProfile <- egocity.emnw %>% select(c(4:10,22)) %>%  mutate(reinf = replace(reinf, is.nan(reinf), 0)) %>% pivot_longer(-CAH)

irec(PlotCAHProfile)
## Recodage de PlotCAHProfile$name en PlotCAHProfile$name_rec
PlotCAHProfile$name_rec <- fct_recode(PlotCAHProfile$name,
               "Nb Projets" = "nbAsso",
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
  labs(y = "", x = "", fill = "Typologies\ndes ego-networks\ndes villes EUCICOP\n(CAH)", 
  caption = "NB : matrice emnw\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 12))

ggsave( filename = "OUT/EUCICOPEGONetwork_CAH5Profiles_BoxPlot_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )
# saveRDS(cah , "DataProd/cah_ego_eucicop_em2nw.rds")
# saveRDS(cah , "DataProd/cah_ego_eucicop_emnw.rds")

## plot profile of CAH with centralities indexes

# Add centralities mesures

CentralCities_emnw <- Centralities[[1]] %>% 
  select(-ends_with("2"))%>% 
  as.data.frame()%>% select(-geometry) %>% 
  select(1:8) %>% rename(name = geonameId)

egocity.emnw <- egocity.emnw %>% left_join(CentralCities_emnw)

PlotCAHProfile <- egocity.emnw %>% select(22, 24,26,28,29) %>% pivot_longer(-CAH)



ggplot(PlotCAHProfile)+ 
  geom_boxplot(aes(x =CAH, fill = CAH, y = value), axes = TRUE, outlier.colour="black", outlier.size=0.5, outlier.alpha = 0.6)+ 
  facet_wrap(~ name, scales = "free_y")+ 
  labs(y = "", x = "", fill = "Typologies\ndes ego-networks\ndes villes EUCICOP\n(CAH)", 
       caption = "Note : ces variables de centralités\nne sont pas prises en compte\ndans la CAH.\nSources : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 12))

ggsave( filename = "OUT/EUCICOPEGONetwork_CAH5Profiles_BoxPlotCentralities_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )
### Explore clustering outcomes in relation to population and administrative status of cities



### Test anova Classe and Pop

# correct population variable
egocity.emnw <- egocity.emnw %>% 
  mutate(population = as.numeric(population)) %>% 
  mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11))

# Perform AnonVa

anovaTab <- AnovaTab(df = egocity.emnw , varx = c("CAH"), 
                     vary = c("PopAdmin11")) 
resume <- ComputeRegression(egocity.emnw, vardep = "PopAdmin11", varindep = "CAH")
resume$TABCOEF ## R2 = 0.08
rm(resume, anovaTab)
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

ggsave(filename = "OUT/HeatMap_PearsonChi2_CAHegoNetandAdmin_EUCICOP_emnw.pdf", device = "pdf",width = 8.3, height = 5.8 , units = "in")
class(resChi2)
testChi2$stdres

library(vcd)

cramer <-assocstats(tableCAH)
cramer  # 0.28

#### ==== MAPPING SOME EGO NETWORKS  =====

#just create  a variable name specific to project such as Acro in ETMUN

NodesInfos <- NodesInfos %>% mutate(Acro = ifelse(str_detect(name, "P"), label, NA))
# create a sample of egonetworks

sampleEgo <- egocity.em2nw  %>% group_by(CAH) %>% sample_n(1, replace = T)

sampleEgo <- egocity.emnw  %>% group_by(CAH) %>% sample_n(1, replace = T)
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
        scale_color_manual( values =c("lightskyblue", "mediumvioletred", "orange"), breaks = c("City", "ego", "Asso"), labels = c("Villes", "Ego", "Projets"))+
        scale_size(range = c(1.5,6))+
        labs(title = paste("Ego Network de ", strInfo$label,", ", strInfo$Country, sep = ""),
             subtitle = paste(strInfo$CAH), color = "Type de sommets",
             caption = "Algorithme de spatialisation : kk")+
        annotate("label",
                 label = paste("Densité : ", round(strInfo$densite,3),
                               "\nNb Villes : ", strInfo$nbVilles,
                               "\nNb Projets : ", strInfo$nbAsso,
                               "\nNb composantes : ", strInfo$nbcomp, 
                               "\nReinforcement : ", round(strInfo$reinf,4), 
                               sep = ""), 
                 x = -1, 
                 y = -0.5,
                 hjust = 0, 
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
egoplot.em2nw <- EgoMapping(df = sampleEgo,egolist = ego.em2nw,NodesInfos =NodesInfos)

egoplot.em2nw[[5]]

egoplot.emnw <- EgoMapping(df = sampleEgo,egolist = ego.emnw,NodesInfos =NodesInfos)
# Grid plot function

Gridplot <- function(myplots, n){
  require(cowplot)
  splitted_plots <- split(myplots, ceiling(seq_along(myplots)/n))
  
  lapply(splitted_plots, function(x) plot_grid(plotlist = x))
  
}


#Try
gEgo <- Gridplot(egoplot.em2nw, 4)

gEgo[[1]]
ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5.pdf", width = 11.7, height = 8.3, units = "in" )

gEgo[[2]]

ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5_UniqueClasse5.pdf", width = 8.3, height = 5.8 , units = "in" )


gEgo <- Gridplot(egoplot.emnw, 4)

gEgo[[1]]
ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5_emnw.pdf", width = 11.7, height = 8.3, units = "in" )

gEgo[[2]]

ggsave( filename = "OUT/EUCICOPEGONetwork_cluster5_UniqueClasse5_emnw.pdf", width = 8.3, height = 5.8 , units = "in" )
