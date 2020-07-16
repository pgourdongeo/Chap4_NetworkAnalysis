############################### MAPPING TOP CITIES CENTRALITIES ###################
#                               
#                          
# DESCRIPTION : Cartographie des TOP50 centralités (STEP2) pour ETMUN et EUCICOP
#
# 
############################################################################## PG juin 2020


### ==== LOAD PACKAGES AND DATA ====
setwd("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/4.2.InterOrgaNet/")

#Packages
library(tidyverse)
library(tidylog)
library(sf)


#####Data (From each STEP 2)

## TOP 50 ETMUN
Top50ETMUN <- readRDS("ETMUN_Inter/DataProd/Top50Cities_CentralitiesETMUN.rds")

## TOP 50 EUCICOP
Top50EUCICOP <- readRDS("EUCICOP_Inter/DataProd/Top50Cities_Centralities_EUCICOP.rds")

## GN info for cities

DbCity <- readRDS("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/CITY_GN/DBCity_LauUmzFua.rds")

## rec


rec <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/rec_3035.geojson")

#EU 
sfEU <- st_read("~/Chap4_NetworkAnalysis/Chap4_NetworkAnalysis/Data/Geometry/fondEuropeLarge.geojson", stringsAsFactors = FALSE,crs = 3035)



### ==== JOINING INFO GN TO THE TOP 50 DATAFRAMES ====


Top50ETMUN <- Top50ETMUN %>% 
  left_join(select(DbCity,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11, geometry))

Top50EUCICOP <- Top50EUCICOP %>% 
  left_join(select(DbCity,geonameId,countryCode,subregion, adminLevel , population, PopAdmin11, geometry))

# Replace NA in pop admin

Top50ETMUN <- Top50ETMUN %>% mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11)) %>% st_as_sf()

Top50EUCICOP <- Top50EUCICOP %>% mutate(PopAdmin11 = ifelse(is.na(PopAdmin11), population, PopAdmin11)) %>% st_as_sf()




### ==== JOINING INFO GN TO THE TOP 50 DATAFRAMES ====


#####  ETMUN ######

### create a simple and pretty scale bar 500km
myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(Top50ETMUN$B)
s[[1]]
bks <- c(s[[1]],s[[2]],s[[3]], s[[5]], s[[6]])
lbs <- c("4000", "35 000","45 000", "110 000", "240 000" )


citiesEtmun2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = Top50ETMUN ,
          mapping = aes(size = B,  colour = NC), show.legend = NA, alpha = 0.9) +
  scale_color_gradient( low= "lightyellow", high ="firebrick3", n.breaks = 5)+
  geom_sf(data = Top50ETMUN  ,
          mapping = aes(size = B), shape = 1, colour = "grey60", show.legend = NA) +
  scale_size(name = "Betweeness Centralités des villes",
             breaks = bks,
             labels = lbs,
             range = c(2, 11))+
  annotate("text", label = "Source : ETMUN 2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = Top50ETMUN, aes(label = asciiName), size = 2.2, color = "#4d4d4d",
               check_overlap = FALSE, nudge_x = +100000) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Closeness Normalisée") +
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7.5))


## display end save
pdf(file = "OUT/ETMUN_top50Map_centralities.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesEtmun2
dev.off()
ggsave( file="OUT/ETMUN_top50Map_centralities.svg", plot=citiesEtmun2, width=8.3, height=5.8, device = "svg")

#####  EUCICOP ######

### create a simple and pretty scale bar 500km
myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(Top50EUCICOP$B)
s[[1]]
bks <- c(s[[1]],s[[2]],s[[3]], s[[5]], s[[6]] )
lbs <- c("1 M", "3 M","4 M", "6 M", "16.5 M" )

citiesEUCICOP2 <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", alpha = 0.9 ) +
  geom_sf(data = Top50EUCICOP ,
          mapping = aes(size = B,  colour = NC), alpha = 0.6) +
  scale_color_gradient( low= "lightyellow", high ="firebrick3", n.breaks = 5)+
  geom_sf(data = Top50EUCICOP  ,
          mapping = aes(size = B), shape = 1, colour = "grey60", show.legend = NA) +
  scale_size(name = "Betweeness Centralités des villes",
             breaks = bks,
             labels = lbs,
             range = c(2, 11)) +
  annotate("text", label = "Source : EUCICOP 2019 ; KEEP Closed Projects 2000-2019 / PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = Top50EUCICOP, aes(label = asciiName), size = 2.2, color = "#4d4d4d",
               check_overlap = FALSE, nudge_x = +100000) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Closeness Normalisée") +
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7.5))


## display end save
pdf(file = "OUT/EUCICOP_top50Map_centralities.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
citiesEUCICOP2
dev.off()
ggsave( file="OUT/EUCICOP_top50Map_centralities.svg", plot=citiesEUCICOP2, width=8.3, height=5.8, device = "svg")

### ==== MAP PRESENCE ABSENCE IN THE TOP 50 ====
# Prepare on sf with cities of the 2 tops
df1 <- Top50ETMUN %>% mutate(ETMUN = 1) %>% select(geonameId, ETMUN) %>%  st_drop_geometry()

df2 <- Top50EUCICOP %>% mutate(EUCICOP = 1) %>% select(geonameId, EUCICOP) %>%st_drop_geometry()

tops <-full_join(df1, df2, by= "geonameId")

tops <- tops %>% mutate(Network = case_when(ETMUN == 1 & EUCICOP == 1 ~ "Top centralité dans les deux réseaux",
                                         ETMUN == 1 & is.na(EUCICOP) ~ "Top centralité dans le réseaux ETMUN",
                                         is.na(ETMUN) & EUCICOP == 1 ~ "Top centralité dans le réseaux EUCICOP"))


tops <- tops %>%  
  left_join(select(DbCity,geonameId, asciiName,countryCode,subregion, adminLevel , population, PopAdmin11, geometry)) %>% st_as_sf()

tops <- tops %>%  mutate(PopAdmin11 = as.numeric(ifelse(is.na(PopAdmin11), population, PopAdmin11)))

### map

myScaleBar <- data.frame(X = c(c(st_bbox(rec)[3]-900000), c(st_bbox(rec)[3]-400000)),
                         Y = c(c(st_bbox(rec)[2]+200000), c(st_bbox(rec)[2]+200000)))

s <-summary(tops$PopAdmin11)
s[[1]]
bks <- c(s[[1]],s[[2]],s[[3]], s[[5]], s[[6]] )
lbs <- c("40 000", "200 000","400 000", "800 000", "9 M" )

topsCitiesNetwork <- ggplot() + 
  geom_sf(data = sfEU, fill = "#bfbfbf", color = "white", size = 0.5) +
  geom_sf(data = tops ,
          mapping = aes(  size = PopAdmin11, colour = Network), alpha = 0.8, show.legend = NA) +
  scale_color_manual(values=c("#ff620080","#3384b0", "#812aa3"))+
  geom_sf(data = tops  ,
          mapping = aes(size = PopAdmin11), shape = 1, colour = "grey80", show.legend = NA)+
  scale_size(name = "Population Administrative 2011",
             breaks = bks,
             labels = lbs,
             range = c(2, 11)) +
  annotate("text", label = "Source : ETMUN 2019, EUCICOP 2019/ PG. 2020",
           size = 2.2, 
           hjust = 1,
           x = c(st_bbox(rec)[3]), y = c(st_bbox(rec)[2]-130000)) +
  labs(x = "", y = "") +
  geom_sf_text(data = tops, aes(label = asciiName), size = 2.2, color = "#4d4d4d",
               check_overlap = FALSE, nudge_x = +100000, nudge_y = -50000) +
  geom_line(data = myScaleBar, aes(x = X, y = Y), size = 0.5, color = "#333333") +
  annotate("text", label = "500 km", size = 2.5, color = "#333333", hjust = 0,
           x = c(st_bbox(rec)[3]-800000), y = c(st_bbox(rec)[2]+280000)) +
  geom_sf(data = rec, fill = NA, color = "ivory4", size = 0.5) +
  coord_sf(crs = 3035, datum = NA,
           xlim = st_bbox(rec)[c(1,3)],
           ylim = st_bbox(rec)[c(2,4)]) + labs(color = "Appartenance aux réseaux") +
  theme_void() +
  theme(legend.position =  c(0.18, 0.60), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 7.5))



## display end save
pdf(file = "OUT/top50Map_centralities_EUCICOP_ETMUN.pdf", width = 8.3, height = 5.8, pagecentre = FALSE)
topsCitiesNetwork
dev.off()

ggsave( file="OUT/top50Map_centralities_EUCICOP_ETMUN.svg", plot=topsCitiesNetwork, width=8.3, height=5.8, device = "svg")

