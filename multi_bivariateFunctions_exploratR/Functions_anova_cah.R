
# compute classification ----
ComputeClassif <- function(df, varquanti, stand, method){
  classifObj <- agnes(x = df[, varquanti], diss = FALSE, metric = "euclidean", 
                      stand = stand, method = method)
  return(classifObj)
}

# plot dendrogram ----
PlotDendro <- function(classifobj){
  dendroPlot <- as.dendrogram(classifobj)
  dendroData <- dendro_data(dendroPlot, type = "rectangle")
  dendroGgplot <- ggplot(segment(dendroData)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous("") + scale_y_continuous("") +
    theme_bw()
  
  return(dendroGgplot)
}

# plot inertia ----
PlotHeight <- function(classifobj){
  sortedHeight <- sort(classifobj$height, decreasing = TRUE)
  relHeigth <- sortedHeight / sum(sortedHeight) * 100
  tabHeight <- data.frame(NODES = factor(1:20),
                          INERTIE = relHeigth[1:20])
  
  heightPlot <- ggplot(tabHeight) +
    geom_bar(aes(x = NODES, y = INERTIE), fill = "grey30", stat = "identity") +
    scale_x_discrete("Nombre de classes") + scale_y_continuous("Niveau") +
    theme_bw()
  
  return(heightPlot)
}

# plot profile ----
PlotProfile <- function(classifobj, nbclus){
  dfOri <- as.data.frame(classifobj$data, stringsAsFactors = FALSE)
  clusId <- cutree(classifobj, k = nbclus)
  dfOri$CLUS <- factor(clusId,
                       levels = 1:nbclus,
                       labels = paste("CLASSE", 1:nbclus))
  clusProfile <- aggregate(dfOri[, 1:ncol(dfOri)-1],
                           by = list(dfOri$CLUS),
                           mean)
  colnames(clusProfile)[1] <- "CLASSE"
  clusLong <- melt(clusProfile, id.vars = "CLASSE")
  
  profilePlot <- ggplot(clusLong) +
    geom_bar(aes(x = variable, y = value), fill = "grey30", position = "identity", stat = "identity") +
    scale_x_discrete("Variable") + scale_y_continuous("Valeur moyenne par classe") +
    facet_wrap(~ CLASSE) + coord_flip() + theme_bw()
  
  return(list(PROFILE = profilePlot, CLUSID = dfOri$CLUS))
}




# ---- ANOVA functions ----
# Commenges H. (2016) ExploratR : outil interactif d'exploration statistique uni- bi- tri- et multi-variée avec R, UMR 8504 Géographie-cités. 
# APPLI : https://analytics.huma-num.fr/geographie-cites/ExploratR/
# CODE : https://zenodo.org/record/155333#.XdZn7dVCfIU

## Anova parameters (1 factor) --
AnovaTab <- function(df, varx, vary){
  groupMean <- round(tapply(df[, vary], df[, varx], mean, na.rm = TRUE), digits = 2)
  groupMedian <- round(tapply(df[, vary], df[, varx], median, na.rm = TRUE), digits = 2)
  groupVar <- round(tapply(df[, vary], df[, varx], var, na.rm = TRUE), digits = 2)
  tabGroup <- data.frame(Modalité = names(groupMean), 
                         Moyenne = groupMean,
                         Médiane = groupMedian,
                         Variance = groupVar, 
                         stringsAsFactors = FALSE)
  tabAll <- data.frame(Modalité = "Ensemble", 
                       Moyenne = round(mean(df[, vary]), digits = 2), 
                       Médiane = round(median(df[, vary]), digits = 2), 
                       Variance = round(var(df[, vary]), digits = 2), 
                       stringsAsFactors = FALSE)
  
  tabVariance <- rbind(tabGroup, tabAll)
  
  return(tabVariance)
}

## Anova plot (1 factor) -- labels not generalized
AnovaPlot <- function(df, varx, vary, tx, ty, source){
  
  xLevels <- sort(unique(df[, varx]))
  df$ID <- df[, varx]
  df$VAR <- df[, vary]
  
  if(length(xLevels) == 2){
    colPal <- brewer.pal(n = 3, name = "Set1")[1:2]
  } else if (length(xLevels) > 2){
    colPal <- brewer.pal(n = length(xLevels), name = "Set1")
  }
  
  # jitter points
  set.seed(99)
  df$JIT <- as.numeric(as.factor(df[, varx])) + sample(x = seq(-0.3, 0.3, 0.01), size = nrow(df), replace = TRUE)
  
  # mean segments
  groupMean <- tapply(df[, vary], df[, varx], mean, na.rm = TRUE)
  avgSegment <- data_frame(ID = names(groupMean), 
                           XMIN = seq(1, length(groupMean), 1) - 0.4,  
                           XMAX = seq(1, length(groupMean), 1) + 0.4, 
                           YMIN = groupMean, 
                           YMAX = groupMean)
  
  # residuals segments
  df <- df %>% left_join(x = ., y = avgSegment[, c(1, 4)], by = "ID")
  
  aovPlot <- ggplot() +
    geom_hline(yintercept = mean(df$VAR, na.rm = TRUE), color = "grey60", size = 1, linetype = 2) +
    geom_segment(data = avgSegment, aes(x = XMIN, xend = XMAX, y = YMIN, yend = YMAX), color = "grey40", size = 2) +
    geom_segment(data = df, aes(x = JIT, xend = JIT, y = YMIN, yend = VAR), color = "grey40", alpha = 0.5) +
    geom_point(data = df, aes(JIT, VAR, color = ID), show.legend = FALSE) +
    scale_color_manual(values = colPal) +
    scale_x_continuous(name = tx, breaks = seq(1, length(groupMean), 1), labels = xLevels) +
    scale_y_continuous(name = ty) +
    labs(x = tx, y = ty) +
    theme_bw() +
    labs(caption = source) +
    theme(axis.text.x = element_text(size = 9, angle = 20, vjust = 0.6),
          plot.caption = element_text(size = 6))
  
  return(aovPlot)
}

### Compute linear model -- 
ComputeRegression <- function(df, vardep, varindep, interact = FALSE){
  if(interact == FALSE){
    linMod <- lm(formula = formula(eval(paste(vardep, "~", paste(varindep, collapse = "+")))), data = df)
    linModSumry <- summary(linMod)
  } else {
    linMod <- lm(formula = formula(eval(paste(vardep, "~", paste(varindep, collapse = "*")))), data = df)
    linModSumry <- summary(linMod)
  }
  coefReg <- round(linModSumry$coefficients, digits = 4)[, 1:2]
  rawR2 <- round(linModSumry$r.squared, digits = 2)
  adjR2 <- round(linModSumry$adj.r.squared, digits = 2)
  
  tabResid <- data.frame(ABSRESID = round(linModSumry$residuals, digits = 3), 
                         RELRESID = round(linModSumry$residuals / (df[, vardep] - linModSumry$residuals), digits = 3))
  
  tabResults <- data.frame(CONCEPT = c("Coef. de détermination",
                                       "Coef. de détermination ajusté",
                                       row.names(coefReg)),
                           VALEUR = c(rawR2, adjR2, coefReg[, 1]),
                           stringsAsFactors = FALSE)
  
  return(list(TABCOEF = tabResults, TABRESID = tabResid, COEF = coefReg))
}
