Packages <- c(
  "ape", "lme4", "ggtern", "nlme", "lmerTest", "openxlsx", 
  "edgeR", "ggplot2", "vegan", "tidyverse", "ggalt", 
  "patchwork", "NST", "ggpmisc", "reshape2", "ggprism", 
  "ggpubr", "ggalluvial", "dplyr", "ggh4x", "car", "NST",
  "sjPlot", "performance", "MuMIn", "mgcv", "emmeans", 
  "multcomp", "stringr", "gamm4", "Matrix", "ggthemes", 
  "tidyr", "glmm.hp", "partR2", "corrplot", "gratia", 
  "ggstar", "ggmagnify", "ggbreak", "picante", "RColorBrewer", 
  "ggforce", "grid", "gridExtra", "linkET", "cowplot"
)
#install.packages(setdiff(Packages, rownames(installed.packages())))
lapply(Packages, library, character.only = TRUE)

setwd("D:/桌面/ASV_analysis/")

###########Common functions################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(length(x[[col]]))))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
##########################Tables######################
###alpha diversity in field survey
set.seed(12345)
diversity <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#diversity$plot <- paste0(diversity$species,"_",diversity$elevation)
diversity$month <- factor(diversity$month,levels = c("Sep","Oct","Dec","Jun","Aug"))

###model of time & elevation
model <- lmer(log(ASVs.of.total) ~ compartment * log(elevation) * month + (1|pair)+(1|species), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 2)
ranova(model)

model <- lmer(ASVs ~ compartment * log(elevation) * month + (1|pair)+(1|species), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 2)
ranova(model)

model <- lmer(read.of.shared ~ compartment * log(elevation) * month + (1|pair)+(1|species), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 2)
ranova(model)

####model of time ////elevation
model <- lmer(log(ASVs.of.total) ~ compartment * species * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(ASVs ~ species * compartment * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(read.of.shared ~ species * compartment * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

###beta diversity in field survey
###PERMANOVA ANALYSIS
total <- readRDS("field_wunifrac.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#group$plot <- paste0(group$species,"_",group$elevation)
with(group, adonis2(total ~ pair + month * species * compartment, data = group, permutations = 999, strata = elevation, by = "term"))

with(group, adonis2(total ~ pair + month * elevation * compartment, data = group, permutations = 999, strata = species, by = "term"))

with(group, adonis2(total ~ pair + elevation * species * compartment, data = group, permutations = 999, strata = month, by = "term"))

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = elevation, by = "term"))
##leaf & environment
group1 <- subset(group, compartment == "L")
data <- total[rownames(total) %in% rownames(group1),]
wunifrac <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(wunifrac ~ species+pH+AT+AM+STN+LTN+SLA+Chl+
                       species:pH+species:AT+species:AM+species:STN+species:LTN+species:SLA+
                       species:Chl,data = group1, permutations = 999,strata = elevation, by = "term"))

##root & environment
group1 <- subset(group, compartment == "R")
data <- total[rownames(total) %in% rownames(group1),]
wunifrac <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(wunifrac ~ species+pH+AT+AM+STN+LTN+SLA+Chl+
                       species:pH+species:AT+species:AM+species:STN+species:LTN+species:SLA+
                       species:Chl,data = group1, permutations = 999,strata = elevation, by = "term"))

###shared & unique analysis
total <- readRDS("field_shared_wunifrac.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#group$plot <- paste0(group$species,"_",group$elevation)

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = elevation, by = "term"))

total <- readRDS("field_unique_wunifrac.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#group$plot <- paste0(group$species,"_",group$elevation)

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = elevation, by = "term"))

############################################################################
#################################graph######################################
############################################################################
#######################Fig1####################
diversity <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

diversity$logASVs.of.total <- log(diversity$ASVs.of.total)
df <- data_summary(diversity, varname = "ASVs", groupnames = c("compartment"))

create_plot <- function(data, response_var, y_label) {
  df <- data_summary(data, varname = response_var, groupnames = c("compartment", "month", "species"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  
  ggplot() +
    geom_line(data = df, mapping = aes_string(x = "month", y = response_var, group = "compartment"), color = "black", size = 0.5) +
    geom_errorbar(data = df, mapping = aes_string(x = "month", y = response_var, ymin = paste0(response_var, " - se"), ymax = paste0(response_var, " + se"), color = "compartment"), width = 0.1, size = 0.5) +
    geom_jitter(data = data, mapping = aes_string(x = "month", y = response_var, shape = "compartment", color = "compartment"), size = 1, width = 0.03, alpha = 0.2) +
    geom_point(data = df, mapping = aes_string(x = "month", y = response_var, color = "compartment", shape = "compartment"), size = 2, width = 0.3) +
    scale_shape_manual(values = c("L" = 16, "R" = 16)) +
    scale_color_manual(values = c("L" = "#108b96", "R" = "#7A8234")) +
    labs(y = y_label) +
    theme_classic() +
    theme(
      axis.title = element_text(color = 'black', size = 12),
      axis.title.x = element_text(colour = 'black', size = 12, vjust = 0),
      axis.title.y = element_text(colour = 'black', size = 12, vjust = 0),
      axis.text = element_text(colour = 'black', size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key = element_blank(),
      legend.position = "none",
      legend.background = element_rect(colour = "black")
    ) +
    facet_wrap(~species) +
    theme(strip.text = element_blank())
}

p1a <- create_plot(diversity, "logASVs.of.total", "No. of total ASVs (log)")
p1b <- create_plot(diversity[-200,], "ASVs", "% ASVs shared")
p1c <- create_plot(diversity, "read.of.shared", "RA of shared ASVs (%)")

layout <- "
AAAAAA
BBBBBB
CCCCCC
"
p1a + p1b + p1c + 
  plot_layout(design = layout)

###Multiple comparisons (based on mixed models)
set.seed(12345)
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
group$ASVs.of.total <- log10(group$ASVs.of.total)
diversity$month <- factor(diversity$month,levels = c("Sep","Oct","Dec","Jun","Aug"))
group$plot <- paste0(group$species,"_",group$elevation)
##ASVs.of.total    |    ASVs     |    read.of.shared
###emmeans packages
###ASVs.of.total
model <- lmer(log(ASVs.of.total) ~ species * compartment * month + (1|pair)+(1|plot), 
              na.action= na.omit, data = group)
car::Anova(model, type = 3)
ranova(model)

emm1 = emmeans(model, specs = pairwise ~ compartment * month|species, type = 'response', adjust = 'tukey')
emm1_multi = multcomp::cld(emm1,alpha=0.05,Letters=letters,adjust="none",decreasing = T)
emm1_multi$.group <- trimws(emm1_multi$.group)
emm1_multi

###ASVs
model <- lmer(ASVs ~ species * compartment * month + (1|pair)+(1|plot), 
              na.action= na.omit, data = group)
car::Anova(model, type = 3)
ranova(model)

emm1 = emmeans(model, specs = pairwise ~ compartment * month|species, type = 'response', adjust = 'tukey')
emm1_multi = multcomp::cld(emm1,alpha=0.05,Letters=letters,adjust="none",decreasing = T)
emm1_multi$.group <- trimws(emm1_multi$.group)
emm1_multi
###read.of.shared
model <- lmer(read.of.shared ~ species * compartment * month + (1|pair)+(1|plot), 
              na.action= na.omit, data = group)
car::Anova(model, type = 3)
ranova(model)

emm1 = emmeans(model, specs = pairwise ~ compartment * month|species, type = 'response', adjust = 'tukey')
emm1_multi = multcomp::cld(emm1,alpha=0.05,Letters=letters,adjust="none",decreasing = T)
emm1_multi$.group <- trimws(emm1_multi$.group)
emm1_multi

################Fig2##############
entire_dis <- readRDS("field_wunifrac.rds")
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
entire_dis <- as.matrix(entire_dis)

perform_pcoa_and_plot <- function(data, group_data, compartment1) {
  compartment_data <- subset(group_data, compartment == compartment1)
  filtered_data <- data[rownames(data) %in% rownames(compartment_data), ]
  leaf_data <- filtered_data[, colnames(filtered_data) %in% rownames(compartment_data)]
  
  pcoa_result <- cmdscale(leaf_data, k = 2, eig = TRUE)
  points_df <- as.data.frame(pcoa_result$points)
  names(points_df)[1:2] <- c("PCoA1", "PCoA2")
  points_df$SampleID <- rownames(points_df)
  
  group_data$SampleID <- rownames(group_data)
  sample_site <- merge(points_df, group_data)
  
  df <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  
  colnames(df)[2:3] <- c("PCoA1_mean", "PCoA2_mean")
  
  sample_site$type <- paste0(sample_site$month)
  df$type <- paste0(df$month)
  sample_site2 <- sample_site %>% left_join(df[, c(2:4)], by = "type")
  
  df_line <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  colnames(df_line)[2:3] <- c("PCoA1_line", "PCoA2_line")
  
  sample_site2$month <- factor(sample_site2$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line$month <- factor(df_line$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line <- df_line %>% arrange(month)

  pcoa_eig <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig), 3)
  
  ggplot() +
    geom_point(sample_site2, mapping = aes(PCoA1, PCoA2, color = month), shape = 21, size = 1.5, alpha = 0.8) +
    geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = month), size = 2.5, alpha = 1, shape = 21) +
    geom_point(df_line, mapping = aes(PCoA1_line, PCoA2_line), color = NA, size = 2.5) +
    scale_color_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    scale_fill_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    geom_encircle(data = sample_site2, mapping = aes(PCoA1, PCoA2, group = month, fill = month), expand = 0, spread = 0.5,
                  s_shape = 1, size = 1, linetype = 1, alpha = 0.1) +
    geom_curve(data = df_line, mapping = aes(x = lag(PCoA1_line), y = lag(PCoA2_line), 
                                             xend = PCoA1_line, yend = PCoA2_line),
               curvature = 0.3, arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    labs(x = paste("PCoA1 (", format(100 * pcoa_eig[1], digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * pcoa_eig[2], digits = 3), "%)", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.001)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.001))
}
p2a <- perform_pcoa_and_plot(entire_dis, group, "L")
p2b <- perform_pcoa_and_plot(entire_dis, group, "R")
p2a + p2b

####p2c
entire_dis <- readRDS("field_wunifrac.rds")
data <- as.matrix(entire_dis)
merge_and_filter <- function(data, group, columns) {
  group1 <- group[, columns]
  group1$group <- rownames(group)
  
  data_3col <- dist.3col(data) %>%
    setNames(c("group", "group2", "distance")) %>%
    left_join(group1, by = "group")
  
  group2 <- group[, columns]
  group2$group2 <- rownames(group)
  
  data_joined <- left_join(data_3col, group2, by = "group2")
  data_joined$species <- ifelse(data_joined$species.x == data_joined$species.y, data_joined$species.x, 0)
  data_filtered <- subset(data_joined, species != 0)
  
  return(data_filtered)
}

columns <- c(2, 3, 5, 6, 7)
data_filtered <- merge_and_filter(data, group, columns)

data_filtered$elevation <- ifelse(data_filtered$elevation.x == data_filtered$elevation.y, data_filtered$elevation.x, 0)
data_filtered <- subset(data_filtered, elevation != 0)

data_filtered$compartment <- ifelse(data_filtered$compartment.x == data_filtered$compartment.y, data_filtered$compartment.x, 0)
data_filtered <- subset(data_filtered, compartment != 0)

data_filtered$time <- abs(data_filtered$time.x - data_filtered$time.y)

leaf <- subset(data_filtered, compartment == "L" & species == "F")
root <- subset(data_filtered, compartment == "R" & species == "F")
###gamm model
fit_leaf <- gamm(distance ~ s(time, bs = 'cr', k = 3), random = list(elevation.x = ~1), data = leaf, family = gaussian, method = "REML")
fit_root <- gamm(distance ~ s(time, bs = 'cr', k = 3), random = list(elevation.x = ~1), data = root, family = gaussian, method = "REML")
summary(lm(leaf$time~leaf$distance))
summary(lm(root$time~root$distance))

draw(fit_leaf$gam) + theme_bw()
summary(fit_leaf$gam)

draw(fit_root$gam) + theme_bw()
summary(fit_root$gam)

###p1/p2/p3
p1 <- ggplot() +
  #geom_jitter(data = root, 
  #           aes(x = time, y = distance), width = 0.2,
  #          height = 0, size = 1, shape = 21, 
  #         stroke = 0.3, color = "#7A8234", alpha = 0.7,
  #        show.legend = FALSE) +
  #geom_jitter(data = leaf, 
  #           aes(x = time, y = distance), width = 0.2,
  #          height = 0, size = 1, shape = 21, 
  #         stroke = 0.3, color = "#108b96", alpha = 0.7,
  #        show.legend = FALSE) +
  stat_smooth(data = leaf,
              aes(x = time, y = distance),
              method = "gam", formula = y ~ s(x, k = 3),se = TRUE, linewidth = 0.5, color = "#108b96", fill = "grey70") +
  stat_smooth(data = root,
              aes(x = time, y = distance),
              method = "gam", formula = y ~ s(x, k = 3),se = TRUE, linewidth = 0.5, color = "#7A8234", fill = "grey70") +
  labs(x = "Time interval (month)", y = "Community dissimilarity \n (weighted Unifrac distance)") +
  theme_bw() +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.001)) +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 12));p3
p1

###AT & Chl
entire_dis <- readRDS("field_wunifrac.rds")
data <- as.matrix(entire_dis)
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
colnames(group)

##functions
process_data <- function(data, group, chl_col_index, variable_name) {
  variable_col_name <- colnames(group)[chl_col_index]
  
  group1 <- group[, c(3, 5, 6, 7, chl_col_index)] %>%
    mutate(group = rownames(group))
  
  data_3col <- dist.3col(data) %>%
    setNames(c("group", "group2", "distance")) %>%
    as.data.frame() %>%
    left_join(group1, by = "group")
  
  group2 <- group[, c(3, 5, 6, 7, chl_col_index)] %>%
    mutate(group2 = rownames(group))
  
  data3 <- left_join(data_3col, group2, by = "group2") %>%
    mutate(
      position = ifelse(compartment.x == compartment.y, compartment.x, 0),
      elevation = ifelse(elevation.x == elevation.y, elevation.x, 0)
    ) %>%
    filter(position != 0, elevation != 0)
  return(data3)
}

data5_temp <- process_data(data, group, 21, "AT")
data5_chl <- process_data(data, group, 24, "Chl")
data5_temp$AT <- abs(data5_temp$AT.x - data5_temp$AT.y)
data5_chl$Chl <- abs(data5_chl$Chl.x - data5_chl$Chl.y)

leaf_temp <- subset(data5_temp, position == "L")
root_temp <- subset(data5_temp, position == "R")

leaf_chl <- subset(data5_chl, position == "L")
root_chl <- subset(data5_chl, position == "R")

# correlation test
cor.test(leaf_temp$AT, leaf_temp$distance)
cor.test(root_temp$AT, root_temp$distance)
cor.test(leaf_chl$Chl, leaf_chl$distance)
cor.test(root_chl$Chl, root_chl$distance)

##graph function
create_plot <- function(data_filtered, x_col, y_col, x_label, y_label, color1, color2) {
  ggplot() +
    geom_jitter(data = filter(data_filtered, position == 'R'), 
                aes_string(x = x_col, y = y_col), width = 0.1,
                height = 0, size = 0.8, shape = 21, 
                stroke = 0.3, color = color2, alpha = 0.2,
                show.legend = FALSE) +
    geom_jitter(data = filter(data_filtered, position == 'L'), 
                aes_string(x = x_col, y = y_col), width = 0.1,
                height = 0, size = 0.8, shape = 21, 
                stroke = 0.3, color = color1, alpha = 0.2,
                show.legend = FALSE) +
    stat_smooth(data = filter(data_filtered, position == 'L'),
                aes_string(x = x_col, y = y_col),
                method = "lm", se = TRUE, linewidth = 0.5, color = color1, fill = "grey70") +
    stat_smooth(data = filter(data_filtered, position == 'R'),
                aes_string(x = x_col, y = y_col),
                method = "lm", se = TRUE, linewidth = 0.5, color = color2, fill = "grey70") +
    labs(x = x_label, y = y_label) +
    theme_bw() +
    #scale_y_continuous(labels = scales::label_comma(accuracy = 0.001)) +
    theme(
      panel.grid.major = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(color = "black", size = 10),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(color = "black", size = 12)
    )
}
p2d <- create_plot(data5_temp, "AT", "distance", "Difference in air temperature", "Compositional dissimilarity", "#108b96", "#7A8234")
p2e <- create_plot(data5_chl, "Chl", "distance", "Difference in Chl", "Compositional dissimilarity", "#108b96", "#7A8234")
p2d+p2e

##ancova analysis
mod <- lm(distance ~ AT * position, data= data5_temp)
summary(mod)
mod <- lm(distance ~ Chl * position, data= data5_chl)
summary(mod)

####dbRDA
total <- readRDS("field_wunifrac.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

# Define time points
time_points <- unique(group$month)
time_points
compartments <- unique(group$compartment)
compartments

# Initialize result storage
r2_results <- data.frame(
  time = character(),
  compartment = character(),
  Type = character(),
  normal = numeric(),
  adj = numeric(),
  stringsAsFactors = FALSE
)


for (i in time_points) {
  for (j in compartments) {
    # Filter the environmental data at the corresponding time point
    env_time <- subset(group, month == i | compartment == j)
    env1 <- data.frame(scale(env_time[c(16:18, 21:24)]))
    # Make sure the row names of the environment data and distance data match
    data <- total[rownames(total) %in% rownames(env_time), ]
    wunifrac <- data[, colnames(data) %in% rownames(env_time)]
    
    if (nrow(env_time) > 1 && ncol(wunifrac) > 1) {
      mod <- dbrda(wunifrac ~ ., env1)
      vif.cca(mod)  # Screening for collinearity
      anova_result <- data.frame(anova(mod, by = "term", permutations = 999))
      anova_result$Term <- row.names(anova_result)
      # Screen out significant environmental factors (p < 0.05)
      env_factors <- anova_result %>%
        filter(Pr..F. < 0.05, Term %in% c("pH", "AT", "AM")) %>%
        pull(Term)
      
      plant_factors <- anova_result %>%
        filter(Pr..F. < 0.05, Term %in% c("STN", "LTN", "SLA")) %>%
        pull(Term)
      
      if (length(env_factors) > 0) {
        env_formula <- as.formula(paste("wunifrac ~", paste(env_factors, collapse = " + ")))
        model_env <- dbrda(env_formula, data = env1)
        r2_env <- RsquareAdj(dbrda(env_formula, env1))
        r2_results <- rbind(r2_results, data.frame(
          time = i,
          compartment = j,
          Type = "env",
          normal = r2_env$r.squared,
          adj = r2_env$adj.r.squared,
          stringsAsFactors = FALSE
        ))
      }
      
      if (length(plant_factors) > 0) {
        plant_formula <- as.formula(paste("wunifrac ~", paste(plant_factors, collapse = " + ")))
        model_plant <- dbrda(plant_formula, data = env1)
        r2_plant <- RsquareAdj(dbrda(plant_formula, env1))
        r2_results <- rbind(r2_results, data.frame(
          time = i,
          compartment = j,
          Type = "plant",
          normal = r2_plant$r.squared,
          adj = r2_plant$adj.r.squared,
          stringsAsFactors = FALSE
        ))
      }
    } 
  }
}
r2_results
#write.table(r2_results, file ="result.csv",sep =",", quote =FALSE)
####################FigS3#######################
entire_dis <- readRDS("field_wunifrac.rds")
entire_dis <- data.frame(as.matrix(entire_dis))
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
leaf <- subset(group, compartment == "L")
soil <- subset(group, compartment == "R")
data <- entire_dis[rownames(entire_dis) %in% rownames(leaf),]
rhizos1 <- data[, colnames(data) %in% rownames(soil)]
diagonal_values <- diag(as.matrix(rhizos1))  
row_names_diag <- rownames(rhizos1) 
col_names_diag <- colnames(rhizos1)  
result <- data.frame(ID = row_names_diag, id = col_names_diag, Value = diagonal_values)
result1 <- left_join(result, leaf[c(1,3,5,6)],by = "ID")
Fdata <- subset(result1, species == "T") 
fit<-gamm(formula(Value~s(time,bs='cr', k = 5)), 
          random=list(elevation=~1),data=Fdata,
          family=gaussian)
summary(fit$gam)
draw(fit$gam)+theme_bw()
##p1/p2/p3
p3 <- ggplot() +
  geom_jitter(data = Fdata, 
              aes(x = time, y = Value), width = 0.2,
              height = 0, size = 2, shape = 21, 
              stroke = 0.3, fill = "black", alpha = 0.3,
              show.legend = FALSE) +
  stat_smooth(data = Fdata,
              aes(x = time, y = Value),
              method = "gam", formula = y ~ s(x, k = 5),se = TRUE, linewidth = 0.5, color = "black", fill = "grey70") +
  theme_bw() +
  labs(y = "Community dissimilarity\n(Paired phyllosphere and rhizosphere)") +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 12))
p1+p2+p3

####################FigS4###############
#Take field survey as an example
total <- read.xlsx("field_survey.xlsx", sheet = "shared_ASV", colNames = TRUE, rowNames = TRUE)
field_group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

#Function: Perform EdgeR analysis
perform_edgeR_analysis <- function(species_name, total_data, field_group_data) {
  group1 <- subset(field_group_data, species == species_name)
  total_core_otu <- total_data[, colnames(total_data) %in% rownames(group1),]
  total_core_otu <- total_core_otu[rowSums(total_core_otu) != 0, ]
  
  group <- factor(group1$compartment, levels = c("L", "R"))
  library_sizes <- colSums(total_core_otu)
  total_core_otu <- total_core_otu[, library_sizes > 0]
  group <- group[library_sizes > 0]  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  
  design <- model.matrix(~group)
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  
  edgeR_diff <- lrt$table[order(lrt$table$FDR, lrt$table$logFC, decreasing = c(FALSE, TRUE)), ]
  edgeR_diff[which(edgeR_diff$logFC >= 1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in root', sep = '')
  edgeR_diff[which(edgeR_diff$logFC <= -1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in leaf', sep = '')
  edgeR_diff[which(abs(edgeR_diff$logFC) <= 1 | edgeR_diff$FDR >= 0.05), 'sig'] <- 'Non-significant'
  
  return(edgeR_diff)
}

# Function: Create a volcano plot
create_volcano_plot <- function(edgeR_diff) {
  edgeR_diff$sig <- as.factor(edgeR_diff$sig)
  
  p <- ggplot(edgeR_diff, aes(logFC, -log(FDR, 10), color = sig)) +
    geom_point(size = 1.5, pch = 21) + 
    scale_color_manual(values = c('Rich in leaf' = "#108b96", 'Non-significant' ='gray40', 'Rich in root' = "#7A8234")) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 12), 
          axis.ticks.y = element_line(), axis.line.y = element_line(), axis.line.x = element_line(color = 'black'),
          axis.text.y = element_text(size = 10, angle = 0, colour = "black"),
          axis.text.x = element_text(size = 10, angle = 0, colour = "black"),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.position = "none",
          legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), color = 'gray', linewidth = 0.25, linetype = 2) + 
    geom_hline(yintercept = -log(0.05, 10), color = 'gray', linewidth = 0.25, linetype = 2) +
    labs(x = expression(Log[2]~(Fold~change)), y = expression(-Log[10]~(P-value)), color = NA)
  
  return(p)
}

# Perform analysis and create graphs
species_list <- c("F", "P", "T")
plots <- list()

for (species in species_list) {
  edgeR_diff <- perform_edgeR_analysis(species, total, field_group)
  plots[[species]] <- create_volcano_plot(edgeR_diff)
  
  filter_up <- subset(edgeR_diff, FDR < 0.05 & logFC >= 1)
  filter_down <- subset(edgeR_diff, FDR < 0.05 & logFC <= -1)
  print(paste('Rich in root for', species, ':', nrow(filter_up)))
  print(paste('Rich in leaf for', species, ':', nrow(filter_down)))
  print(paste('Total for', species, ':', nrow(edgeR_diff)))
}

p4 <- plots[["F"]]
p5 <- plots[["P"]]
p6 <- plots[["T"]]

(p4+p5+p6)

#############FigS5#############
entire_dis <- readRDS("field_unique_wunifrac.rds")
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
entire_dis <- as.matrix(entire_dis)
perform_pcoa_and_plot <- function(data, group_data, compartment1) {
  
  compartment_data <- subset(group_data, compartment == compartment1)
  filtered_data <- data[rownames(data) %in% rownames(compartment_data), ]
  leaf_data <- filtered_data[, colnames(filtered_data) %in% rownames(compartment_data)]
  
  pcoa_result <- cmdscale(leaf_data, k = 2, eig = TRUE)
  points_df <- as.data.frame(pcoa_result$points)
  names(points_df)[1:2] <- c("PCoA1", "PCoA2")
  points_df$SampleID <- rownames(points_df)
  
  
  group_data$SampleID <- rownames(group_data)
  sample_site <- merge(points_df, group_data)
  
  df <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  
  colnames(df)[2:3] <- c("PCoA1_mean", "PCoA2_mean")
  
  sample_site$type <- paste0(sample_site$month)
  df$type <- paste0(df$month)
  sample_site2 <- sample_site %>% left_join(df[, c(2:4)], by = "type")
  
  df_line <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  colnames(df_line)[2:3] <- c("PCoA1_line", "PCoA2_line")
  
  sample_site2$month <- factor(sample_site2$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line$month <- factor(df_line$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line <- df_line %>% arrange(month)

  pcoa_eig <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig), 3)
  
  ggplot() +
    geom_point(sample_site2, mapping = aes(PCoA1, PCoA2, color = month), shape = 21, size = 1.5, alpha = 0.8) +
    geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = month), size = 2.5, alpha = 1, shape = 21) +
    geom_point(df_line, mapping = aes(PCoA1_line, PCoA2_line), color = NA, size = 2.5) +
    scale_color_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    scale_fill_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    geom_encircle(data = sample_site2, mapping = aes(PCoA1, PCoA2, group = month, fill = month), expand = 0, spread = 0.5,
                  s_shape = 1, size = 1, linetype = 1, alpha = 0.1) +
    geom_curve(data = df_line, mapping = aes(x = lag(PCoA1_line), y = lag(PCoA2_line), 
                                             xend = PCoA1_line, yend = PCoA2_line),
               curvature = 0.3, arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    labs(x = paste("PCoA1 (", format(100 * pcoa_eig[1], digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * pcoa_eig[2], digits = 3), "%)", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.0001)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.0001))
}
p2a <- perform_pcoa_and_plot(entire_dis, group, "L")
p2b <- perform_pcoa_and_plot(entire_dis, group, "R")
p2c <- perform_pcoa_and_plot(entire_dis, group, "L")
p2d <- perform_pcoa_and_plot(entire_dis, group, "R")
(p2a + p2b)/(p2c + p2d)

