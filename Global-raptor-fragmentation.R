
library(ggplot2)
library(DescTools)
library(terra)
library(dplyr)
library(car)
library(multcompView)
library(agricolae)
library(mgcv)
library(showtext)
library(devEMF)
library(scales)
library(dunn.test)
library(rstatix)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(ggcorrplot)
library(randomForest)
library(caret)
library(sf)
library(stringr)
library(exactextractr)


####### FFI2000、FFI2020####

data2020 <- read.csv("G:/*/FFI-2020.csv")
data2000 <- read.csv( "G:/*/FFI-2000.csv")

get_limits <- function(x){
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  
  IQR_value <- Q3 - Q1
  upper <- Q3 + 1.5 * IQR_value
  lower <- Q1 - 1.5 * IQR_value
  
  return(list(
    Q1 = Q1,
    Q3 = Q3,
    IQR = IQR_value,
    upper = upper,
    lower = lower
  ))
}

clip_with_limits <- function(x, lower, upper){
  
  x[x > upper] <- upper
  x[x < lower] <- lower
  
  return(x)
}

MPA_limits <- get_limits(data2020$MPA)
ED_limits <- get_limits(data2020$ED)
PD_limits <- get_limits(data2020$PD)

data2020$MPA_clean <- clip_with_limits(data2020$MPA,MPA_limits$lower,MPA_limits$upper)
data2020$ED_clean <- clip_with_limits(data2020$ED,ED_limits$lower,ED_limits$upper)
data2020$PD_clean <- clip_with_limits(data2020$PD,PD_limits$lower,PD_limits$upper)

data2000$MPA_clean <- clip_with_limits(data2000$MPA,MPA_limits$lower,MPA_limits$upper)
data2000$ED_clean <- clip_with_limits(data2000$ED,ED_limits$lower,ED_limits$upper)
data2000$PD_clean <- clip_with_limits(data2000$PD,PD_limits$lower,PD_limits$upper)

MPA_min <- min(data2020$MPA_clean, na.rm = TRUE)
MPA_max <- max(data2020$MPA_clean, na.rm = TRUE)

ED_min <- min(data2020$ED_clean, na.rm = TRUE)
ED_max <- max(data2020$ED_clean, na.rm = TRUE)

PD_min <- min(data2020$PD_clean, na.rm = TRUE)
PD_max <- max(data2020$PD_clean, na.rm = TRUE)


cat("\nMPA:\n")
cat("min =", MPA_min, "\n")
cat("max =", MPA_max, "\n")

cat("\nED:\n")
cat("min =", ED_min, "\n")
cat("max =", ED_max, "\n")

cat("\nPD:\n")
cat("min =", PD_min, "\n")
cat("max =", PD_max, "\n")

normalize_fixed <- function(x, xmin, xmax){
  (x - xmin) / (xmax - xmin)
}

data2020$MPA_nor <- normalize_fixed(data2020$MPA_clean, MPA_min,MPA_max)
data2020$ED_nor <- normalize_fixed(data2020$ED_clean,ED_min,ED_max)
data2020$PD_nor <- normalize_fixed(data2020$PD_clean,PD_min,PD_max)

data2000$MPA_nor <- normalize_fixed(data2000$MPA_clean,MPA_min,MPA_max)
data2000$ED_nor <- normalize_fixed(data2000$ED_clean,ED_min,ED_max)
data2000$PD_nor <- normalize_fixed(data2000$PD_clean,PD_min,PD_max)

data2020$FFI <- (data2020$ED_nor +data2020$PD_nor +(1 - data2020$MPA_nor)) / 3
data2000$FFI <- (data2000$ED_nor +data2000$PD_nor +(1 - data2000$MPA_nor)) / 3


write.csv(data2020,
  "G:/*/fragmentation_2020.csv",
  row.names = FALSE)

write.csv(data2000,
  "G:/*/fragmentation_2000.csv",
  row.names = FALSE)


####T-test_Forest_Dependency####


data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)

t_model <- t.test(FFI_2020 ~ Forest_Dependency, data = data,var.equal = FALSE)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`FFI_2020` ~ Forest_Dependency, data=data, mean)
aggregate(`FFI_2020` ~ Forest_Dependency, data=data, sd)


boxplot(`FFI_2020` ~ Forest_Dependency,
        data=data,
        xlab="Forest dependency",
        ylab="Forest Fragmentation Index (FFI)",
        main="Fragmentation differences between forest dependency classes")

####T-test_Body_Size####

data <- read.csv("G:*/Species_information.csv")
data$Body_Size <- as.factor(data$Body_Size)

t_model <- t.test(FFI_2020 ~ Body_Size, data = data)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")

sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`FFI_2020` ~ Body_Size, data=data, mean)
aggregate(`FFI_2020` ~ Body_Size, data=data, sd)

boxplot(`FFI_2020` ~ Body_Size,
        data=data,
        xlab="Body size class",
        ylab="Forest Fragmentation Index (FFI)",
        main="Fragmentation differences among body-size classes")


##### GAM_Range_Size_vs_FFI####

data <- read.csv("G:/*/Species_information.csv")

model_gam <- gam(FFI_2020 ~ s(Species_Area, k = 10), data = data)
summary(model_gam)

new_data <- data.frame(Species_Area = seq(min(data$Species_Area), max(data$Species_Area), length.out = 500))
pred <- predict(model_gam, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

new_data <- new_data %>%
  mutate(
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 

sci_label <- function(x) {
  lab <- scientific_format()(x)
  lab <- gsub("e\\+", "e", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
windowsFonts(Arial = windowsFont("Arial"))
ggplot(new_data, aes(x = Species_Area, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#fdbb84",
              alpha = 0.5) +
  geom_line(color = "#e34a33",
            linewidth = 1.2) +labs(x = "Species Area",
                                   y = "2020 FFI") +
  coord_cartesian(ylim = c(0.4, 1))+
  scale_x_continuous(
    labels = sci_label,
    breaks = c(0, 1e7, 2e7,3e7,4e7)
  )+
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line = element_blank(), 
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    panel.background = element_blank(),
    plot.background = element_blank()
  )


summary(model_gam)

ggsave("GAM.pdata",width = 6,height = 6,device = cairo_pdata)
getwd()

##### Boxplot_IUCN_vs_FFI####

data <- read.csv("G:/*/Species_information.csv")

data$IUCN_status <- factor(data$IUCN_status,
                         levels = c("LC","NT","VU","EN","CR"))

data$IUCN_status <- factor(data$IUCN_status,
                         levels = c("LC","NT","VU","EN","CR"))
anova_model <- aov(FFI_2020 ~ IUCN_status, data = data)
summary(anova_model)
anova_result <- summary(anova_model)

tukey <- TukeyHSD(anova_model)
tukey

tukey_letters <- multcompLetters4(anova_model, tukey)

letters_data <- data.frame(
  IUCN_status = names(tukey_letters$IUCN_status$Letters),
  Letters = tukey_letters$IUCN_status$Letters
)

letters_data


F_value <- anova_result[[1]]$`F value`[1]
p_value <- anova_result[[1]]$`Pr(>F)`[1]

cat("F value =", F_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(FFI_2020 ~ IUCN_status, data=data, mean)
aggregate(FFI_2020 ~ IUCN_status, data=data, sd)


font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 
par(family = "Arial")

axis_font_size <- 1.3  
title_font_size <- 1.4 
label_font_size <- 1.4 

boxplot(FFI_2020 ~ IUCN_status,
        data = data,
        xlab = "IUCN status",
        ylab = "FFI_2020",
        main = "Fragmentation differences among IUCN status classes",
        ylim = c(0, 1.1),                 
        yaxt = "n",                     
        family = "Arial",
        cex.axis = axis_font_size,       
        cex.lab = label_font_size,       
        cex.main = title_font_size      
)

axis(2,
     at = pretty(c(0, 1.3), n = 4),
     las = 1,
     cex.axis = axis_font_size)


ggsave("GAM.pdata",width = 6,height = 6,device = cairo_pdata)


##### Correlation_Forest_Loss_vs_FFI ####
data <- read.csv("G:/*/Species_information.csv")

vars <- data[, c(
               "FFI_2020",
               "D_FFI",
               "loss_rate_2020",
               "ENN_2020",
               "loss_2020",
               "loss",
               "D_ENN",
               "LossRateSlope"
               )]
names(data)
cor_matrix <- cor(vars,
                  method = "spearman",
                  use = "complete.obs")
cor_matrix

cor_test <- rcorr(as.matrix(vars), type = "spearman")

cor_matrix <- cor_test$r
p_matrix <- cor_test$P
p_matrix

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()
par(family = "Arial")

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()
ggcorrplot(cor_matrix,
           type = "upper",
           lab = TRUE,
           p.mat = p_matrix,
           sig.level = 0.05,
           insig = "blank",
           colors = c("#2166ac", "white", "#b2182b"))+
           theme(text = element_text(family = "Arial"))


##### Quadrant_Analysis####

data$IUCN_status <- factor(data$IUCN_status,
                           levels = c("LC","NT","VU","EN","CR"),
                           ordered = TRUE)
ggplot(data, aes(x = loss_rate_2020,
                 y = FFI_2020,
                 color = IUCN_status)) +
  geom_point(size = 3, alpha = 0.8) +
  
  scale_color_manual(values = c(
    "LC" = "#4CAF50",
    "NT" = "#8BC34A",
    "VU" = "#FFC107",
    "EN" = "#FF7043",
    "CR" = "#D32F2F"
  )) +
  theme_minimal(base_size = 16) +
  
  labs(x = "Forest loss rate (2020)",
       y = "Forest Fragmentation Intensity (FFI)",
       color = "IUCN Status",
       title = "Species exposure to forest loss and fragmentation") +
  
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


#### ENN &FFI GAM####

data <- read.csv("G:/*/Species_information.csv")

model_gam <- gam(FFI_2020 ~ s(ENN_2020, k = 10), data = data)
summary(model_gam)

new_data <- data.frame(ENN_2020 = seq(min(data$ENN_2020), max(data$ENN_2020), length.out = 500))
pred <- predict(model_gam, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

new_data <- new_data %>%
  mutate(
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 

sci_label <- function(x) {
  lab <- scientific_format()(x)
  lab <- gsub("e\\+", "e", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
windowsFonts(Arial = windowsFont("Arial"))
ggplot() +
  geom_point(data = data, aes(x = ENN_2020, y = FFI_2020),
             size = 4, alpha = 0.4, shape = 16, color = "grey50") +
  geom_ribbon(data = new_data, aes(x = ENN_2020, ymin = lower, ymax = upper),
              fill = "#2b8cbe",
              alpha = 0.3) +
  geom_line(data = new_data, aes(x = ENN_2020, y = fit),
            color = "#2b83ba",
            linewidth = 2) +
  coord_cartesian(ylim = c(-0.25, 1))+
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.5))+
  scale_x_continuous()+
  labs(
    x = expression(ENN[2020]),
    y = expression(FFI[2020])
  )+
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line = element_blank(),  
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.title.x = element_text(size = 25, face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = ggplot2::margin(r = 10)),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
coord_fixed(ratio = 1)

summary(model_gam)


#### loss &FFI GAM####

data <- read.csv("G:/*/Species_information.csv")

model_gam <- gam(FFI_2020 ~ s(loss_rate_2020_log, k = 10), data = data)
summary(model_gam)

new_data <- data.frame(loss_rate_2020_log = seq(min(data$loss_rate_2020_log), max(data$loss_rate_2020_log), length.out = 500))
pred <- predict(model_gam, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

new_data <- new_data %>%
  mutate(
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 

sci_label <- function(x) {
  lab <- scientific_format()(x)
  lab <- gsub("e\\+", "e", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
windowsFonts(Arial = windowsFont("Arial"))
ggplot() +
  geom_point(data = data, aes(x = loss_rate_2020_log, y = FFI_2020),
             size = 4, alpha = 0.4, color = "grey50") +
  geom_ribbon(data = new_data, aes(x = loss_rate_2020_log, ymin = lower, ymax = upper),
              fill = "#FF0000",alpha = 0.25) +
  geom_line(data = new_data, aes(x = loss_rate_2020_log, y = fit),
            color = "#FF0000",
            linewidth = 2)  +
  coord_cartesian(ylim = c(-0.5, 1.2))+
  scale_x_continuous( )+
  labs(
    #x = expression(Loss[2020]),
    x = expression(Log(loss~rate[2020])),
    #y = expression(FFI[2020])
  )+
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line = element_blank(), 
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.title.x = element_text(size = 25,face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = element_blank(),
    #axis.title.y = element_text(size = 25,face = "bold", margin = ggplot2::margin(r = 10)),
    panel.background = element_blank(),
    plot.background = element_blank()
  )


summary(model_gam)


###### Body_Size_and_Forest_Dependency ###########
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)
data$Body_Size <- as.factor(data$Body_Size)

table(data$Forest_Dependency, data$Body_Size)

prop_table <- prop.table(table(data$Forest_Dependency, data$Body_Size), margin = 1)
print(prop_table)

chisq_test <- chisq.test(data$Forest_Dependency, data$Body_Size)
print(chisq_test)


plot_data <- data %>%
  count(Forest_Dependency, Body_Size) %>%
  group_by(Forest_Dependency) %>%
  mutate(proportion = n / sum(n))

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()

ggplot(plot_data, aes(x = Forest_Dependency, y = proportion, fill = Body_Size)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Large_bodied" = "#E41A1C", "Small_medium" = "#377EB8"),
                    name = "Body size") +
  labs(x = "Forest dependency", y = "Proportion") +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top"
  )

ggsave("ForestDependency_BodySize_proportions.tiff", 
       width = 6, height = 5, units = "cm", dpi = 1200, compression = "lzw")


#####T-test_Forest_Dependency_vs_Range_Size####
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)

t_model <- t.test(Species_Area ~ Forest_Dependency, data = data,var.equal = FALSE)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`Species_Area` ~ Forest_Dependency, data=data, mean)
aggregate(`Species_Area` ~ Forest_Dependency, data=data, sd)


boxplot(`Species_Area` ~ Forest_Dependency,
        data=data,
        xlab="Forest dependency",
        ylab="Species_Area",
        main="Species_Area between forest dependency classes")

#####T-test_Forest_Dependency_vs_Forest_Area####
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)

t_model <- t.test(Forest_coverage ~ Forest_Dependency, data = data,var.equal = FALSE)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`Forest_coverage` ~ Forest_Dependency, data=data, mean)
aggregate(`Species_Area` ~ Forest_Dependency, data=data, sd)



boxplot(`Forest_coverage` ~ Forest_Dependency,
        data=data,
        xlab="Forest dependency",
        ylab="Forest_coverage",
        main="Forest_coverage between forest dependency classes")

#####T-test_Forest_Dependency_vs_Risk####
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)

t_model <- t.test(Risk ~ Forest_Dependency, data = data,var.equal = FALSE)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`Risk` ~ Forest_Dependency, data=data, mean)
aggregate(`Risk` ~ Forest_Dependency, data=data, sd)



boxplot(`Risk` ~ Forest_Dependency,
        data=data,
        xlab="Forest dependency",
        ylab="Risk",
        main="Risk between forest dependency classes")

#####T-test_Forest_Dependency_vs_Fragmentation_Change####
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- as.factor(data$Forest_Dependency)
t_model <- t.test(D_FFI ~ Forest_Dependency, data = data,var.equal = FALSE)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")


sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`D_FFI` ~ Forest_Dependency, data=data, mean)
aggregate(`D_FFI` ~ Forest_Dependency, data=data, sd)

boxplot(`D_FFI` ~ Forest_Dependency,
        data=data,
        xlab="Forest dependency",
        ylab="Change_FFI",
        main="Change_FFI between forest dependency classes")

#####T-test_Body_Size_and_Forest_Dependency_vs_Fragmentation_Change####
data <- read.csv("G:/*/Species_information.csv")

data$Body_Size <- as.factor(data$Body_Size)          
data$Forest_Dependency <- as.factor(data$Forest_Dependency)  

data$Group <- with(data, interaction(Body_Size, Forest_Dependency, sep = "_"))
data$Group <- factor(data$Group, levels = c("Large_bodied_High", "Large_bodied_Medium",
                                            "Small_medium_High", "Small_medium_Medium"))

cat("\n========== Large-bodied: High vs Medium dependency ==========\n")
large_data <- subset(data, Body_Size == "Large_bodied")
t_large <- t.test(D_FFI ~ Forest_Dependency, data = large_data)
print(t_large)

cat("\n========== Small-medium: High vs Medium dependency ==========\n")
small_data <- subset(data, Body_Size == "Small_medium")
t_small <- t.test(D_FFI ~ Forest_Dependency, data = small_data)
print(t_small)

t_large_val <- round(t_large$statistic, 2)
p_large_val <- t_large$p.value
t_small_val <- round(t_small$statistic, 2)
p_small_val <- t_small$p.value

sig_label <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

annotation_df <- data.frame(
  Body_Size = c("Large_bodied", "Small_medium"),
  label = c(paste0("t = ", t_large_val, ", p = ", formatC(p_large_val, format = "e", digits = 2)),
            paste0("t = ", t_small_val, ", p = ", formatC(p_small_val, format = "e", digits = 2))),
  y = c(max(data$D_FFI, na.rm = TRUE) * 0.9, max(data$D_FFI, na.rm = TRUE) * 0.9)
)

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()

dependency_colors <- c("High" = "#2E8B57", "Medium" = "#4169E1")

p <- ggplot(data, aes(x = Body_Size, y = D_FFI, fill = Forest_Dependency)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.size = 0.8) +
  scale_fill_manual(values = dependency_colors, name = "Forest dependency") +
  labs(x = "Body size class", y = expression(Delta * FFI)) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )


p <- p + geom_text(data = annotation_df, aes(x = Body_Size, y = y, label = label),
                   inherit.aes = FALSE, size = 4, vjust = 0, family = "Arial")

print(p)

ggsave("BodySize_Dependency_D_FFI.tiff", plot = p, width = 6, height = 5, 
       units = "cm", dpi = 1200, compression = "lzw")

#####T-test_Body_Size_vs_Forest_and_Range_Area####

data <- read.csv("G:/*/Species_information.csv")

data$Body_Size <- as.factor(data$Body_Size)

t_model <- t.test(Species_Area ~ Body_Size, data = data)

t_value <- t_model$statistic
p_value <- t_model$p.value

cat("t value =", t_value, "\n")
cat("p value =", p_value, "\n")

sig <- ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01, "**",
                     ifelse(p_value < 0.05, "*", "ns")))

cat("Significance:", sig, "\n")

aggregate(`Species_Area` ~ Body_Size, data=data, mean)
aggregate(`Species_Area` ~ Body_Size, data=data, sd)

boxplot(`Species_Area` ~ Body_Size,
        data=data,
        xlab="Body size class",
        ylab="Species_Area",
        main="Fragmentation Change among body-size classes")
#####GAM_Range_Size_vs_Fragmentation_Change#####

data <- read.csv("G:/*/Species_information.csv")
data <- data[is.finite(data$Species_Area) & is.finite(data$D_FFI), ]

data$Log_Area <- log10(data$Species_Area)

data$Forest_Dependency <- factor(data$Forest_Dependency,
                                 levels = c("High", "Medium"),
                                 labels = c("Specialist", "Generalist"))

model_all <- gam(D_FFI ~ s(Log_Area, k = 5), data = data)
summary(model_all)

model_int <- gam(D_FFI ~ Forest_Dependency + s(Log_Area, by = Forest_Dependency, k = 5),
                 data = data)
summary(model_int)

new_x <- seq(min(data$Log_Area), max(data$Log_Area), length.out = 200)

new_all <- data.frame(Log_Area = new_x)
pred_all <- predict(model_all, newdata = new_all, se.fit = TRUE)
new_all$fit <- pred_all$fit
new_all$upper <- new_all$fit + 1.96 * pred_all$se.fit
new_all$lower <- new_all$fit - 1.96 * pred_all$se.fit
new_all$Forest_Dependency <- "All species"


new_data <- expand.grid(
  Log_Area = new_x,
  Forest_Dependency = levels(data$Forest_Dependency)
)
pred_int <- predict(model_int, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred_int$fit
new_data$se <- pred_int$se.fit
new_data$upper <- new_data$fit + 1.96 * new_data$se
new_data$lower <- new_data$fit - 1.96 * new_data$se

new_all <- new_all[, c("Log_Area", "fit", "upper", "lower", "Forest_Dependency")]
new_data <- new_data[, c("Log_Area", "fit", "upper", "lower", "Forest_Dependency")]
pred_all <- bind_rows(new_all, new_data)
pred_all$Forest_Dependency <- factor(pred_all$Forest_Dependency, 
                                     levels = c("All species", "Specialist", "Generalist"))


font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()

ggplot() +

  geom_point(data = data,
             aes(x = Log_Area, y = D_FFI),
             color = "grey70", alpha = 0.3, size = 3) +

  geom_point(data = data[data$Forest_Dependency == "Specialist", ],
             aes(x = Log_Area, y = D_FFI),
             color = "#006F00", alpha = 0.3, size = 3) +

  geom_point(data = data[data$Forest_Dependency == "Generalist", ],
             aes(x = Log_Area, y = D_FFI),
             color = "#CBCC00", alpha = 0.3, size = 3) +

  geom_ribbon(data = pred_all,
              aes(x = Log_Area, ymin = lower, ymax = upper, fill = Forest_Dependency),
              alpha = 0.15, color = NA) +
  geom_line(data = pred_all[pred_all$Forest_Dependency == "All species", ],
            aes(x = Log_Area, y = fit, color = Forest_Dependency),
            linewidth = 2.5, linetype = "dashed") +
  geom_line(data = pred_all[pred_all$Forest_Dependency != "All species", ],
            aes(x = Log_Area, y = fit, color = Forest_Dependency),
            linewidth = 2.5) +
  scale_color_manual(
    values = c("All species" = "grey50",
               "Specialist" = "#006F00",
               "Generalist" = "#CBCC00"),
    breaks = c("All species", "Specialist", "Generalist")
  ) +
  scale_fill_manual(
    values = c("All species" = "grey50",
               "Specialist" = "#006F00",
               "Generalist" = "#CBCC00"),
    breaks = c("All species", "Specialist", "Generalist")
  ) +
  scale_y_continuous(
    breaks = seq(-0.1, 0.1, by = 0.1)) +
  labs(x = expression("Log"[10] * " (Species range size)"),
       y = expression(Delta * FFI),
       color = NULL, fill = NULL) +
  coord_cartesian(
    xlim = c(2, 7.5),
    ylim = c(-0.1, 0.1)) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text = element_text(size = 30, color = "black"),
    axis.title = element_text(size = 30),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(size = 20),   
    legend.title = element_text(size = 22), 
    legend.position = c(0.18, 0.88),
    legend.background = element_rect(fill = NA, color = NA, linewidth = 0.3),
    legend.key = element_blank()
  )

summary(model_gam)

ggsave("Fig3e_Interaction_GAM_log_v2.tiff", width = 8, height = 8, units = "cm", dpi = 1200, compression = "lzw")

getwd()

##### Boxplot_IUCN_vs_Fragmentation_Change####

data <- read.csv("G:/*/Species_information.csv")

data <- data %>% filter(IUCN_status %in% c("LC","NT","VU","EN","CR"))
data$IUCN_status <- factor(data$IUCN_status, levels = c("LC","NT","VU","EN","CR"))

run_anova <- function(data_sub, group_name) {
  cat("\n========== ", group_name, " ==========\n")
  
  anova_model <- aov(D_FFI ~ IUCN_status, data = data_sub)
  anova_result <- summary(anova_model)
  
  F_value <- anova_result[[1]]$`F value`[1]
  p_value <- anova_result[[1]]$`Pr(>F)`[1]
  cat("F value =", F_value, "\n")
  cat("p value =", p_value, "\n")
  sig <- ifelse(p_value < 0.001, "***",
                ifelse(p_value < 0.01, "**",
                       ifelse(p_value < 0.05, "*", "ns")))
  cat("Significance:", sig, "\n")
  
  cat("\nMean D_FFI by IUCN:\n")
  print(aggregate(D_FFI ~ IUCN_status, data = data_sub, mean))
  cat("\nSD D_FFI by IUCN:\n")
  print(aggregate(D_FFI ~ IUCN_status, data = data_sub, sd))

  if (p_value < 0.05) {
    tukey <- TukeyHSD(anova_model)
    tukey_letters <- multcompLetters4(anova_model, tukey)
    letters_data <- data.frame(
      IUCN_status = names(tukey_letters$IUCN_status$Letters),
      Letters = tukey_letters$IUCN_status$Letters
    )
    cat("\nTukey letters:\n")
    print(letters_data)
  } else {
    cat("\nTukey test skipped (ANOVA not significant).\n")
  }
}

run_anova(data, "All species")

data_high <- data[data$Forest_Dependency == "High", ]
run_anova(data_high, "Forest specialist")

data_med <- data[data$Forest_Dependency == "Medium", ]
run_anova(data_med, "Forest generalist")

data$Group <- "All species"
data_high$Group <- "Forest specialist"
data_med$Group <- "Forest generalist"
plot_data <- bind_rows(data, data_high, data_med)
plot_data$Group <- factor(plot_data$Group, levels = c("All species", "Forest specialist", "Forest generalist"))

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto()

group_colors <- c("All species" = "#999999", "Forest specialist" = "#2E8B57", "Forest generalist" = "#4169E1")

ggplot(plot_data, aes(x = IUCN_status, y = D_FFI, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual(values = group_colors, guide = "none") +
  labs(x = "IUCN status", y = "ΔFFI") +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 11)
  )

ggsave("GAM_three_lines.tiff", width = 8, height = 4, units = "cm", dpi = 1200, compression = "lzw")

##### GAM_ENN_Change_vs_Fragmentation_Change ####

data <- read.csv("G:/*/Species_information.csv")
data <- data[is.finite(data$D_ENN) & is.finite(data$D_FFI), ] 
model_gam <- gam(D_FFI ~ s(D_ENN, k = 8), data = data)
summary(model_gam)

new_data <- data.frame(D_ENN = seq(min(data$D_ENN), max(data$D_ENN), length.out = 500))
pred <- predict(model_gam, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

new_data <- new_data %>%
  mutate(
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 

sci_label <- function(x) {
  lab <- scientific_format()(x)
  lab <- gsub("e\\+", "e", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
windowsFonts(Arial = windowsFont("Arial"))
ggplot() +
  geom_point(data = data, aes(x = D_ENN, y = D_FFI),
             size = 4, alpha = 0.4, shape = 16, color = "grey50") +
  geom_ribbon(data = new_data, aes(x = D_ENN, ymin = lower, ymax = upper),
              fill = "#2b8cbe",
              alpha = 0.3) +
  geom_line(data = new_data, aes(x = D_ENN, y = fit),
            color = "#2b83ba",
            linewidth = 2) +
  coord_cartesian(ylim = c(-0.2, 0.2))+
  scale_y_continuous(breaks = seq(-0.2, 0.2, by = 0.2))+
  scale_x_continuous()+
  labs(
    x = expression(ΔENN),
    y = expression(ΔFFI)
  )+
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line = element_blank(), 
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.title.x = element_text(size = 25, face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = ggplot2::margin(r = 8)),
    panel.background = element_blank(),
    plot.background = element_blank()
  )
coord_fixed(ratio = 1)

summary(model_gam)


##### GAM_Forest_Loss_Trend_vs_Fragmentation_Change ####

data <- read.csv("G:/*/Species_information.csv")
data <- data[is.finite(data$LossRateSlope) & is.finite(data$D_FFI), ]

model_gam <- gam(D_FFI ~ s(LossRateSlope, k = 2), data = data)
summary(model_gam)

new_data <- data.frame(LossRateSlope = seq(min(data$LossRateSlope), max(data$LossRateSlope), length.out = 500))
pred <- predict(model_gam, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

new_data <- new_data %>%
  mutate(
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se
  )

font_add("Arial", "C:/Windows/Fonts/arial.ttf") 
showtext_auto() 

sci_label <- function(x) {
  lab <- scientific_format()(x)
  lab <- gsub("e\\+", "e", lab)
  parse(text = gsub("e", " %*% 10^", lab))
}
windowsFonts(Arial = windowsFont("Arial"))
ggplot() +
  geom_point(data = data, aes(x = LossRateSlope, y = D_FFI),
             size = 4, alpha = 0.4, color = "grey50") +
  geom_ribbon(data = new_data, aes(x = LossRateSlope, ymin = lower, ymax = upper),
              fill = "#FF0000",alpha = 0.25) +
  geom_line(data = new_data, aes(x = LossRateSlope, y = fit),
            color = "#FF0000",
            linewidth = 2)  +
  coord_cartesian(ylim = c(-0.2, 0.2))+
  scale_y_continuous(breaks = seq(-0.2, 0.2, by = 0.2))+
  scale_x_continuous( )+
  labs(
    x = "Loss Rate Slope",
   y = expression(Delta * FFI)
  )+
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line = element_blank(),  
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    
    axis.title.x = element_text(size = 25,face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 25,face = "bold", margin = ggplot2::margin(r = 6)),
    panel.background = element_blank(),
    plot.background = element_blank()
  )


summary(model_gam)

#### Z-score_Transformation #####

data <- read.csv("G:/*/Species_information.csv")

data$FFI_z <- scale(data$D_FFI)
data$Loss_z <- scale(data$LossRateSlope)

data$Risk <- (data$FFI_z + data$Loss_z) / 2
getwd()

write.csv(data, "species_risk_index-6-24.csv", row.names = TRUE)

###### Risk_vs_Forest_Dependency ####
data <- read.csv("G:/*/Species_information.csv")
data <- data[!is.na(data$LossRateSlope), ]

data$Forest_Dependency <- factor(data$Forest_Dependency,
                                 levels = c("High", "Medium"),
                                 labels = c("Specialist", "Generalist"))

 data$FFI_z <- scale(data$D_FFI)
 data$Loss_z <- scale(data$LossRateSlope)
 data$Risk <- (data$FFI_z + data$Loss_z) / 2

t_risk <- t.test(Risk ~ Forest_Dependency, data = data)
print(t_risk)

aggregate(Risk ~ Forest_Dependency, data = data, mean)
aggregate(Risk ~ Forest_Dependency, data = data, sd)

p1 <- ggplot(data, aes(x = Forest_Dependency, y = Risk, fill = Forest_Dependency)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("Specialist" = "#2E8B57", "Generalist" = "#4169E1")) +
  labs(x = "Forest dependency", y = "Risk index", 
       title = paste0("Risk difference (t = ", round(t_risk$statistic, 2), 
                      ", p = ", formatC(t_risk$p.value, format = "e", digits = 2), ")")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

t_loss <- t.test(LossRateSlope ~ Forest_Dependency, data = data)
print(t_loss)

aggregate(LossRateSlope ~ Forest_Dependency, data = data, mean)
aggregate(LossRateSlope ~ Forest_Dependency, data = data, sd)

p2 <- ggplot(data, aes(x = Forest_Dependency, y = LossRateSlope, fill = Forest_Dependency)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("Specialist" = "#2E8B57", "Generalist" = "#4169E1")) +
  labs(x = "Forest dependency", y = "Loss rate slope", 
       title = paste0("Loss trend difference (t = ", round(t_loss$statistic, 2), 
                      ", p = ", formatC(t_loss$p.value, format = "e", digits = 2), ")")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p1)
print(p2)

ggsave("Risk_by_dependency.tiff", plot = p1, width = 6, height = 5, units = "cm", dpi = 1200, compression = "lzw")
ggsave("LossRateSlope_by_dependency.tiff", plot = p2, width = 6, height = 5, units = "cm", dpi = 1200, compression = "lzw")

##### Quadrant_Proportions_1 #####
data <- read.csv("G:/*/Species_information.csv")

data$Forest_Dependency <- factor(data$Forest_Dependency,
                                 levels = c("High", "Medium"),
                                 labels = c("Specialist", "Generalist"))
data$IUCN_status <- factor(data$IUCN_status, levels = c("LC", "NT", "VU", "EN", "CR"))

high_data <- data[data$Forest_Dependency == "Specialist", ]

high_data <- high_data[!is.na(high_data$D_FFI) & !is.na(high_data$LossRateSlope), ]

high_data$Q1 <- (high_data$D_FFI > 0) & (high_data$LossRateSlope > 0)

quadrant_summary <- high_data %>%
  group_by(IUCN_status) %>%
  summarise(
    total_n = n(),
    in_Q1_n = sum(Q1),
    prop_Q1 = mean(Q1) * 100
  ) %>%
  ungroup()

print(quadrant_summary)

contingency_table <- table(high_data$IUCN_status, high_data$Q1)
print(contingency_table)

chi_test <- chisq.test(contingency_table)
print(chi_test)

fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 2000)
print(fisher_test)

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()
windowsFonts(Arial = windowsFont("Arial"))

p <- ggplot(quadrant_summary, aes(x = IUCN_status, y = prop_Q1, fill = IUCN_status)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(prop_Q1, 1), "%\n(", in_Q1_n, "/", total_n, ")")), 
            vjust = -0.3, size = 4, family = "Arial") +
  scale_fill_manual(values = c("LC" = "#66C2A5", "NT" = "#8DA0CB", "VU" = "#E78AC3", 
                               "EN" = "#FC8D62", "CR" = "#D73027")) +
  labs(x = "IUCN status", y = "Percentage in Quadrant I (%)",
       subtitle = paste0("χ² p = ", formatC(chi_test$p.value, format = "e", digits = 2))) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

print(p)

ggsave("Fig_Quadrant1_IUCN_raw.tiff", plot = p, width = 8, height = 6, units = "cm", dpi = 1200, compression = "lzw")

##### Quadrant_Proportions_2 #####
data <- read.csv("G:/*/Species_information.csv")

high_data <- data[data$Forest_Dependency == "High", ]

high_data <- high_data[!is.na(high_data$D_FFI) & !is.na(high_data$LossRateSlope), ]

high_data$Quadrant <- with(high_data, {
  case_when(
    D_FFI > 0 & LossRateSlope > 0   ~ "Q1 (ΔFFI↑, Loss↑)",   
    D_FFI > 0 & LossRateSlope <= 0  ~ "Q2 (ΔFFI↑, Loss↓)",   
    D_FFI <= 0 & LossRateSlope <= 0 ~ "Q3 (ΔFFI↓, Loss↓)",   
    D_FFI <= 0 & LossRateSlope > 0  ~ "Q4 (ΔFFI↓, Loss↑)"    
  )
})

high_data$IUCN_status <- factor(high_data$IUCN_status, levels = c("LC", "NT", "VU", "EN", "CR"))

quadrant_percent <- high_data %>%
  group_by(IUCN_status, Quadrant) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(IUCN_status) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()


print(quadrant_percent)

font_add("Arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()
windowsFonts(Arial = windowsFont("Arial"))

quadrant_colors <- c(
  "Q1 (ΔFFI↑, Loss↑)" = "#D73027",
  "Q2 (ΔFFI↑, Loss↓)" = "#4575B4",
  "Q3 (ΔFFI↓, Loss↓)" = "#999999",
  "Q4 (ΔFFI↓, Loss↑)" = "#FDB863"
)

p <- ggplot(quadrant_percent, aes(x = IUCN_status, y = percent, fill = Quadrant)) +
  geom_col(position = "fill", width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = quadrant_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = "IUCN status",
    y = "Percentage of species",
    fill = "Quadrant",
    title = "Distribution of forest specialists among landscape risk quadrants"
  ) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p)
write.csv(quadrant_percent, "high_dependency_quadrant_percent.csv", row.names = FALSE)

ggsave("Fig_S2_Quadrant_StackedBar.tiff", plot = p, width = 10, height = 7, units = "cm", dpi = 1200, compression = "lzw")
