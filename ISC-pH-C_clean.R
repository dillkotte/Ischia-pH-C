

##### R script for analysis of carbon cycling in Posidonia oceanica meadows in CO2 vent systems in Ischia, Italy ####
##  Author: Theodor Kindeberg
##  Affiliation: Laboratoire d'Océanographie de Villefranche, Sorbonne Université-CNRS
##  Contact: theo.kindeberg@gmail.com
##  Study: "Ocean acidification enhances carbon burial in seagrass meadows: new insights from CO2 vents"
##  Study authors: Theodor Kindeberg, Núria Teixido, Steeve Comeau, Jean-Pierre Gattuso, Beat Gasser, Alice Mirasole, Samir Alliouane, Ioannis Kalaitzakis, Denisa Berbece, Christopher Cornwall, Pere Masque
##  

rm(list = ls())
#### 01 Load libraries ####
library(ggplot2) #for plotting
library(dplyr) #for data handling, formatting etc
library(tidyr) #for tidy data
library(zoo) #for interpolation
library(patchwork) #for attaching multiple panels into a figure
library(tidypaleo) #For tidy core profile plots
library(ggpubr)
library(lme4) #For linear mixed effects models
library(lmerTest) #For ANOVA test of LMMs
library(nlme) #for GLS models
library(MuMIn) # for model selection
library(emmeans) #for post hoc of lmer and nlme models
library(sjPlot) #for output lmer table
library(simmr) #for stable isotope mixing model
library(car) #for leveneTest
library(rstatix) #For games howell post hoc test
library(DHARMa) #to check assumptions of models
library(sf) #to create map
library(ggmap) #to create map
library(rnaturalearth) #to create map
library(rnaturalearthdata) #to create map
library(ggspatial) #to create map
library(osmdata) #to create map
library(xlsx) # to export xlsx tables

#Clear environment
rm(list = ls())

#### 02 Set orders, colors etc. for figures ####
pH.colors <- c("#03A1D7", "#E09017", "#C00918")
pH.order <- as.factor(c("Ambient", "Low", "Extreme Low"))
theme_set(theme_pubclean(14))

#### 03 Load data ####
core_prof <- read.csv("data/share/1_Sediment.csv", header=T) #Sediment data
dating <- read.csv("data/share/2_Dating.csv", header=T) #210Pb dating
sources <- read.csv("data/share/3_Stable_isotope_sources.csv", header=T) #Stable isotope source data
sites <- read.csv("data/share/4_Sites.csv", header=T) #Coordinates of sampling sites
core_info <- read.csv("data/share/5_Cores.csv", header=T) #Core info

## Merge core_info with core_prof
core_prof <- core_prof %>%
  left_join(core_info, by = "Core")

## Add core surface area
core_prof$core_sa <- pi*(6.3/2)^2 #Core surface area (diameter 6.3 cm)

## Calculate compression corrected core depths
core_prof <- core_prof %>%
  mutate(
    Compaction = ((Corer.length - External.core.depth) -
      (Corer.length - Internal.core.depth)) / (Corer.length - External.core.depth),
    Depth.top.correct = Depth.top *
      (Corer.length - External.core.depth) /
      (Corer.length - Internal.core.depth),
    Depth.bottom.correct = Depth.bottom * 
      (Corer.length - External.core.depth) /
      (Corer.length - Internal.core.depth)
  )

## Calculate inorganic carbon content
# assuming C makes up 27% of CO2 (mw_C = 12.01 g/mol; mw_O2 = 15.99*2 = 31.98 g/mol)
core_prof$IC <- core_prof$X.Acid*12.01/(12.01+(15.99*2))

## Create custom function to interpolate core profiles of OC, IC, d13C and d15N
#Add new column "Depth.mid" equaling the middle section of each core slice
core_prof <- core_prof %>%
  mutate(Depth = (Depth.top + Depth.bottom) / 2)
#Sort by Core and Depth
# interpolate missing values using nearest true values in each core
interpolate_nearest <- function(df, core_col = "Core", input_col, output_col) {
  # Copy original column
  df[[output_col]] <- df[[input_col]]
  
  # Loop through rows
  for (i in seq_len(nrow(df))) {
    if (is.na(df[[output_col]][i])) {
      current_core <- df[[core_col]][i]  # Current core ID
      
      # Find previous non-NA value within the same core
      prev_values <- df[[input_col]][1:(i-1)][
        df[[core_col]][1:(i-1)] == current_core & !is.na(df[[input_col]][1:(i-1)])
      ]
      prev_value <- ifelse(length(prev_values) > 0, tail(prev_values, 1), NA)
      
      # Find next non-NA value within the same core
      next_value <- NA
      if (i < nrow(df)) {
        next_values <- df[[input_col]][(i+1):nrow(df)][
          df[[core_col]][(i+1):nrow(df)] == current_core & !is.na(df[[input_col]][(i+1):nrow(df)])
        ]
        next_value <- ifelse(length(next_values) > 0, head(next_values, 1), NA)
      }
      
      # Interpolate between true values
      if (!is.na(prev_value) & !is.na(next_value)) {
        df[[output_col]][i] <- (prev_value + next_value) / 2
      } else if (!is.na(prev_value)) {
        df[[output_col]][i] <- prev_value
      } else if (!is.na(next_value)) {
        df[[output_col]][i] <- next_value
      }
    }
  }
  
  return(df)
}

## Apply the function to create a new interpolated column for each Core
core_prof <- interpolate_nearest(core_prof, "Core","OC", "OC_Interpolated")
core_prof <- interpolate_nearest(core_prof, "Core","IC", "IC_Interpolated")
core_prof <- interpolate_nearest(core_prof, "Core","TN", "TN_Interpolated")
core_prof <- interpolate_nearest(core_prof, "Core","d13C", "d13C_Interpolated")
core_prof <- interpolate_nearest(core_prof, "Core","d15N", "d15N_Interpolated")


#### 04 Stocks and Rates #### 
# Join MAR ± Error from dating into core_prof by Core
core_prof <- core_prof %>%
  left_join(
    dating %>% dplyr::select(Core, MAR, Error_MAR),
    by = "Core"
  )

#Replace "Mixing" with "NA" when no MAR is available due to mixed core
core_prof$MAR[core_prof$MAR == "Mixing"] <- NA
core_prof$MAR <- as.numeric(core_prof$MAR) #Convert MAR to numeric

#Calculate carbon and nitrogen stocks for each slice (g/cm2)
core_prof$OC_stock <- (core_prof$OC_Interpolated/100) * core_prof$Dry.mass/core_prof$core_sa * 10000 # 
core_prof$IC_stock <- (core_prof$IC_Interpolated/100) * core_prof$Dry.mass/core_prof$core_sa * 10000 # 
core_prof$TN_stock <- (core_prof$TN_Interpolated/100) * core_prof$Dry.mass/core_prof$core_sa * 10000 # 

## Calculate stocks integrated from surface (g/cm2)
core_prof <- core_prof %>%
  group_by(Core) %>%
  mutate(
    OC_stock_sfc = cumsum((OC_Interpolated / 100) * Dry.mass / core_sa * 10000),
    IC_stock_sfc = cumsum((IC_Interpolated / 100) * Dry.mass / core_sa * 10000),
    TN_stock_sfc = cumsum((TN_Interpolated / 100) * Dry.mass / core_sa * 10000),
  ) %>%
  ungroup()  

## Burial rates:  Burial rate (g/m2/yr) = Content (%) * MAR (g/cm2/yr)

core_prof$OC_burial_rate <- (core_prof$OC_Interpolated/100)*core_prof$MAR*10000 
core_prof$Error_OC_burial_rate <- (core_prof$OC_Interpolated/100)*core_prof$Error_MAR*10000
core_prof$TN_burial_rate <- (core_prof$TN_Interpolated/100)*core_prof$MAR*10000
core_prof$Error_TN_burial_rate <- (core_prof$TN_Interpolated/100)*core_prof$Error_MAR*10000
core_prof$IC_burial_rate <- (core_prof$IC_Interpolated/100)*core_prof$MAR*10000
core_prof$Error_IC_burial_rate <- (core_prof$IC_Interpolated)/100*core_prof$Error_MAR*10000

#Calculate OC, IC and TN mass
core_prof$OC_mass <- (core_prof$OC_Interpolated/100)*core_prof$Dry.mass
core_prof$IC_mass <- (core_prof$IC_Interpolated/100)*core_prof$Dry.mass
core_prof$TN_mass <- (core_prof$TN_Interpolated/100)*core_prof$Dry.mass
#Carbon:nitrogen ratio (OC:TN)
core_prof$C_N_mass <- core_prof$OC_mass/core_prof$TN_mass
#Carbon to nitrogen molar mass ratio is 12.01 g/mol / 14.01 g/mol = 0.857
core_prof$C_N_molar <- core_prof$C_N_mass/0.857

#IC:OC ratio
core_prof$OC_IC_mass <- core_prof$OC_mass/core_prof$IC_mass #assuming 12% IC in CaCO3
#core_prof$OC_IC_mass[is.infinite(core_prof$OC_IC_mass)] <- NA #Replace Inf with NA. Inf occurs because dry mass = 0 in some slices.
core_prof$OC_IC_molar <- core_prof$OC_Interpolated/core_prof$IC_Interpolated

#CO2eq removal rate (Mg CO2eq ha-1 yr-1) #based only on OC burial. Each gram of carbon corresponds to 44.01/12.01 = 3.664 g CO2 / g C
core_prof$CO2eqRR_OC <- (core_prof$OC_burial_rate * 3.664) /100
#IC burial release rate (g CO2eq m-2 yr-1) 
core_prof$CO2eqRR_IC <- (core_prof$IC_burial_rate * 3.664 * 0.6) /100
core_prof$CO2eqRR_net <- (core_prof$OC_burial_rate - 
                            core_prof$IC_burial_rate * 0.6) * (44.01/12.01) 
core_prof$gC_RR_net <- core_prof$CO2eqRR_net * (12.01 / 44.01) * 100 #in g C/m2/yr



#Subset only cores with reliable dating
core_prof_subset <- subset(core_prof, Core!="ISC-12" &
                             Core!="ISC-14" &
                             Core!="ISC-11" & 
                             Core!="ISC-17" & 
                             Core!="ISC-3" & 
                             Core!="ISC-4") 

#Calculate common oldest date
group_by(core_prof_subset, Core) %>%
  summarise(
    count = n(),
    min = min(Date, na.rm = TRUE), 
    max = max(Date, na.rm = TRUE))

#Subset dataset to common oldest date (1954)
core_prof_common <- subset(core_prof_subset, Date>1953)


#### 05 Calculate descriptive statistics for the whole dataset grouped by pH regime #### 

#Water depth
core_prof %>% 
  group_by(pH.regime) %>%  
  summarise(
    n = n(),
    mean = mean(Water.depth, na.rm = TRUE),
    sd = sd(Water.depth, na.rm = TRUE),
    se = sd / sqrt(n),
    min = min(Water.depth, na.rm = TRUE), 
    max = max(Water.depth, na.rm = TRUE),
    median = median(Water.depth, na.rm = TRUE),
    .groups = 'drop'
  )

#Depth bottom core
core_prof %>% 
  group_by(Site, Core) %>%  # Group by pH.regime and Site to get max Depth.bottom per core
  summarise(max_depth = max(Depth.bottom.correct, na.rm = TRUE), .groups = 'drop') %>% 
  group_by(Site) %>% 
  summarise(
    n = n(),
    mean = mean(max_depth, na.rm = TRUE),
    sd = sd(max_depth, na.rm = TRUE),
    se = sd / sqrt(n),
    min = min(max_depth, na.rm = TRUE), 
    max = max(max_depth, na.rm = TRUE),
    median = median(max_depth, na.rm = TRUE),
    .groups = 'drop'
  )

#Calculate OC stock of deepest core parts
#Subset only deep layers (beyond remineralization, c.f. Johannessen 2023)
core_prof_deep <- subset(core_prof, Depth.bottom >= 15 ) #stable profiles on average below 15 cm depth based on visual inspection of core profiles


#### 06 Group means ####
SE <- function(x) sd(x, na.rm=TRUE)/(sqrt(sum(!is.na(x))))
SD <- function(x) sd(x, na.rm=TRUE)


#Average±SD 1954 to today (oldest common age) per core and pH regime
avg_common <- core_prof_common %>%
  group_by(Core, pH.regime) %>%
  summarise(
    across(
      .cols = c(
        OC_burial_rate, OC_stock, OC_stock_sfc,
        TN_burial_rate, TN_stock,
        d13C, d15N, C_N_molar,
        IC_burial_rate, IC_stock, IC_stock_sfc,
        OC_IC_molar, CO2eqRR_net
      ),
      .fns = list(mean = ~mean(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )


#Average deep parts of core
avg_deep <- core_prof_deep %>%
  group_by(Core, pH.regime) %>%
  summarise(
    OC_stock_mean = mean(OC_stock, na.rm = TRUE),
    OC_stock_sd   = sd(OC_stock, na.rm = TRUE),
    OC_burial_mean = mean(OC_burial_rate, na.rm = TRUE),
    OC_burial_sd   = sd(OC_burial_rate, na.rm = TRUE),
    IC_stock_mean = mean(IC_stock, na.rm = TRUE),
    IC_stock_sd   = sd(IC_stock, na.rm = TRUE),
    IC_burial_mean = mean(IC_burial_rate, na.rm = TRUE),
    IC_burial_sd   = sd(IC_burial_rate, na.rm = TRUE),
    .groups = "drop"
  )

##Core mean average C stock in top 10 cm (depth corrected for compaction)
OC_stock_10cm_calc <- core_prof_subset %>%
  # Select depth closest to 10 for each core
  group_by(Core) %>%
  slice_min(abs(Depth.bottom.correct - 10), n = 1, with_ties = FALSE) %>% 
  mutate(
    distance_from_10 = Depth.bottom.correct - 10,
    abs_distance = abs(distance_from_10)
  ) %>%
  ungroup() %>%
  mutate(sd_distance = sd(distance_from_10, na.rm = TRUE)) %>%
  # Caclulate 10 cm OC stock by pH regime
  group_by(factor(pH.regime, levels = pH.order)) %>%
  summarise(
    n = sum(!is.na(OC_stock_sfc)),
    mean = mean(OC_stock_sfc, na.rm = TRUE),
    sd = sd(OC_stock_sfc, na.rm = TRUE),
    se = sd / sqrt(n()),
    median = median(OC_stock_sfc, na.rm = TRUE),
    iqr = IQR(OC_stock_sfc, na.rm = TRUE),
    sd_distance_from_10 = first(sd_distance))
OC_stock_10cm_calc

## Core average C stock in top 20 cm (depth corrected for compaction)
OC_stock_20cm_calc <- core_prof_subset %>%
# Select depth closest to 20 for each core
  group_by(Core) %>%
  slice_min(abs(Depth.bottom.correct - 20), n = 1, with_ties = FALSE) %>% 
  mutate(
    distance_from_20 = Depth.bottom.correct - 20,
    abs_distance = abs(distance_from_20)
  ) %>%
  ungroup() %>%
  mutate(sd_distance = sd(distance_from_20, na.rm = TRUE)) %>%
# Caclulate 20 cm OC stock by pH regime
  group_by(factor(pH.regime, levels = pH.order)) %>%
  summarise(
    n = sum(!is.na(OC_stock_sfc)),
    mean = mean(OC_stock_sfc, na.rm = TRUE),
    sd = sd(OC_stock_sfc, na.rm = TRUE),
    se = sd / sqrt(n()),
    median = median(OC_stock_sfc, na.rm = TRUE),
    iqr = IQR(OC_stock_sfc, na.rm = TRUE),
    sd_distance_from_20 = first(sd_distance))
OC_stock_20cm_calc


#### 07 Response across pH regimes ####
#Calculate CORE MEANS (n=8) of downcore means to avoid pseudoreplication within each core

# Function to calculate group means and propagated standard deviations

# Function to calculate summary statistics including propagated SD and SE
calculate_summary <- function(data, variable, sd_variable, order = NULL) {
  if (!is.null(order)) {
    data <- data %>%
      mutate(pH.regime = factor(pH.regime, levels = order))
  }
  
  data %>%
    group_by(pH.regime) %>%
    summarise(
      n = n(),
      mean = mean({{ variable }}, na.rm = TRUE),
      sd = sd({{ variable }}, na.rm = TRUE),
      se = sd / sqrt(n),
      median = median({{ variable }}, na.rm = TRUE),
      min = min({{ variable }}, na.rm = TRUE),
      max = max({{ variable }}, na.rm = TRUE),
      IQR = IQR({{ variable }}, na.rm = TRUE),
      cv = round(sd / mean, 2),
      
      # Calculate propagated standard deviation and standard error
      sd_propagated = sqrt(sum(({{ sd_variable }})^2, na.rm = TRUE)),
      se_propagated = sd_propagated / n
    )
}

# Applying the function to different datasets and variables
#Core mean burial rate
mean_OC_burial_pH_common <- calculate_summary(avg_common, OC_burial_rate_mean, OC_burial_rate_sd, pH.order)
mean_IC_burial_pH_common <- calculate_summary(avg_common, IC_burial_rate_mean, IC_burial_rate_sd, order=pH.order) #Inorganic C

#Core mean OC stock from surface
mean_OC_stock_pH_common <- calculate_summary(avg_common, OC_stock_sfc_mean, OC_stock_sfc_sd, order=pH.order)
#Core mean IC stock <- 
mean_IC_stock_pH_common <- calculate_summary(avg_common, IC_stock_sfc_mean, IC_stock_sfc_sd, order=pH.order)
#Core mean OC:IC ratio
mean_OC_IC_pH_common <- calculate_summary(avg_common, OC_IC_molar_mean, OC_IC_molar_sd, order=pH.order)
#Core mean net CO2eq flux
mean_CO2eq_pH_common <- calculate_summary(avg_common, CO2eqRR_net_mean, CO2eqRR_net_sd, order=pH.order)

#Core mean deep stock
mean_OC_stock_pH_deep <- calculate_summary(avg_deep, OC_stock_mean, OC_stock_sd, order=pH.order)
mean_OC_stock_pH_deep

#Core mean deep stock
mean_OC_burial_pH_deep <- calculate_summary(avg_deep, OC_burial_mean, OC_burial_sd, order=pH.order)
mean_OC_burial_pH_deep


#### 08 Figure 1 ####
#OC burial rate since 1954 until today (calculated from surface)
p.OC_burial_since1954 <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                           y=OC_burial_rate, fill=factor(pH.regime, levels=pH.order))) +
  geom_violin(alpha=0.2, color="white") +
  geom_point(data=mean_OC_burial_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=4) +
  geom_errorbar(data=mean_OC_burial_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  geom_point(data=mean_OC_burial_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=median, color=factor(pH.regime, levels=pH.order)),
             size=2, shape=5) +
  ylab(bquote(OC ~burial~ rate ~ (g~ m^-2~yr^-1))) +
  xlab("pH regime") +
  ylim(0,23) +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")

#OC stock from surface since 1954 until today

p.OC_stock_since1954 <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                          y=OC_stock_sfc, fill=factor(pH.regime, levels=pH.order))) +
  #geom_violin(alpha=0.2, color="white") +
  geom_point(data=mean_OC_stock_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=4) +
  geom_violin(alpha=0.2, color="white") +
  geom_errorbar(data=mean_OC_stock_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  geom_point(data=mean_OC_stock_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=median, color=factor(pH.regime, levels=pH.order)),
             size=2, shape=5) +
  ylab(bquote(OC ~stock ~ (g~ m^-2))) +
  xlab("pH regime") +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")
p.OC_stock_since1954 + p.OC_burial_since1954

#Core mean stable isotopes
#mean_d13C_pH <- calculate_summary(avg_pH, d13C, SD.d13C,  order=pH.order)
#mean_d15N_pH <- calculate_summary(avg_pH, d15N, SD.d15N, order=pH.order)
mean_d13C_pH_common <- calculate_summary(avg_common, d13C_mean, d13C_sd, order=pH.order)
mean_d15N_pH_common <- calculate_summary(avg_common, d15N_mean, d15N_sd, order=pH.order)
#Core mean C:N ratio (molar)
mean_C_N_pH_common <- calculate_summary(avg_common, C_N_molar_mean, C_N_molar_sd, order=pH.order)



#### 09 Linear models ####
#Linear regression of OC and IC
#all data (n=14)
fit.OC_IC_amb <- lm(OC_mass~IC_mass, data=subset(core_prof, pH.regime =="Ambient"))
summary(fit.OC_IC_amb)
fit.OC_IC_low <- lm(OC_mass~IC_mass, data=subset(core_prof, pH.regime =="Low"))
summary(fit.OC_IC_low)
fit.OC_IC_ext <- lm(OC_mass~IC_mass, data=subset(core_prof, pH.regime =="Extreme Low"))
summary(fit.OC_IC_ext)

#since 1954 (n=8 cores)
fit.OC_IC_amb <- lm(OC_mass~IC_mass, data=subset(core_prof_common, pH.regime =="Ambient"))
summary(fit.OC_IC_amb)
fit.OC_IC_low <- lm(OC_mass~IC_mass, data=subset(core_prof_common, pH.regime =="Low"))
summary(fit.OC_IC_low)
fit.OC_IC_ext <- lm(OC_mass~IC_mass, data=subset(core_prof_common, pH.regime =="Extreme Low"))
summary(fit.OC_IC_ext)


#### 10 Figure 2 #### 
#a
#OC:IC ratio since 1954
p.a <- ggplot(core_prof_common, aes(x=IC_mass, y=OC_mass, color=factor(pH.regime, level=pH.order))) +
  geom_point(aes(color=factor(pH.regime, level = pH.order))) +
  stat_smooth(method = "lm", formula = y ~ x, se=T) + 
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_x_continuous(limits=c(0.0, 0.6), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_y_continuous(limits=c(0.0, 0.6), breaks=c(0, 0.2, 0.4, 0.6)) +
  theme_classic(14) 
  #theme(legend.position = "none")
p.a
#b OC:IC  all data
p.b <- ggplot(core_prof, aes(x=IC_mass, y=OC_mass, color=factor(pH.regime, level=pH.order))) +
  geom_point(aes(color=factor(pH.regime, level = pH.order))) +
  stat_smooth(method = "lm", formula = y ~ x, se=T) + 
  scale_color_manual(values=pH.colors, name="pH regime") +
  xlab("") + ylab("") +
  scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1.0)) +
  scale_y_continuous(limits=c(0.0,1), breaks=c(0, 0.5, 1.0)) +
  theme_classic(14) +
  theme(legend.position = "none")
  
p.b

#c 
#Net CO2 sequestered since 1954
p.c <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                  y=CO2eqRR_net, fill=factor(pH.regime, levels=pH.order))) +
  geom_violin(alpha=0.2, color="white") +
  geom_point(data=mean_CO2eq_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=median, color=factor(pH.regime, levels=pH.order)),
             size=4) +
  geom_errorbar(data=mean_CO2eq_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=median, ymin=median-IQR, ymax=median+IQR,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  
  xlab("pH regime") +
  scale_y_continuous(name=bquote(Net~CO[2]~sequestered~(g~CO[2]~eq~yr^-1)), limits=c(-80, 80), breaks=c(-80, -40, 0, 40, 80)) +
  geom_abline(slope=0, intercept=0, linetype="dashed") +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none",
        text = element_text(size = 12, color = "black"),  # Sets all text elements
        axis.title = element_text(size = 12, color = "black"),  # Axis labels
        axis.text = element_text(size = 12, color = "black"),   # Axis tick labels
        legend.text = element_text(size = 12, color = "black"), # Legend text
        plot.title = element_text(size = 12, color = "black")   # Plot title
  ) +
  theme_classic(14) +
  theme(legend.position = "none")
p.c
p.a + inset_element(p.b, 0.4, 0.4, 1, 1) | p.c


#### 11 LMM and GLS ####

### Dated cores 1954-2021 (dataset=core_prof_common) ###

## OC burial rate since 1954

#Proceed with generalized linear mixed model (package nmle) because of heteroscedascity (use Maximum likelihood (ML) for model selection)
mod.OC.burial_gls_full <- gls(OC_burial_rate ~ pH.regime+Water.depth+Site, data = core_prof_common,
                              correlation = corCompSymm(form = ~1 | Core),
                              weights = varIdent(form = ~1 | Core), 
                              na.action = na.fail,
                              method="ML") #Apply Maximum likelihood for model selection
#Construct all subsequent models and perform model selection
OC_burial_dredge <- dredge(mod.OC.burial_gls_full)
print(OC_burial_dredge)
OC_burial_modsel <- as.data.frame(OC_burial_dredge)
write.xlsx(OC_burial_modsel, "tables/OC_burial_modsel.xlsx")

#Best model (lowest AICc): pH regime + depth
mod.OC.burial_gls_selected <- gls(OC_burial_rate ~ pH.regime+Water.depth, data = core_prof_common,
                                   correlation = corCompSymm(form = ~1 | Core),
                                   weights = varIdent(form = ~1 | Core), 
                                  na.action = na.fail,
                                  method="REML") #apply REML for final model
summary(mod.OC.burial_gls_selected)
mod.OC.burial_table <-  summary(mod.OC.burial_gls_selected)$tTable
write.xlsx(mod.OC.burial_table, "tables/mod.OC.burial_table.xlsx")

#OC stock since 1954
mod.OC_stock <- lmer(OC_stock_sfc~pH.regime+(1|Core), 
                      data=core_prof_common)
simulateResiduals(mod.OC_stock, plot=T) #Unequal variances
ggqqplot(resid(mod.OC_stock))
shapiro.test(resid(mod.OC_stock)) #normal residuals

mod.OC.stock_gls_full <- gls(OC_stock_sfc ~ pH.regime + Water.depth + Site,
                             data = core_prof_common,
                             correlation = corCompSymm(form = ~1 | Core), #drop varIdent so model can converge (add it again below)
                             method="ML")
#Construct all subsequent models and perform model selection
OC_stock_dredge <- dredge(mod.OC.stock_gls_full)
print(OC_stock_dredge)
OC_stock_modsel <- as.data.frame(OC_stock_dredge)
write.xlsx(OC_stock_modsel, "tables/OC_stock_modsel.xlsx")
#Best model (lowest AICc): pH regime + Site
mod.OC.stock_gls_selected <- gls(OC_stock_sfc ~ pH.regime+Site, data = core_prof_common,
                                  correlation = corCompSymm(form = ~1 | Core),
                                  weights = varIdent(form = ~1 | Core), 
                                  na.action = na.fail,
                                  method="REML") #apply REML
summary(mod.OC.stock_gls_selected)
mod.OC.stock_table <-  summary(mod.OC.stock_gls_selected)$tTable
write.xlsx(mod.OC.stock_table, "tables/mod.OC.stock_table.xlsx")


### Inorganic carbon (IC) ###
mod.IC.burial_gls_full <- gls(IC_burial_rate ~ pH.regime+Water.depth+Site, data = core_prof_common,
                              correlation = corCompSymm(form = ~1 | Core),
                              weights = varIdent(form = ~1 | Core), 
                              na.action = na.fail,
                              method="ML")
#Construct all subsequent models and perform model selection
IC_burial_dredge <- dredge(mod.IC.burial_gls_full)
print(IC_burial_dredge)
IC_burial_modsel <- as.data.frame(IC_burial_dredge)
write.xlsx(IC_burial_modsel, "tables/IC_burial_modsel.xlsx")
#Best model (lowest AICc): null model
mod.IC.burial_gls_selected <- gls(IC_burial_rate ~ 1, data = core_prof_common,
                                  correlation = corCompSymm(form = ~1 | Core),
                                  weights = varIdent(form = ~1 | Core), 
                                  na.action = na.fail,
                                  method="REML") #apply REML
summary(mod.IC.burial_gls_selected)


#IC stock since 1954
mod.IC_stock <- lmer(IC_stock_sfc~pH.regime+(1|Core), 
                     data=core_prof_common)
anova(mod.IC_stock)
simulateResiduals(mod.IC_stock, plot=T) #Unequal variances
ggqqplot(resid(mod.IC_stock))
shapiro.test(resid(mod.IC_stock)) #normal residuals

mod.IC.stock_gls_full <- gls(IC_stock_sfc ~ pH.regime + Water.depth + Site,
                             data = core_prof_common,
                             correlation = corCompSymm(form = ~1 | Core), #drop varIdent so model can converge (add it again below)
                             method="ML")
#Construct all subsequent models and perform model selection
IC_stock_dredge <- dredge(mod.IC.stock_gls_full)
print(IC_stock_dredge)
IC_stock_modsel <- as.data.frame(IC_stock_dredge)
write.xlsx(IC_stock_modsel, "tables/IC_stock_modsel.xlsx")
#Best model (lowest AICc): pH regime + Site
mod.IC.stock_gls_selected <- gls(IC_stock_sfc ~ 1, data = core_prof_common,
                                 correlation = corCompSymm(form = ~1 | Core),
                                 weights = varIdent(form = ~1 | Core), 
                                 na.action = na.fail,
                                 method="REML") #apply REML
summary(mod.IC.stock_gls_selected)


### All cores (dataset=core_prof) ###

#d13C (whole dataset)
mod.d13C_full <- lmer(d13C_Interpolated ~pH.regime + Water.depth + Site + (1|Core), data=core_prof,
                      na.action = na.fail, REML=F)
mod.d13C_full
d13C_full_dredge <- dredge(mod.d13C_full)
print(d13C_full_dredge)
d13C_full_modsel <- as.data.frame(d13C_full_dredge)
write.xlsx(d13C_full_modsel, "tables/d13C_full_modsel.xlsx")

ggqqplot(resid(mod.d13C_full))
anova(mod.d13C_full)

#Best model (lowest AICc): pH regime + Water depth
mod.d13C_selected <- lmer(d13C_Interpolated ~ pH.regime + Water.depth + (1|Core), data=core_prof,
                               na.action = na.fail, REML=TRUE)#apply REML
summary(mod.d13C_selected)
simulateResiduals(mod.d13C_selected, plot=T) #Equal variances
anova(mod.d13C_selected)


#d13C common age (1954-2021)
mod.d13C_common <- gls(d13C_Interpolated ~ pH.regime + Water.depth + Site,
                     data = core_prof_common,
                     correlation = corCompSymm(form = ~1 | Core), 
                     weights = varIdent(form = ~1 | Core),
                     na.action = na.fail,
                     method="ML")
d13C_common_dredge <- dredge(mod.d13C_common)
print(d13C_common_dredge)
d13C_common_modsel <- as.data.frame(d13C_common_dredge)
write.xlsx(d13C_common_modsel, "tables/d13C_common_modsel.xlsx")

#Best model (lowest AICc): pH regime + Water depth
mod.d13C_common_gls_selected <- gls(IC_stock_sfc ~ pH.regime+Site, data = core_prof_common,
                                  correlation = corCompSymm(form = ~1 | Core),
                                  weights = varIdent(form = ~1 | Core), 
                                  na.action = na.fail,
                                  method="REML") #apply REML
summary(mod.d13C_common_gls_selected)
mod.d13C_common_gls_table <-  summary(mod.d13C_common_gls_selected)$tTable
write.xlsx(mod.d13C_common_gls_table, "tables/mod.d13C_common_gls_table.xlsx")
anova(mod.d13C_common_gls_selected)

#d15N 1954-2021
mod.d15N <- lmer(d15N_Interpolated~pH.regime+(1|Core), 
                 data=core_prof_common)
summary(mod.d15N)

simulateResiduals(mod.d15N, plot=T) #Equal variances
ggqqplot(resid(mod.d15N))
shapiro.test(resid(mod.d15N)) #normal residuals
anova(mod.d15N)


#C:N ratio 1954-2021
mod.CN <- lmer(C_N_molar~pH.regime+(1|Core), 
                 data=core_prof_common)
summary(mod.CN)

simulateResiduals(mod.CN, plot=T) #Equal variances
ggqqplot(resid(mod.CN))
shapiro.test(resid(mod.CN)) #normal residuals
anova(mod.CN)

#d13C common age (1954-2021)
mod.d13C_common <- gls(d13C_Interpolated ~ pH.regime + Water.depth + Site,
                       data = core_prof_common,
                       correlation = corCompSymm(form = ~1 | Core), 
                       weights = varIdent(form = ~1 | Core),
                       na.action = na.fail,
                       method="ML")
d13C_common_dredge <- dredge(mod.d13C_common)
print(d13C_common_dredge)
d13C_common_modsel <- as.data.frame(d13C_common_dredge)
write.xlsx(d13C_common_modsel, "tables/d13C_common_modsel.xlsx")



#### 12 Figure 3 ####

p.d13C_common <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                   y=d13C, fill=factor(pH.regime, levels=pH.order))) +
  #geom_violin(alpha=0.2, color="white") +
  geom_point(data=mean_d13C_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=3) +
  geom_errorbar(data=mean_d13C_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  ylab(bquote(delta^13*C~("\u2030"))) +
  xlab("") +
  ylim(-20,-16.5) +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")

#d15N
p.d15N_common <-ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                  y=d15N.Interpolated, fill=factor(pH.regime, levels=pH.order))) +
  geom_point(data=mean_d15N_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=3) +
  geom_errorbar(data=mean_d15N_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  ylab(bquote(delta^15*N~("\u2030"))) +
  ylim(2,6) +
  xlab("") +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")

#C:N ratio
p.C_N_common <-ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                 y=C_N_molar, fill=factor(pH.regime, levels=pH.order))) +
  geom_point(data=mean_C_N_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=3) +
  geom_errorbar(data=mean_C_N_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  ylab("C:N ratio") +
  xlab("") +
  ylim(10,30) +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")

p.source.d13C <- ggplot(sources, aes(x=factor(pH.regime, levels=pH.order), y=mean_d13C, 
                                     color=factor(pH.regime, levels=pH.order), 
                                     shape=Source.group,
                                     group=Source.group)) +
  geom_line(color="black", alpha=0.3, linetype="dashed") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_d13C-SE_d13C, ymax=mean_d13C+SE_d13C), width=0) + 
  scale_color_manual(values=c("#03A1D7", "#E09017", "#C00918"), name="pH regime") +
  scale_y_continuous(name=bquote(delta^13*C~("\u2030"))) +
  xlab("") +
  scale_shape_manual(values=c("triangle","diamond", "square")) +
  theme(legend.position = "none")
p.source.d15N <- ggplot(sources, aes(x=factor(pH.regime, pH.order), y=mean_d15N, 
                                     color=factor(pH.regime, levels=pH.order), 
                                     shape=Source.group, group=Source.group)) +
  geom_line(color="black", alpha=0.3, linetype="dashed") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_d15N-SE_d15N, ymax=mean_d15N+SE_d15N), width=0) +
  scale_y_continuous(name=bquote(delta^15*N~("\u2030"))) +
  scale_color_manual(values=c("#03A1D7", "#E09017", "#C00918"), name="pH regime") +
  scale_shape_manual(values=c("triangle","diamond", "square")) +
  xlab("") +
  theme(legend.position = "none")
p.source.CN <- ggplot(sources, aes(x=factor(pH.regime, pH.order), y=mean_CN, color=factor(pH.regime, levels=pH.order), 
                                   shape=Source.group, group=Source.group)) +
  geom_line(color="black", alpha=0.3, linetype="dashed") +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=mean_CN-SE_CN, ymax=mean_CN+SE_CN), width=0) +
  scale_y_continuous(name="C:N ratio") +
  scale_color_manual(values=c("#03A1D7", "#E09017", "#C00918"), name="pH regime") +
  scale_shape_manual(values=c("triangle","diamond", "square")) +
  xlab("") +
  theme(legend.position = "none")
(p.source.d13C + p.source.d15N + p.source.CN) /
  (p.d13C_common + p.d15N_common + p.C_N_common)



#### 13 Stable isotopes mixing model ####
#Calculate group means of sources
sources <- sources %>%
  group_by(Source.group) %>%
  summarise(
    mean_d13C = mean(mean_d13C, na.rm = TRUE), 
    SE_d13C = sqrt(sum(SE_d13C^2, na.rm = TRUE)) / n(), #pooled SE
    mean_d15N = mean(mean_d15N, na.rm = TRUE),
    SE_d15N = sqrt(sum(SE_d15N^2, na.rm = TRUE)) / n(), #pooled SE
    mean_CN = mean(mean_CN, na.rm = TRUE),
    SE_CN = sqrt(sum(SE_CN^2, na.rm = TRUE)) / n() #pooled SE
  )
#Extract sediment core d13C and d15N values 
compare_all <-  core_prof[c("Site", "Core", "pH.regime", "d13C_Interpolated","d15N_Interpolated","d13C","d15N")]
compare_all$pH.regime <- factor(compare_all$pH.regime, levels=pH.order)
compare_common <- core_prof_common[c("Site", "Core", "pH.regime", "d13C_Interpolated","d15N_Interpolated","d13C","d15N")]
compare_common$pH.regime <- factor(compare_common$pH.regime, levels=pH.order)
mixtures_all <- matrix(c(compare_all$d13C_Interpolated, compare_all$d15N_Interpolated), ncol=2)
colnames(mixtures_all) <- c("d13C_Interpolated", "d15N_Interpolated")
mixtures_common <- matrix(c(compare_common$d13C_Interpolated, compare_common$d15N_Interpolated), ncol=2)
colnames(mixtures_common) <- c("d13C_Interpolated", "d15N_Interpolated")

#Define sources
source_names <- sources$Source.group
source_means <- matrix(c(sources$mean_d13C, sources$mean_d15N), ncol=2)
source_sds <- matrix(c(sources$SE_d13C, sources$SE_d15N), ncol=2)
#Define group
group <- compare_all$pH.regime

### Run mixing model on all data (n=14 cores, data=core_prof) ###
simmr_in_all <- simmr_load(
  mixtures = mixtures_all,
  source_names = source_names,
  source_means = source_means,
  source_sds = source_sds,
  group=group
)

simmr_out_all <- simmr_mcmc(simmr_in_all)
summary(simmr_out_all, type = "diagnostics", group=c(1:3))
post_pred <- posterior_predictive(simmr_out_all)
prior_viz(simmr_out_all)
summary(simmr_out_all, type = "statistics", group=1:3)

### Check posterior probability
#Ambient pH
posterior_ambient <- simmr_out_all$output$Ambient$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_ambient) <- simmr_out_all$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_ambient <- mean(
  apply(posterior_ambient, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_ambient

#Low pH
posterior_low <- simmr_out_all$output$Low$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_low) <- simmr_out_all$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_low <- mean(
  apply(posterior_low, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_low

#Extreme Low pH
posterior_extlow <- simmr_out_all$output$`Extreme Low`$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_extlow) <- simmr_out_all$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_extlow <- mean(
  apply(posterior_extlow, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_extlow


### Run mixing model since 1954 (n=8 cores, data=core_prof_common) ###

group <- compare_common$pH.regime
#Run mixing model on common data (n=14 cores)
simmr_common <- simmr_load(
  mixtures = mixtures_common,
  source_names = source_names,
  source_means = source_means,
  source_sds = source_sds,
  group=group
)

simmr_out_common <- simmr_mcmc(simmr_common)
summary(simmr_out_common, type = "diagnostics", group=c(1:3))
post_pred <- posterior_predictive(simmr_out_common)
prior_viz(simmr_out_common)
summary(simmr_out_common, type = "statistics", group=1:3)

### Check posterior probability
#Ambient pH
posterior_ambient <- simmr_out_common$output$Ambient$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_ambient) <- simmr_out_common$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_ambient <- mean(
  apply(posterior_ambient, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_ambient

#Low pH
posterior_low <- simmr_out_common$output$Low$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_low) <- simmr_out_common$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_low <- mean(
  apply(posterior_low, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_low

#Extreme Low pH
posterior_extlow <- simmr_out_common$output$`Extreme Low`$BUGSoutput$sims.list$p
# Add source names
colnames(posterior_extlow) <- simmr_out_common$input$source_names

# Check probability that Posidonia is the top source at Ambient pH
posterior_prob_posidonia_extlow <- mean(
  apply(posterior_extlow, 1, function(row) row["Posidonia"] == max(row))
)
posterior_prob_posidonia_extlow



#### 14 Figure 4 ####
#All data (n=14 cores)
p.PO_all <- compare_groups(simmr_out_all,
                       source = "Posidonia",
                       groups = 1:3)
p.PO_all <- p.PO_all$plot + ggtitle("P. oceanica") + ylim(0,1) + theme_classic() + theme(legend.position="none")
p.MA_all <- compare_groups(simmr_out_all,
                       source = "Macroalgae",
                       groups = 1:3)
p.MA_all <- p.MA_all$plot + ylim(0,1) + ggtitle("Macroalgae") + ylim(0,1) +theme_classic() + theme(legend.position="none")
p.EP_all <- compare_groups(simmr_out_all,
                       source = "Epiphytes",
                       groups = 1:3)

p.EP_all <- p.EP_all$plot + ylim(0,1) + ggtitle("Epiphytes") + ylim(0,1) +theme_classic() + theme(legend.position="none")

#Since 1954 (n=8 cores)
p.PO_common <- compare_groups(simmr_out_common,
                           source = "Posidonia",
                           groups = 1:3)
p.PO_common <- p.PO_common$plot + ggtitle("P. oceanica") + ylim(0,1) + theme_classic() + theme(legend.position="none")
p.MA_common <- compare_groups(simmr_out_common,
                           source = "Macroalgae",
                           groups = 1:3)
p.MA_common <- p.MA_common$plot + ylim(0,1) + ggtitle("Macroalgae") + ylim(0,1) + theme_classic() + theme(legend.position="none")
p.EP_common <- compare_groups(simmr_out_common,
                           source = "Epiphytes",
                           groups = 1:3)

p.EP_common <- p.EP_common$plot + ylim(0,1) + ggtitle("Epiphytes") + ylim(0,1) +theme_classic() + theme(legend.position="none")

#Combine plots into one
(p.PO_all + p.MA_all + p.EP_all) /
  (p.PO_common + p.MA_common + p.EP_common)




          #### SUPPLEMENTARY MATERIAL ####

#### 15 Figure S1 ####
p.IC_stock_since1954 <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                          y=IC_stock_sfc, fill=factor(pH.regime, levels=pH.order))) +
  geom_violin(alpha=0.2, color="white") +
  geom_point(data=mean_IC_stock_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=4) +
  geom_errorbar(data=mean_IC_stock_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  geom_point(data=mean_IC_stock_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=median, color=factor(pH.regime, levels=pH.order)),
             size=3, shape=5) +
  ylab(bquote(IC ~stock ~ (g~ m^-2))) +
  xlab("pH regime") +
  ylim(0,4000) +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  theme(legend.position="none")
p.IC_burial_since1954 <- ggplot(data=core_prof_common, aes(x=factor(pH.regime, levels=pH.order), 
                                                           y=IC_burial_rate, fill=factor(pH.regime, levels=pH.order))) +
  geom_violin(alpha=0.3, color="white") +
  geom_point(data=mean_IC_burial_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=mean, color=factor(pH.regime, levels=pH.order)),
             size=4) +
  geom_errorbar(data=mean_IC_burial_pH_common, 
                aes(x=factor(pH.regime, levels=pH.order), y=mean, ymin=mean-se_propagated, ymax=mean+se_propagated,
                    color=factor(pH.regime, levels=pH.order)), width=0)+
  geom_point(data=mean_IC_burial_pH_common, aes(x=factor(pH.regime, levels=pH.order), y=median, color=factor(pH.regime, levels=pH.order)),
             size=3, shape=5) +
  ylab(bquote(IC ~burial~ rate ~ (g~ m^-2~yr^-1))) +
  xlab("pH regime") +
  ylim(0,80) +
  scale_color_manual(values=pH.colors, name="pH regime") +
  scale_fill_manual(values=pH.colors, name="pH regime") +
  #ggtitle("Between 1954 - today") +
  #geom_text(aes(
  #  y = 21,  # Adjust the y-value as needed
  #  label = paste0("n slices = ", n, "\n n cores = ", n_cores)
  #), data = sample_counts_common,
  #vjust = -0.5, size = 5) +
  theme(legend.position="none")
p.IC_stock_since1954 + p.IC_burial_since1954



#### 16 Figure S2 ####
#MAP of sampling sites using OSM

# Define bounding box for Ischia (Full map)
bbox_ischia <- c(13.85, 40.68, 14.00, 40.80)

# Fetch coastline data for full Ischia
ischia_coast <- opq(bbox = bbox_ischia) %>%
  add_osm_feature(key = "natural", value = "coastline") %>%
  osmdata_sf()

# Convert coastline lines into a polygon
ischia_poly <- st_polygonize(st_union(ischia_coast$osm_lines))
ischia_sf <- st_sf(geometry = ischia_poly)

# Define zoomed-in bounding box
bbox_zoom <- c(13.95, 40.70, 13.98, 40.75)

#Full Ischia map
main_map <- ggplot() +
  geom_rect(aes(xmin = 13.85, xmax = 14.00, ymin = 40.68, ymax = 40.77), fill = "lightblue") +  # Water
  geom_sf(data = ischia_sf, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(13.85, 14.00), ylim = c(40.68, 40.77)) +  # Adjusted zoom-in extent
  geom_point(data = sites, aes(x = Lon, y = Lat), color = "red", size = 3) +  # Sites
  geom_text(data = sites, aes(x = Lon, y = Lat, label = Site), hjust=-0.2, vjust=0, size=4) +  # Labels
  #geom_rect(aes(xmin = bbox_zoom[1], xmax = bbox_zoom[3], ymin = bbox_zoom[2], ymax = bbox_zoom[4]), 
  #          color = "black", fill = NA, linetype = "dashed", size = 1) +  # Zoomed-in box
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.4, "in"), pad_y = unit(2.7, "in"),
                         style = north_arrow_fancy_orienteering) 
main_map

# Convert bounding box to an sf polygon for cropping
zoom_bbox_sf <- st_as_sfc(st_bbox(c(xmin = bbox_zoom[1], xmax = bbox_zoom[3], 
                                    ymin = bbox_zoom[2], ymax = bbox_zoom[4]), 
                                  crs = st_crs(ischia_sf)))

# Crop Ischia polygon to only include the zoomed-in area
ischia_zoom_sf <- st_intersection(ischia_sf, zoom_bbox_sf)

zoom_map <- ggplot() +
  geom_rect(aes(xmin = bbox_zoom[1]-0.01, xmax = bbox_zoom[3]+0.01, 
                ymin = bbox_zoom[2]-0.01, ymax = bbox_zoom[4]+0.01), 
            fill = "lightblue") +  # Water
  geom_sf(data = ischia_zoom_sf, fill = "lightgray", color = "black") +  # Cropped land
  geom_point(data = sites, aes(x = Lon, y = Lat), color = "red", size = 4) +  # Sites
  geom_text(data = sites, aes(x = Lon, y = Lat, label = Site), hjust=-0.2, vjust=0, size=5) +
  annotation_scale(location = "br", width_hint = 0.5)


#### 17 Figure S3 ####
#Core profiles of Excess 210Pb
p.210Pb <- ggplot(core_prof, aes(x=Excess.210Pb..bulk., y=Accumulated.mass, Group=Core, color=factor(pH.regime, level = pH.order))) +
  geom_point(size=0.5) +
  scale_x_continuous(limits = c(0,300), breaks = c(0,150,300)) +
  geom_errorbar(aes(xmin=Excess.210Pb..bulk.-Error_Excess.210Pb..bulk., xmax=Excess.210Pb..bulk.+Error_Excess.210Pb..bulk.)) +
  geom_ribbon(data = core_prof %>% filter(Mixing == 1),
              aes(xmin=-Inf, xmax=Inf), fill = "grey", col=NA, alpha = 0.2) +
  scale_y_reverse(name="Accumulated mass (g cm-2)") +
  facet_wrap(.~factor(pH.regime,levels=pH.order)*Site*Core, ncol=7, scales="free_y") +
  scale_color_manual(values=pH.colors) +
  theme_classic() +
  theme(strip.background=element_rect(colour="white"), 
        strip.text.x = element_text(size = 0), legend.position = "none",
        axis.text = element_text(size=11, color="black"))
p.210Pb







############ END OF SCRIPT ##############
