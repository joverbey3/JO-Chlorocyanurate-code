## ----setup, include=FALSE, echo=FALSE------------------------------------------------------------------------

#creation of libraries and setting the seed

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(purrr)

set.seed(1)


## ----impoft and plot, echo=FALSE-----------------------------------------------------------------------------
#in this section we will be importing the .csv files for the planktonic and the biofilm counts

#planktonic data
planktonic_data <- read.csv("Plate_counts.csv")

cat("This is how the planktonic culturable data is displayed.\n",
    "The data is organized by disinfectant type, temperature,and replicate number \n\n")
head(planktonic_data)


#uploading planktonic dpcr data
planktonic_dpcr <- read.csv("Planktonic_dPCR_Legionella.csv")

#biofilm data 

cat("This is how the biofilm culturable data alongside the dPCR data is displayed.\n",
    "The data is organized by disinfectant type, temperature,and replicate number \n\n")
biofilm_data <- read.csv("Biofilm_counts.csv")
head(biofilm_data)



## ----planktonic graphs, echo=FALSE, warning=FALSE------------------------------------------------------------
#setting up the summary for the planktonic graphs

#setting lod value
lod <- 2e-7


#replacing all 0 plate counts with LOD

planktonic_data$Concentration <- ifelse(planktonic_data$Concentration == 0, lod, planktonic_data$Concentration)

#need to convert concentration to log concentration

log_Concentration <- log(planktonic_data$Concentration)
planktonic_data$log_concentration <- log_Concentration

#head(planktonic_data)


#Calculating log removal

log_removal_data <- planktonic_data %>%
  group_by(Disinfectant.Type, Temperature, Replicate)%>%
  mutate(
    N0 = log_concentration[Time==0],
    log_removal = log_concentration-N0
    )%>%
  ungroup()



#organizing the data


log_removal_data1 <- log_removal_data %>%
  mutate(
    is_LOD = Concentration == lod,
    log_removal_num = ifelse(is_LOD, NA, log_removal),
    Temperature = factor(Temperature)
  )



##want to regress through all of the points so using lm function

group_models <- log_removal_data1 %>%
  filter(!is.na(log_removal_num)) %>%   # ignore LOD points
  group_by(Disinfectant.Type, Temperature) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log_removal_num ~ Time, data = .x)),
    tidied = map(model,tidy),
    glanced = map(model,glance),
    pred = map2(model, data, ~augment(.x, newdata = .y))
  ) 
  

group_models_pred <- group_models %>%
  select(Disinfectant.Type,Temperature,pred) %>%
  unnest(pred)

#creating a list of regression results



regression_summary <- group_models %>%
  mutate(
    slope = map_dbl(tidied, ~ .x$estimate[.x$term == "Time"]),
    #intercept = map_dbl(tidied, ~ .x$estimate[.x$term == "(Intercept)"]),
    r2 = map_dbl(glanced, "r.squared")
  ) %>%
  select(Disinfectant.Type, Temperature, slope, r2) %>%
  mutate(
    slope = round(slope, 3),
   # intercept = round(intercept, 3),
    r2 = round(r2, 3)
  ) %>%
  arrange(Disinfectant.Type, Temperature)


# what happens if i regress through zero instead of the points, how does that change the regression and fit of the model?



group_models2 <- log_removal_data1 %>%
  filter(!is.na(log_removal_num)) %>%   # ignore LOD points
  group_by(Disinfectant.Type, Temperature) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log_removal_num ~ Time+0, data = .x)),
    tidied = map(model,tidy),
    glanced = map(model,glance),
    pred = map2(model, data, ~augment(.x, newdata = .y))
  ) 
  

group_models_pred2 <- group_models2 %>%
  select(Disinfectant.Type,Temperature,pred) %>%
  unnest(pred)

#creating a list of regression results



regression_summary2 <- group_models2 %>%
  mutate(
    slope = map_dbl(tidied, ~ .x$estimate[.x$term == "Time"]),
    #intercept = map_dbl(tidied, ~ .x$estimate[.x$term == "(Intercept)"]),
    r2 = map_dbl(glanced, "r.squared")
  ) %>%
  select(Disinfectant.Type, Temperature, slope, r2) %>%
  mutate(
    slope = round(slope, 3),
   # intercept = round(intercept, 3),
    r2 = round(r2, 3)
  ) %>%
  arrange(Disinfectant.Type, Temperature)













## ----setup2, echo = FALSE,fig.cap= "Planktonic Legionella pneumophila counts. Plots are separated by the four disinfectant conditions: No Disinfectant, Monochlorocyanurate, Trichlorocyanurate, and Free Chlorine. All replicates are plotted, with LoD values in asterisks. Red is at 25 °C and blue is at 42 °C."----

plot1 <- ggplot() +
  
  # Non-LOD points
  geom_point(
    data = filter(log_removal_data1, !is_LOD),
    aes(x = Time, y = log_removal_num, color = as.factor(Temperature)),
    size = 2, alpha = 0.8
  ) +
  
  # LOD points as asterisks
  geom_text(
    data = filter(log_removal_data1, is_LOD),
    aes(x = Time, y = log_removal, label = "*", color = as.factor(Temperature)),
    size = 6,
    show.legend = FALSE
  ) +
  
  # Regression lines
  geom_line(
    data = group_models_pred,
    aes(x = Time, y = .fitted, color = as.factor(Temperature)),
    linewidth = 1
  ) +
  
  # Axis labels
  labs(
    y = expression(log(N[t]/N[0])),
    x = "Time (minutes)",
    color = "Temperature (°C)"
  ) +
  
  # Facets with nicer names
  facet_wrap(
    ~Disinfectant.Type,
    labeller = labeller(
      Disinfectant.Type = c(
        "Free_Chlorine" = "Free Chlorine",
        "Monochlor" = "Monochlorocyanurate",
        "No_Disinfectant" = "No Disinfectant",
        "Trichlor" = "Trichlorocyanurate"
      )
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

plot1
print(regression_summary)





## ----plot2, echo=FALSE,fig.cap= "Same as above but the points are forced through zero instead of regressed through the points"----

plot2 <- ggplot() +
  
  # Non-LOD points
  geom_point(
    data = filter(log_removal_data1, !is_LOD),
    aes(x = Time, y = log_removal_num, color = as.factor(Temperature)),
    size = 2, alpha = 0.8
  ) +
  
  # LOD points as asterisks
  geom_text(
    data = filter(log_removal_data1, is_LOD),
    aes(x = Time, y = log_removal, label = "*", color = as.factor(Temperature)),
    size = 6,
    show.legend = FALSE
  ) +
  
  # Regression lines
  geom_line(
    data = group_models_pred2,
    aes(x = Time, y = .fitted, color = as.factor(Temperature)),
    linewidth = 1
  ) +
  
  # Axis labels
  labs(
    y = expression(log(N[t]/N[0])),
    x = "Time (minutes)",
    color = "Temperature (°C)"
  ) +
  
  # Facets with nicer names
  facet_wrap(
    ~Disinfectant.Type,
    labeller = labeller(
      Disinfectant.Type = c(
        "Free_Chlorine" = "Free Chlorine",
        "Monochlor" = "Monochlorocyanurate",
        "No_Disinfectant" = "No Disinfectant",
        "Trichlor" = "Trichlorocyanurate"
      )
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

plot2
print(regression_summary2)


## ----plot 3, echo=FALSE, fig.cap="Same planktonic data but organized and faceted by temperature"-------------


log_removal_data1$Disinfectant.Type <- factor(log_removal_data1$Disinfectant.Type)
group_models_pred2$Disinfectant.Type <- factor(group_models_pred2$Disinfectant.Type)


plot3 <- ggplot() +
  
  # Non-LOD points
  geom_point(
    data = filter(log_removal_data1, !is_LOD),
    aes(x = Time, y = log_removal_num, color = Disinfectant.Type),
    size = 2, alpha = 0.8
  ) +
  
  # LOD points as asterisks
  geom_text(
    data = filter(log_removal_data1, is_LOD),
    aes(x = Time, y = log_removal, label = "*", color = Disinfectant.Type),
    size = 6,
    show.legend = FALSE
  ) +
  
  # Regression lines
  geom_line(
    data = group_models_pred2,
    aes(x = Time, y = .fitted, color = Disinfectant.Type),
    linewidth = 1
  ) +
  
  # Axis labels
  labs(
    y = expression(log(N[t]/N[0])),
    x = "Time (minutes)",
    color = "Disinfectant Type"
  ) +
  scale_color_brewer(palette = "Set2",
    labels = c(
    "No_Disinfectant" = "No Disinfectant",
    "Free_Chlorine"   = "Free Chlorine",
    "Monochlor"   = "Monochlorocyanurate",
    "Trichlor"  = "Trichlorocyanurate"
  )) +
  
  facet_wrap(~Temperature) +
  
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

plot3+theme(aspect.ratio=1)
#print(regression_summary2)


## ----ANOVA_planktonic, echo=FALSE----------------------------------------------------------------------------

###this computes the slopes per replicate rather than just one slope per condition

slopes_per_replicate <- log_removal_data1 %>%
  filter(!is.na(log_removal_num)) %>%
  group_by(Disinfectant.Type, Temperature, Replicate) %>%
  summarise(
    slope = coef(lm(log_removal_num ~ Time, data = cur_data_all()))[2],
    .groups = "drop"
  )


slopes_per_replicate$Disinfectant.Type <- factor(slopes_per_replicate$Disinfectant.Type)
slopes_per_replicate$Temperature <- factor(slopes_per_replicate$Temperature)

aov_slopes <- aov(slope ~ Disinfectant.Type * Temperature, data = slopes_per_replicate)

# ANOVA table
anova_table <- anova(aov_slopes)

anova_table

# Tukey post-hoc (by Temperature)
tukey_by_temp <- slopes_per_replicate %>%
  group_by(Temperature) %>%
  group_map(~ {
    model <- aov(slope ~ Disinfectant.Type, data = .x)
    broom::tidy(TukeyHSD(model)) %>%
      mutate(Temperature = unique(.x$Temperature))
  }, .keep = TRUE)

# Combine Tukey results
tukey_by_temp_df <- bind_rows(tukey_by_temp)



## ----tukey_table, echo= FALSE--------------------------------------------------------------------------------
tukey_by_temp_df


## ----ttest, echo=FALSE, warning=FALSE------------------------------------------------------------------------

### calculating the anova test which will compare the disinfectants against eachother at different temperatures


anova_by_disinfectant <- slopes_per_replicate %>%
  group_by(Disinfectant.Type) %>%
  group_modify(~ {
    # Run ANOVA for this disinfectant: slope ~ Temperature
    model <- aov(slope ~ Temperature, data = .x)
    broom::tidy(model)
  }) %>%
  ungroup()

anova_by_disinfectant


## ----tukey_by_disinfectant, echo=FALSE, warning=FALSE--------------------------------------------------------

###organizing the tukeys test by disinfectant

tukey_by_disinfectant <- slopes_per_replicate %>%
  group_by(Disinfectant.Type) %>%
  group_modify(~ {
    model <- aov(slope ~ Temperature, data = .x)
    broom::tidy(TukeyHSD(model))
  }) %>%
  ungroup()

tukey_by_disinfectant




## ----dpcr, echo=FALSE, fig.cap= "This plot shows the relationship between time and log reduction of DNA in gc/uL. All raw data points are plotted.The line shown is the line that is regressed through both replicates."----

extraction_factor <- 100 / 5000 #100 uL elutant / 5000 uL of input in the filter

planktonic_dpcr$conc_with_extraction <- planktonic_dpcr$Conc_undiluted_with_dilution_factor*extraction_factor

head(planktonic_dpcr)



# log concentration
planktonic_dpcr$log_conc <- log(planktonic_dpcr$conc_with_extraction)

# compute log removal 
dpcr_logremoval <- planktonic_dpcr %>%
  group_by(Disinfectant.Type, Temperature, Replicate) %>%
  mutate(
    N0 = log_conc[Time == 0],
    log_removal = log_conc - N0,
  ) %>%
  ungroup()

# convert factors
dpcr_logremoval$Temperature <- factor(dpcr_logremoval$Temperature)
dpcr_logremoval$Disinfectant.Type <- factor(dpcr_logremoval$Disinfectant.Type)



dpcr_models <- dpcr_logremoval %>%
  group_by(Disinfectant.Type, Temperature) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log_removal ~ Time, data = .x)),   # simple linear regression
    tidied = map(model, tidy),
    glanced = map(model, glance),
    pred = map2(model, data, ~ augment(.x, newdata = .y))
  )


regression_summary_dpcr <- dpcr_models %>%
  mutate(
    slope = map_dbl(tidied, ~ .x$estimate[.x$term == "Time"]),
    r2 = map_dbl(glanced, "r.squared")
  ) %>%
  select(Disinfectant.Type, Temperature, slope, r2) %>%
  mutate(
    slope = round(slope, 3),
    r2 = round(r2, 3)
  ) %>%
  arrange(Disinfectant.Type, Temperature)



dpcr_pred <- dpcr_models %>% unnest(pred)


### --- Plot dPCR log-removal curves --- ###

dpcr_plot <- ggplot() +
  
  # raw points
  geom_point(
    data = filter(dpcr_logremoval),
    aes(x = Time, y = log_removal, color = Temperature),
    size = 2, alpha = 0.8
  ) +
 
  geom_line(
    data = dpcr_pred,
    aes(x = Time, y = .fitted, color = Temperature),
    linewidth = 1
  ) +
  
  labs(
    title = "Planktonic dPCR Log Removal",
    x = "Time (minutes)",
    y = expression(log(N[t]/N[0])),
    color = "Temperature (°C)"
  ) +

  facet_wrap(
    ~ Disinfectant.Type,
    labeller = labeller(
      Disinfectant.Type = c(
        "Free_Chlorine" = "Free Chlorine",
        "Monochlor" = "Monochlorocyanurate",
        "No_Disinfectant" = "No Disinfectant",
        "Trichlor" = "Trichlorocyanurate"
      )
    )
  ) +

  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

dpcr_plot

regression_summary_dpcr



## ----pcr_data, echo=FALSE, fig.cap=" This is the same thing but forced through zero."------------------------


dpcr_models1 <- dpcr_logremoval %>%
  group_by(Disinfectant.Type, Temperature) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log_removal ~ Time+0, data = .x)),   # simple linear regression
    tidied = map(model, tidy),
    glanced = map(model, glance),
    pred = map2(model, data, ~ augment(.x, newdata = .y))
  )


regression_summary_dpcr1 <- dpcr_models1 %>%
  mutate(
    slope = map_dbl(tidied, ~ .x$estimate[.x$term == "Time"]),
    r2 = map_dbl(glanced, "r.squared")
  ) %>%
  select(Disinfectant.Type, Temperature, slope, r2) %>%
  mutate(
    slope = round(slope, 3),
    r2 = round(r2, 3)
  ) %>%
  arrange(Disinfectant.Type, Temperature)



dpcr_pred <- dpcr_models1 %>% unnest(pred)


### --- Plot dPCR log-removal curves --- ###

dpcr_plot1 <- ggplot() +
  
  # raw points
  geom_point(
    data = filter(dpcr_logremoval),
    aes(x = Time, y = log_removal, color = Temperature),
    size = 2, alpha = 0.8
  ) +
 
  geom_line(
    data = dpcr_pred,
    aes(x = Time, y = .fitted, color = Temperature),
    linewidth = 1
  ) +
  
  labs(
    title = "Planktonic dPCR Log Removal",
    x = "Time (minutes)",
    y = expression(log(N[t]/N[0])),
    color = "Temperature (°C)"
  ) +

  facet_wrap(
    ~ Disinfectant.Type,
    labeller = labeller(
      Disinfectant.Type = c(
        "Free_Chlorine" = "Free Chlorine",
        "Monochlor" = "Monochlorocyanurate",
        "No_Disinfectant" = "No Disinfectant",
        "Trichlor" = "Trichlorocyanurate"
      )
    )
  ) +

  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

dpcr_plot1

regression_summary_dpcr1



## ----plot6, echo=FALSE, fig.cap="Same data but faceted by temperature. Points are colored by disinfectant type. Regression lines are forced through zero."----

# Ensure factors
dpcr_logremoval$Disinfectant.Type <- factor(dpcr_logremoval$Disinfectant.Type)
dpcr_logremoval$Temperature <- factor(dpcr_logremoval$Temperature)
dpcr_pred$Disinfectant.Type <- factor(dpcr_pred$Disinfectant.Type)
dpcr_pred$Temperature <- factor(dpcr_pred$Temperature)


plot6 <- ggplot() +
  geom_point(
    data = dpcr_logremoval,
    aes(x = Time, y = log_removal, color = Disinfectant.Type),
    size = 2, alpha = 0.8
  ) +

  geom_line(
    data = dpcr_pred,
    aes(x = Time, y = .fitted, color = Disinfectant.Type),
    linewidth = 1
  ) +
  labs(
    y = expression(log(N[t]/N[0])),
    x = "Time (minutes)",
    color = "Disinfectant Type"
  ) +

  scale_color_brewer(palette = "Set2",
    labels = c(
      "No_Disinfectant" = "No Disinfectant",
      "Free_Chlorine"   = "Free Chlorine",
      "Monochlor"       = "Monochlorocyanurate",
      "Trichlor"        = "Trichlorocyanurate"
    )) +

  facet_wrap(~Temperature) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

plot6 + theme(aspect.ratio = 1)



## ----pcr stat analysis---------------------------------------------------------------------------------------

dpcr_logremoval$Temperature <- factor(dpcr_logremoval$Temperature)
dpcr_logremoval$Disinfectant.Type <- factor(dpcr_logremoval$Disinfectant.Type)

##### doing the slopes per replicate
slopes_per_replicate_dpcr <- dpcr_logremoval %>%
  group_by(Disinfectant.Type, Temperature, Replicate) %>%
  summarise(
    slope = coef(lm(log_removal ~ Time, data = cur_data_all()))[2],
    .groups = "drop"
  )

###running an anova on this
aov_slopes_dpcr <- aov(slope ~ Disinfectant.Type * Temperature, data = slopes_per_replicate_dpcr)
anova_table_dpcr <- anova(aov_slopes_dpcr)
anova_table_dpcr

# Tukey post-hoc comparisons by temperature
tukey_by_temp_dpcr <- slopes_per_replicate_dpcr %>%
  group_by(Temperature) %>%
  group_map(~ {
    model <- aov(slope ~ Disinfectant.Type, data = .x)
    broom::tidy(TukeyHSD(model)) %>%
      mutate(Temperature = unique(.x$Temperature))
  }, .keep = TRUE) %>%
  bind_rows()

tukey_by_temp_dpcr

tukey_by_disinfectant_dpcr <- slopes_per_replicate_dpcr %>%
  group_by(Disinfectant.Type) %>%
  group_map(~ {
    model <- aov(slope ~ Temperature, data = .x)
    
    broom::tidy(TukeyHSD(model)) %>%
      mutate(Disinfectant.Type = unique(.x$Disinfectant.Type))   # <- add disinfectant label
  }, .keep = TRUE) %>%
  bind_rows()

tukey_by_disinfectant_dpcr



