# Gc 26.11.24
#Process model outputs - e.g. extract model abundance estimates per habitat type and age.

library(brms)
# library(cmdstanr)
library(tidyverse)
library(stringdist)
library(fuzzyjoin)
library(tidybayes)
library(data.table)
library(posterior)
library(bayesplot)
library(data.table)

# Turn off scientific notation globally
options(scipen = 999)

#---- NR2 config ----
source(file.path("Scripts", "Nature_Revision2", "00_config.R"))

#Inputs ####

# read in raw abundance data  

# With singletons and doubletons included
#beetles <- read.csv(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withSingletons.csv"))

# Without singletons and doubletons 
beetles <- read.csv(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv"))

#read in model Zero-inflated negative binom ouput 
rmodel<- readRDS(file.path(nr2_models_dir, "DB_zi_full_Nov24.rds"))

summary(rmodel)
formula(rmodel)

#----plot the raw data ----
beetles %>%   mutate(habitat = case_when(
  habitat == "primary" ~ "P",
  habitat == "restored" ~ "SP",
  habitat == "twice-logged" ~ "2L",
  habitat == "once-logged" ~ "1L",
  habitat == "eucalyptus" ~ "EC",
  habitat == "albizia" ~ "A",
  TRUE ~ habitat  # Keep other values unchanged
)) %>% 
  ggplot(aes(habitat, sum_count, colour = habitat)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~spp, scales = "free") +
  theme(strip.background = element_blank())

beetles %>%
  group_by(habitat) %>%
  slice_max(time_since_intervention) %>% 
  select(habitat, time_since_intervention) %>% unique()


## Create a data frame with all combinations of species, habitat type, and habitat age I want to predict
time_since_intervention_prediction_points = seq(0, 60, by = 1) # years we want to predict for
time_since_intervention_ctr = (time_since_intervention_prediction_points-mean(beetles$time_since_intervention))/sd(beetles$time_since_intervention)
timesJoin <- cbind(time_since_intervention_prediction_points,time_since_intervention_ctr) %>% as.data.frame()

#for habitat types that don't change in time-since-intervention (primary and twice-logged)
#set to mean0val
mean0val <- beetles %>% filter(habitat=="twice-logged") %>% select(time_since_intervention_ctr) %>% unique() %>% pull

#note - we need to provide the scaled values for when we are making predictions
predict_hab <- expand.grid(
  spp = unique(beetles$spp),
  #habitat = c("albizia","once-logged","eucalyptus","restored","primary","twice-logged"),
  #this puts time_since_intervention_prediction_points onto the correct scale as model that was build on scaled data
  time_since_intervention_ctr = time_since_intervention_ctr,
  # for intervention 
  #trait_index = unique(beetles$trait_index),
  # site = unique(beetles$site),
  total_effort = 4) %>%     # assume all data was sampled over 4 days 
  # join a true age column
  left_join(timesJoin)  

#bring predict hab to the correct structure for extracting data from model 
predict_hab2 <- rbind(data.frame(predict_hab, albizia = 1, eucalyptus = 0,
                                 once_logged= 0, twice_logged= 0, primary= 0, restored = 0),
                      data.frame(predict_hab, albizia = 0, eucalyptus = 1,
                                 once_logged= 0, twice_logged= 0, primary= 0, restored = 0),
                      data.frame(predict_hab, albizia = 0, eucalyptus = 0,
                                 once_logged= 1, twice_logged= 0, primary= 0, restored = 0),
                      data.frame(predict_hab, albizia = 0, eucalyptus = 0,
                                 once_logged= 0, twice_logged= 1, primary= 0, restored = 0),
                      data.frame(predict_hab, albizia = 0, eucalyptus = 0,
                                 once_logged= 0, twice_logged= 0, primary= 1, restored = 0),
                      data.frame(predict_hab, albizia = 0, eucalyptus = 0,
                                 once_logged= 0, twice_logged= 0, primary= 0, restored = 1)) %>%  
  #filter data to remove ages we don't care about
  filter(
    !(eucalyptus == 1 & time_since_intervention_prediction_points > 6) &
      !(albizia ==1 & time_since_intervention_prediction_points > 12)) %>%
  
  #predict primary and twice-logged with 0 variation in time 
  mutate(time_since_intervention_ctr = case_when(
    primary == 1 | twice_logged == 1 ~ mean0val, TRUE~
      time_since_intervention_ctr))

#get 500 draws from the model predictions for a given spp, habitat type, and time since intervenioj 
formula(rmodel)

modelDraws <- add_epred_draws(rmodel, newdata = predict_hab2, re_formula = ~(
  0 + albizia + primary + twice_logged + once_logged + eucalyptus + restored + albizia:time_since_intervention_ctr + once_logged:time_since_intervention_ctr + eucalyptus:time_since_intervention_ctr + 
    restored:time_since_intervention_ctr + time_since_intervention_ctr | spp), se.fit = TRUE, type = "response",ndraws = 500)

head(modelDraws)

#bring model to the correct format for  analysis
modelDraws <- modelDraws %>% mutate(habitat = case_when(albizia == 1 ~ "albizia",
                                                        primary == 1 ~ "primary",
                                                        twice_logged == 1 ~ "twice-logged",
                                                        once_logged == 1 ~ "once-logged",
                                                        restored == 1 ~ "restored",
                                                        eucalyptus == 1 ~ "eucalyptus")) %>%
  ungroup() %>%
  select(-c(albizia, primary, twice_logged, once_logged, restored, eucalyptus))

#organise column names
modelDraws <- modelDraws %>% ungroup() %>%  select(spp, habitat, time_since_intervention_prediction_points, .draw, .epred) %>%
  rename(functionalhabAge = time_since_intervention_prediction_points,
         abundance = .epred, iteration = .draw,
         species = spp)

head(modelDraws) 

modelDraws <- as.data.table(modelDraws)



# 
# # ---- visualise the mean model predictions per spp ------------
# 
# unique_spp <- beetles %>%ungroup %>%  select(spp) %>% unique() %>% slice(1:50) %>% pull()
# test <- modelDraws %>% filter(species %in% unique_spp) %>%
#   filter(!habitat == "restored") %>% 
#   group_by(species, functionalhabAge,
#            habitat) %>% summarise(mean = mean(abundance),
#                                   lower_percentile = quantile(abundance, 0.05),
#                                   upper_percentile = quantile(abundance, 0.95))# %>% filter(!habitat == "restored")
# ggplot(test, aes(functionalhabAge, mean, colour = habitat)) +
#   geom_line() + 
#   #  geom_ribbon(aes(ymin = lower_percentile, ymax = upper_percentile),fill = NA) +
#   facet_wrap(~species, scales = "free")


#------------ process draws to provide full data and platue unobserved years  -----------------------


#---- process beetles data ----

# #define missing field data years where we want to replace model-derived estimates
# ages <- data.frame(functionalhabAge = rep(seq(0,60))) 
# #early years of once-logged and restored 
# missing1L_R_rs <- data.frame(functionalhabAge = rep(seq(0,20))) 
# #late years of restored
# missing20yrsR <- data.frame(functionalhabAge = rep(seq(41,60))) 
# 
# names(modelDraws)
# 
# 
process_beetle_data <- function(x) {
  
  # Define age grids
  ages            <- data.frame(functionalhabAge = seq(0, 60))
  missing1L_R_rs  <- data.frame(functionalhabAge = seq(0, 20))
  missing20yrsR   <- data.frame(functionalhabAge = seq(40, 60))
  
  modelDraws <- x  # keep everything inside function
  
  #--------------------------------------------------------------
  # Step 1: Normalise habitat names
  #--------------------------------------------------------------
  modelDraws <- modelDraws %>%  
    mutate(habitat = case_when(
      habitat == "eucalyptus" ~ "eucalyptus_current",
      TRUE ~ habitat
    ))
  
  #--------------------------------------------------------------
  # Step 2: Plateau early years (0–20) for once-logged & restored
  #--------------------------------------------------------------
  beetles_1L_19 <- modelDraws %>%
    filter(habitat == "once-logged", functionalhabAge == 20) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing1L_R_rs)
  
  beetles_R_19 <- modelDraws %>%
    filter(habitat == "restored", functionalhabAge == 20) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing1L_R_rs)
  
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("once-logged", "restored") &
               functionalhabAge < 20)) %>%
    bind_rows(beetles_1L_19, beetles_R_19)
  
  #--------------------------------------------------------------
  # Step 3: Plateau late restored years (40–60)
  #--------------------------------------------------------------
  beetles_R_40 <- modelDraws %>%
    filter(habitat == "restored", functionalhabAge == 40) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing20yrsR)
  
  modelDraws <- modelDraws %>%
    filter(!(habitat == "restored" & functionalhabAge > 40)) %>%
    bind_rows(beetles_R_40)
  
  #--------------------------------------------------------------
  # Step 4: Remove eucalyptus values >6 years
  #--------------------------------------------------------------
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("eucalyptus_current", "eucalyptus_improved") &
               functionalhabAge > 6))
  
  #--------------------------------------------------------------
  # Step 5: Define deforested (mean of plantation yrs 0–2)
  #--------------------------------------------------------------
  beetles_deforested <- modelDraws %>%
    filter(habitat %in% c("eucalyptus_current", "albizia_current"),
           functionalhabAge %in% c(0, 1, 2)) %>%
    select(-functionalhabAge, -habitat) %>%
    group_by(species, iteration) %>%
    mutate(abundance = mean(abundance)) %>%
    crossing(ages) %>%
    mutate(habitat = "deforested")
  
  modelDraws <- modelDraws %>%
    bind_rows(beetles_deforested) %>%
    filter(functionalhabAge <= 60)
  
  return(modelDraws)
}


# allow prediction up to yr 60 (ie beyond the data range (yr 40) for restored)

process_beetle_data_linear_20only <- function(x) {
  
  ages            <- data.frame(functionalhabAge = seq(0, 60))
  missing1L_R_rs  <- data.frame(functionalhabAge = seq(0, 20))
  
  modelDraws <- x
  
  # Normalize habitat names
  modelDraws <- modelDraws %>%
    mutate(habitat = case_when(
      habitat == "eucalyptus" ~ "eucalyptus_current",
      TRUE ~ habitat
    ))
  
  #--------------------------------------------------------------
  # Early plateau for 0–20
  #--------------------------------------------------------------
  beetles_1L_19 <- modelDraws %>%
    filter(habitat == "once-logged", functionalhabAge == 20) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing1L_R_rs)
  
  beetles_R_19 <- modelDraws %>%
    filter(habitat == "restored", functionalhabAge == 20) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing1L_R_rs)
  
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("once-logged", "restored") &
               functionalhabAge < 20)) %>%
    bind_rows(beetles_1L_19, beetles_R_19)
  
  #--------------------------------------------------------------
  # NO plateau for ages >40 (KEEP model estimates!)
  #--------------------------------------------------------------
  
  #--------------------------------------------------------------
  # Eucalyptus >6
  #--------------------------------------------------------------
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("eucalyptus_current", "eucalyptus_improved") &
               functionalhabAge > 6))
  
  #--------------------------------------------------------------
  # Deforested = mean of 0–2
  #--------------------------------------------------------------
  beetles_deforested <- modelDraws %>%
    filter(habitat %in% c("eucalyptus_current", "albizia_current"),
           functionalhabAge %in% c(0, 1, 2)) %>%
    select(-functionalhabAge, -habitat) %>%
    group_by(species, iteration) %>%
    mutate(abundance = mean(abundance)) %>%
    crossing(ages) %>%
    mutate(habitat = "deforested")
  
  modelDraws <- modelDraws %>%
    bind_rows(beetles_deforested) %>%
    filter(functionalhabAge <= 60)
  
  return(modelDraws)
}


processed_beetles     <- process_beetle_data(modelDraws)
processed_beetles_lin <- process_beetle_data_linear_20only(modelDraws)

saveRDS(processed_beetles, file.path(nr2_rds_dir, "processedOccBeetles_plateau.rds"))
saveRDS(processed_beetles_lin, file.path(nr2_rds_dir, "processedOccBeetles_no40yrplateau.rds"))


#PLOT PROCESSED DATA 

# ---- Summarise posterior draws ----
beetle_summ <- processed_beetles %>%
  group_by(species, habitat, functionalhabAge) %>%
  summarise(
    mid = mean(abundance),
    lwr = quantile(abundance, 0.1),
    upr = quantile(abundance, 0.9),
    .groups = "drop"
  )

# ---- Clean habitat labels to match bird figure ----
unique(beetle_summ$habitat)
beetle_summ <- beetle_summ %>%
  mutate(
    habitat = case_when(
      habitat == "primary" ~ "Primary",
      habitat == "once-logged" ~ "Once logged",
      habitat == "twice-logged" ~ "Twice logged",
      habitat == "restored" ~ "Restored",
      habitat == "eucalyptus_current" ~ "Eucalyptus pellita",
      habitat == "albizia" ~ "Albizia falcataria",
      TRUE ~ habitat
    )
  )

#remove unrealistically predictions for 	
#Un_matched _ nr. Caccobius sp. A_FE and Onthophagus.quasijohkii in restored forest 
beetle_summ <- beetle_summ %>%  filter(!species %in% c("Un_matched _ nr. Caccobius sp. A_FE",
                                                       "Un_matched _ nr. Caccobius bawangensis_FE",
                                                       "Onthophagus.quasijohkii"))

unique(beetle_summ$species)
# ---- Plantation dataframe ----
plantation_pred_df <- beetle_summ %>%
  filter(habitat %in% c("Eucalyptus pellita", "Albizia falcataria")) %>%
  rename(plantation_age = functionalhabAge)

# ---- Logging dataframe ----
logging_pred_df <- beetle_summ %>%
  filter(habitat %in% c("Once logged", "Restored")) %>%
  rename(time_since_logging = functionalhabAge)

# ---- Primary baseline points ----
primary_points <- beetle_summ %>%
  filter(habitat == "Primary") %>%
  select(species, mid) %>%
  distinct() %>%
  mutate(
    habitat = "Twice logged",
    time_since_logging = 0
  )

# ---- Twice logged trajectories ----
twice_log_seq <- seq(5, 60)

twice_logged <- beetle_summ %>%
  filter(habitat == "Twice logged") %>%
  select(species, mid) %>%
  crossing(time_since_logging = twice_log_seq) %>%
  mutate(habitat = "Twice logged")

twice_logged_full <- bind_rows(primary_points, twice_logged)

# ---- Combine logging data ----
plot_df <- bind_rows(logging_pred_df, twice_logged_full)

# ---- Plantation figure ----
plantation_fig <- plantation_pred_df %>%  
  filter(plantation_age >= 0) %>%
  ggplot(aes(plantation_age, mid, group = species)) +
  geom_line(alpha = .5, col = 'grey0') +
  geom_point(data = plantation_pred_df %>% filter(plantation_age == 0),
             alpha = .5, col = 'grey0') +
  geom_line(data = plantation_pred_df %>% filter(plantation_age <= 0),
            lty = 'longdash', alpha = .5, col = 'grey0') +
  facet_wrap(~habitat, nrow = 1) +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Relative abundance", x = "Plantation age") +
  scale_x_continuous(breaks = c(0,5,10),
                     labels = c("Primary",5,10))

# ---- Logging figure ----
logging_fig <- plot_df %>%  
  filter(time_since_logging >= 19 | habitat == "Twice logged") %>%
  ggplot(aes(time_since_logging, mid, group = species)) +
  geom_line(alpha = .5, col = 'grey0') +
  geom_point(data = plot_df %>% filter(time_since_logging == 0),
             alpha = .5, col = 'grey0') +
  geom_line(data = plot_df %>% filter(time_since_logging <= 21),
            lty = 'longdash', alpha = .3, col = 'grey0') +
  facet_wrap(~habitat, nrow = 1) +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Abundance", x = "Time since logging") +
  ylim(0,25)+
  # scale_x_continuous(breaks = c(0,40,60),
  #                    labels = c("Primary",20,40,60))
  scale_x_continuous(limits = c(0,40),
                    breaks = c(0,10,20,30,40),
                   labels = c("Primary",10,20,30, 40))

# ---- Combine panels ----
all_beetle_curves <- cowplot::plot_grid(logging_fig, plantation_fig,
                                        ncol = 1,
                                        rel_heights = c(1.2, 0.8))

ggsave(file.path(nr2_figures_dir, "all_beetle_curves_clipped_40yrs.png"),
       all_beetle_curves,
       units = "mm",
       height = 297,
       width = 210)


#----plot processed data -----
restored <- processed_beetles %>% filter(habitat == "restored")   #restored data seems a bit high; come back to
oneL <- processed_beetles %>% filter(habitat == "once-logged")

fun_plot <- function(x){
  x %>%  
    filter(species %in% unique(species)[1:30]) %>%
    
    ggplot( aes(functionalhabAge, abundance, group = iteration)) +
    geom_line(alpha = 0.2) + 
    #  geom_ribbon(aes(ymin = lower_percentile, ymax = upper_percentile),fill = NA) +
    facet_wrap(~species, scales = "free")
}

fun_plot(restored)
fun_plot(oneL)

#COME BACK HERE - I APPEAR TO BE GETTING SOME NON-LINEAR RELATIONSHIPS! 
#-----calculate species categories --------

#remove improved improved yields
# Filter out rows where 'habitat' contains the string "improved" (case-insensitive)
processed_beetles <- as.data.table(processed_beetles)
#processed_beetles <- processed_beetles[!grepl("improved", habitat, ignore.case = TRUE)]

#----summarise across posterior draws (for visualisation only) ----

#calculate 90% confidence intervals 
#take median because otherwise uncertainty around restored habitat swamps mean values 
dbMeans <- processed_beetles[, .(abundance = median(abundance), 
                                 ab_lwr = quantile(abundance, 0.2), 
                                 ab_upr = quantile(abundance, 0.8)),
                             by = .(species, habitat, functionalhabAge)]

#-----calculate species categories ----- 
CR <- dbMeans %>% filter(species == 'Catharsius.renaudpauliani')

losers <- dbMeans %>% 
  filter(functionalhabAge < 30) %>%  
  group_by(species)  %>%  
  filter(abundance == max(abundance)) %>% 
  mutate(spp_category = case_when(habitat =="primary" ~"loser", TRUE ~ NA_character_)) %>% 
  filter(spp_category == "loser") %>% select(species,spp_category) %>% unique()

losersCR <- CR %>% 
  filter(functionalhabAge < 30) %>%  
  filter(abundance == max(abundance))

intermediates1L <- dbMeans %>%
  filter(functionalhabAge < 30) %>%  
  group_by(species) %>% filter(abundance == max(abundance)) %>% 
  mutate(spp_category = case_when(habitat =="once-logged" ~"intermediate1L",
                                  habitat == "restored" ~ "intermediate1L",
                                  TRUE ~ NA_character_)) %>% 
  filter(spp_category == "intermediate1L") %>% select(species,spp_category) %>% unique

intermediates2L <- dbMeans %>%
  filter(functionalhabAge < 30) %>%  
  group_by(species) %>% filter(abundance == max(abundance)) %>% 
  mutate(spp_category = case_when(habitat =="twice-logged" ~"intermediate2L",
                                  TRUE ~ NA_character_)) %>% 
  filter(spp_category == "intermediate2L") %>% select(species,spp_category) %>% unique

winners <- dbMeans %>%
  filter(functionalhabAge < 30) %>%  
  group_by(species) %>%
  filter(abundance == max(abundance)) %>%
  mutate(spp_category = case_when(
    habitat == "albizia_current" ~ "winner",
    habitat == "albizia" ~ "winner",
    habitat == "eucalyptus_current" ~ "winner",
    habitat == "eucalyptus" ~ "winner",
    
    TRUE ~ NA_character_
  )) %>% 
  filter(spp_category == "winner") %>% select(species,spp_category) %>% unique()

#spp categories 
sppCategories <- rbind(losers,intermediates1L,intermediates2L,winners) %>% ungroup

#check all good 
all_species <- unique(dbMeans$species)

categorized_species <- unique(sppCategories$species)
species_not_categorized <- setdiff(all_species, categorized_species)
species_not_categorized

#----save outtputs ----
saveRDS(processed_beetles, file.path(nr2_rds_dir, "DBs_abundance_by_habAge_iterations.rds"))
saveRDS(sppCategories, file.path(nr2_rds_dir, "DBsppCategories.rds"))



