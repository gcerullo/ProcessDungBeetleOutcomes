# Gc 26.11.24
#Process model outputs - e.g. extract model abundance estiamtes per habitat type and age.

library(brms)
library(cmdstanr)
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

#Inputs ####

# read in raw abundance data  

# With singletons and doubletons included
beetles <- read.csv("Outputs/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv")

# Without singletons and doubletons 
#read.csv("Outputs/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv")

#read in model Zero-inflated negative binom ouput 
#rmodel<- readRDS("Models/DB_zi_full4.rds")
rmodel<- readRDS("Models/DB_zi_full_Nov24.rds")

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


# ---- visualise the model predictions per spp ------------

unique_spp <- beetles %>%ungroup %>%  select(spp) %>% unique() %>% slice(1:50) %>% pull()
test <- modelDraws %>% filter(species %in% unique_spp) %>%
  filter(!habitat == "restored") %>% 
  group_by(species, functionalhabAge,
           habitat) %>% summarise(mean = mean(abundance),
                                  lower_percentile = quantile(abundance, 0.05),
                                  upper_percentile = quantile(abundance, 0.95))# %>% filter(!habitat == "restored")
ggplot(test, aes(functionalhabAge, mean, colour = habitat)) +
  geom_line() + 
  #  geom_ribbon(aes(ymin = lower_percentile, ymax = upper_percentile),fill = NA) +
  facet_wrap(~species, scales = "free")

#------------ process draws to provide full data and platue unobserved years  -----------------------


#---- process beetles data ----

#define missing field data years where we want to replace model-derived estimates
ages <- data.frame(functionalhabAge = rep(seq(0,60))) 
#early years of once-logged and restored 
missing1L_R_rs <- data.frame(functionalhabAge = rep(seq(0,20))) 
#late years of restored
missing20yrsR <- data.frame(functionalhabAge = rep(seq(30,60))) 

names(modelDraws)

#function to process raw beetle data across posterior distributions
process_beetle_data <- function(x) {
  
  
  
  # Step 1: Add "improved" habitat (obsolete)
  modelDraws <- modelDraws %>%  
    mutate(habitat = case_when(
    #  habitat == "albizia" ~ "albizia_current",
      habitat == "eucalyptus" ~ "eucalyptus_current",
      TRUE ~ habitat))
  
  # improved <- modelDraws %>%
  #   filter(habitat %in% c("eucalyptus_current", "albizia_current")) %>%
  #   mutate(
  #     habitat = case_when(
  #       habitat == "eucalyptus_current" ~ "eucalyptus_improved",
  #       habitat == "albizia_current" ~ "albizia_improved",
  #       TRUE ~ habitat
  #     )
  #   )
  
 # modelDraws <- bind_rows(modelDraws, improved)
  
  # Step 3: Add missing years and assume 19 years of recovery for "once-logged" and "restored"
  
  #assume each year before functionalhabAge 19 has the same occupancy as yr 19 (i.e immediate recovery to yr 19 levels) 
  beetles_1L_19 <- modelDraws %>% 
    filter(habitat == "once-logged" & functionalhabAge == 19) %>% 
    select(-functionalhabAge) %>%
    group_by(species,habitat,iteration) %>% 
    crossing(missing1L_R_rs) 
  
  beetles_R_19 <- modelDraws %>% 
    filter(habitat == "restored" & functionalhabAge == 19) %>% 
    select(-functionalhabAge) %>%
    group_by(species,habitat,iteration) %>% 
    crossing(missing1L_R_rs) 
  
  #remove 1L and restored model-interpolated yrs and add yr 19 levels
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("once-logged", "restored") & functionalhabAge < 19 )) %>%
    bind_rows(beetles_1L_19) %>%  
    bind_rows(beetles_R_19)
  
  
  # Step 4: Assume no further recovery for "restored" habitat beyond ~40 years
  beetles_R_40 <- modelDraws %>% 
    filter(habitat == "restored" & functionalhabAge == 30) %>% 
    select(-functionalhabAge) %>%
    group_by(species,habitat,iteration) %>% 
    crossing(missing20yrsR) 
  
  #remove restored beyond 40 yrs and add in plateud data 
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("restored") & functionalhabAge > 30 )) %>%
    bind_rows(beetles_R_40) 
  
  # Step 5: Remove data for "eucalyptus" habitat beyond 6 years
  modelDraws <- modelDraws %>%
    filter(!(habitat %in% c("eucalyptus_current", "eucalyptus_improved") & functionalhabAge > 6))
  
  # Step 6: Interpolate deforested data from early plantation years
  beetles_deforested <-  modelDraws %>%
    filter((habitat %in% c("eucalyptus_current", "albizia_current")) & (functionalhabAge %in% c(0, 1, 2))) %>% 
    select(-c(functionalhabAge,habitat)) %>%
    group_by(species,iteration) %>% 
    mutate(abundance = mean(abundance)) %>%  
    crossing(ages) %>% 
    cbind(habitat = "deforested")
  
  modelDraws<- modelDraws %>% rbind(beetles_deforested) %>% 
    filter(functionalhabAge < 62)
  
  return(modelDraws)
}

# Call the function to process your 'beetles' data frame
processed_beetles <- process_beetle_data(modelDraws)

#----plot processed data -----
restored <- processed_beetles %>% filter(habitat == "restored")   #restored data seems a bit high; come back to
oneL <- processed_beetles %>% filter(habitat == "once-logged")

fun_plot <- function(x){
  x %>%  ggplot( aes(functionalhabAge, abundance, group = iteration)) +
    geom_line(alpha = 0.2) + 
    #  geom_ribbon(aes(ymin = lower_percentile, ymax = upper_percentile),fill = NA) +
    facet_wrap(~species, scales = "free")
}

fun_plot(restored)
fun_plot(oneL)
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
saveRDS(processed_beetles,"Outputs/DBs_abundance_by_habAge_iterations.rds")
saveRDS(sppCategories,"Outputs/DBsppCategories.rds")



