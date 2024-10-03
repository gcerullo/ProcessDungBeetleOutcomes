#GC 17/06/24 

#Process abundances ready for model fitting 

rm(list = ls())

#library(Rtools43)
library(brms)
library(cmdstanr)
library(tidyverse)
library(stringdist)
library(fuzzyjoin)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

#read in dung beetle data
df <- read.csv("Outputs/dungBeetlesForAbundanceAnalysis.csv")

#quick summaries of field sampling effort
total_trap_days <- df %>% select(site,trap,day) %>% unique %>%  group_by(site,trap) %>%  
  count() %>% mutate(traphrs = 24 *n)
sum(total_trap_days$traphrs)
sum(total_trap_days$n)


#quick additional summaries 
#num spp
df %>% select(spp) %>%  unique %>% count
#num individuals 
sum(df$abundance)
#num sites 
df %>% select(site) %>% unique %>% count
#num traps 
df %>%  select(trap,site, habitat) %>% unique() %>%  count
#traps per hab
df %>%  select(trap,site, habitat) %>% unique() %>%  group_by(habitat) %>% count() %>% ungroup
#colnames
names(df)

#numtraps collected using shorter sampling effort per trap  
df %>% filter(sampler %in% c("Slade_et al_2011", "SAFE_2018", "SAFE_2011", "SAFE_2015")) %>%  
  dplyr::select(trap) %>% unique %>% count

# Obsolete: add trait data  ####
# traits <- read.csv("Inputs/sppTraits.csv")  
# df <- df %>% left_join(traits, by = "spp")
# #body size (body length small: <10 mm; medium: 10–20 mm; and large: >20 mm)) and diel (diurnal/nocturnal)
# 

#format data ####

# (Slade data) has day = NA; basically it's collected after 48hrs, and only once 
df <- df %>% mutate(day_numeric = ifelse(is.na(day), 2, as.numeric(substr(day, 2,2))), #make eleanor data day 2 and make day numeric
                    time_since_intervention = coalesce(plantation_age,time_since_logging)) %>% # combine time since logging and plantation age as since column (note; it no longer therefore makes sense to include time_since logging as a fixed effectst - only in interaction with habitat type;)
  mutate(time_since_intervention = tidyr::replace_na(time_since_intervention, 0))

#get sum count 

#summarise trap count 
sum_df <- df %>% group_by(spp,site,transect, trap) %>% summarise(sum_count = sum(abundance),total_effort = max(day_numeric),
                                                                 habitat = unique(habitat),
                                                                 time_since_intervention = unique(time_since_intervention)                                                             sample_year = unique(sample_year)) %>% ungroup



# SAFE LFE doesn't have a time since logging; set as 20 years (roughly correct)
sum_df <- sum_df %>% 
  mutate(time_since_intervention = case_when(
    site == "SAFE LFE" ~ 20,
    TRUE ~ time_since_intervention
  ))


#deal with singletons and doubletons #### 

##54 species are either singletons or doubletons across all surveying effort! 
singleton_doubleton_spp <-  sum_df %>%
  filter(sum_count > 0) %>%
  group_by(spp) %>% 
  count() %>%
  arrange(n) %>%  
  filter(n < 3) %>%
  select(spp)

#what habitat were these species found in? 
singleton_doubleton_spp_HAB <- singleton_doubleton_spp %>%
  left_join(sum_df, by= "spp") %>%  
  filter(sum_count > 0 ) %>% 
  group_by(habitat) %>% count() 

#total sampling days per hab? 
trap_days <- sum_df %>%
  select(site, trap, total_effort, habitat) %>%
  group_by(site,trap, total_effort, habitat) %>% 
  distinct() %>% 
  group_by(habitat) %>%
  summarise(trap_days = sum(total_effort))

#Scale data ####

# We need to set the time_since_intervention for primary and twice-logged to 0 (otherwise the scale function spits up an NA, and stops the model running)
#The scale() function scales the variable by subtracting the mean and dividing by the standard deviation. However, if the standard deviation of the variable is 0 (which can happen when all values are the same, like in the case of 0), it will result in division by zero and produce NaN values.

sum_df <- sum_df %>%
  mutate(
    time_since_intervention_ctr = scale(time_since_intervention, scale = TRUE))

sum_df <- sum_df %>%  mutate(time_since_intervention_ctr = case_when(
  is.na(time_since_intervention_ctr) ~ 0,
  TRUE ~ time_since_intervention_ctr))

#make sample year categorical AND UNFROUP
sum_df$sample_year<- as.factor(sum_df$sample_year)

sum_df$spp %>% unique


#Export outputs #### 

# With singletons and doubletons included
write.csv(sum_df, "Outputs/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv")

#without singletons and doubletons 
df_filt <- sum_df %>% anti_join(singleton_doubleton_spp, by = "spp") 

write.csv(sum_df, "Outputs/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv")
