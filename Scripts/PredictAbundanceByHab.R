#Process model outputs - e.g. extract model abundance estiamtes per habitat type and age.

rm(list = ls())

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install_cmdstan() 
#install.packages("Rtools42")
#install_cmdstan() 
#install.packages("Rtools42")

#library(Rtools43)
library(brms)
library(cmdstanr)
library(tidyverse)
library(stringdist)
library(fuzzyjoin)

#check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
#check_cmdstan_toolchain()
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
#install_cmdstan(cores = 2)

setwd("C:/Users/Gianluca Cerullo/OneDrive - University of Cambridge/PhD/Chapter_4_Borneo/CompleteFolder/R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis")

#read in dung beetle data
df <- read.csv("Inputs/dungBeetlesForAbundanceAnalysis.csv")

total_trap_days <- df %>% select(site,trap,day) %>% unique %>%  group_by(site,trap) %>%  
  count() %>% mutate(traphrs = 24 *n)
sum(total_trap_days$traphrs)
sum(total_trap_days$n)
sum(df$abundance)
Slade <- df %>%filter(!(sampler %in% c("GC", "TL", "FE", "CF"))) 
unique(Slade$sampler)
Slade %>% select(site, trap) %>%  unique() %>% count()

#quick summaries 
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

#num slade traps 
df %>% filter(sampler %in% c("Slade_et al_2011", "SAFE_2018", "SAFE_2011", "SAFE_2015")) %>%  
  dplyr::select(trap) %>% unique %>% count


##================================================
#ADD TRAIT DATA 
##================================================
traits <- read.csv("Inputs/sppForTraits.csv")  
df <- df %>% left_join(traits, by = "spp")
#make a unique identifier for each unique combination of nesting habit (dweller, roller, tunnerler, arboreal tunneler), 
#body size (body length small: <10 mm; medium: 10–20 mm; and large: >20 mm)) and diel (diurnal/nocturnal)
#FOR NOW, WE DON'T CONSIDER DIEL. JUST NESTING HABIT AND BODY LENGTH 
df$trait_index <- paste(df$size, df$nest, sep = "_")

#================================================
#BOILER CODE STRUCTURE
#================================================

#20268 (Slade data) has day = NA; basically it's collected after 48hrs, and only once 
df <- df %>% mutate(day_numeric = ifelse(is.na(day), 2, as.numeric(substr(day, 2,2))), #make eleanor data day 2 and make day numeric
                    time_since_intervention = coalesce(plantation_age,time_since_logging)) %>% # combine time since logging and plantation age as since column (note; it no longer therefore makes sense to include time_since logging as a fixed effectst - only in interaction with habitat type;)
  mutate(time_since_intervention = tidyr::replace_na(time_since_intervention, 0))

#get sum count 

#summarise trap count 
sum_df <- df %>% group_by(spp,site,transect, trap) %>% summarise(sum_count = sum(abundance),total_effort = max(day_numeric),
                                                                 habitat = unique(habitat),
                                                                 time_since_intervention = unique(time_since_intervention),
                                                                 sample_year = unique(sample_year),
                                                                 trait_index = unique(trait_index)) %>% ungroup


#for twice-logged and primary, need to constrain to time-since-intervention 
#to between 3 SD of the normal distribution to stop infinit model runaway from 
#try to fit slope from one data point (0 years) 
p <- sum_df %>% filter(habitat == "primary")


#CENTRE TIME_SINCE_INTERVENTION by habitat
# SAFE LFE doesn't have a time since logging; set as 20 years (roughly correct)
sum_df <- sum_df %>% 
  mutate(time_since_intervention = case_when(
    site == "SAFE LFE" ~ 20,
    TRUE ~ time_since_intervention
  ))


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#remove rare species - e.g. singletons and doubletons
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

#plot number of sites that each species is detected at
sp_occ <-  sum_df %>% filter(sum_count > 0) %>% group_by(spp) %>% count() %>% arrange(n) %>% group_by(n) %>% count() %>%  
  rename(num_site = n, 
         detections_across_site = nn) 

sp_occ %>% 
  ggplot(aes(x = num_site, y = detections_across_site))+
  geom_bar(stat = "identity")

##54 species are either singletons or doubletons across all surveying effort! 
singleton_doubleton_spp <- sp_occ <-  sum_df %>% filter(sum_count > 0) %>% group_by(spp) %>% count() %>% arrange(n) %>%  
  filter(n < 3) %>% select(spp)

#what habitat were these species found in? 
singleton_doubleton_spp_HAB <- singleton_doubleton_spp %>% left_join(sum_df, by= "spp") %>%  
  filter(sum_count > 0 ) %>% 
  group_by(habitat) %>% count() 

#total sampling days per hab? 
trap_days <- sum_df %>% select(site, trap, total_effort, habitat) %>% group_by(site,trap, total_effort, habitat) %>% 
  distinct() %>% 
  group_by(habitat) %>% summarise(trap_days = sum(total_effort))

#effort corrected abundance 
singleton_doubleton_spp_HAB %>% left_join(trap_days, by = "habitat") %>%  mutate(effort_corrected_ab = n / trap_days *1000)


#REMOVE 54 SINGLETON AND DOUBLETON SPECIES
sum_df <- sum_df %>% anti_join(singleton_doubleton_spp, by = "spp") 
write.csv(sum_df, "raw_data_minus_singletonAndDoubletons")

#CENTRE SCALE DATA 


# We need to set the time_since_intervention for primary and twice-logged to the mean, which is 
#equivalent to setting it for zero (otherwise the scale function spits up an NA, and stops the model running)
#The scale() function scales the variable by subtracting the mean and dividing by the standard deviation. However, if the standard deviation of the variable is 0 (which can happen when all values are the same, like in the case of 0), it will result in division by zero and produce NaN values.

sum_df <- sum_df %>%
  mutate(
    time_since_intervention_ctr = scale(time_since_intervention, scale = TRUE))

mean_replacement_value <- sum_df$time_since_intervention_ctr %>% na.omit %>% mean
options(scipen = 999)

sum_df <- sum_df %>%  mutate(time_since_intervention_ctr = case_when(
  is.na(time_since_intervention_ctr) ~ mean_replacement_value,
  TRUE ~ time_since_intervention_ctr))


#make sample year categorical AND UNFROUP
sum_df$sample_year<- as.factor(sum_df$sample_year)

sum_df$spp %>% unique

#write.csv(sum_df, "Inputs/full_DB_dataframeFor_BRMS_analysis.csv")
write.csv(sum_df, "Inputs/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv")
dbb <- read.csv("Inputs/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv")

# Store the scaling information as attributes
attr(sum_df$time_since_intervention_ctr, "scaled:center") <- mean(sum_df$time_since_intervention, na.rm = TRUE)
attr(sum_df$time_since_intervention_ctr, "scaled:scale") <- sd(sum_df$time_since_intervention, na.rm = TRUE)


# Access the stored scaling information
mean_used_for_scaling <- attr(sum_df$time_since_intervention_ctr, "scaled:center")
sd_used_for_scaling <- attr(sum_df$time_since_intervention_ctr, "scaled:scale")




#MOST SIMPLE MODEL [CONVERGES]

#file.edit("~/.R/Makevars.win") 
#file.edit("~/.R/Makevars")
simple_test_model = brm(
  sum_count ~ total_effort + habitat:time_since_intervention_ctr + trait_index,
  prior = set_prior('normal(0,10)', class = 'b'),#this sets a NOrmal(0, sd=10) prior on all fixed effects, which is what's called a weakly regularizing/informative    data = sum_df,
  family = poisson,
  chains = 2,
  iter = 2000,
  backend = 'cmdstanr',
  cores = 2,
  data=sum_df)

saveRDS(simple_test_model, "test_brms_model_output1")
rmodel<- readRDS("test_brms_model_output1")
print(rmodel)
summary(rmodel)



#2ND MOST SIMPLE MODEL  [removes all SITE information] - performs better, nearly converging. But now we are assuming that the 
#points within sites are independent, which they are not, so I am not happy with this model ecologically. 

simple_test_model3 = brm(
  sum_count ~ total_effort + habitat:time_since_intervention_ctr + trait_index + 
    
    #allow abundance to be informed by trait index 
    habitat:trait_index+
    # species random intercepts for each habitat type
    (1+habitat|spp),
  
  prior = set_prior('normal(0,10)', class = 'b'),#this sets a NOrmal(0, sd=10) prior on all fixed effects, which is what's called a weakly regularizing/informative    data = sum_df,
  family = poisson,
  chains = 2,
  iter = 2000,
  backend = 'cmdstanr',
  cores = 2,
  data=sum_df)

saveRDS(simple_test_model3, "test_brms_model_output3")
rmodel<- readRDS("test_brms_model_output3")
print(rmodel)
summary(rmodel)


#3rd MOST SIMPLE MODEL (incorporates site as a random effect, runs more iterations and chains)
runtime <- system.time({
  simple_test_model5 = brm(
    sum_count ~ total_effort + habitat:time_since_intervention_ctr +
      
      #allow abundance to be informed by trait index 
      habitat:trait_index+
      # species random intercepts for each habitat type
      (1+habitat:trait_index|spp)+
      # site random effect 
      (1|site),
    prior = set_prior('normal(0,10)', class = 'b'),#this sets a NOrmal(0, sd=10) prior on all fixed effects, which is what's called a weakly regularizing/informative    data = sum_df,
    family = poisson,
    chains = 3,
    iter = 3000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 3,
    data=sum_df)
})
print(runtime)
saveRDS(simple_test_model5, "test_brms_model_output5")


#SIMON'S SUGGESTIONS 
#4TH MOST SIMPLE MODEL (incorporates site as a random effect, runs more iterations and chains)
runtime <- system.time({
  simple_test_model6 = brm(
    sum_count ~ total_effort + habitat:time_since_intervention_ctr +
      habitat +
      # species random intercepts for each habitat type
      (1+habitat|spp)+ 
      #partial pooling of variance 
      (1+habitat|trait_index) +
      # site random effect 
      (1|site),
    prior = set_prior('normal(0,3)', class = 'b'),
    family = poisson,
    chains = 3,
    iter = 3000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 3,
    data=sum_df)
})
print(runtime)
saveRDS(simple_test_model6, "test_brms_model_output6")
rmodel <- readRDS("test_brms_model_output6")
summary(rmodel)


# ADVICE FROM SIMON: Running it longer won't help with divergences. Divergences indicate that there is a bit of parameter space that the model is struggling to sample: more draws just results in more divergences. What I would do is get rid of the correlation terms in the random effects part of the model. To do this you write (habitat||spp) rather than (habitat|spp). Do this for both traits and spp. I'd also only bother with 1000 post-warmup draws. You might also try bumping adapt_delta to .9, but I would probably do this if you still have divergences after dropping the correlation terms


# 6th most complex model 
runtime <- system.time({
  simple_test_model7 = brm(
    sum_count ~ total_effort + habitat:time_since_intervention_ctr +
      habitat +
      # species random intercepts for each habitat type
      (1+habitat||spp)+ 
      #partial pooling of variance 
      (1+habitat||trait_index) +
      # site random effect 
      (1|site),
    prior = set_prior('normal(0,3)', class = 'b'),
    family = poisson,
    chains = 4,
    iter = 4000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 4,
    data=sum_df)
})

print(runtime)
#saveRDS(simple_test_model7, "test_brms_model_output7")
rmodel <- readRDS("test_brms_model_output7")
summary(rmodel)

#RERUN CONVERGIN MODEL WHILE REMOVING SINGLETON AND DOUBLETONS AND SCALING ACROSS HABITATS FOR TIME-SINE INTERVENTION (
#PREVIOUSLY I GROUPED_BY habitat BEFORE SCALING 
runtime <- system.time({
  simple_test_model8 = brm(
    sum_count ~ total_effort + habitat:time_since_intervention_ctr +
      habitat +
      # species random intercepts for each habitat type
      (1+habitat||spp)+ 
      # (1+habitat+habitat:time_since_intervention_ctr||spp)+ #option 2
      # (1+habitat+habitat:time_since_intervention_ctr||trait_index)+ #option 1 - 
      
      #partial pooling of variance 
      (1+habitat||trait_index) + #OPTION 1 REPLACES THIS LINE
      # site random effect 
      (1|site),
    prior = set_prior('normal(0,3)', class = 'b'),
    family = poisson,
    chains = 4,
    iter = 4000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 4,
    data=sum_df)
})

#allow each trait index group to vary in their post-logging recovery 
runtime <- system.time({
  simple_test_model9 = brm(
    sum_count ~ total_effort + habitat:time_since_intervention_ctr +
      habitat +
      # species random intercepts for each habitat type
      (1+habitat||spp)+ 
      # (1+habitat+habitat:time_since_intervention_ctr||spp)+ #option 2 - allows each spp age recovery trajectory to vary independently 
      # (1+habitat+habitat:time_since_intervention_ctr||trait_index)+ #option 1 - 
      
      #partial pooling of variance 
      (1+ habitat + habitat:time_since_intervention_ctr||trait_index)+
      # site random effect 
      (1|site),
    prior = set_prior('normal(0,3)', class = 'b'),
    family = poisson,
    chains = 4,
    iter = 4000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 4,
    data=sum_df)
})

saveRDS(simple_test_model9, "test_brms_model_output9")

#apply a zero-inflated negative-binomial model instead of a poisson distribution

# FULL zero-inflated NEGATIVE BINOMIal runs VERY SLOW
runtime <- system.time({
  simple_test_model11 = brm(
    bf(
      sum_count ~ total_effort + habitat:time_since_intervention_ctr +
        habitat +
        # sample_year + (1|sample_year:spp) +
        # species random intercepts for each habitat type
        (1+habitat||spp)+ 
        #partial pooling of variance 
        (1+ habitat + habitat:time_since_intervention_ctr||trait_index)+
        # site random effect 
        (1|site),
      #(1|site:spp),
      
      zi ~ 1 + total_effort + habitat:time_since_intervention_ctr +
        habitat +
        #sample_year + (1|sample_year:spp) +
        # species random intercepts for each habitat type
        (1+habitat||spp)+ 
        #partial pooling of variance 
        (1+ habitat + habitat:time_since_intervention_ctr||trait_index)+
        # site random effect 
        (1|site)),
    #    (1|site:spp)),
    prior = c(set_prior('normal(0,3)', class = 'b'),
              set_prior('normal(0,3)', class = 'Intercept'),
              set_prior('normal(0,3)', class = 'sd')),
    # family = poisson,
    #family = negbinomial()
    #family = zero_inflated_poisson(),
    family = zero_inflated_negbinomial(),
    chains = 4,
    iter = 000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 4,
    data=sum_df)
})

saveRDS(simple_test_model11, "test_brms_model_output11")

## Run a hopefully faster zero-inflated NEGATIVE BINOMIal -STILL SUPER SLOW
runtime <- system.time({
  simple_test_model13 = brm(
    bf(
      sum_count ~ total_effort + habitat:time_since_intervention_ctr +
        habitat +
        # sample_year + (1|sample_year:spp) +
        # species random intercepts for each habitat type
        (1+habitat||spp)+ 
        #partial pooling of variance 
        (1+ habitat + habitat:time_since_intervention_ctr||trait_index)+
        # site random effect 
        (1|site),
      #(1|site:spp),
      
      zi ~ 1 + time_since_intervention_ctr), 
    prior = c(set_prior('normal(0,3)', class = 'b'),
              set_prior('normal(0,3)', class = 'Intercept'),
              set_prior('normal(0,3)', class = 'sd')),
    # family = poisson,
    #family = negbinomial()
    #family = zero_inflated_poisson(),
    family = zero_inflated_negbinomial(),
    chains = 4,
    iter = 2000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 8,
    threads = 2,
    data=sum_df)
})

saveRDS(simple_test_model13, "test_brms_model_output13")


## RuN zero-inflated NEGATIVE BINOMIal, without any temporal effects  
runtime <- system.time({
  simple_test_model13 = brm(
    bf(
      sum_count ~ total_effort + 
        habitat +
        (1+habitat||spp)+ 
        #partial pooling of variance 
        (1+ habitat||trait_index)+
        (1|site),
      #(1|site:spp),
      
      zi ~ 1 + (1+habitat||spp)), 
    prior = c(set_prior('normal(0,3)', class = 'b'),
              set_prior('normal(0,3)', class = 'Intercept'),
              set_prior('normal(0,3)', class = 'sd')),
    # family = poisson,
    #family = negbinomial()
    #family = zero_inflated_poisson(),
    family = zero_inflated_negbinomial(),
    chains = 4,
    iter = 2000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 8,
    threads = 2,
    data=sum_df)
})

saveRDS(simple_test_model13, "test_brms_model_output13")


#run simple zero inflation model 
#check run-time with just a simple negative binomial 
runtime <- system.time({
  simple_test_model12 = brm(
    bf(
      sum_count ~ total_effort + habitat:time_since_intervention_ctr +
        habitat +
        # species random intercepts for each habitat type
        (1+habitat||spp)+
        #partial pooling of variance
        (1+ habitat + habitat:time_since_intervention_ctr|trait_index) +
        (1+ habitat + habitat:time_since_intervention_ctr|spp:trait_index) +
        # site random effect
        (1|site)),
    prior = c(set_prior('normal(0,3)', class = 'b'),
              set_prior('normal(0,3)', class = 'Intercept'),
              set_prior('normal(0,3)', class = 'sd')),
    family = negbinomial(),
    chains = 4,
    iter = 2000,
    warmup =1000,
    backend = 'cmdstanr',
    cores = 8,
    threads = threading(2),
    data=sum_df)
})



# 
# #allow each trait index group to vary in their post-logging recovery 
# runtime <- system.time({
#   simple_test_model10 = brm(
#     bf(
#     sum_count ~ total_effort + habitat +
#       # species random intercepts for each habitat type
#       (1+habitat||spp)+ 
#       #partial pooling of variance 
#       (1+ habitat||trait_index)+
#       # site random effect 
#       (1|site),
#     
#     zi ~ 1 + total_effort + habitat + # species random intercepts for each habitat type
#       (1+habitat||spp)+ 
#       #partial pooling of variance 
#       (1+ habitat||trait_index)+
#       # site random effect 
#       (1|site)),
#     prior = c(set_prior('normal(0,3)', class = 'b'),
#                   set_prior('normal(0,3)', class = 'Intercept'),
#                   set_prior('normal(0,3)', class = 'sd')),
#    # family = poisson,
#     #family = negbinomial()
#     #family = zero_inflated_poisson(),
#     family = zero_inflated_negbinomial(),
#     chains = 4,
#     iter = 4000,
#     warmup =1000,
#     backend = 'cmdstanr',
#     cores = 4,
#     data=sum_df)
# })


saveRDS(simple_test_model10, "test_brms_model_output10")
summary(simple_test_model10)

print(runtime)
saveRDS(simple_test_model8, "test_brms_model_output8")
rmodel <- readRDS("test_brms_model_output8")
summary(rmodel)


#MOST COMPLEX MODEL [DOESN'T CONVERGE]
runtime <- system.time({
  
  test_mothdel = brm(sum_count ~ total_effort + habitat + habitat:time_since_intervention +
                       sample_year + (1|sample_year:spp) + 
                       # site plus site x species intercepts
                       site + (1|site:spp) + 
                       #allow abundance to be informed by trait index 
                       trait_index + habitat:trait_index+
                       # species random intercepts for each habitat type
                       (1+habitat|spp),
                     data = sum_df, 
                     prior = set_prior('normal(0,10)', class = 'b'),
                     family = poisson,  chains = 2, backend = 'cmdstanr', iter = 2000, cores = 3)
})
print(runtime)

saveRDS(test_model, "test_brms_model_output2")
rmodel2 <- readRDS("test_brms_model_output2")
summary(rmodel2)
print(rmodel2)
#make slade data 2 days 



lm(sum_count = ~ sum_duaration, primary + # intercept is primary forest
     twice-logged + # limited variation in time since second logging, so not including this in model structure 
     once-logged + once-logged:time_since_logging + 
     restored + restored:time_since_logging +
     eucalyptus + eucalyptus:plantation_age + 
     albizia + albizia:plantation_age +
     # year plus year x species intercepts
     year + (1|sample_year:species) + 
     # site plus site x species intercepts
     site + (1|site:species) + 
     # species random intercepts for each habitat type
     (1 + twice-logged + once-logged + logged-restored + eucalyptus + 
        albizia|species) 
)


