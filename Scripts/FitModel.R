#GC 17/06/24
#Fit zero-inflated negative binomial model to dung beetle data


library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)

#Read in data ####


#Read in processed data (where time since intervention has already been scaled)
#This is the output from FormatDBforNegBinom.R 

sum_df<- data.table::fread("Outputs/full_DB_dataframeFor_BRMS_analysis_withSingletons.csv")

#NB 17.06.24 - NEED TO RERUN THIS WITH SINGLETONS AND DOUBLETONS REMOVED!!!
#sum_df<- data.table::fread("Outputs/full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv")

#last data formatting ####
sum_df <- sum_df %>%
  mutate(albizia = ifelse(habitat == "albizia", 1, 0),
         primary = ifelse(habitat == "primary", 1, 0),
         twice_logged = ifelse(habitat == "twice-logged", 1, 0),
         once_logged = ifelse(habitat == "once-logged", 1, 0),
         eucalyptus = ifelse(habitat == "eucalyptus", 1, 0),
         restored = ifelse(habitat == "restored", 1, 0),
         Year_factor = as.factor(sample_year))


# final model  ####

#NB can take >3 days to run on 8 cores
runtime <- system.time({
  full_model = brm(
    bf(
      sum_count ~ 0 +total_effort + albizia + primary + twice_logged + once_logged +
        eucalyptus + restored + albizia:time_since_intervention_ctr + 
        once_logged:time_since_intervention_ctr +
        eucalyptus:time_since_intervention_ctr + restored:time_since_intervention_ctr +
        time_since_intervention_ctr  +
        (1|Year_factor) + (1|site) +
        (0 + albizia + primary + twice_logged + once_logged +
           eucalyptus + restored + albizia:time_since_intervention_ctr + 
           once_logged:time_since_intervention_ctr +
           eucalyptus:time_since_intervention_ctr + restored:time_since_intervention_ctr +
           time_since_intervention_ctr|spp),
      zi ~ 0 + total_effort + albizia + primary + twice_logged + once_logged +
        eucalyptus + restored + albizia:time_since_intervention_ctr + 
        once_logged:time_since_intervention_ctr +
        eucalyptus:time_since_intervention_ctr + restored:time_since_intervention_ctr +
        time_since_intervention_ctr + total_effort +
        (1|Year_factor) + (1|site) +
        (0 + albizia + primary + twice_logged + once_logged +
           eucalyptus + restored + albizia:time_since_intervention_ctr + 
           once_logged:time_since_intervention_ctr +
           eucalyptus:time_since_intervention_ctr + restored:time_since_intervention_ctr +
           time_since_intervention_ctr|spp)),
    prior = c(set_prior('normal(0,1)', class = 'b'),
              set_prior('normal(0,1)', class = 'sd'),
              set_prior('normal(0,1)', class = 'b', dpar = "zi"),
              set_prior('normal(0,1)', class = 'sd', dpar = "zi")),
    data=sum_df,
    file = "Models/DB_zi_full4.rds",
    family = zero_inflated_negbinomial(),
    chains = 4, iter = 2500, warmup = 1250, 
    control = list(adapt_delta = 0.9),
    backend = 'cmdstanr', threads = threading(2),
    cores = 8, silent = 0
  )
})

# Model checks ####

#read
full_model <- readRDS("Models/DB_zi_full4.rds")

## Posterior predictive distribution check (dark and blue lines should look similar)
pp_check(full_model)
pp_check(full_model) + coord_cartesian(xlim = c(0, 100))

## Distributions mean check dark blue line is the mean of your data, historgram is the mean of
## the predicted distr (should be focsed around the dark blue line)
pp_check(full_model, type = "stat")

## Proportion of zeroes in real and simulated data - super important for these data 
## types where there is a bunch of zeros
prop_zero <- function(x) {
  sum(x == 0)/length(x)
  }

ppc_stat(y = sum_df$sum_count, yrep = posterior_predict(full_model, draws = 100), stat="prop_zero")

## Residual distribution (might take some time to run...) should just be centred on zero with 
## random uncertainty (+/-) that shows no trends or patterns.
sum_df %>%
  add_residual_draws(full_model) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval(alpha = .4, colour = "lightsteelblue3") + theme_classic() 


#### Average species count ####

## Get data for the range of data we have per hab and predict the counts
newdat <- sum_df %>% 
  group_by(habitat, time_since_intervention_ctr, time_since_intervention) %>%
  tally()


#come back to - doesnt currently run like this 
Ave_sp <- add_epred_draws(newdat, full_model, re_formula = NA, ndraws = 250) %>%
  group_by(habitat, time_since_intervention) %>%
  median_hdci(.epred, .width = .9)


ggplot(Ave_sp, aes(time_since_intervention, .epred, colour = habitat)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA) +
  facet_wrap(~habitat) +
  theme_minimal()

##  nb we don't actual predict across the full age range (ie. we dont extrapolate beyond the data; see ProcessScenarioOutcomes.R)
newdat2 <- expand.grid(habitat = unique(sum_df$habitat), time_since_intervention_ctr = 
                         unique(sum_df$time_since_intervention_ctr))

Ave_sp2 <- add_epred_draws(newdat2, simple_test_model13, re_formula = NA, ndraws = 250) %>%
  group_by(habitat, time_since_intervention_ctr) %>%
  median_hdci(.epred, .width = .9)

## It largely works, because we're extrapolating to never observed states
## the uncertainty is pretty big (but see above point)
ggplot(Ave_sp2, aes(time_since_intervention_ctr, .epred, colour = habitat)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA) +
  facet_wrap(~habitat, scales = "free") +
  theme_minimal()


#### Species-level ####

newdat3 <- sum_df %>% group_by(habitat, time_since_intervention_ctr, time_since_intervention, spp) %>% tally()

Ave_sp3 <- add_epred_draws(newdat3, full_model, re_formula = NULL, ndraws = 250) %>%
  group_by(spp, habitat, time_since_intervention_ctr, time_since_intervention) %>%
  median_hdci(.epred, .width = .9)

ggplot(Ave_sp3, aes(time_since_intervention, .epred, colour = habitat)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = NA) +
  facet_wrap(~spp, scales = "free") +
  theme_minimal()

#previous model runs ####

# #previous model 1 ####
# runtime <- system.time({
#   simple_test_model13 = brm(
#     bf(
#       sum_count ~ habitat*time_since_intervention_ctr +
#         (1+ habitat + habitat:time_since_intervention_ctr|spp),
#       # sample_year + (1|sample_year:spp) +
#       # species random intercepts for each habitat type
#       #(1+habitat||spp)+
#       #partial pooling of variance
#       #(1+ habitat + habitat:time_since_intervention_ctr||trait_index)+
#       # site random effect
#       #(1|site:spp),
#       zi ~ habitat*time_since_intervention_ctr +
#         (1+ habitat + habitat:time_since_intervention_ctr|spp)),
#     prior = c(set_prior('normal(0,1)', class = 'b'),
#               set_prior('normal(0,1)', class = 'Intercept'),
#               set_prior('normal(0,1)', class = 'sd'),
#               set_prior('normal(0,1)', class = 'b', dpar = "zi"),
#               set_prior('normal(0,1)', class = 'Intercept', dpar = "zi"),
#               set_prior('normal(0,1)', class = 'sd', dpar = "zi")),
#     data=sum_df,
#     family = zero_inflated_negbinomial(),
#     chains = 4, iter = 1500, warmup = 750, 
#     backend = 'cmdstanr', 
#     cores = 4, silent = 0
#     #threads = 2,
#   )
# })
# 
# saveRDS(simple_test_model13, "Models/DB_zi1.rds")

