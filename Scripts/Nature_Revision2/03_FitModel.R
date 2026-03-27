#GC 17/06/24
#Fit zero-inflated negative binomial model to dung beetle data

library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)

#---- NR2 path configuration ----
# NR2 INPUT: If running this script outside the project root, set this path explicitly.
project_root <- "."
nr2_root <- file.path(project_root, "Outputs", "NR2")
nr2_figures_dir <- file.path(nr2_root, "figures")
nr2_rds_dir <- file.path(nr2_root, "rds")
nr2_models_dir <- file.path(nr2_root, "models")

dir.create(nr2_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nr2_models_dir, recursive = TRUE, showWarnings = FALSE)

#Read in data ####


#Read in processed data (where time since intervention has already been scaled)
#This is the output from FormatDBforNegBinom.R 

#sum_df<- data.table::fread(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withSingletons.csv"))

sum_df<- data.table::fread(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv"))

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
    file = file.path(nr2_models_dir, "DB_zi_full.rds"),
    family = zero_inflated_negbinomial(),
    chains = 4, iter = 2500, warmup = 1250, 
    control = list(adapt_delta = 0.9),
    backend = 'cmdstanr', threads = threading(2),
    cores = 8, silent = 0
  )
})

# Model checks ####

#read
full_model <- readRDS(file.path(nr2_models_dir, "DB_zi_full.rds"))

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
