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
beetles <- read.csv(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withSingletons.csv"))

# Without singletons and doubletons 
#beetles <- read.csv(file.path(nr2_rds_dir, "full_DB_dataframeFor_BRMS_analysis_withoutSingletonsAndDoubletons.csv"))

#read in model Zero-inflated negative binom ouput 
rmodel<- readRDS(file.path(nr2_models_dir, "DB_zi_full.rds"))

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

# Build a version with no plateau assumptions, then compare species-curve shapes side by side.
process_beetle_data_no_plateau <- function(x) {
  ages <- data.frame(functionalhabAge = seq(0, 60))

  md <- x %>%
    mutate(habitat = case_when(
      habitat == "eucalyptus" ~ "eucalyptus_current",
      TRUE ~ habitat
    ))

  # Keep original model trajectory (no early/late plateaus), but keep original domain rules.
  md <- md %>%
    filter(!(habitat %in% c("eucalyptus_current", "eucalyptus_improved") & functionalhabAge > 6))

  # Keep the original deforested construction so comparison isolates plateau choices.
  beetles_deforested <- md %>%
    filter(habitat %in% c("eucalyptus_current", "albizia_current"),
           functionalhabAge %in% c(0, 1, 2)) %>%
    select(-functionalhabAge, -habitat) %>%
    group_by(species, iteration) %>%
    mutate(abundance = mean(abundance)) %>%
    crossing(ages) %>%
    mutate(habitat = "deforested")

  md %>%
    bind_rows(beetles_deforested) %>%
    filter(functionalhabAge <= 60)
}

#for restored habitat for which we don't have data beyond 40 yrs, so plateau recovery here
process_beetle_data_late_only <- function(x) {
  ages <- data.frame(functionalhabAge = seq(0, 60))
  missing20yrsR <- data.frame(functionalhabAge = seq(40, 60))

  md <- x %>%
    mutate(habitat = case_when(
      habitat == "eucalyptus" ~ "eucalyptus_current",
      TRUE ~ habitat
    ))

  # Late plateau only for restored (no early plateau for once-logged/restored).
  beetles_R_40 <- md %>%
    filter(habitat == "restored", functionalhabAge == 40) %>%
    select(-functionalhabAge) %>%
    group_by(species, habitat, iteration) %>%
    crossing(missing20yrsR)

  md <- md %>%
    filter(!(habitat == "restored" & functionalhabAge > 40)) %>%
    bind_rows(beetles_R_40)

  md <- md %>%
    filter(!(habitat %in% c("eucalyptus_current", "eucalyptus_improved") & functionalhabAge > 6))

  beetles_deforested <- md %>%
    filter(habitat %in% c("eucalyptus_current", "albizia_current"),
           functionalhabAge %in% c(0, 1, 2)) %>%
    select(-functionalhabAge, -habitat) %>%
    group_by(species, iteration) %>%
    mutate(abundance = mean(abundance)) %>%
    crossing(ages) %>%
    mutate(habitat = "deforested")

  md %>%
    bind_rows(beetles_deforested) %>%
    filter(functionalhabAge <= 60)
}

processed_beetles_late_only <- process_beetle_data_late_only(modelDraws)
processed_beetles_no_plateau <- process_beetle_data_no_plateau(modelDraws)

# Use late-only plateau as the main downstream dataset in this script.
processed_beetles <- processed_beetles_late_only

saveRDS(processed_beetles_late_only, file.path(nr2_rds_dir, "processedOccBeetles_lateOnlyPlateau.rds"))
saveRDS(processed_beetles_no_plateau, file.path(nr2_rds_dir, "processedOccBeetles_noPlateau.rds"))
# 
# 
# #compare side by side plots of plateus
# summarise_curves_for_plot <- function(dat, method_label) {
#   dat %>%
#     filter(habitat %in% c("once-logged", "restored", "twice-logged", "eucalyptus_current", "albizia")) %>%
#     mutate(
#       habitat = case_when(
#         habitat == "once-logged" ~ "Once logged",
#         habitat == "restored" ~ "Restored",
#         habitat == "twice-logged" ~ "Twice logged",
#         habitat == "eucalyptus_current" ~ "Eucalyptus pellita",
#         habitat == "albizia" ~ "Albizia falcataria",
#         TRUE ~ habitat
#       )
#     ) %>%
#     group_by(species, habitat, functionalhabAge) %>%
#     summarise(mid = median(abundance), .groups = "drop") %>%
#     mutate(method = method_label) %>%
#     filter(functionalhabAge <= 60)
# }
# 
# primary_points_for_plot <- function(dat, method_label) {
#   primary_vals <- dat %>%
#     filter(habitat == "primary") %>%
#     group_by(species) %>%
#     summarise(mid = median(abundance), .groups = "drop")
# 
#   primary_vals %>%
#     crossing(habitat = c("Once logged", "Restored", "Twice logged", "Eucalyptus pellita", "Albizia falcataria")) %>%
#     mutate(
#       functionalhabAge = 0,
#       method = method_label
#     )
# }
# 
# build_transition_segments <- function(curve_df, primary_df) {
#   target_points <- curve_df %>%
#     filter(functionalhabAge > 0) %>%
#     group_by(species, habitat) %>%
#     summarise(
#       # Requested rule:
#       # - Albizia connectors go to year 2.
#       # - All other habitats connect to year 10.
#       target_age_requested = if_else(habitat == "Albizia falcataria", 2, 10),
#       target_age_actual = {
#         ages <- sort(unique(functionalhabAge))
#         requested <- unique(if_else(habitat == "Albizia falcataria", 2, 10))
#         if (requested %in% ages) {
#           requested
#         } else if (any(ages <= requested)) {
#           max(ages[ages <= requested])
#         } else {
#           min(ages)
#         }
#       },
#       .groups = "drop"
#     ) %>%
#     left_join(curve_df, by = c("species", "habitat", "target_age_actual" = "functionalhabAge")) %>%
#     rename(target_mid = mid)
# 
#   transition_seed <- primary_df %>%
#     select(species, habitat, primary_mid = mid) %>%
#     left_join(target_points %>% select(species, habitat, target_age_actual, target_mid),
#               by = c("species", "habitat")) %>%
#     filter(!is.na(target_mid))
# 
#   bind_rows(
#     transition_seed %>%
#       transmute(species, habitat, functionalhabAge = 0, mid = primary_mid),
#     transition_seed %>%
#       transmute(species, habitat, functionalhabAge = target_age_actual, mid = target_mid)
#   ) %>%
#     arrange(species, habitat, functionalhabAge)
# }
# 
# plot_curves_method <- function(curve_df, primary_df, transition_df, title_text) {
#   ggplot(curve_df %>% filter(functionalhabAge >= 10), aes(functionalhabAge, mid, group = species)) +
#     geom_line(alpha = 0.45, col = "grey15") +
#     geom_point(data = primary_df, aes(functionalhabAge, mid), alpha = 0.55, size = 0.7, col = "grey15") +
#     geom_line(data = transition_df, aes(functionalhabAge, mid, group = species),
#               linetype = "dotted", alpha = 0.45, col = "grey15") +
#     facet_wrap(~habitat, scales = "fixed", nrow = 2) +
#     scale_x_continuous(limits = c(0, 60), breaks = c(0, 10, 20, 30, 40, 50, 60)) +
#     theme_bw() +
#     theme(
#       strip.text = element_text(face = "bold"),
#       strip.background = element_blank(),
#       axis.text = element_text(colour = "black"),
#       panel.grid = element_blank(),
#       plot.title = element_text(face = "bold")
#     ) +
#     labs(title = title_text, x = "Functional habitat age", y = "Relative abundance")
# }
# 
# curves_no_plateau <- summarise_curves_for_plot(processed_beetles_no_plateau, "No plateau")
# curves_plateau_late <- summarise_curves_for_plot(processed_beetles_late_only, "Plateau (late only)")
# 
# primary_no_plateau <- primary_points_for_plot(processed_beetles_no_plateau, "No plateau")
# primary_plateau_late <- primary_points_for_plot(processed_beetles_late_only, "Plateau (late only)")
# 
# transition_no_plateau <- build_transition_segments(curves_no_plateau, primary_no_plateau)
# transition_plateau_late <- build_transition_segments(curves_plateau_late, primary_plateau_late)
# 
# p_no_plateau <- plot_curves_method(curves_no_plateau, primary_no_plateau, transition_no_plateau, "No plateau")
# p_plateau_late <- plot_curves_method(curves_plateau_late, primary_plateau_late, transition_plateau_late, "Plateau (late only)")
# 
# curve_shape_comparison <- cowplot::plot_grid(
#   p_no_plateau, p_plateau_late,
#   ncol = 2,
#   labels = c("A", "B"),
#   label_size = 12
# )
# 
# ggsave(
#   file.path(nr2_figures_dir, "beetle_curves_no_plateau.png"),
#   p_no_plateau,
#   width = 210, height = 210, units = "mm"
# )
# 
# ggsave(
#   file.path(nr2_figures_dir, "beetle_curves_late_only_plateau.png"),
#   p_plateau_late,
#   width = 210, height = 210, units = "mm"
# )
# 
# ggsave(
#   file.path(nr2_figures_dir, "beetle_curve_shape_comparison_noPlateau_vs_lateOnlyPlateau.png"),
#   curve_shape_comparison,
#   width = 420, height = 210, units = "mm"
# )


#PLOT curves for ms figure 

# ---- Summarise posterior draws ----
beetle_summ <- processed_beetles %>%
  group_by(species, habitat, functionalhabAge) %>%
  summarise(
    mid = mean(abundance),
    lwr = quantile(abundance, 0.1),
    upr = quantile(abundance, 0.9),
    .groups = "drop"
  )

#remove unrealistic predictions for a few taxa
beetle_summ <- beetle_summ %>%
  filter(!species %in% c("Un_matched _ nr. Caccobius sp. A_FE",
                         "Un_matched _ nr. Caccobius bawangensis_FE",
                         "Onthophagus.quasijohkii")) %>%
  filter(habitat %in% c("primary", "once-logged", "twice-logged", "restored", "eucalyptus_current", "albizia"))

# Build species-specific primary baseline once, and force the same baseline across all non-primary habitats.
primary_vals <- beetle_summ %>%
  filter(habitat == "primary") %>%
  select(species, primary_mid = mid) %>%
  distinct()

logging_pred_df <- beetle_summ %>%
  filter(habitat %in% c("once-logged", "restored", "twice-logged")) %>%
  rename(time_since_logging = functionalhabAge) %>%
  left_join(primary_vals, by = "species") %>%
  mutate(mid = if_else(time_since_logging == 0, primary_mid, mid)) %>%
  mutate(
    habitat = case_when(
      habitat == "once-logged" ~ "Once logged",
      habitat == "restored" ~ "Restored",
      habitat == "twice-logged" ~ "Twice logged",
      TRUE ~ habitat
    )
  )

plantation_pred_df <- beetle_summ %>%
  filter(habitat %in% c("eucalyptus_current", "albizia")) %>%
  rename(plantation_age = functionalhabAge) %>%
  left_join(primary_vals, by = "species") %>%
  mutate(mid = if_else(plantation_age == 0, primary_mid, mid)) %>%
  mutate(
    habitat = case_when(
      habitat == "eucalyptus_current" ~ "Eucalyptus pellita",
      habitat == "albizia" ~ "Albizia falcataria",
      TRUE ~ habitat
    )
  )

logging_primary_points <- primary_vals %>%
  crossing(habitat = c("Once logged", "Restored", "Twice logged")) %>%
  mutate(time_since_logging = 0, mid = primary_mid)

plantation_primary_points <- primary_vals %>%
  crossing(habitat = c("Eucalyptus pellita", "Albizia falcataria")) %>%
  mutate(plantation_age = 0, mid = primary_mid)

plantation_transition_seed <- plantation_pred_df %>%
  filter(plantation_age == 1) %>%
  select(species, habitat, primary_mid, first_age = plantation_age, first_mid = mid) %>%
  distinct()

plantation_transition <- bind_rows(
  plantation_transition_seed %>%
    transmute(species, habitat, plantation_age = 0, mid = primary_mid),
  plantation_transition_seed %>%
    transmute(species, habitat, plantation_age = first_age, mid = first_mid)
) %>%
  arrange(species, habitat, plantation_age)

logging_transition_seed <- logging_pred_df %>%
  filter(habitat %in% c("Once logged", "Restored"), time_since_logging == 20) %>%
  select(species, habitat, primary_mid, first_age = time_since_logging, first_mid = mid) %>%
  distinct()

logging_transition <- bind_rows(
  logging_transition_seed %>%
    transmute(species, habitat, time_since_logging = 0, mid = primary_mid),
  logging_transition_seed %>%
    transmute(species, habitat, time_since_logging = first_age, mid = first_mid)
) %>%
  arrange(species, habitat, time_since_logging)

# ---- Plantation figure ----
plantation_fig <- plantation_pred_df %>%
  filter(plantation_age >= 1) %>%
  ggplot(aes(plantation_age, mid, group = species)) +
  geom_line(alpha = .5, col = 'grey0') +
  geom_point(data = plantation_primary_points,
             alpha = .5, col = 'grey0') +
  geom_line(data = plantation_transition,
            lty = 'longdash', alpha = .5, col = 'grey0') +
  facet_wrap(~habitat, nrow = 1) +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Relative abundance", x = "Plantation age") +
  ylim(0,30)+
  
  scale_x_continuous(breaks = c(0,5,10),
                     labels = c("Primary",5,10))

# ---- Logging figure ----
logging_fig <- logging_pred_df %>%
  filter((habitat %in% c("Once logged", "Restored") & time_since_logging >= 20) |
           (habitat == "Twice logged" & time_since_logging > 10)) %>%
  ggplot(aes(time_since_logging, mid, group = species)) +
  geom_line(alpha = .5, col = 'grey0') +
  geom_point(data = logging_primary_points,
             alpha = .5, col = 'grey0') +
  geom_line(data = logging_transition,
            lty = 'longdash', alpha = .3, col = 'grey0') +
  geom_line(data = logging_pred_df %>%
              filter(habitat == "Twice logged", time_since_logging <= 10),
            lty = 'longdash', alpha = .3, col = 'grey0') +
  facet_wrap(~habitat, nrow = 1) +
  theme_bw() +
  theme(strip.text = element_text(hjust = 0, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Abundance", x = "Time since logging") +
  ylim(0,30)+
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

# 
# #----plot processed data -----
# restored <- processed_beetles %>% filter(habitat == "restored")   #restored data seems a bit high; come back to
# oneL <- processed_beetles %>% filter(habitat == "once-logged")
# 
# fun_plot <- function(x){
#   x %>%  
#     filter(species %in% unique(species)[1:30]) %>%
#     
#     ggplot( aes(functionalhabAge, abundance, group = iteration)) +
#     geom_line(alpha = 0.2) + 
#     #  geom_ribbon(aes(ymin = lower_percentile, ymax = upper_percentile),fill = NA) +
#     facet_wrap(~species, scales = "free")
# }
# 
# fun_plot(restored)
# fun_plot(oneL)
# 
# #COME BACK HERE - I APPEAR TO BE GETTING SOME NON-LINEAR RELATIONSHIPS! 
# #-----calculate species categories --------
# 
# #remove improved improved yields
# # Filter out rows where 'habitat' contains the string "improved" (case-insensitive)
# processed_beetles <- as.data.table(processed_beetles)
# #processed_beetles <- processed_beetles[!grepl("improved", habitat, ignore.case = TRUE)]
# 
# #----summarise across posterior draws (for visualisation only) ----
# 
# #calculate 90% confidence intervals 
dbMeans <- processed_beetles[, .(abundance = median(abundance),
                                 ab_lwr = quantile(abundance, 0.2),
                                 ab_upr = quantile(abundance, 0.8)),
                             by = .(species, habitat, functionalhabAge)]

# # # Top 5 rows per species by abundance across habitat and age.
# top5_abundance_rows_by_species <- dbMeans %>%
#   group_by(species) %>%
#   slice_max(order_by = abundance, n = 5, with_ties = FALSE) %>%
#   ungroup()
# 

#-----calculate species categories ----- 
dbMeans_for_categories <- dbMeans %>%
  filter(
    #only consider the for restored and once-logged the places where we have sufficient data when defining species habitat affiliations
    habitat != "restored"  | dplyr::between(functionalhabAge, 20, 40),
    habitat != "once-logged" | dplyr::between(functionalhabAge, 20, 60)
  )
unique(dbMeans_for_categories$species)
best_rows <- dbMeans_for_categories %>%
  group_by(species) %>%
  filter(abundance == max(abundance, na.rm = TRUE)) %>%
  ungroup()

losers <- best_rows %>%
  filter(habitat == "primary") %>%
  transmute(species, spp_category = "loser") %>%
  distinct()

intermediates1L <- best_rows %>%
  filter(habitat == "once-logged"|habitat == "restored") %>%
  transmute(species, spp_category = "intermediate1L") %>%
  distinct()

intermediates2L <- best_rows %>%
  filter(habitat == "twice-logged") %>%
  transmute(species, spp_category = "intermediate2L") %>%
  distinct()

winners <- best_rows %>%
  filter(habitat %in% c("albizia", "albizia_current", "eucalyptus", "eucalyptus_current")) %>%
  transmute(species, spp_category = "winner") %>%
  distinct()


# 
# # Restrict age windows for loser classification:
# # - once-logged: 20-60 (ie don't use period where we don't have data)
# # - don't use restored to quantify losers; it occupies too small a habitat type at landscape scale
# # - all other habitats: all available ages
# dbMeans_for_loser <- dbMeans %>%
#   filter(
#     habitat != "restored",
#     habitat != "once-logged" | dplyr::between(functionalhabAge, 30, 60)
#   )
# 
# # # Original loser definition (all ages across habitats).
# # losers <- dbMeans_for_loser %>%
# #   group_by(species)  %>%
# #   filter(abundance == max(abundance, na.rm = TRUE)) %>%
# #   mutate(spp_category = case_when(habitat =="primary" ~"loser", TRUE ~ NA_character_)) %>%
# #   filter(spp_category == "loser") %>% select(species,spp_category) %>% unique()
# 
# intermediates1L <- dbMeans %>%
#   filter(functionalhabAge < 30) %>%  
#   group_by(species) %>% filter(abundance == max(abundance)) %>% 
#   mutate(spp_category = case_when(habitat =="once-logged" ~"intermediate1L",
#                                   habitat == "restored" ~ "intermediate1L",
#                                   TRUE ~ NA_character_)) %>% 
#   filter(spp_category == "intermediate1L") %>% select(species,spp_category) %>% unique
# 
# intermediates2L <- dbMeans %>%
#   filter(functionalhabAge < 30) %>%  
#   group_by(species) %>% filter(abundance == max(abundance)) %>% 
#   mutate(spp_category = case_when(habitat =="twice-logged" ~"intermediate2L",
#                                   TRUE ~ NA_character_)) %>% 
#   filter(spp_category == "intermediate2L") %>% select(species,spp_category) %>% unique
# 
# winners <- dbMeans %>%
#   filter(functionalhabAge < 30) %>%  
#   group_by(species) %>%
#   filter(abundance == max(abundance)) %>%
#   mutate(spp_category = case_when(
#     habitat == "albizia_current" ~ "winner",
#     habitat == "albizia" ~ "winner",
#     habitat == "eucalyptus_current" ~ "winner",
#     habitat == "eucalyptus" ~ "winner",
#     
#     TRUE ~ NA_character_
#   )) %>% 
#   filter(spp_category == "winner") %>% select(species,spp_category) %>% unique()
# 
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



