#GC 17.06.24

#Assess the dung beetle outcomes of different scenarios, where each scenario is disaggregated by age

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)
library(foreach)
library(doParallel)
library(purrr)
library(profvis)
# install.packages("bayestestR", repos = "https://easystats.r-universe.dev")
library(bayestestR)

#------------------------------------------
# #DBspp correction PUSH UPSTREAM!! #####
# #Notes (these are all singletons in the dataset)
# #CaccobiusSp1_GC & Caccobius_sp1_GC should have been removed; not dung beetles 
# #Onthophagus_sp1_GC should be Onthophagus trituber - remove for now 
# #Onthophagus_sp3_GC should be Caccobius.bawangensis - rmove for now
# #Onthophagus_sp8_GC should be Onthophagus fujii Ochi & Kon - remove for now 
# #Onticellus_sp1_GC ; not a dung beetle - remove
# # DBs_incorrect <- DBs %>%
# #   filter(grepl("_GC", species, ignore.case = FALSE)) %>%  select(species)  %>% unique()
# 
# # Define the function to remove the specified species names
# remove_specific_species <- function(df) {
#   species_to_remove <- c("CaccobiusSp1_GC", "Caccobius_sp1_GC", 
#                          "Onthophagus_sp1_GC", "Onthophagus_sp3_GC", 
#                          "Onthophagus_sp8_GC", "Onticellus_sp1_GC")
#   
#   df %>%
#     filter(!species %in% species_to_remove)
# }
#------------------------------------------------

#Define inputs ####

#read in the scenario parametres containing conversion factors for converting from point to parcel/entire landscape  
source("Inputs/ScenarioParams.R")


#read in DB trap-level abundance for each posterior draw 
#this is calculated in PredictAbundanceByHab.R
DBs <- readRDS("Outputs/DBs_abundance_by_habAge_iterations.rds")
DBs <- remove_specific_species(DBs)


#-----read in scenarios without delays to get scenario composition -------
scenarios <- readRDS("Inputs/MasterAllScenarios.rds")
scenario_composition <- rbindlist(scenarios, use.names=TRUE) # get scenario composition
rm(scenarios)


#read in scenarios WITH delays, where every scenarioType is a single csv
#NB this is the same as MasterAllScenarios_withDelays.rds, except that all list elements are csvs

#get the csv file name for each scenario 
scenario_folder <- "Inputs/ScenariosWithDelaysCSVs"
csv_files <- list.files(scenario_folder, pattern = "*.csv", full.names = TRUE)


#define folder for saving, per scenario, final dung beetle outputs (time-averaged-abundance per scenario and species) 
final_output_folder <- "Outputs/FinalDBScenarioAbundSept24"

#remove eucalpyus and albizia improved if we don't need it 
processed_beetles <- DBs %>%
  filter(!(habitat == "eucalyptus_improved"|habitat == "albizia_improved")) %>%  
  #for the sake of analysis, call abundance as occ
  rename(occ = abundance)


#define all starting landscapes: 
all_start_landscape

#define total points of a landscape 
total_landscape_pts <- DB_CF*1000


# #=============  CALCULATE starting LANDSCAPE Abundance THRU TIME UNCERTAINTY ===========
#calculate the error about the summing of landscape_occ thru to give occ_60yrs across bootstraps 

species <- processed_beetles %>% select(species) %>%  unique() %>% pull()
all_start_landscape_scaled <- all_start_landscape %>%
  mutate(num_points = num_parcels*DB_CF) %>% select(-num_parcels) %>% 
  as.data.table()

unique_SL <- all_start_landscape_scaled %>% select(scenarioStart) %>% unique() %>% as.vector()
unique_spp <- processed_beetles %>% select(species) %>%  unique() %>% as.vector()
combinations_SL <- as.data.table(expand.grid(scenarioStart = unique_SL$scenarioStart, species = unique_spp$species))

function_SL_60yr_uncertainty <- function(single_scenario_i, processed_beetles_i) {
  
  beetle_join <- processed_beetles_i[single_scenario_i, on = .(habitat), allow.cartesian = TRUE]
 
   #calculate hab_occ
  result <- beetle_join[, hab_occ := occ * num_points]
  
  #[Across habitat type transitions (e.g for ALL hab_parcel transitions) in a scenario, calculate occupancy for a given year]
  result <- beetle_join[, .(landscape_occ = sum(hab_occ)/total_landscape_pts), 
                        by = .(species, scenarioStart, iteration, functionalhabAge)]
  
  # Step 4: Calculate occ_60yr for each iteration and species
  #[calculate occ60 for each iteration and species]
  result <- result[, .(occ_60yr = sum(landscape_occ)), 
                   by = .(species, scenarioStart, iteration)]
}

execute_SL_fun <-function(zeta) {
  single_scenario_name <- combinations_SL$scenarioStart[zeta]
  processed_beetles_name <- combinations_SL$species[zeta]
  
  single_scenario_i <- all_start_landscape_scaled[scenarioStart == single_scenario_name]
  processed_beetles_i <- processed_beetles[species == processed_beetles_name]
  
  result <- function_SL_60yr_uncertainty(single_scenario_i, processed_beetles_i)
  # Print the current iteration
  cat("Running iteration", zeta, "\n")
  
  return(result)
}

#apply the starting landscape function for each starting landscape and species in combinations SL
#do this to calculate occ_60yr for each species, scenarioStart, and posterior draw iteration 
result_list_SL <- lapply(1:nrow(combinations_SL), execute_SL_fun) 

saveRDS(result_list_SL,"Outputs/SL_occ60yr_perIterationJan25.rds")



# ---- Calculate the SCENARIO OCCUPANCY THRU TIME UNCERTAINTY ------


# # Define a function for calculating 0-60yr occupancy per species and scenario for each posterior draw

##....... Faster data-table function ........................

function_scenario_60yr_uncertainty <- function(single_scenario_i, processed_beetles_i) {
  # Step 1: Perform a left join and calculate hab_occ 
  # [join birds and scenarios and get num point for a given habitat type and stagger]
  
  
  beetle_join <- merge(single_scenario_i, processed_beetles_i, by.x = c("functional_habitat", "functionalhabAge"), 
                       by.y = c("habitat", "functionalhabAge"), all.x = TRUE, allow.cartesian = TRUE)
  beetle_join[, hab_occ := occ * num_points]
  beetle_join[, parcel_occ_stag := hab_occ / harvest_window]
  
  # Step 2: Group and summarize occupancy for specific year and habitat transition
  #[for each true year and habitat transition, calculate occupancy combined across the staggered
  # harvesting schedule (i.e. the occupancy in a given habitat transition for a given year)]
  
  result <- beetle_join[, .(occ_hab_year = sum(parcel_occ_stag)), 
                        by = .(species, index, iteration, production_target, true_year, original_habitat, habitat)]
  
  # Step 3: Calculate landscape occupancy
  #[Across habitat type transitions (e.g for ALL hab_parcel transitions) in a scenario, calculate occupancy for a given year]
  result <- result[, .(landscape_occ = sum(occ_hab_year)/total_landscape_pts), 
                   by = .(species, index, iteration, production_target, true_year)]
  
  # Step 4: Calculate occ_60yr for each iteration and species
  #[calculate occ60 for each iteration and species]
  result <- result[, .(occ_60yr = sum(landscape_occ)), 
                   by = .(species, index, iteration, production_target)]
  
  # # Step 5: Summarise across posterior draws
  # result <- result [, .(mean_60yr = mean(occ_60yr),
  #                       sd_60yr_error = sd(occ_60yr),
  #                       se_60yr_error = sd(occ_60yr) / sqrt(.N)),
  #                       by = .(index, species, production_target)]
  
  return(result)
}

#-----read in scenario group and define harvest delay ------
#get the csv file name for each scenario 
csv_files <- list.files(scenario_folder, pattern = "*.csv", full.names = TRUE)

#DEFINE DELAY FILTER ####
#(we have to subset only a few delay schedules to improve computational efficiency)
#delay filter (availalbe 0-29 in 1 year increments basically allow each scenario to be delayed in its first conversion by the delay filter)
delayFilters <- c("delay 0", "delay 29")
#delayFilters <- c("delay 0")
harvest_window <-  length(delayFilters)##how many harvest delays?

# Load necessary libraries for parallel processing
library(parallel)
library(doParallel)


#set a folder for saving outputs, showing for each species and scenario and iteration, occ_60 for lanscape
rds_folder <- "Outputs/Ab60perScenarioIterationJan25"

# Detect cores and set cluster (use detectCores()-1 to leave one core for system processes)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary libraries to cluster
clusterEvalQ(cl, {
  library(tidyr)
  library(data.table)
  library(dplyr)
  library(purrr)
})

for (k in seq_along(csv_files)){
  
  #read in single scenario type and add total points   
  scenario_group  <- read.csv(csv_files[[k]]) %>% as.data.table() %>% 
    mutate(num_points = num_parcels*DB_CF) %>%  
    filter(harvest_delay == delayFilters)  # only filter subset of delay to enhance efficieny 
  
  # Generate a unique rds file name
  rds_file_name <- sub(".csv$", "", unique(scenario_group$scenarioName))
  rds_file_name <- paste(unique(rds_file_name), "ab60.rds", sep = "_")
  
  # Combine folder path and file name to create full file path
  rds_file_path <- file.path(rds_folder, rds_file_name)   
  
  # Get unique scenarios and species
  unique_index <- scenario_group %>% select(index) %>% unique() %>% as.vector()
  unique_spp <- processed_beetles %>% select(species) %>% unique() %>% pull()
  combinations <- as.data.table(expand.grid(index = unique_index$index, species = unique_spp))
  
  # Export the scenario-specific variables to the cluster
  clusterExport(cl, c("combinations", "scenario_group", "processed_beetles", 
                      "function_scenario_60yr_uncertainty", "harvest_window", "total_bird_pts"))
  
  # Execute in parallel
  result_list <- foreach(zeta = 1:nrow(combinations), .packages = c('data.table', 'dplyr')) %dopar% {
    single_scenario_name <- combinations$index[zeta]
    processed_beetles_name <- combinations$species[zeta]
    
    single_scenario_i <- scenario_group[index == single_scenario_name]
    processed_beetles_i <- processed_beetles[species == processed_beetles_name]
    
    # Print the current iteration
    cat("Running iteration", zeta, "\n")
    
    result <- function_scenario_60yr_uncertainty(single_scenario_i, processed_beetles_i)
    return(result)
  }
  
  # Save the output to an rds file
  saveRDS(result_list, file = rds_file_path)
  
  # Clean up the result list
  rm(result_list)
  
  # Print the completion message for the current CSV file
  cat("Completed scenario", k, "\n")
}

# Stop the cluster after finishing all tasks
stopCluster(cl)


# #!!! change folder name here when running on desktop!!!
# #set a folder for saving outputs, showing for each species and scenario and iteration, occ_60 for lanscape
# rds_folder <- "Outputs/Ab60perScenarioIterationSept24"
# 
# for (k in seq_along(csv_files)){
#   
#   #read in single scenario type and add total points   
#   scenario_group  <- read.csv(csv_files[[k]]) %>% as.data.table() %>% 
#     mutate(num_points = num_parcels*DB_CF) %>%  
#     filter(harvest_delay == delayFilters)  # only filter subset of delay to enhance efficieny 
#   
#   # Generate a unique rds file name
#   rds_file_name <- sub(".csv$", "", unique(scenario_group$scenarioName))
#   rds_file_name <- paste(unique(rds_file_name), "ab60.rds", sep = "_")
#   
#   # Combine folder path and file name to create full file path
#   rds_file_path <- file.path(rds_folder, rds_file_name)   
#   
#   #=====
#   # # Choose a smaller number of indices and species for testing
#   #.............................................................
#   # # For example, let's use the first 5 indices and first 3 species
#   # selected_indices <- scenario_group %>% select(index) %>% unique %>%  slice(1:30) %>%  pull()
#   # selected_species <- processed_beetles %>% select(species) %>% unique %>%  slice(1:5) %>%  pull()
#   # # Create subsets of scenario_group and processed_beetles
#   # subset_scenario_group <- scenario_group %>% filter(index %in% selected_indices)
#   # subset_processed_beetles <- processed_beetles %>% filter(species %in% selected_species)
#   # 
#   # single_scenario_i <- as.data.table(subset_scenario_group)
#   # processed_beetles_i <- as.data.table(subset_processed_beetles)
#   # names(single_scenario2)
#   # names(processed_beetles2)
#   #=====
#   # 
#   
#   #----apply the function using parrellised approach  -----
#   
#   #calculate full set of combinatations of species and scenario
#   #get unique scenarios asn species 
#   unique_index <- scenario_group %>% select(index) %>% unique() %>% as.vector()
#   unique_spp <- processed_beetles %>% select(species) %>%  unique() %>% pull()
#   combinations <- as.data.table(expand.grid(index = unique_index$index, species = unique_spp))
#   
#   # Pre-allocate a list for results
#   result_list <- vector("list", nrow(combinations))
#   
#   
#   # Register a parallel backend
#   cl <- makeCluster(detectCores())
#   cl <- makeCluster(8)
#   registerDoParallel(cl)
#   clusterEvalQ(cl, c(library(tidyr),library(data.table),library(dplyr),library(purrr),library(profvis)))
#   clusterExport(cl,c("combinations","scenario_group","processed_beetles","function_scenario_60yr_uncertainty","harvest_window","total_bird_pts"))
#   
#   
#   # Use lapply to apply the function to each combination
#   #function works by taking the first row (defined by zeta) of combinations, filtering that given species and scenario, and then applying the function_scenario_60yr_uncertainty to these  
#   
#   execute_uncertainty_fun <-function(zeta) {
#     single_scenario_name <- combinations$index[zeta]
#     processed_beetles_name <- combinations$species[zeta]
#     
#     single_scenario_i <- scenario_group[index == single_scenario_name]
#     processed_beetles_i <- processed_beetles[species == processed_beetles_name]
#     
#     # Print the current iteration
#     cat("Running iteration", zeta, "\n")
#     
#     result <- function_scenario_60yr_uncertainty(single_scenario_i, processed_beetles_i)
#     
#     return(result)
#   }
#   
#   
#   timing <- system.time({
#     # result_list <- lapply(1:nrow(combinations), execute_uncertainty_fun)
#     
#     # result_list <- parLapplyLB(cl, 1:1000, execute_uncertainty_fun)
#     result_list <- parLapplyLB(cl, 1:nrow(combinations), execute_uncertainty_fun)
#   })
#   
#   cat("Elapsed time: ", timing[3], " seconds\n")
#   # Clean up parallel backend
#   stopCluster(cl)
#   
#   
#   #save the output to an rds folder 
#   saveRDS(result_list, file = rds_file_path)
#   
#   #saveRDS(result_list, "scenario60yrUncertainty1.rds")
#   
#   rm(result_list)
#   
#   # Print the execution time
#   cat("Elapsed time: ", timing[3], " seconds\n")
# }


#-----------------calculate median for each  iteration and species category ----
#cap <- 1.5 # don't allow scenario occ to be more than 1.5 starting landscape occ [only used if calculating geometric mean]
sppCategories <- readRDS("Outputs/DBsppCategories.rds")
sppCategories<- as.data.table(sppCategories)
abundance_folder <- "Outputs/Ab60perScenarioIterationJan25"
ab60_files <- list.files(abundance_folder, pattern = "*.rds", full.names = TRUE)


#scenario start abundance
#SL_60yrOcc 
SL_occ60 <- readRDS("Outputs/SL_occ60yr_perIterationJan25.rds")
SL_occ60_dt <- rbindlist(SL_occ60) %>%
  rename(SL_occ_60yr = occ_60yr)

#EXTRACT ONLY THE BASELINE ALL_PRIMARY SL
SL_all_primary_dt<- SL_occ60_dt %>% filter(scenarioStart == "all_primary") 

# Allocate folder for summarised results (geometric means and relative abundance )
relative_ab_folder <- "Outputs/RelativeAbundancePerIterationJan25"


for (w in seq_along(ab60_files)){
  occ60 <- readRDS(ab60_files[[w]])
  occ60_dt <- rbindlist(occ60)
  
  # Generate a unique  file name
  rds_file_name <- paste("OGbaseline_", basename(ab60_files[[w]]), sep = "")
  
  #rds_file_name <- paste("geometricMean_", basename(ab60_files[[w]]), sep = "")
  #relOcc_file_name <- paste("relOcc_", basename(occ60_files[[w]]), sep = "")
  
  
  # Combine folder path and file name to create full file path
  geom_file_path <- file.path(relative_ab_folder, rds_file_name)   
  #relOcc_file_path <- file.path(raw_rel_occ_folder, relOcc_file_name) 
  
  #add starting landscape to scenarios 
  scenarioStart <- occ60_dt %>% select(index, production_target) %>%
    unique() %>% left_join(scenario_composition, by = c("index", "production_target")) %>%  
    select(scenarioStart) %>% unique() %>% drop_na()
  occ60_dt[, scenarioStart := scenarioStart]
  
  #join SL and scenarios for each iteration
  #NB1: this conveys each species starting landscape occupancy
  #  occ_comb <- occ60_dt[SL_occ60_dt, on = .(species, scenarioStart), nomatch = 0]
  
  #NB2: this instead denotes each species old-growth baseline ouccupancy 
  occ_comb <- occ60_dt[SL_all_primary_dt, on = .(species, iteration), nomatch = 0]
  
  #calculate rel_occ; if rel_occ is > cap, replace with cap, to ensure scenario landscape cannot be more than 1.5 of starting landscape
  #occ_comb[, rel_occ_capped := pmin((occ_60yr / SL_occ_60yr), cap)]
  
  #calculate rel_occ -uncapped, for calculating median 
  
  occ_comb <- occ_comb[, rel_occ := occ_60yr / SL_occ_60yr]
  
  #export raw relative occupancy values 
  #saveRDS(occ_comb, file = "relOcc_file_path")
  
  #add in species categories 
  occ_comb <- sppCategories[occ_comb, on = "species"]
  
  #calculate geometric mean of each iteration of posterior draw. Thus we will end up with 500 geometric means 
  #per spp category,and scenario 
  
  # geom_means <- occ_comb %>% group_by(iteration, spp_category, index, production_target) %>%  
  #   summarise(geometric_mean = exp(mean(log(rel_occ_capped),na.rm = TRUE)), 
  #             medRelOcc = median(rel_occ,na.rm = TRUE)) %>% as.data.table()
  
  #summarise group median rel occ across 500 iterations 
  
  # Summarize species-level median rel occ across iterations
  rel_abs <- occ_comb %>%  group_by(species, index, production_target) %>%  
    summarize(
      medianRelativeOccupancy = median(rel_occ, na.rm = TRUE),
      meanRelativeOccupancy = mean(rel_occ, na.rm = TRUE),
      p1_medianRelativeOccupancy = quantile(rel_occ, 0.1, na.rm = TRUE),
      p9_medianRelativeOccupancy = quantile(rel_occ, 0.9, na.rm = TRUE),
    )
  
  
  # #add back in species categories if not involved in the grouping variable above  
  # rel_abs <- sppCategories[rel_abs, on = "species"]
  # 
  # Print the current iteration
  cat("Running iteration", rds_file_name, "\n")
  
  #save the output to an rds folder 
  saveRDS(rel_abs, file = geom_file_path)
}

#---

###################################################
#FOR UNCERTAINTY -calculate proportion of scenarios where logging is better than plantations
#set older for storing best scenario (logging or plantation) for each production target
best_scenario_folder <- "Outputs/BestScenarioUncertainty"

#folder storing abudances across 60yrs
ab60 <- "Outputs/Ab60perScenarioIterationJan25"
ab60files <- list.files(ab60, pattern = "*.rds", full.names = TRUE)

#scenario start abundance
#SL_60yrOcc 
SL_occ60 <- readRDS("Outputs/SL_occ60yr_perIterationJan25.rds")
SL_occ60_dt <- rbindlist(SL_occ60) %>%
  rename(SL_occ_60yr = occ_60yr)

#EXTRACT ONLY THE BASELINE ALL_PRIMARY SL
SL_all_primary_dt<- SL_occ60_dt %>% filter(scenarioStart == "all_primary") 

#how often are plantation scenarios better than logging scenarios 
logging_or_plantation_scenarios <- scenario_composition %>%  
  group_by(index, production_target) %>%  
  # Add information on proportion of plantation
  mutate(propPlant = sum(num_parcels[habitat %in% c("eucalyptus_current", "albizia_current", "albizia_future", "eucalyptus_future")]) / 1000) %>%  
  select(index, production_target, propPlant, scenarioStart) %>% 
  mutate(treatment_strategy = case_when(
    propPlant > 0 ~ "plantation",
    propPlant == 0 ~ "logging"
  )) %>%    select(-propPlant) %>%  unique() %>%  
  as.data.table()

for (w in seq_along(ab60files)){
  occ60 <- readRDS(ab60files[[w]])
  occ60_dt <- rbindlist(occ60)
  rds_file_name <- paste("BestScenario_", basename(ab60files[[w]]), sep = "")
  best_scenario_file_path <- file.path(best_scenario_folder, rds_file_name)   
  
  #add starting landscape to scenarios 
  scenarioStart <- occ60_dt %>% select(index, production_target) %>%
    unique() %>% left_join(scenario_composition, by = c("index", "production_target")) %>%  
    select(scenarioStart) %>% unique() %>% drop_na()
  occ60_dt[, scenarioStart := scenarioStart]
  
  #NB this conveys each species starting landscape occupancy
  occ_comb <- occ60_dt[SL_all_primary_dt, on = .(species, iteration), nomatch = 0]
  
  #calculate rel_occ
  occ_comb <- occ_comb[, rel_occ := occ_60yr / SL_occ_60yr]
  
  #add in scenario composition to highlight logging vs plantation scenarios
  occ_comb <- occ_comb %>% left_join(logging_or_plantation_scenarios)
  
  #for each production target and species, find the proportion of iterations where the best scenario 
  #(ie with the highest relOcc) is plantation-dominated. 
  #Do this for PAIRED scenario draws (ie iterations of the model)
  best_scenario <- occ_comb %>%  group_by(species, production_target, iteration) %>% 
    filter(rel_occ == max(rel_occ)) %>%  
    select(species, iteration, production_target, treatment_strategy, rel_occ) %>%  
    rename(max_rel_occ = rel_occ) %>% unique()
  
  #PUSH UPSTREAM (remove incorrect singleton species)
  best_scenario  <- best_scenario %>% remove_specific_species()
  
  #save the output to an rds folder 
  saveRDS(best_scenario, file = best_scenario_file_path)
}

#This would be an example of paired scenario draws (ie based on the same model params)
PairedExample <- occ_comb %>% filter(species == "Anoctus.laevis" & iteration == '459' & production_target == 0.39)

x %>% filter(rel_occ == max(rel_occ))

###################################################

#--------  read in relative abundance ----------------------
#read in data that has been baselined the fully old-growth starting landscape 
# Allocate folder for summarised results 
relative_ab_folder <- "Outputs/RelativeAbundancePerIterationJan25"
relAb_files <- list.files(relative_ab_folder, pattern = "^OGbaseline.*\\.rds$", full.names = TRUE)
rel_abs <- lapply(relAb_files, readRDS)
rel_occ_df <- rbindlist(rel_abs) %>%  
  left_join(sppCategories)

#----- summarise relative abundance across groups of species for which have median relative occupancy  -----

#summarised across data that has already been medianed per spp from 500 draws


#for species groupings (winner, loser, intermediate)
summarise_across_posterior_fun <- function(x){
  x %>% 
    remove_specific_species() %>% #move ustream 
    group_by(spp_category, index, production_target) %>%  
    summarise(medianRelativeOccupancy = median(medianRelativeOccupancy),
              p1_medianRelativeOccupancy = quantile(medianRelativeOccupancy, 0.1),
              p9_medianRelativeOccupancy = quantile(medianRelativeOccupancy, 0.9))
  
  
}
#---output of summarised statistics across posterior draws ----
final_relOcc <- summarise_across_posterior_fun(rel_occ_df)


#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()


output <- final_relOcc %>% select(index, production_target,
                                 medianRelativeOccupancy, p1_medianRelativeOccupancy, p9_medianRelativeOccupancy,
                                 spp_category) %>% cbind(outcome = "dungBeetles ")

#add back in scenarioNames
scenarioNames <- scenario_composition %>% select(index,scenarioStart,production_target, scenarioName) %>% unique()
output <- output %>% left_join(scenarioNames)

saveRDS(output, "FinalPerformanceOutput/MasterDBPerformance.rds")

