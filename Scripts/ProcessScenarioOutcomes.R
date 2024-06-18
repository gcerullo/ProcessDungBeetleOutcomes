#GC 17.06.24

#Assess the dung beetle outcomes of different scenarios, where each scenario is disaggregated by age

rm(list = ls())

library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(stringr) 
library(cowplot)
library(foreach)
library(doParallel)


#read in the scenario parametres containing conversion factors for converting from point to parcel/entire landscape  
source("Inputs/ScenarioParams.R")


#set a cap for relative occupancy so that the scenario occupancy of a species cannot be more than CAPtimes the occupancy of the starting landscape
cap <- 1  

#read in DB trap-level abundance for each posterior draw 
#this is calculated in PredictAbundanceByHab.R
DBs <- readRDS("Outputs/DBs_abundance_by_habAge_iterations.rds")

names(DBs)
#remove eucalpyus and albizia improved if we don't need it 
processed_beetles <- DBs %>%
  filter(!(habitat == "eucalyptus_improved"|habitat == "albizia_improved")) %>%  
  #for the sake of analysis, call abundance as occ
  rename(occ = abundance)

#read in folder where scenarios are stored as seperate CSVs
scenario_folder <- "R_code/AssessBiodiversityOutcomes/Outputs/scenariosForBirdsToBatch"

#define folder for saving, per scenario, final dung beetle outputs (time-averaged-abundance per scenario and species) 
final_output_folder <- "R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/FinalDBScenarioAbund"

#get scenario composition by reading in scenarios as an rds
#(not processed with code above to be temporal, #but allows us to get scenario composition efficiently 
scenarios <- readRDS("R_code/BuildScenarios/BuildingHarvestingScenarios/allScenariosStaggered.rds")
scenario_composition <- rbindlist(scenarios, use.names=TRUE) # get scenario composition
rm(scenarios) 

#get the csv file name for each scenario 
csv_files <- list.files(scenario_folder, pattern = "*.csv", full.names = TRUE)

#define all starting landscapes: 
all_start_landscape

#define total points of a landscape 
total_landscape_pts <- DB_CF*1000




# #=============  CALCULATE starting LANDSCAPE OCCUPANCY THRU TIME UNCERTAINTY ===========
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

saveRDS(result_list_SL,"R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/SL_occ60yr_perIteration.rds")



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

#!!! change folder name here when running on desktop!!!
#set a folder for saving outputs, showing for each species and scenario and iteration, occ_60 for lanscape
rds_folder <- "R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/Ab60perScenarioIteration"

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
  
  #=====
  # # Choose a smaller number of indices and species for testing
  #.............................................................
  # # For example, let's use the first 5 indices and first 3 species
  # selected_indices <- scenario_group %>% select(index) %>% unique %>%  slice(1:30) %>%  pull()
  # selected_species <- processed_beetles %>% select(species) %>% unique %>%  slice(1:5) %>%  pull()
  # # Create subsets of scenario_group and processed_beetles
  # subset_scenario_group <- scenario_group %>% filter(index %in% selected_indices)
  # subset_processed_beetles <- processed_beetles %>% filter(species %in% selected_species)
  # 
  # single_scenario_i <- as.data.table(subset_scenario_group)
  # processed_beetles_i <- as.data.table(subset_processed_beetles)
  # names(single_scenario2)
  # names(processed_beetles2)
  #=====
  # 
  
  #----apply the function using parrellised approach  -----
  
  #calculate full set of combinatations of species and scenario
  #get unique scenarios asn species 
  unique_index <- scenario_group %>% select(index) %>% unique() %>% as.vector()
  unique_spp <- processed_beetles %>% select(species) %>%  unique() %>% pull()
  combinations <- as.data.table(expand.grid(index = unique_index$index, species = unique_spp))
  
  # Pre-allocate a list for results
  result_list <- vector("list", nrow(combinations))
  
  
  # Register a parallel backend
  cl <- makeCluster(detectCores())
  cl <- makeCluster(8)
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(tidyr),library(data.table),library(dplyr),library(purrr),library(profvis)))
  clusterExport(cl,c("combinations","scenario_group","processed_beetles","function_scenario_60yr_uncertainty","harvest_window","total_bird_pts"))
  
  
  # Use lapply to apply the function to each combination
  #function works by taking the first row (defined by zeta) of combinations, filtering that given species and scenario, and then applying the function_scenario_60yr_uncertainty to these  
  
  execute_uncertainty_fun <-function(zeta) {
    single_scenario_name <- combinations$index[zeta]
    processed_beetles_name <- combinations$species[zeta]
    
    single_scenario_i <- scenario_group[index == single_scenario_name]
    processed_beetles_i <- processed_beetles[species == processed_beetles_name]
    
    # Print the current iteration
    cat("Running iteration", zeta, "\n")
    
    result <- function_scenario_60yr_uncertainty(single_scenario_i, processed_beetles_i)
    
    return(result)
  }
  
  
  timing <- system.time({
    # result_list <- lapply(1:nrow(combinations), execute_uncertainty_fun)
    
    # result_list <- parLapplyLB(cl, 1:1000, execute_uncertainty_fun)
    result_list <- parLapplyLB(cl, 1:nrow(combinations), execute_uncertainty_fun)
  })
  
  cat("Elapsed time: ", timing[3], " seconds\n")
  # Clean up parallel backend
  stopCluster(cl)
  
  
  #save the output to an rds folder 
  saveRDS(result_list, file = rds_file_path)
  
  #saveRDS(result_list, "scenario60yrUncertainty1.rds")
  
  rm(result_list)
  
  # Print the execution time
  cat("Elapsed time: ", timing[3], " seconds\n")
}



#-----------------calculate geometric for each  iteration and species category ----
cap <- 1.5 # don't allow scenario occ to be more than 1.5 starting landscape occ [only used if calculating geometric mean]
sppCategories <- readRDS("R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/DBsppCategories.rds")
sppCategories<- as.data.table(sppCategories)
rds_folder <- "R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/Ab60perScenarioIteration"
ab60_files <- list.files(rds_folder, pattern = "*.rds", full.names = TRUE)


#scenario start abundance
#SL_60yrOcc 
SL_occ60 <- readRDS("R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/SL_occ60yr_perIteration.rds")
SL_occ60_dt <- rbindlist(SL_occ60) %>%
  rename(SL_occ_60yr = occ_60yr)

#EXTRACT ONLY THE BASELINE ALL_PRIMARY SL
SL_all_primary_dt<- SL_occ60_dt %>% filter(scenarioStart == "all_primary") 

# Allocate folder for summarised results (geometric means and relative abundance )
geom_result_folder <- "R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/GeometricMeansPerIteration"


for (w in seq_along(ab60_files)){
  occ60 <- readRDS(ab60_files[[w]])
  occ60_dt <- rbindlist(occ60)
  
  # Generate a unique  file name
  rds_file_name <- paste("OGbaseline_", basename(ab60_files[[w]]), sep = "")
  
  #rds_file_name <- paste("geometricMean_", basename(ab60_files[[w]]), sep = "")
  #relOcc_file_name <- paste("relOcc_", basename(occ60_files[[w]]), sep = "")
  
  
  # Combine folder path and file name to create full file path
  geom_file_path <- file.path(geom_result_folder, rds_file_name)   
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
  
  #summarise species-level median rel occ across 500 iterations 
  
  # Summarize species-level median rel occ across iterations
  geom_means <- occ_comb[, .(SppMedRelOcc = median(rel_occ, na.rm = TRUE)), 
                         by = .(species, index, production_target)]
  
  #add back in species categories if not involved in the grouping variable above  
  geom_means <- sppCategories[geom_means, on = "species"]
  
  # Print the current iteration
  cat("Running iteration", rds_file_name, "\n")
  
  #save the output to an rds folder 
  saveRDS(geom_means, file = geom_file_path)
}

#---


#--------  read in geom means ----------------------
#geomMean_files <- list.files(geom_result_folder, pattern = "*.rds", full.names = TRUE)
#read in data that has been baselined the fully old-growth starting landscape  
geomMean_files <- list.files(geom_result_folder, pattern = "^OGbaseline.*\\.rds$", full.names = TRUE)
geomMeans <- lapply(geomMean_files, readRDS)

#----- summarise geom means and across posterior draws  -----

#check; are values normally distributed across posterior draws - if not then take median of median rel occ
#or median of geometric means 
subGeom <- geomMeans[[1]] %>% filter(spp_category == "loser") %>% select(index) %>% 
  unique() %>% slice(1:80) %>% pull()
checkDist <- geomMeans[[1]] %>% filter(index %in% subGeom)

checkDist %>% ggplot(aes(x = SppMedRelOcc)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.6) +
  facet_wrap(~index, scales = "free") +
  labs(
    title = "Histogram of medianRelOcc by Index",
    x = "medianRelOcc",
    y = "Frequency"
  )+
  xlim(0, 1)

checkDist %>% filter(index == "all_primary_CY_D.csv 6") %>%
  ggplot(aes(x = SppMedRelOcc)) +
  geom_histogram(binwidth = 0.0001, color = "black", fill = "blue", alpha = 0.6) +
  facet_wrap(~index, scales = "free") +
  labs(
    title = "Histogram of medianRelOcc by Index",
    x = "medianRelOcc",
    y = "Frequency"
  )+
  xlim(0, 1)

subGeom <- geomMeans[[1]]
subGeom %>% filter(index == "all_primary_CY_D.csv 6")

names(subGeom)  
# summarise_across_posterior_fun <- function(x){
#   x %>% group_by(spp_category, index, production_target) %>%  
#     summarise(geom_mean =median(geometric_mean), 
#               medianRelativeOccupancy = median(medRelOcc),
#               p5_medianRelativeOccupancy = quantile(medianRelativeOccupancy, 0.05),
#               p95_medianRelativeOccupancy = quantile(medianRelativeOccupancy, 0.95), 
#               IQR)
#               
# }

#summarised across data that has already been medianed per spp from 500 draws

summarise_across_posterior_fun <- function(x){
  x %>% #left_join(sppCategories, by = "species") %>%
    group_by(spp_category, index, production_target) %>%  
    summarise(medianRelativeOccupancy = median(SppMedRelOcc),
              p5_medianRelativeOccupancy = quantile(SppMedRelOcc, 0.05),
              p95_medianRelativeOccupancy = quantile(SppMedRelOcc, 0.95), 
              IQR = IQR(SppMedRelOcc))
  
}
#---output of summarised statistics across posterior draws ----
final_geoms <- lapply(geomMeans, summarise_across_posterior_fun)
final_geoms <- rbindlist(final_geoms)

#add back in key information 
final_geoms <- final_geoms %>%  left_join(scenario_composition, by = c("index", "production_target"), relationship = "many-to-many")# %>% 



#-----EXPORT OUTCOME PERFORMANCE for consolidated figure of all outcomes -----
getwd()
names(final_geoms)
output <- final_geoms %>% select(index, production_target, scenarioName,scenarioStart,
                                 medianRelativeOccupancy, p5_medianRelativeOccupancy, p95_medianRelativeOccupancy,
                                 spp_category) %>% cbind(outcome = "dungBeetles ")
saveRDS(output, "R_code/AllOutcomesFigure/Data/OG_baseline_dungBeetles.rds")

#--------------








#---- Calculate Prop OG in scenario ----

#get the amount of hab in each starting landscape 
primaryInStart <- all_start_landscape %>% filter(habitat == "primary") %>%  
  rename(SL_primary_parcels = num_parcels) %>% dplyr::select(-habitat)
habInStart <- all_start_landscape %>% select(scenarioStart) %>% unique() %>% 
  mutate(originalOG = c(1,0.2,0.2,0.8,0.2,0.2), 
         original1L = c(0,0.8,0,0,0.6,0), 
         original2L = c(0,0,0.8,0,0,0.6))


#build a function that calculates proportion of remaining habitat 
#in each scenario 

prop_OG_fun <- function(x){
  
  #proportion of TOTAL landscape [1000 parcels] in different habitat type 
  x %>% group_by(index, production_target) %>% 
    #total OG
    mutate(propOG = sum(num_parcels[habitat == "primary"])/1000,
           propPlant = sum(num_parcels[habitat %in% c("eucalyptus_current", "albizia_current", "albizia_future","eucalyptus_future")])/1000,   
           #prop-1L in the scenario landscape
           prop1L = sum(num_parcels[habitat == "once-logged"])/1000,
           #proportion of 2-L in the scenario landscape
           prop2L = sum(num_parcels[habitat == "twice-logged"])/1000) %>%  
    
    #get starting landscape
    mutate(scenarioStart = scenarioName) %>% 
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_D.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_D.csv")) %>% 
    ungroup %>% 
    
    #get total amount of each habitat in STARTING landscape for a scenario
    left_join(habInStart, by = "scenarioStart") %>% 
    
    #calculate PROPORTION of REMAINING original habitat type 
    #(nb there can actually be more once-logged or twice-logged forest in scenario than scenarioStart, if primary forest is logged)
    mutate(remainingOG = propOG/originalOG, 
           remaining1L = prop1L/original1L, 
           remaining2L = prop2L/original2L) %>%  
    #correct for INF values for if dividing by 0
    mutate_at(vars(remainingOG, remaining1L, remaining2L), ~ ifelse(is.infinite(.) | is.nan(.), 0, .)) %>%
    
    select(index, production_target, scenarioName,scenarioStart,
           propOG, propPlant,prop1L,prop2L,
           remainingOG,remaining1L,remaining2L) %>% unique()
  
}

propOGcomp <- prop_OG_fun(scenario_composition) %>% ungroup

#for each scenario, add the proportion starting landscapes
propOGcomp_dt <- as.data.table(propOGcomp)
geom_results <- as.data.table(final_geoms)

#if index is numeric make character
geom_results <- geom_results[, index := as.character(index)]
geom_results <- geom_results[, production_target := as.numeric(production_target)]
propOGcomp_dt <- propOGcomp_dt[, index := as.character(index)]
propOGcomp_dt <- propOGcomp_dt[, production_target := as.numeric(production_target)]

geom_results_df <- propOGcomp_dt[geom_results, on = .(index, production_target)] 


#---------- #BIVARIATE PLOTTING PARAMETRES --------------------

library(biscale)
COL <- "DkBlue2" # define colour pallete
COL <- "BlueOr"
#get colours for bivariate plotting
biv_pallete <- bi_pal(COL, dim =4 ) # for plotting
cols <- data.frame(bi_pal(COL, dim = 4, preview = FALSE))
colnames(cols) <- c("hex")
cols <- cols %>% mutate(bi_class = rownames(.))

textSize  <- 15

#make bivar legend
primary_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                            xlab = "Proportion old-growth", 
                            ylab = "Proportion once-logged", size = textSize)

onceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                          xlab = "Proportion remaining old-growth", 
                          ylab = "Proportion remainng once-logged",size = textSize)

twiceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                           xlab = "Proportion remaining old-growth", 
                           ylab = "Proportion remainng twice-logged",size = textSize)

all_legend <- plot_grid(primary_legend,onceL_legend,twiceL_legend, ncol =3)

#assign scenarios the colours from the bivariate plot for primary start
bivariate_colours_PRIM <- function(X){
  X %>%  bi_class(x = propOG, y = prop1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_PRIM(geom_results_df) %>% rename(hexP = hex)

#assign scenarios the colours from the bivariate plot for mostly 1L start
bivariate_colours_1L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_1L(geom_results_df) %>% rename(hex1L = hex)

#assign scenarios the colours from the bivariate plot for mostly 2L start
bivariate_colours_2L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining2L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_2L(geom_results_df) %>% rename(hex2L = hex)

#hex shows colours for 1L vs primary.
#hex_2L shows colours for 2L vs primary


# final_carbon_df_4DR <- bi_class(final_carbon_df_4DR, x = propOriginalOG, y = prop1L, dim = 4, style = "equal") %>%  
#   left_join(cols, by = "bi_class") # add hex colours

#================= build some summary plts ==========================================================

#filter by category
legend_plot <-  geom_results_df %>% filter(spp_category == "loser" & scenarioName == "all_primary_CY_D.csv") 
losers <- geom_results_df %>% filter(spp_category == "loser") 
intermediate1L <- geom_results_df %>% filter(spp_category == "intermediate1L") 
intermediate2L <- geom_results_df %>% filter(spp_category == "intermediate2L") 
winners <-  geom_results_df %>% filter(spp_category == "winner") 



scenario_filters <- c("all_primary_CY_D.csv")#, "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")



#build PLOTTING FUNCTION #### 
plot_fun <- function(x){
  
  x <- x %>%
    filter(scenarioName %in% scenario_filters)
  
  #if scenario contains plantation add cross 
  x <- x %>%
    mutate(is_cross = ifelse(propPlant > 0, "Cross", "Point"))
  
  
  max_propPlant <- max(x$propPlant, na.rm = TRUE)
  
  
  #reorder facet order 
  
  x$scenarioName <- factor(x$scenarioName, levels = c(
    "all_primary_CY_D.csv"))#,
  #"all_primary_CY_ND.csv","all_primary_IY_D.csv", "all_primary_IY_ND.csv",
  ##"mostly_1L_CY_D.csv",
  #"mostly_1L_CY_ND.csv", "mostly_1L_IY_D.csv", "mostly_1L_IY_ND.csv",
  #"mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv", "mostly_1L_deforested_IY_D.csv", "mostly_1L_deforested_IY_ND.csv", 
  #"mostly_2L_CY_D.csv"))
  #"mostly_2L_CY_ND.csv", "mostly_2L_IY_D.csv", "mostly_2L_IY_ND.csv",
  #"mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv", "mostly_2L_deforested_IY_D.csv","mostly_2L_deforested_IY_ND.csv",
  #"primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv", "primary_deforested_IY_D.csv","primary_deforested_IY_ND.csv"))
  
  
  x %>%  ggplot(aes(x = production_target, y = medianRelativeOccupancy))+
    
    # conditionally colour so that if we plot bivariate between proportion of primary and proportion of least logged (either 1L or 2L depending on starting landscape) in the scenario    
    geom_point(aes(
      x = production_target,
      # y = geom_mean,
      y = medianRelativeOccupancy,
      colour = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        #  scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        #  scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L
      ),
      shape = is_cross, 
    ), position = position_jitter(width = 0.05, height = -0.03)) +
    scale_colour_identity()+
    #GIVE CORSSES TO PLANTATION CONTAINING SCENAIOS ####
  scale_shape_manual(values = c("Point" = 19, "Cross" = 3)) + # Define shape mapping
    #scale_size_continuous(range = c(3, 8), breaks = seq(0, max_propPlant, by = 0.05)) + # Adjust size range and breaks
    
    xlim(0, 1)+
    xlab("Production target")+
    ylab(   "Median relative abundance 
  (averaged over posterior draws)")+
    
    #labs(colour = "Proportion of remaining old-growth forest spared")+
    # labs(colour = "Proportion of plantation in remaining landscape")+
    
    facet_wrap(~scenarioName, ncol = 4)+
    #   geom_hline(aes(yintercept = SL_geom_mean))+
    theme_bw(base_size = textSize)+
    theme(legend.position = "none")
  
}


#get legend
legend_plot <- plot_fun(legend_plot)
legend <- get_legend(legend_plot + theme(legend.position = "bottom",         # c(0.5, 0.15),
                                         legend.spacing.x = unit(1.0, 'cm'),
                                         legend.title  = element_text(size  = 30, face = "bold"))) 
#plot figures (without legends)
losers <- plot_fun(losers)
intermediate1L <- plot_fun(intermediate1L)
intermediate2L <- plot_fun(intermediate2L)
winners <- plot_fun(winners)

#add legend function
add_legend <-  function(x){
  plot_grid(x, legend, 
            nrow =2 , ncol = 1,
            rel_heights = c(1, 0.1))
} 

#plot final figs for each above-defined category 
add_legend(losers)
plot_grid(losers, all_legend, nrow =2)
add_legend(intermediate1L)
plot_grid(intermediate1L, all_legend, nrow =2)
add_legend(intermediate2L)
plot_grid(intermediate2L, all_legend, nrow =2)
add_legend(winners)
plot_grid(winners, all_legend, nrow =2)


#---UNUSED CODE----






















































































# ---------------RUN ONCE -------------------------------
#CALCULATE TIME-AVERAGED OCCUPANCY FOR EACH SCENARIO ####
#takes 30 mins for 1-delay period 

#code takes a csv file for a scenario, joins bird data and processes time-averaged biodiversity,
#which it then stores as a csv file - RUN ONCE 
#nb - need to rerun for different delay periods
for (csv_file_path  in csv_files) {
  
  # Get a single  scenario 
  scenario <- read.csv(csv_file_path)
  #make scenario data table format
  scenario <- as.data.table(scenario)
  
  #filter subset of harvest delays
  scenario <- filtDelay(scenario)
  
  # #set the join keys to match on for merge 
  # setkeyv(birds_10km2, c("habitat", "functionalhabAge"))
  # setkeyv(scenario, c("functional_habitat", "functionalhabAge"))
  
  
  
  # # Join that scenarios (based on above set join keys) to bird occupancy data
  scen_bio <- scenario[DBs_10km2,
                       on = .(functional_habitat == habitat,
                              functionalhabAge ==  functionalhabAge),
                       nomatch = NA,
                       allow.cartesian=TRUE]
  
  # on = c("functional_habitat" = "habitat", "functionalhabAge" = "functionalhabAge"), 
  #  nomatch = NA]#, #if not exact match, don't join
  #  allow.cartesian=TRUE] #allows "many-to-many" join
  
  # Reset keys (to remove grouping)
  setkey(scen_bio, NULL)
  
  # Apply the landscape_occ function to calculate the landscape-wide occupancy of that scenario 
  #this shows occupancy between 0-1 over the entire scenario landscape, for each species and year (0-60)
  scen_bio_processed <- landscape_occ(scen_bio)
  
  #remove chunky scen_bio
  rm(scen_bio)
  
  #-------------------------- apply the timed_occ function-----------------------------
  # apply the timed_occ function to get occupancy for each scenario and species summed across 60years 
  time_scen_bio_processed <- timed_occ(scen_bio_processed)   #this bit sums occupancy across 60 years and returns one row per scenario for time-averaged occupancy
  # Store the processed data in the outcomes list            #REMVOVE if you want annual occupancy per species per year per scenario 
  
  #............... save the processed scenario as a csv  ........................
  
  # Define a new file name for the processed data
  processed_csv_file_name <- paste("processed_", unique(time_scen_bio_processed$scenarioName), sep = "")
  #make file path by combinning folder name and file name
  processed_csv_file_path <- file.path(final_output_folder, processed_csv_file_name)
  
  # Save the processed data as a new CSV file
  # ----UNCOMMENT ----
  write.csv(time_scen_bio_processed, file = processed_csv_file_path, row.names = FALSE)
  
  # Remove the scenario data from memory
  rm(scenario_data)
  rm(processed_scenario_data)
  
  
  # Print a progress message
  cat("Processed CSV file", processed_csv_file_name, "saved.\n")
  
}

#------------ read back in occupancy calculated through time for each scenario (currently) ----------------------------------------------

# Get a list of all CSV files in the folder
csv_files <- list.files(final_output_folder, pattern = "\\.csv$", full.names = TRUE)

#make one master dataframe of all DB results
combined_df <- data.table()

# Loop through each CSV file and rbind into the combined dataframe
for (csv_file_path in csv_files) {
  # Read the CSV file
  csv_data <- fread(csv_file_path)
  
  # rbind the data to the combined dataframe
  combined_df <- rbind(combined_df, csv_data)
}

#====================== calculate starting landscape  time-averaged abundance ==================================================================

hab_by_year <- read.csv("Tables/HabByYears.csv", strip.white = TRUE) %>%  
  rename(true_year = year, 
         functionalhabAge = functional_habAge, 
         habitat = transition_habitat) %>% select(-X)

#add time
time_forSL <- hab_by_year %>% 
  select(original_habitat, true_year) %>% unique() %>%  
  rename(habitat = original_habitat, 
         functionalhabAge = true_year)
all_start_landscape <- all_start_landscape %>%
  left_join(time_forSL, by = "habitat", relationship = "many-to-many")


#add biodiversity
all_start_landscape_DBs <- all_start_landscape %>% 
  left_join(DBs_10km2, by = c("habitat","functionalhabAge"), relationship = "many-to-many")

#calculate time averaged biodiversity 
all_start_landscape_DBs  <- all_start_landscape_DBs %>%   #add columns to match function format
  mutate(index = scenarioStart, 
         production_target = 0,
         original_habitat = habitat, 
         true_year = functionalhabAge, 
         harvest_delay = "delay 0", 
         scenarioName = scenarioStart ) %>% 
  as.data.table() %>% 
  #calculate starting landscape occupancy per year and species
  landscape_occ() %>%  
  #calculate time-averaged occupancy 
  timed_occ() %>%  
  #rename to make clear this is occupancy for starting landscape 
  mutate(SLocc_60yr = occ_60yr, 
         SLocc_60yr_lwr = occ_60yr_lwr, 
         SLocc_60yr_upr = occ_60yr_upr) %>% 
  #remove unwanted columns 
  select(-c(index, scenarioName,occ_60yr,occ_60yr_lwr,occ_60yr_upr, production_target))  

#Join each scenario to the relevant starting landscape occupancy ####
combined_df <- combined_df[all_start_landscape_DBs, on = .(species, scenarioStart)] %>% 
  #remove scenarios with producion target of zero 
  filter(!production_target == 0)

# ------ save output -------------- 
write.csv(combined_df, "R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/masterDungBeeetleOutcomes_df.csv")



#---------- CAN START AGAIN HERE -------------------- 
setwd("C:/Users/Gianluca Cerullo/OneDrive - University of Cambridge/PhD/Chapter_4_Borneo/CompleteFolder")

#read in all dung beetle outcomes, calculated for each scenario
combined_df <- read.csv("R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/masterDungBeeetleOutcomes_df.csv") 

#read in 10km2 abundance outcomes 
DBs_10km2 <- read.csv("R_code/CleaningHistoricDungBeetleData/AbundanceAnalysis/Outputs/DBs_10km2.csv")


# ------ Summarise winner, loser, intermediate sppp ----------
#based abudance in different habitat tpyes 

#at the moment we are selecting losers so that if occupancy is higher at ANY age of logging, then 
#the species is not a loser. 
losers <- DBs_10km2 %>%  group_by(species)  %>%  
  #filter(true_age < 20 ) %>% 
  filter(parcel_occ == max(parcel_occ)) %>% 
  mutate(spp_category = case_when(habitat =="primary" ~"loser", TRUE ~ NA_character_)) %>% 
  filter(spp_category == "loser") %>% select(species,spp_category) %>% unique()


intermediates <- DBs_10km2 %>%  group_by(species) %>% 
  #filter(true_age < 20 ) %>% 
  filter(parcel_occ == max(parcel_occ)) %>% 
  mutate(spp_category = case_when(habitat =="once-logged" ~"intermediate",
                                  habitat == "twice-logged" ~ "intermediate", 
                                  habitat == "restored" ~ "intermediate",
                                  TRUE ~ NA_character_)) %>% 
  filter(spp_category == "intermediate") %>% select(species,spp_category) %>% unique


winners <- DBs_10km2 %>%
  group_by(species) %>%
  filter(parcel_occ == max(parcel_occ)) %>%
  mutate(spp_category = case_when(
    habitat == "albizia_current" ~ "winner",
    habitat == "eucalyptus_current" ~ "winner",
    TRUE ~ NA_character_
  )) %>% 
  filter(spp_category == "winner") %>% select(species,spp_category) %>% unique()

#spp categories 
winner_loser <- rbind(losers,intermediates,winners) %>% ungroup

#===========   join species category information for each scenario      ============================================================= 

add_spp_inf <- function(x){
  x %>%  left_join(winner_loser, by = "species")
}

combined_df <- add_spp_inf(combined_df)

#============================  Calculate the  geometric mean with errors   =============================================
# Calculate the  geometric mean with errors



#1. calculate relative occupancy (occupancy in scenario landscape versus in starting landscape
#NB calculating delta_occupancy(the difference in time-averaged occupancy between starting landscape and scenario landscape) results in many negative values (as higher occ in scenrio landscape) which means you can't take geometric mean (as you can't log a negative value)
rel_occ <- function(x){
  x %>% 
    #calculate 
    mutate(rel_occ  = occ_60yr/ SLocc_60yr
    ) %>%  
    #APPLY CAP ON RELATIVE ABUNDANCE SO THAT IT CAN'T BE > THAN CAP TIMES STARTING LANDSCAPE OCCUPANCY
    mutate(rel_occ = ifelse(rel_occ > cap, cap, rel_occ)) %>%  
    
    #CALCULATE ERROR
    
    #THIS is the individual error for a specific species comparing starting landscape versus scenario landscape
    mutate(errorSL = sqrt(
      ( 
        #calcaulate total error for starting landscape for a given species 
        abs(SLocc_60yr - SLocc_60yr_upr)^2 + abs(SLocc_60yr - SLocc_60yr_lwr)^2) 
    )/2 
    ) %>%  
    #calcaulate total error for scenario landscape for a given species 
    mutate(errorScen = sqrt(
      (
        abs(occ_60yr - occ_60yr_upr)^2 + abs(occ_60yr - occ_60yr_lwr)^2)
    )/2
    )
  
}

combined_df <- rel_occ(combined_df)

#The geometric mean function incorporates two types of error:
#1. Species-level error, which is the error associated with calculating occupancy for the starting and scenario landscape for a given species
#2. Geometric mean error, which is the error from taking a mean over diferent species with different errors

#note; we use weighted geometric means as this means that we weight species where we are less sure of species-level occupancy less strongly
geom_means <- function(x){
  x %>% group_by(index, spp_category)  %>%  
    mutate(
      
      #calculate the geometric mean of relative occupancy (the difference in time-averaged occupancy between starting landscape and scenario landscape)
      #scenario geom mean and error 
      geom_mean = exp(mean(log(rel_occ))),
      
      #calculate nrows in geom_mean calculation (equivalent to num species in the subset)
      nrows_geom_mean = n(),
      
      #calculate the error of calculatinng spp error between starting landscpape and scenario landscape 
      sppError = sqrt(
        (errorScen^2+ errorSL^2)/2
      ),
      
      #calcluate arithmetic mean of relative occupancy 
      arith_mean = mean(rel_occ)
      
      #calculate the erorr from calculating geometric mean across species with different errors   
    ) %>% 
    mutate(
      # multSppsd = sd(rel_occ),
      geomError = geom_mean *   (
        (prod(sppError / geom_mean))
        ^(1/nrows_geom_mean)     )
    ) %>% 
    
    mutate(
      # Calculate the weighted geometric mean using the reciprocal of the square of the species-specific errors as weights
      weighted_geom_mean = exp(sum(log(rel_occ) / sppError^2) / sum(1 / sppError^2)),
      
      # Calculate the error from calculating weighted geometric mean across species
      weighted_geomError = weighted_geom_mean * (
        prod(sppError / weighted_geom_mean) ^ (1 / nrows_geom_mean)
      )
    ) %>% 
    
    #NEED TO INPUT HOW TO
    slice(1) %>%  
    dplyr::select(index, spp_category, geom_mean,arith_mean,scenarioName,production_target,
                  sppError,geomError, nrows_geom_mean,weighted_geom_mean,weighted_geomError)
}
geom_results <- geom_means(combined_df)

#=========================    Summarise the composition of scenarios (proportion of forest) ==================================================

#get the amount of hab in each starting landscape 
primaryInStart <- all_start_landscape %>% filter(habitat == "primary") %>%  
  rename(SL_primary_parcels = num_parcels) %>% dplyr::select(-habitat)
habInStart <- all_start_landscape %>% select(scenarioStart) %>% unique() %>% 
  mutate(originalOG = c(1,0.2,0.2,0.8,0.2,0.2), 
         original1L = c(0,0.8,0,0,0.6,0), 
         original2L = c(0,0,0.8,0,0,0.6))


prop_OG_fun <- function(x){
  
  #proportion of TOTAL landscape [1000 parcels] in different habitat type 
  x %>% group_by(index, production_target) %>% 
    #total OG
    mutate(propOG = sum(num_parcels[habitat == "primary"])/1000,
           propPlant = sum(num_parcels[habitat %in% c("eucalyptus_current", "albizia_current", "albizia_future","eucalyptus_future")])/1000,   
           #prop-1L in the scenario landscape
           prop1L = sum(num_parcels[habitat == "once-logged"])/1000,
           #proportion of 2-L in the scenario landscape
           prop2L = sum(num_parcels[habitat == "twice-logged"])/1000) %>%  
    
    #get starting landscape
    mutate(scenarioStart = scenarioName) %>% 
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_ND.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_IY_D.csv")) %>%
    mutate(scenarioStart = str_remove(scenarioStart, "_CY_D.csv")) %>% 
    ungroup %>% 
    
    #get total amount of each habitat in STARTING landscape for a scenario
    left_join(habInStart, by = "scenarioStart") %>% 
    
    #calculate PROPORTION of REMAINING original habitat type 
    #(nb there can actually be more once-logged or twice-logged forest in scenario than scenarioStart, if primary forest is logged)
    mutate(remainingOG = propOG/originalOG, 
           remaining1L = prop1L/original1L, 
           remaining2L = prop2L/original2L) %>%  
    #correct for INF values for if dividing by 0
    mutate_at(vars(remainingOG, remaining1L, remaining2L), ~ ifelse(is.infinite(.) | is.nan(.), 0, .)) %>%
    
    select(index, production_target, scenarioName,scenarioStart,
           propOG, propPlant,prop1L,prop2L,
           remainingOG,remaining1L,remaining2L) %>% unique()
  
}

propOGcomp <- prop_OG_fun(scenario_composition) %>% ungroup

#for each scenario, add the proportion starting landscapes
propOGcomp_dt <- as.data.table(propOGcomp)
geom_results <- as.data.table(geom_results)
geom_results_df <- propOGcomp_dt[geom_results, on = .(index, production_target,scenarioName)] 



#----------    BIVARIATE PLOTTING PARAMETRES--------------------

library(biscale)
COL <- "DkBlue2" # define colour pallete
COL <- "BlueOr"
#get colours for bivariate plotting
biv_pallete <- bi_pal(COL, dim =4 ) # for plotting
cols <- data.frame(bi_pal(COL, dim = 4, preview = FALSE))
colnames(cols) <- c("hex")
cols <- cols %>% mutate(bi_class = rownames(.))

textSize  <- 15

#make bivar legend
primary_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                            xlab = "Proportion old-growth", 
                            ylab = "Proportion once-logged", size = textSize)

onceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                          xlab = "Proportion remaining old-growth", 
                          ylab = "Proportion remainng once-logged",size = textSize)

twiceL_legend <- bi_legend(pal = "BlueOr", dim = 4, 
                           xlab = "Proportion remaining old-growth", 
                           ylab = "Proportion remainng twice-logged",size = textSize)

all_legend <- plot_grid(primary_legend,onceL_legend,twiceL_legend, ncol =3)

#assign scenarios the colours from the bivariate plot for primary start
bivariate_colours_PRIM <- function(X){
  X %>%  bi_class(x = propOG, y = prop1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_PRIM(geom_results_df) %>% rename(hexP = hex)

#assign scenarios the colours from the bivariate plot for mostly 1L start
bivariate_colours_1L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining1L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_1L(geom_results_df) %>% rename(hex1L = hex)

#assign scenarios the colours from the bivariate plot for mostly 2L start
bivariate_colours_2L <- function(X){
  X %>%  bi_class(x = remainingOG, y = remaining2L, dim = 4, style = "equal") %>%  
    left_join(cols, by = "bi_class") # add hex colours
}
geom_results_df<- bivariate_colours_2L(geom_results_df) %>% rename(hex2L = hex)
#hex shows colours for 1L vs primary.
#hex_2L shows colours for 2L vs primary


# final_carbon_df_4DR <- bi_class(final_carbon_df_4DR, x = propOriginalOG, y = prop1L, dim = 4, style = "equal") %>%  
#   left_join(cols, by = "bi_class") # add hex colours

#======================= #build some summary plts: ====================================================


#filter by category
legend_plot <-  geom_results_df %>% filter(spp_category == "loser" & scenarioName == "all_primary_CY_D.csv") 
losers <- geom_results_df %>% filter(spp_category == "loser") 
intermediate <- geom_results_df %>% filter(spp_category == "intermediate") 
winners <-  geom_results_df %>% filter(spp_category == "winner") 

#==========================#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
scenario_filters <- c("all_primary_CY_D.csv", "mostly_1L_CY_D.csv", "mostly_2L_CY_D.csv")

#build plotting function
plot_fun <- function(x){
  
  x <- x %>%
    filter(scenarioName %in% scenario_filters)
  
  
  
  #reorder facet order 
  
  x$scenarioName <- factor(x$scenarioName, levels = c(
    "all_primary_CY_D.csv",
    #"all_primary_CY_ND.csv","all_primary_IY_D.csv", "all_primary_IY_ND.csv",
    "mostly_1L_CY_D.csv",
    #"mostly_1L_CY_ND.csv", "mostly_1L_IY_D.csv", "mostly_1L_IY_ND.csv",
    #"mostly_1L_deforested_CY_D.csv", "mostly_1L_deforested_CY_ND.csv", "mostly_1L_deforested_IY_D.csv", "mostly_1L_deforested_IY_ND.csv", 
    "mostly_2L_CY_D.csv"))
  #"mostly_2L_CY_ND.csv", "mostly_2L_IY_D.csv", "mostly_2L_IY_ND.csv",
  #"mostly_2L_deforested_CY_D.csv", "mostly_2L_deforested_CY_ND.csv", "mostly_2L_deforested_IY_D.csv","mostly_2L_deforested_IY_ND.csv",
  #"primary_deforested_CY_D.csv", "primary_deforested_CY_ND.csv", "primary_deforested_IY_D.csv","primary_deforested_IY_ND.csv"))
  
  
  x %>%  ggplot(aes(x = production_target, y = weighted_geom_mean))+
    # conditionally colour so that if we plot bivariate between proportion of primary and proportion of least logged (either 1L or 2L depending on starting landscape) in the scenario    
    geom_point(aes(
      x = production_target,
      y = weighted_geom_mean,
      colour = case_when(
        scenarioStart %in% c("all_primary", "primary_deforested") ~ hexP,
        scenarioStart %in% c("mostly_1L", "mostly_1L_deforested") ~ hex1L,
        scenarioStart %in% c("mostly_2L", "mostly_2L_deforested") ~ hex2L
      )
    ), position = position_jitter(width = 0.05, height = -0.03)) +
    scale_colour_identity()+
    
    xlim(0, 1)+
    xlab("Production target")+
    ylab("Weighted geometric mean change")+
    
    #labs(colour = "Proportion of remaining old-growth forest spared")+
    # labs(colour = "Proportion of plantation in remaining landscape")+
    
    facet_wrap(~scenarioName, ncol = 4)+
    #   geom_hline(aes(yintercept = SL_geom_mean))+
    theme_bw(base_size = textSize)+
    theme(legend.position = "none")
  
}

#get legend
legend_plot <- plot_fun(legend_plot)
legend <- get_legend(legend_plot + theme(legend.position = "bottom",         # c(0.5, 0.15),
                                         legend.spacing.x = unit(1.0, 'cm'),
                                         legend.title  = element_text(size  = 30, face = "bold"))) 
#plot figures (without legends)
losers <- plot_fun(losers)
intermediate <- plot_fun(intermediate)
winners <- plot_fun(winners)

#add legend function
add_legend <-  function(x){
  plot_grid(x, legend, 
            nrow =2 , ncol = 1,
            rel_heights = c(1, 0.1))
} 

#-------------------- PLOT FINAL FIGURES ----------------------------------------
#plot final figs for each above-defined category 
add_legend(losers)
plot_grid(losers, all_legend, nrow =2)
add_legend(intermediate)
plot_grid(intermediate, all_legend, nrow =2)
add_legend(winners)
plot_grid(winners, all_legend, nrow =2)







