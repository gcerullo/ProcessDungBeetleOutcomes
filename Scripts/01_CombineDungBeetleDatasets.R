#31.05.2024
#MASTER CODE FOR COMBINING ALL DUNG BEETLES DATASETS 

#this code is for cleaning all historic dung beetle for Sabah, including 
#1.  Gianluca Master's 2017 dung beetle data (by day)
#2.  Finlayson 2019 data (by day)
#3.  Felicity Data (2011,2014,2017 data) (pooled across 4 days)
#4.  Trond Larsen 2009 Data (pooled from Master's 2017 data - if needed I can do by day; I have the data)
#5   Gianluca 2022 plantation and logged forest data 
#6.  Slade 2005 Brl and Danum Data. 
#7   Slade SAFE sites (once-logged)
#Combine all with matching species names into one dataframe


#IMPORTANT TAXONOMIC NOTES:' 
#dung beetle taxonomy has improved through time, leading to species that were identified in earlier datasetes 
#being split to a greater number of species in subsequent datasets. 

#THUS: Catherine 2017 seperates: 
#Ochicanthon dysticoides from O masomotoi 
# Onthophagus negrobscurior from O. cervicapra 
#Onthophagus paviobscurior from 0. obscurior 

#THUS: Gianluca (2017) seperates: 
#C.rend from C.dayacus 
#Onthophagus aff delinesis from Onthophagus hidiaki


library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(forcats)
library(geosphere)
library(stringr)


#read in functions
source("Functions/DB_functions.R")

#READ IN INPUTS AND DATA: 

#READ IN DATA ####

#DUNG BEETLE DATA ####

#Gianluca 2017 data ####

Primary <- read.csv("RawData/MastersPrimaryForestSites.csv") %>% cbind(habitat = "primary") %>% rename(spp = 1) 
Logged <- read.csv("RawData/LoggedMastersSites.csv") %>% cbind(habitat = "once-logged")  %>% rename(spp = 1)
Restored <- read.csv("RawData/MastersRestoredSites.csv") %>% cbind(habitat = "restored") %>% rename(spp = 1)
#combine into a list
mastersList <- list(Primary,Logged,Restored)

#Catherine's data ####
#NAs on day 1,2,3 = not collected (means these beetles were just collected the subsequent day, so fine to pool)
#NAS on day 4 means dung beetles were not collected AT ALL = a bias when pooling across days - corrected below  

catherine <- read.csv("RawData/DB_Data_2019_Finlayson_2L_and_Primary.csv") %>% rename(spp =1) %>%  slice(-1:-5)


#Felicity's Data ####
#Nb - this is pooled across 4 days and must have Gianluca Master's and Trond sites removed 

trond_felicity <- read.csv("RawData/FelicityTrondGianData.csv") 

#Trond's Data ####
#read in Trond 2wice logged data (split by day)
trond_2l <- read.csv("RawData/TwiceLoggedTrondData.csv") %>% slice(-1:-5)

#Eleanor Slade's data ####
#Read in data from Eleanor Slades 2005 BRL and Danum PhD
slade <- read.csv("RawData/Slade_2005_DanumBrL.csv") %>%  rename(traptype = trap, 
                                                                 trap = trapID, 
                                                                 sampler = author.study) %>% 
  filter(traptype == "BPT") %>% # take only BPT Trap  
  filter(habitat == "Old-growth Forest" | habitat == "Logged Forest") # Take only OG and LF data 

#Gianluca 2022 data ####
#Primary, logged forest, restored forest and timber plantation
gian2022 <- read.csv("RawData/DungBeetle2022DataRawDanum_SSBPlantation.csv") %>% mutate_all(as.character)


#OTHER INPUTS - Read in Necessary inputs  #### 

#sppName backbone ####
#Import species-name backbone to match names for species across datasets
allsppNameBackbone <- read.csv("Inputs/AllNamesBackbone.csv") %>% dplyr::select(-X)

#GPS locations ####
#NB need to come back to this as lots of points are certain points are missing GPS data 

HistoricGPS <- read.csv("Inputs/allHistoricSites.csv") %>% dplyr::select(site,trap, Long, Lat)

#add gian2022 forest point count data 
GC2022_gps <- read.csv("Inputs/allForest2022PointLocations.csv") %>% rename(Long = X, Lat = Y) %>% 
  dplyr::select(Long,Lat, trap, site)

#add Eleanor's GPS points (I don't have SAFE or selected Danum points, and some traps are inferred manually by names, so there may be slight inconsitencies in location)
slade_GPS <- read.csv("Inputs/EleanorSladeDungBeetle_GPS_points.csv") %>% select(lat, lon, GC_ID) %>% 
  rename(trap = GC_ID, 
         Long = lon, 
         Lat = lat)


#Plantation age info ####
plantation <- read.csv("Inputs/Plantation_Habitat_Structure.csv") %>% dplyr::select(Site,Point, Age) %>% 
  mutate(Point = paste0("t", Point)) %>% 
  rename(site = Site, 
         trap = Point,
         plantation_age = Age) 

#add plantation GPS sites to Historic GPS 
plantationGPS <- read.csv("Inputs/DaveSimonDanielPoints.csv") %>%
  filter(Habitat == "Albizia_falcataria"| Habitat ==  "Eucalyptus_pellita") %>%  
  mutate(point = str_replace(point, "^P", "t")) %>%  
  rename(site = Site, 
         trap = point, 
         Lat = Lat, 
         Long = Long, 
         habitat = Habitat) %>% 
  mutate(habitat = case_when(
    habitat == "Albizia_falcataria" ~ "albizia",
    habitat == "Eucalyptus_pellita" ~ "eucalyptus",
    TRUE ~ habitat  # Keep other values unchanged
  )) %>% 
  select(site,trap, Lat, Long)


HistoricGPS <- HistoricGPS %>% rbind(plantationGPS)

#add in CrossWalk for correcting 35 species crosschecked for ID at Singapore Uni
#corect the 35 individuals re-analyses @ singapore University by Slade et al. 

correctIDsSingapore <- read.csv("RawData/SingaporeIDCheckCrossWalk.csv") %>%  
  rename() %>% 
  mutate(abundance = 1) %>%  
  rename(day = Day, 
         trap = Trap, 
         site = Site)

#Clean data and fill in gaps ####
#Gianluca's 2017 data 

#apply pivot function to the list then rbind list
masters <- lapply(mastersList, pivot_fun) %>% rbindlist()
masters <- masters %>% 
  mutate_at('abundance', ~replace_na(.,0)) %>% 
  cbind(sample_year = 2017) %>% cbind(sampler = "GC") %>% 
  filter(!spp == "Species")

#make data with no data for a trap/day combo NA & make site = transect 
masters <- zeroToNa(masters) %>% ungroup %>%  mutate(transect = site)

#View missing trap/day combinations 
NoSample<- masters %>% group_by(trap,site,day) %>% summarise(sum_ab = sum(abundance)) %>%  
  filter(sum_ab == 0 | is.na(sum_ab))

#simulate data for these missing trap/day combos 
masters <- simulateMissingDays(masters)


sum(masters$abundance)

#Catherine's data ####

#pivot the data to long format 

catherine <-   catherine %>% pivot_longer(cols = !c(spp),
                                          names_to = c("location","transect", "trap","day"),  # seperate "site_trap_day" into three columns..  
                                          names_sep = "_",  # 'split by presence of underscore
                                          values_to = "abundance")   #and assigning values to a column called "abundance" 

#add habitat info
catherine <- catherine %>% 
  mutate(habitat = case_when(location == "Malua" ~ "twice-logged",
                             TRUE ~ "primary")) %>% 
  mutate(site = transect) %>% 
  dplyr::select(-location)

#add year and sampler 
catherine <- catherine %>% cbind(sample_year = 2019) %>% cbind(sampler = "CF")
catherine$abundance <- as.numeric(catherine$abundance)

#correct minor errors 
#1. Block five in trap 74 has a random "." in it 
catherine <- catherine %>% 
  mutate(site = gsub("\\.", "", site))


#view all zero day/trap combinations
CathNoSamplling<- catherine %>% group_by(trap,site,day) %>% summarise(sum_ab = sum(abundance)) %>%  
  filter(sum_ab == 0 | is.na(sum_ab)) %>% ungroup()

#make data with all zeros nas 
catherine <- zeroToNa(catherine)

#simulate data for all traps not sampled
catherine <- simulateMissingDays(catherine)

sum(catherine$abundance)

#Felicity's Data ####
#pivot data

trond_felicity <-   trond_felicity %>% pivot_longer(cols = !c(spp),
                                                    names_to = c("habitat","site", "trap"),  # seperate "site_trap_day" into three columns..  
                                                    names_sep = "_",  # 'split by presence of underscore
                                                    values_to = "abundance") %>%    #and assigning values to a column called "abundance" 
  mutate(habitat = case_when(
    habitat == "UL" ~ "primary", 
    habitat == "L" ~ "once-logged", 
    habitat == "R" ~ "restored"))


#dplyr::select only sites sampled by Felicity (the rest of the data we have elsewhere, split by day)
f_primary <- trond_felicity %>% filter(site %in% c("tb1","tb2","t15","t2a","rh1","rh2")) 
f_restored <- trond_felicity %>% filter(site  %in% c("c5", "c1", "c7")) 
f_logged <- trond_felicity %>% filter(site  %in% c("westus1", "westus2", "c17","l2.1","l2.2")) 

#combine
felicity <- f_primary %>% rbind(f_logged) %>% rbind(f_restored) %>% cbind(sampler = "FE") 

#add sampling year for each 
felicity <- felicity %>% mutate(sample_year = case_when(
  site == "tb1" ~ 2011,
  site == "tb2" ~ 2011,
  site == "t15" ~ 2014,
  site == "t2a" ~ 2014,#primary
  site == "rh1" ~ 2014,
  site == "rh2" ~ 2014,
  
  site == "c5" ~ 2014,
  site == "c1" ~ 2014,
  site == "c7" ~ 2014,
  
  site == "westus1" ~ 2011,
  site == "westus2" ~ 2011,  
  site == "c17" ~ 2014,
  site == "l2.1" ~ 2014,
  site == "l2.2" ~ 2014
))


#view all zero day/trap combinations (none - no missing days)
FelicityNoSamplling<- felicity %>% group_by(trap,site) %>% summarise(sum_ab = sum(abundance)) %>%  
  filter(sum_ab == 0 | is.na(sum_ab))

#split Felicity's pooled data into 4 days, allocated abundance equally for each spp*trap*site
#to each day sampled & add site = transect
felicity <- splitPooledData(felicity) %>% mutate(transect = site)

sum(felicity$abundance)


#Trond's data #####

#READ in Trond Data 1L,Primary, restored //note I am shortcutting here and using my Master's data for Trond's restored and primary sites
#which encompassed Trond's data, but pooled;
#{if needed I have the data seperated by day}

t_primary <- trond_felicity %>% filter(site %in% c("er1","er2","tr1","tr2","brl1","brl2")) 
t_restored <- trond_felicity %>% filter(site  %in% c("danum1", "danum2")) 
t_logged <- trond_felicity %>% filter(site  %in% c("t45", "t55", "mbp1"," mbp2")) 

trond <- t_primary %>% rbind(t_restored) %>% rbind(t_logged) %>% 
  cbind(sampler = "TL") %>% cbind(sample_year = 2009) %>%
  mutate(transect = site)

#check which days didn't have sampling and assign trap-level means (none)
TrondNoSamplling<- trond %>% group_by(trap,site) %>% summarise(sum_ab = sum(abundance)) %>%  
  filter(sum_ab == 0 | is.na(sum_ab))

#split pooled data by day, assigning species abundances to  
trond <- splitPooledData(trond)

#pivot trond 2L data
trond_2l<- trond_2l %>% pivot_longer(cols = !c(spp),
                                     names_to = c("site","transect", "trap","day"),  # seperate "site_trap_day" into three columns..  
                                     names_sep = "_",  # 'split by presence of underscore
                                     values_to = "abundance") %>%   #and assigning values to a column called "abundance" %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  #add details 
  cbind(sample_year = 2009) %>% cbind(sampler = "TL") %>% cbind(habitat = "twice-logged") %>% 
  mutate_at('abundance', ~replace_na(.,0))


# #check which days didn't have sampling an assign trap-level means (none)
# TrondNoSamplling<- trond_2l %>% group_by(trap,site) %>% summarise(sum_ab = sum(abundance)) %>%  
#   filter(sum_ab == 0 | is.na(sum_ab))

#combine trond twice_logged and other trond data 
trond <- trond %>% rbind(trond_2l)

sum(trond$abundance)

#Eleanor Slade's data ####

slade <- slade %>%  pivot_longer(
  cols = !c(trap, sampler,site,month,year,traptype,habitat),
  names_to = "spp", 
  values_to = "abundance")

sladeSPP <- slade %>% dplyr::select(spp) %>% unique()


#remove data from a fragment (e.g. where site contains the word fragment)
slade <- slade %>% filter(!str_detect(site, "Fragment"))

#reorder and rename #
slade <- slade %>% dplyr::select(spp,habitat,abundance, site,trap,year, month,sampler) %>% 
  rename(sample_year = year) %>% 
  mutate(habitat = case_when(
    habitat == "Old-growth Forest" ~ "primary",
    habitat =="Logged Forest" ~ "once-logged")) %>%  
  filter(!sampler == "Gray_et al_2014_2017") %>%  
  mutate(transect = site)

names(slade)
#GPS points to ask Eleanor for
# TrapsforGPS <- slade %>% dplyr::select(site, trap,sample_year, sampler) %>% unique()
# SitesforGPS <- slade %>% dplyr::select(site) %>% unique()
# slade %>% dplyr::select(habitat,site,trap) %>% unique() %>%  group_by(habitat) %>%  count()
# write.csv(TrapsforGPS, "TrapsToAskEleanorForGPS.csv")

sum(slade$abundance)

#Gianluca's 2022 data  
gian2022 <- gian2022 %>% pivot_longer(cols = !c(spp),
                                      names_to = c("site", "trap","day"),  # seperate "site_trap_day" into three columns..  
                                      names_sep = "_",  # 'split by presence of underscore
                                      values_to = "abundance")   #and assigning values to a column called "abundance" 


#add a transect column that's same as site
gian2022 <- gian2022 %>% mutate(transect = site) %>%  
  #replace nas as zero 
  mutate(abundance = as.numeric(abundance)) %>% 
  mutate(abundance = replace(abundance, is.na(abundance), 0)) %>% 
  #make sure NA for trap;days where I did not sample 
  mutate(abundance = ifelse(site == "GCP" & trap %in% c("t6","t7","t8","t9","t10","t11","t12") & day == "d4", NA, abundance)) %>%   
  mutate(abundance = ifelse(site == "GC1L2" & day == "d4", NA, abundance))  %>% 
  cbind(sample_year = 2022) %>% 
  cbind(sampler = "GC")

#add habitat info 
gian2022 <- gian2022 %>% mutate(habitat = case_when(
  site == "EC1" ~ "eucalyptus",
  site == "EC2" ~ "eucalyptus",
  site == "EC3" ~ "eucalyptus",
  site == "AL1" ~ "albizia",
  site == "AL2" ~ "albizia",
  site == "AL3" ~ "albizia",
  site == "GC1L" ~ "once-logged",
  site == "GC1L2" ~ "once-logged",
  site == "GCR1" ~ "restored",
  site == "GCR2" ~ "restored",
  site == "GC2L" ~ "twice-logged",
  site == "GCP" ~ "primary"))


#simulate  data for traps with missing Nas on the 4th day 
gian2022 <- simulate4thDay(gian2022)

sum(gian2022$abundance, na.rm = TRUE)

names(gian2022)

#now we need to correct the IDs of species that were checked over at Singapore university
# Join the full_df and CorrectedSpID dataframes

gian2022 <- gian2022 %>%
  full_join(correctIDsSingapore, 
            by = c("spp" = "GcnameOldName", "site" = "site", "trap" = "trap", "day" = "day")) %>%
  mutate(
    # Update spp where CorrectedSpName is not NA, otherwise keep original spp
    spp = if_else(!is.na(CorrectedSpName), CorrectedSpName, spp), 
   
     # Fill in habitat, sampler, sample_year, and transect from correctIDsSingapore where needed
    habitat = coalesce(habitat.x, habitat.y),
    sampler = coalesce(sampler.x, sampler.y),
    sample_year = coalesce(sample_year.x, sample_year.y),
    transect = coalesce(transect.x, transect.y), 
    abundance = coalesce(abundance.x, abundance.y))  %>%  
  # Remove the CorrectedSpName column as it was only needed for the conditional replacement
  select(-c(CorrectedSpName,SpeciesID_Singapore, Label.codes,
            sample_year.x,sampler.x,habitat.x,transect.x,
            sample_year.y,sampler.y,habitat.y,transect.y,
            abundance.x, abundance.y))


#Clean and combine data sets #### 

#remove species that never occur in the dataset
masters <- removeAbsentSpp(masters)
catherine <- removeAbsentSpp(catherine)
felicity <- removeAbsentSpp(felicity)
trond <- removeAbsentSpp(trond)
slade <- removeAbsentSpp(slade)
gian2022 <- removeAbsentSpp(gian2022)

#make a list of spp found in each dataset
sppMaster <- uniqueSpp(masters) #%>% rename(sppMaster =1)
sppCatherine <- uniqueSpp(catherine) #%>% rename(sppCatherine =1)
sppFelicity <- uniqueSpp(felicity) #%>%  rename(sppFelicity =1)
sppTrond <-  uniqueSpp(trond)# %>% rename(sppTrond =1)
sppSlade <- uniqueSpp(slade) #%>%  rename(sppSlade =1)
sppGian <- uniqueSpp(gian2022) #%>% rename(sppGian =1)

#All species names
allsppName <- rbind(sppMaster,sppCatherine,sppFelicity,sppTrond,sppSlade,sppGian) %>% unique()

#Match all spp names to Eleanor backbone names 
matchingNames <- sppSlade %>% left_join(allsppName) %>% mutate(MatchedNames = spp) 
nonMatchingNames <- allsppName %>% anti_join(sppSlade) %>% mutate(MatchedNames = NA )
allsppNamecurrent <- rbind(matchingNames,nonMatchingNames)

#export names to make AllNamesCrossOver with matching names, manually
#write.csv(allsppName, "Outputs/ManuallyEditAllNamesBackbone.csv")


#COMBINE DATA 
full_df <- rbind(masters,catherine,felicity,trond,slade,gian2022, fill = TRUE) %>%
  ungroup %>% 
  left_join(allsppNameBackbone) %>%
  dplyr::select(-spp) %>%
  rename(spp = MatchedNames)

#add missing month data using custom function
full_df <- addMonthData(full_df)


#Add gps locations ####
#NB need to come back to this as lots of points are certain points are missing GPS data 

full_df <- full_df %>% ungroup() %>% 
  #now combine the new information with the db data, get rid of old columns and rename everything correctly 
 
  #add historic gps points 
   left_join(HistoricGPS, by = c("site", "trap"))%>%
  #add GC 2022 gps points 
  left_join(GC2022_gps, by = c("trap", "site")) %>% 
  mutate(Long = coalesce(Long.x, Long.y),
         Lat = coalesce(Lat.x, Lat.y)) %>% 
  dplyr::select(!c(Long.x,Long.y,Lat.x, Lat.y)) %>%
  
  #add eleanor gps points 
  left_join(slade_GPS, by = "trap", relationship = "many-to-many") %>%  
  mutate(Long = coalesce(Long.x, Long.y),
         Lat = coalesce(Lat.x, Lat.y)) %>% 
  dplyr::select(!c(Long.x,Long.y,Lat.x, Lat.y)) 

names(full_df)

#which traps are missing GPS data? - need to come back to and fill once I get these
#points from other folks 
missingGPS <- full_df %>% dplyr::select(site,transect, trap,Long,Lat) %>% unique 

#all spp
uniqueSpp(full_df)

# Unnest the dataframe and fill missing values with 0
#create full dataframe where all sites are also summarised as zeros if no individual is found 
full_df <-   full_df %>% ungroup() %>% 
  dplyr::select(spp) %>% 
  unique() %>% 
  crossing(full_df %>% ungroup() %>% dplyr::select(site, trap, habitat, day,sample_year,sampler, transect, month) %>% unique()) %>%
  left_join(full_df, by = c("spp", "site", "trap", "day", "transect", "habitat", "sample_year", "sampler", "month")) %>%
  mutate(abundance = ifelse(is.na(abundance), 0, abundance)) 


sum(full_df$abundance)

# Add time since logging and plantation age info ####

#add plantation age info
full_df <- full_df %>% ungroup %>%left_join(plantation, by = c("trap", "site")) 

#Logging age info 

full_df <- full_df %>%mutate(logging_year = case_when(site == "GC1L" ~ 1960,
                                                      site == "GC1L2" ~ 1981,
                                                      site == "GCR1" ~ 1981,
                                                      site == "GCR2" ~ 1991,
                                                      
                                                      # if NOT SURE exact year WHEN LOGGING HAPPENED - ASSUME LOGGING HAPPENED IN 1989
                                                      
                                                      site == "c1" ~ 1988, 
                                                      site == "c17" ~ 1991, 
                                                      site == "c5" ~ 1989, 
                                                      site == "c7" ~ 1987, 
                                                      site == "Coupe81" ~ 1981, 
                                                      site == "Coupe88" ~ 1988, 
                                                      site == "danum1" ~ 1989, 
                                                      site == "danum2" ~ 1989, 
                                                      site == "l2.1" ~ 1988, 
                                                      site == "l2.2" ~ 1988, 
                                                      site == "Malua1" ~ 1989, 
                                                      site == "Malua2" ~ 1989, 
                                                      site == "mbp1" ~ 1985, 
                                                      site == "SAFE LF1" ~ 1989, 
                                                      site == "SAFE LF2" ~ 1989, 
                                                      site == "SAFE LF3" ~ 1989, 
                                                      site == "SAFE VJR" ~ 1989, 
                                                      site == "T21" ~ 1989, 
                                                      site == "t45" ~ 1989, 
                                                      site == "t55" ~ 1989, 
                                                      site == "T81.18" ~ 1981, 
                                                      site == "TNR1" ~ 1989, 
                                                      site == "TNR2" ~ 1989, 
                                                      site == "westus1" ~ 1987, 
                                                      site == "westus2" ~ 1987, 
                                                      site == "X78C" ~ 1989, 
                                                      site == "X88E" ~ 1988, 
                                                      site == "X88H" ~ 1988))

#calculate time since logging 
full_df <- full_df %>% mutate(time_since_logging = sample_year - logging_year) %>%  
  filter(!is.na(spp))


#EXPORT OUTPUTS ####
#save the master Dung beetle data 
write.csv(full_df, "Outputs/dungBeetlesForAbundanceAnalysis.csv")

#Non-manuscript exports #####----------
# #write master copy of Gcerullo 2022 Sampling (need to add plantation GPS points) 
# gian2022ForestData <- full_df %>%
#   filter(sampler == "GC" & sample_year == 2022) %>% 
#   filter(abundance > 0) %>%  
# # filter(!(habitat %in% c("albizia", "eucalyptus"))) %>%  
#   select(-plantation_age) %>% 
#   mutate(SamplingMethod = "HBPT") %>%  
#   left_join(HistoricGPS)
# 
# #export species for Singapore export license
# gianSPP <- gian2022ForestData %>% select(spp, abundance) %>%  
#   group_by(spp) %>% 
#   summarise(count = sum(abundance))


# write.csv(gian2022ForestData,"Outputs/2022_GianlucaForestData_Clean_Upload_v2.csv")
# gian2022ForestData %>%  select(site) %>% unique()
# write.csv(gianSPP,"Outputs/2025_SingaporeExportSpp.csv")
# 
# #gian2022ForestData <- full_df %>%
# gian2022allData <- full_df %>%  
#   filter(sampler == "GC" & sample_year == 2022) %>% 
#   filter(abundance > 0)

#EXTRA UNRELATED##### 
#!!!!TO DO
# #ADD  SAFE GPS POINTS 
# #ADD GIAN2022 PLANTATIONS GPS POINTS [need to correct plantation data]
# 
# #extract data (gps) for crossChecked in Singapore spp  
# 
# correctIDsSingapore <- correctIDsSingapore %>%left_join(full_df, by= c( "CorrectedSpName" = "spp", "site", "trap", 
#                                                       "sample_year", "transect", "abundance", "day", "sampler", "habitat")) %>% 
#   filter(abundance >0) %>%  select(-c(Long, Lat,month))
# 
# gps_filt <- full_df %>% select(site,trap, Long,Lat, transect) %>% unique() %>%  
#   filter(!is.na(Long)) %>%  
#   filter(!is.na(Lat))
# 
# correctIDsSingapore <- correctIDsSingapore %>% left_join(gps_filt) %>% unique() %>% 
#   select(-c(plantation_age,time_since_logging, logging_year))
# 
# write.csv(correctIDsSingapore, "Outputs/SingaporeFullyLabelledSpecimens.csv")
         