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
library(data.table)
library(geosphere)
library(stringr)


# read in functions ####


#this function checks the number of beetles sampled at a trap ina day. If zero individuals total were sampled at a trap then 
#this suggests the trap was never sampled and we assign the abundance as NA  

zeroToNa <- function(x) {
  noSampling <- x %>% group_by(trap,site,day) %>% summarise(sum_ab = sum(abundance)) %>%  
    filter(sum_ab == 0 | is.na(sum_ab)) %>% dplyr::select(site,trap, day) %>% unique %>% cbind(abundanceNA =100000) %>% ungroup()
  x <-   x %>% ungroup() %>% 
    left_join(noSampling, by = c("site", "trap", "day")) %>%
    mutate(abundanceNUM = ifelse(is.na(abundanceNA), 100001, abundanceNA)) %>%
    mutate(abundance = ifelse(abundanceNUM== 100000, NA, abundance)) %>% 
    dplyr::select(-c(abundanceNA, abundanceNUM))
  
  
}


#this function simulates data for a dataset where only day 4 data is missing 
#(use this on Gian2022 data, where plantation data contains genuine zero samples (e.g. where there were no beetles in trap))
simulate4thDay <- function(x){
  
  #REPLACE missing DAY 4 DATA 
  #filter sites with 4th day NAs 
  naday4sites <- x %>% filter(is.na(abundance) & day == "d4") %>% dplyr::select(trap, site) %>% unique()  #calculate mean spp abundance for other 3 days in this site
  simD4 <- x %>% filter(trap %in% naday4sites$trap & site %in% naday4sites$site ) %>%
    mutate_at('abundance', ~replace_na(.,0)) %>% 
    group_by(spp,trap,site) %>% summarise(sim_abundance = mean(abundance)) %>% 
    mutate(sim_abundance = ceiling(sim_abundance)) %>% # round up simulated abundance 
    cbind(day = "d4") %>% 
    ungroup()
  #replace NAs with simulated abundances in correct place
  x <- x %>% 
    left_join(simD4, by = c("spp", "trap", "day","site")) %>% 
    mutate(abundance = ifelse(is.na(abundance), sim_abundance, abundance)) %>% 
    dplyr::select(-sim_abundance)
  x
}

#this function simulates data for all missing days, by assigning mean abundance  from other 3 days sampled to species
simulateMissingDays <- function(x){
  
  #REPLACE missing DAY 4 DATA 
  #filter sites with 4th day NAs 
  naday4sites <- x %>% filter(is.na(abundance) & day == "d4") %>% dplyr::select(trap, site) %>% unique()  #calculate mean spp abundance for other 3 days in this site
  simD4 <- x %>% filter(trap %in% naday4sites$trap & site %in% naday4sites$site ) %>%
    mutate_at('abundance', ~replace_na(.,0)) %>% 
    group_by(spp,trap,site) %>% summarise(sim_abundance = mean(abundance)) %>% 
    mutate(sim_abundance = ceiling(sim_abundance)) %>% # round up simulated abundance 
    cbind(day = "d4") %>% 
    ungroup()
  #replace NAs with simulated abundances in correct place
  x <- x %>% 
    left_join(simD4, by = c("spp", "trap", "day","site")) %>% 
    mutate(abundance = ifelse(is.na(abundance), sim_abundance, abundance)) %>% 
    dplyr::select(-sim_abundance)
  
  #REPLACE missing DAY 3 DATA 
  
  naday3sites <- x %>% filter(is.na(abundance) & day == "d3") %>% dplyr::select(trap, site) %>% unique()  #calculate mean spp abundance for other 3 days in this site
  simD3 <- x %>% filter(trap %in% naday3sites$trap & site %in% naday3sites$site ) %>%
    mutate_at('abundance', ~replace_na(.,0)) %>% 
    group_by(spp,trap,site) %>% summarise(sim_abundance = mean(abundance)) %>% 
    mutate(sim_abundance = ceiling(sim_abundance)) %>% # round up simulated abundance 
    cbind(day = "d3") %>% 
    ungroup()
  #replace NAs with simulated abundances in correct place
  x <- x %>% 
    left_join(simD3, by = c("spp", "trap", "day","site")) %>% 
    mutate(abundance = ifelse(is.na(abundance), sim_abundance, abundance)) %>% 
    dplyr::select(-sim_abundance)
  
  
  #REPLACE missing DAY 2 DATA 
  
  naday2sites <- x %>% filter(is.na(abundance) & day == "d2") %>% dplyr::select(trap, site) %>% unique()  #calculate mean spp abundance for other 3 days in this site
  simD2 <- x %>% filter(trap %in% naday2sites$trap & site %in% naday2sites$site ) %>%
    mutate_at('abundance', ~replace_na(.,0)) %>% 
    group_by(spp,trap,site) %>% summarise(sim_abundance = mean(abundance)) %>% 
    mutate(sim_abundance = ceiling(sim_abundance)) %>% # round up simulated abundance 
    cbind(day = "d2") %>% 
    ungroup()
  #replace NAs with simulated abundances in correct place
  x <- x %>% 
    left_join(simD2, by = c("spp", "trap", "day","site")) %>% 
    mutate(abundance = ifelse(is.na(abundance), sim_abundance, abundance)) %>% 
    dplyr::select(-sim_abundance)
  
  
  #REPLACE missing DAY 1 DATA 
  
  naday1sites <- x %>% filter(is.na(abundance) & day == "d1") %>% dplyr::select(trap, site) %>% unique()  #calculate mean spp abundance for other 3 days in this site
  simD1 <- x %>% filter(trap %in% naday1sites$trap & site %in% naday1sites$site ) %>%
    mutate_at('abundance', ~replace_na(.,0)) %>% 
    group_by(spp,trap,site) %>% summarise(sim_abundance = mean(abundance)) %>% 
    mutate(sim_abundance = ceiling(sim_abundance)) %>% # round up simulated abundance 
    cbind(day = "d1") %>% 
    ungroup()
  #replace NAs with simulated abundances in correct place
  x <- x %>% 
    left_join(simD1, by = c("spp", "trap", "day","site")) %>% 
    mutate(abundance = ifelse(is.na(abundance), sim_abundance, abundance)) %>% 
    dplyr::select(-sim_abundance)
  
  
  x
}

#this function takes data that is pooled and splits into 4 days, where each day gets one quarter of data (this is to enable 
#the combination of data from Eleanor (48hrs) and felicty (pooled across 4 days, as she lost the by-day data)) 
#for species with remaining individuals after dividing by 4, randomly allocate individuals across the 4 days  

splitPooledData <- function(data) {
  allocated_df <- data %>%
    group_by(trap, site, spp) %>%
    mutate(
      Allocated_Abundance = floor(abundance/ 4),  # Allocate abundance by dividing by 4 and taking the floor
      Remaining_Abundance = abundance - Allocated_Abundance*4)   # Calculate the remaining abundance after initial allocation
  
  
  #filter out where remaining abundance to allocate is 1,2 or 3 - these are a special case 
  uneven_df <- allocated_df %>% filter(Remaining_Abundance < 4 & Remaining_Abundance >0)
  
  #for species*traps where abundance is divisible by four, give names for d1,d2,d3,d4
  even_df <- allocated_df %>% 
    #repeat each row 4 times
    slice(rep(row_number(), each = 4)) 
  
  #get days column & combine
  even_days <- data.frame(day = rep(c('d1', 'd2', 'd3', 'd4'), times = nrow(even_df)/4 )) 
  even_df <- even_df %>% cbind(even_days)  %>% 
    dplyr::select(-c(Remaining_Abundance, abundance)) %>% 
    rename(abundance = Allocated_Abundance)
  
  #now we allocate the remaining individuals (3,2,1 ind per spp, for non-divisible traps)
  
  # Add a new column 'day' with randomly assigned values
  uneven_df_days<- data.frame(day = rep(c('d1', 'd2', 'd3', 'd4'), times = nrow(uneven_df)))
  i <- nrow(uneven_df_days) - nrow(uneven_df) #next few lines cbind with diff lenth cols
  JoinNas <-  data.frame(spp = rep(NA, i))
  uneven_df <- uneven_df %>% rbind(JoinNas)
  uneven_df <- uneven_df %>% cbind(uneven_df_days) %>% na.omit()
  uneven_df <- uneven_df %>% dplyr::select(-c(Allocated_Abundance, abundance)) %>% 
    rename(abundance = Remaining_Abundance) 
  
  #combine even and uneven df to give fill disaggregate by day data 
  df_full_split_by_day <- rbind(uneven_df, even_df) 
  
  #check the abundance in the original data is the same as that split by day
  sum(df_full_split_by_day$abundance)
  sum(data$abundance)
  
  #return full data frame 
  df_full_split_by_day
}

#remove species that never occur in the dataset
removeAbsentSpp <- function(x){
  AbsentSpp <- x %>% ungroup() %>% group_by(spp) %>% summarise(sppTot = sum(abundance)) %>% dplyr::select(spp,sppTot)
  x <- x %>% left_join(AbsentSpp, by = "spp") %>% 
    filter(!sppTot == 0 ) %>%  
    dplyr::select(-sppTot)
}

#create dataframe of species present in each dataset
uniqueSpp <- function(x){
  x %>% ungroup()  %>% dplyr::select(spp) %>% unique()
}

#for full dataset after all data is combined, this fills in extra month information
addMonthData <- function(x){
  
  full_df <- x %>% mutate(month = case_when(
    site == "W5" ~ "June",
    site =="W10" ~ "June",
    site == "TNR" ~ "May",
    site == "T21" ~ "May",
    site == "T81.18" ~ "June",
    site == "X88H" ~ "June",
    site == "X78C" ~ "June",
    site == "X88E" ~ "July",
    site == "NorthBenta" ~ "April",
    site =="Benta" ~ "March",
    site == "MuluaFar" ~ "April",
    site == "ElephantRidge" ~ "August",
    site == "tr1" ~ "September",
    site == "tr2" ~ "September",
    site == "brl1" ~ "October",
    site == "brl2" ~ "October",
    site == "danum1" ~ "August",
    site =="danum2" ~ "August",
    site == "INFAPRO" ~ "June",
    site == "mbp1" ~ "Junes",
    site == "tb1" ~ "July",
    site == "tb2" ~ "July",
    site == "coffintrail" ~ "June",
    site == "NorthDanum" ~ "June",
    site == "RhinoRidge" ~ "July",
    site =="westus1" ~ "August",
    site == "westus2" ~ "August",
    site == "c17" ~ "July",
    site == "l2.1" ~ "August",
    site == "l2.2" ~ "August",
    site == "c5" ~ "July",
    site == "c1" ~ "July",
    site == "c7" ~ "July", 
    site == "S1" ~ "July", 
    site == "TNR1" ~ "June", 
    site == "TNR2" ~ "June", 
    site == "Block1" ~ "July", 
    site == "Block2" ~ "June", 
    site == "Block3" ~ "June", 
    site == "Block4" ~ "August", 
    site == "Block5" ~ "October", 
    site == "CatherineW1" ~ "May", 
    site == "CatherineW2" ~ "May", 
    site == "CatherineT1" ~ "May", 
    site == "er1" ~ "August",
    site == "er2" ~ "August", 
    site == "t15" ~ "September",
    site == "t2a" ~ "September",
    site == "rh1" ~ "July", 
    site == "rh2" ~ "July", 
    site == "t45" ~ "September", 
    site == "t55" ~ "September", 
    site == "BentaNew" ~ "October", 
    site == "EC1" ~ "July", 
    site == "EC2" ~ "August", 
    site == "EC3" ~ "August", 
    site == "AL1" ~ "July", 
    site == "AL2" ~ "July", 
    site == "AL3" ~ "August", 
    site == "GC1L" ~ "June", 
    site == "GC1L2" ~ "June",
    site == "GCR1" ~ "July", 
    site == "GCR2" ~ "August", 
    site == "GC2L" ~ "July", 
    site == "GCP" ~ "July",
    TRUE ~  month
  ))
  full_df
}

#pivot the data wider
pivot_fun <- function(x){
  x %>% pivot_longer(cols = !c(spp, habitat),
                     names_to = c("site", "trap","day"),  # seperate "site_trap_day" into three columns..  
                     names_sep = "_",  # 'split by presence of underscore
                     values_to = "abundance")   #and assigning values to a column called "abundance" 
} 

#READ IN DATA ####

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
trond_2l <- read.csv("TwiceLoggedTrondData.csv") %>% slice(-1:-5)

#Eleanor Slade's data ####
#Read in data from Eleanor Slades 2005 BRL and Danum PhD####
slade <- read.csv("Slade_2005_DanumBrL.csv") %>%  rename(traptype = trap, 
                                                         trap = trapID, 
                                                         sampler = author.study) %>% 
  filter(traptype == "BPT") %>% # take only BPT Trap  
  filter(habitat == "Old-growth Forest" | habitat == "Logged Forest") # Take only OG and LF data 











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


#pivot data

trond_felicity <-   trond_felicity %>% pivot_longer(cols = !c(spp),
                                                    names_to = c("habitat","site", "trap"),  # seperate "site_trap_day" into three columns..  
                                                    names_sep = "_",  # 'split by presence of underscore
                                                    values_to = "abundance") %>%    #and assigning values to a column called "abundance" 
  mutate(habitat = case_when(
    habitat == "UL" ~ "primary", 
    habitat == "L" ~ "once-logged", 
    habitat == "R" ~ "restored"))

#now we only dplyr::select Felicity's data
trond_felicity %>% dplyr::select(site) %>%  unique()

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


#view all zero day/trap combinations (none)
FelicityNoSamplling<- felicity %>% group_by(trap,site) %>% summarise(sum_ab = sum(abundance)) %>%  
  filter(sum_ab == 0 | is.na(sum_ab))

#split Felicity's pooled data into 4 days, allocated abundance equally for each spp*trap*site
#to each day sampled & add site = transect
felicity <- splitPooledData(felicity) %>% mutate(transect = site)

sum(felicity$abundance)
#--------------------------------------------------------------------------------------------
#TROND DATA#####
#--------------------------------------------------------------------------------------------

#read in Trond 2wice logged data (split by day)
trond_2l <- read.csv("TwiceLoggedTrondData.csv") %>% slice(-1:-5)

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


#------------------------------------------------------------
# ELEANOR SLADE DATA ####
#------------------------------------------------------------

#Read in data from Eleanor Slades 2005 BRL and Danum PhD####
slade <- read.csv("Slade_2005_DanumBrL.csv") %>%  rename(traptype = trap, 
                                                         trap = trapID, 
                                                         sampler = author.study) %>% 
  filter(traptype == "BPT") %>% # take only BPT Trap  
  filter(habitat == "Old-growth Forest" | habitat == "Logged Forest") # Take only OG and LF data 


slade <- slade %>%  pivot_longer(
  cols = !c(trap, sampler,site,month,year,traptype,habitat),
  names_to = "spp", 
  values_to = "abundance")

sladeSPP <- slade %>% dplyr::select(spp) %>% unique()


#remove data from a fragment (e.g. where site contains the word fragment)
slade <- slade %>% filter(!str_detect(site, "Fragment"))

#reorder and rename #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHECK LOGGED SITES ARE ONCE-LOGGED AND NOT SALVAGE LOGGING 
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

#------------------------------------------------------------------------------------------------
# READ IN GIANLUCA 2022 DUNG BEETLE DATA 
#------------------------------------------------------------------------------------------------
gian2022 <- read.csv("DungBeetle2022DataRawDanum_SSBPlantation.csv") %>% mutate_all(as.character)


gian2022 <- gian2022 %>% pivot_longer(cols = !c(spp),
                                      names_to = c("site", "trap","day"),  # seperate "site_trap_day" into three columns..  
                                      names_sep = "_",  # 'split by presence of underscore
                                      values_to = "abundance")   #and assigning values to a column called "abundance" 


#add a transect column that's same a site
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
#-----------------------------------------------------------------------------------------------
# COMBINE and organise DATASET ####
#-----------------------------------------------------------------------------------------------

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

#Match all spp names to Eleanor names 
matchingNames <- sppSlade %>% left_join(allsppName) %>% mutate(MatchedNames = spp) 
nonMatchingNames <- allsppName %>% anti_join(sppSlade) %>% mutate(MatchedNames = NA )
allsppName <- rbind(matchingNames,nonMatchingNames)

#export names to make AllNamesCrossOver with matching names
#write.csv(allsppName, "AllNamesCrossOver.csv")

#Import species crossover to match names for all species
allsppName <- read.csv("AllNamesCrossOver.csv") %>% dplyr::select(-X)

#COMBINE DATA 
full_df <- rbind(masters,catherine,felicity,trond,slade,gian2022, fill = TRUE) %>% ungroup %>% 
  left_join(allsppName) %>% dplyr::select(-spp) %>% rename(spp = MatchedNames)


#add missing month data using custom function
full_df <- addMonthData(full_df)






#-----------------------------------------------------------------------
#ADD gps locations ####
#-----------------------------------------------------------------------
newInfo <- read.csv("allHistoricSites.csv") %>% dplyr::select(site,trap, Long, Lat)

full_df <- full_df %>% ungroup() %>% 
  #now combine the new information with the db data, get rid of old columns and rename everything correctly 
  left_join(newInfo, by = c("site", "trap"))

#add gian2022 forest point count data 
gianForest2022 <- read.csv("allForest2022PointLocations.csv") %>% rename(Long = X, Lat = Y) %>% 
  dplyr::select(Long,Lat, trap, site)

full_df <- full_df %>% left_join(gianForest2022, by = c("trap", "site")) %>% 
  mutate(Long = coalesce(Long.x, Long.y),
         Lat = coalesce(Lat.x, Lat.y)) %>% 
  dplyr::select(!c(Long.x,Long.y,Lat.x, Lat.y))

names(full_df)
#which traps are missing GPS data?
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

#==========================
##ADD in: 
#1 time since first logging 
#2 plantation age
#===========================
plantation <- read.csv("Plantation_Habitat_Structure.csv") %>% dplyr::select(Site,Point, Age) %>% 
  mutate(Point = paste0("t", Point)) %>% 
  rename(site = Site, 
         trap = Point,
         plantation_age = Age) 

full_df <- full_df %>% ungroup %>%left_join(plantation, by = c("trap", "site")) 



Y <- full_df %>% dplyr::select(trap, site) %>% unique() 
Z <- plantation %>% dplyr::select(trap, site,plantation_age) %>% unique()
NN <- Z %>% left_join(Y, by = c("site", "trap"))

x<- full_df %>% filter(habitat=="once-logged"|habitat=="restored") %>% dplyr::select(site) %>% unique
#ADD IN LOGGING YEAR - NEED TO GET HOLD OF REAL DATA FOR THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!![THIS IS MAE UP ESTIMATES]
full_df <- full_df %>%mutate(logging_year = case_when(site == "GC1L" ~ 1960,
                                                      site == "GC1L2" ~ 1989,
                                                      site == "GCR1" ~ 1989,
                                                      site == "GCR2" ~ 1991,
                                                      
                                                      # NOT SURE WHEN LOGGING HAPPENED - ASSUME LOGGING HAPPENED IN 1989
                                                      
                                                      site == "c1" ~ 1989, 
                                                      site == "c17" ~ 1989, 
                                                      site == "c5" ~ 1989, 
                                                      site == "c7" ~ 1989, 
                                                      site == "Coupe81" ~ 1989, 
                                                      site == "Coupe88" ~ 1989, 
                                                      site == "danum1" ~ 1989, 
                                                      site == "danum2" ~ 1989, 
                                                      site == "l2.1" ~ 1989, 
                                                      site == "l2.2" ~ 1989, 
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
                                                      site == "T81.18" ~ 1989, 
                                                      site == "TNR1" ~ 1989, 
                                                      site == "TNR2" ~ 1989, 
                                                      site == "westus1" ~ 1989, 
                                                      site == "westus2" ~ 1989, 
                                                      site == "X78C" ~ 1989, 
                                                      site == "X88E" ~ 1989, 
                                                      site == "X88H" ~ 1989))

# CALCULATE TIME SINCE LOGGING...AND TIME SINCE RESTORATION...
full_df <- full_df %>% mutate(time_since_logging = sample_year - logging_year)

#save the master Dung beetle data 
write.csv(full_df, "dungBeetlesForAbundanceAnalysis.csv")
#==========================================================================
#plot for each speces how abundance changes 24,48,72,&96 hrs after sampling [ to show Phillip]  
#==========================================================================
#add sampl 
plot_df <- full_df %>% 
  mutate(trap_hours = case_when(
    is.na(day)~ 48, 
    TRUE~ 24))

#plot data without Eleanor's dung beetles 
#plot_df <- full_df 

#call all Slade's data day two for plotting purposes (as cumulatively collected over 48hrs)
plot_df <- plot_df %>%  
  mutate(day = case_when(trap_hours == 48~ "d2", 
                         TRUE~day))

#summarise cumalative abundance
d1 <- plot_df %>% filter(day =="d1") %>%   group_by(spp) %>% mutate(cum_abundance = sum(abundance)) %>% slice(1)  %>% ungroup %>% dplyr::select(spp,day,cum_abundance)
d2 <- plot_df %>% filter(day =="d1"|day =="d2") %>%  group_by(spp) %>% mutate(cum_abundance = sum(abundance)) %>% slice(1) %>% mutate(day = "d2") %>% ungroup %>% dplyr::select(spp,day,cum_abundance)
d3 <- plot_df %>% filter(day =="d1"|day =="d2"|day =="d3") %>%  group_by(spp) %>% mutate(cum_abundance = sum(abundance)) %>% slice(1)  %>% mutate(day = "d3")%>% ungroup %>% dplyr::select(spp,day,cum_abundance)
d4 <- plot_df %>% filter(day =="d1"|day =="d2"|day =="d3"|day =="d4") %>% group_by(spp) %>% mutate(cum_abundance = sum(abundance)) %>% slice(1) %>% mutate(day = "d4")%>% ungroup %>% dplyr::select(spp,day,cum_abundance)

plot_df <- rbind(d1,d2,d3,d4) 

#filter species that we detected >100 times across the dataset 
reg_spp <- plot_df %>% ungroup() %>% group_by(spp)  %>% filter(cum_abundance > 100) %>%dplyr::select(spp) %>% unique
plot_df_filt <- reg_spp %>% left_join(plot_df, by = "spp", relationship = "many-to-many") %>% ungroup %>% mutate(day = factor(day, levels = c("d1", "d2", "d3", "d4")))

#summary stats
plot_df_filt <- plot_df_filt %>% group_by(spp, day) %>% 
  summarize(
    meanAb = mean(cum_abundance),
    se = sd(cum_abundance)/sqrt(n()),
    sd = sd(cum_abundance), 
    minAb = meanAb - se, 
    maxProp = meanAb + se)

plot_df <- reg_spp %>%  left_join(plot_df, by ="spp") %>% group_by(spp, day) %>%  
  ungroup %>% mutate(day = factor(day, levels = c("d1", "d2", "d3", "d4")))

ggplot() +
  #plt raw data
  # geom_point(data = plot_df, aes(x = day, y = abundance), alpha = 0.1, colour = "grey90", fill = "transparent",
  #             position = position_jitter(width = 0.1, height = 0.1),  size =1) + 
  #plt std error
  geom_errorbar(data = plot_df_filt, aes(x = day, y = meanAb, ymin = minAb, ymax = minAb, colour = day), width = 0.2) +
  #plt mean
  geom_point(data = plot_df_filt, aes(day, meanAb, colour= day), size =2)+
  #  ylim(0,1)+
  facet_wrap(~ spp, scales = "free_y") +
  scale_fill_viridis_d(alpha = 0.9)+
  labs(x = "Habitat", y = "Mean cumulative abundance") +
  theme_bw(base_size = 13)+
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom")




#!!!!TO DO
#ADD SLADE GPS POINTS 

#ADD GIAN2022 PLANTATIONS GPS POINTS [need to correct plantation data]


#----------------------------------------------------------------------
#CARRY OUT DATA SUMMARIES ####
#----------------------------------------------------------------------
#how many traps per site? 
siteSummary <- full_df %>% dplyr::select(trap,site, habitat) %>% group_by(habitat,site) %>% dplyr::select(trap) %>% unique() %>% count()

#how many traps per habitat 
siteSummary2 <- full_df %>% dplyr::select(trap,site, habitat) %>% unique %>%  group_by(habitat) %>% count()

sites <- full_df %>% dplyr::select(site) %>% unique()

#for each species, abundance/trap hours 

#spp abundandance by habitat
trap_hrs <- full_df %>% 
  mutate(trap_hours = case_when(
    is.na(day)~ 48, 
    TRUE~ 24
  )) %>%  dplyr::select(habitat,site,trap,day,trap_hours) %>% unique() %>% 
  group_by(habitat) %>% 
  summarise(hab_trap_hours = sum(trap_hours)) %>% ungroup()

hab_abund <- full_df %>% group_by(spp,habitat) %>% summarise(hab_abund = sum(abundance))

hab_abund <- hab_abund %>% left_join(trap_hrs, by = "habitat") %>% 
  mutate(abundance_1hr = hab_abund/hab_trap_hours) %>% 
  mutate(standardised_abundance = abundance_1hr * 100000) # gives us the the abundance we would expect if we sampled 100000 hours in a habitat type

primaryLovers <-  hab_abund %>% group_by(spp) %>% filter(standardised_abundance ==max(standardised_abundance)) %>% filter(habitat == "primary")

#reorder factor order and rename
hab_abund<- hab_abund %>% mutate(habitat = case_when(
  habitat == "albizia" ~ "AL", 
  habitat == "eucalyptus" ~ "EC",
  habitat == "restored" ~ "R", 
  habitat == "once-logged" ~ "1L", 
  habitat == "twice-logged" ~ "2L", 
  habitat == "primary" ~ "P")) %>% 
  mutate(habitat = fct_relevel(habitat, 
                               "EC", "AL", "2L", 
                               "R", "1L", "P")) 

plot_facets <- function(x,i,j){
  x %>% ungroup() %>% 
    slice(i:j) %>% 
    ggplot(aes(x = habitat, y = standardised_abundance,fill = habitat))+
    geom_bar(stat = "identity", position=position_dodge())+
    facet_wrap(~spp, scales = "free_y")+
    ylab("effort_standardised_abundance")+
    scale_fill_manual(values=c('#d73027', '#fc8d59','#fee08b', '#d9ef8b','#91cf60','#1a9850'))+
    theme_pubr(base_size = 9)+
    theme(axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 16))
}


p1 <- plot_facets(hab_abund, 1, 240)
p2 <- plot_facets(hab_abund, 241, 480)
p3 <- plot_facets(hab_abund, 481, 750)


wide_stand <- hab_abund %>% dplyr::select(spp,habitat,standardised_abundance) %>% 
  pivot_wider(names_from = habitat, values_from = standardised_abundance) %>%  
  dplyr::select(spp, eucalyptus, albizia, 7,restored,4,primary) %>% 
  rename(EC = 2, 
         AL = 3, 
         L2 = 4, 
         R = 5, 
         L1 = 6, 
         P = 7)



