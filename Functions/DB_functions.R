
# FUNCTIONS for dung beetle processing ####

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(forcats)
library(data.table)
library(geosphere)
library(stringr)

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


