library(tidyverse)
library(summarytools)

#Data for Juliette MEP project test scenario

df<- read.csv("Outputs/dungBeetlesForAbundanceAnalysis.csv")

df_gc <- df %>% filter(sampler %in% c("GC", "TL", "FE")) %>% 
  select(-sampler)

df_fe <- df %>% filter(sampler %in% c("FE")) %>% 
  select(-sampler) %>% 
  filter(abundance >0)

#write.csv(df_gc, "Outputs/DBs_MEPproject.csv")


#explore data 
summary(df_gc)
view(dfSummary(df_gc))

write.csv(df_fe, 'Outputs/MEP_Edwards_Dung_Beetles.csv')

summary(df_fe)
view(dfSummary(df_fe))
# Step 1: Define column names and types
column_types <- list(
  species_name = "character",
  abundance = "numeric",
  forest_type = "factor",
  sampling_year = "integer",
  latitude = "numeric",
  longitude = "numeric",
  site = "character",
  transect = "character",
  point = "character"
)

# Step 2: Create a template dataframe
template_df <- tibble(
  species_name = character(),
  abundance = numeric(),
  forest_type = factor(levels = c("primary", "once-logged", "restored", "twice-logged")),
  sampling_year = integer(),
  latitude = numeric(),
  longitude = numeric(),
  site = character(),
  transect = character(),
  point = character()
)

# Step 3: Generalized standardization function
standardize_data <- function(data, rename_map, column_types, template_df) {
  # Debug: Check initial column names
  message("Initial column names in dataset: ", paste(names(data), collapse = ", "))
  
  # Rename columns based on the provided mapping
  data <- data %>%
    rename_with(~ rename_map[.x], .cols = names(rename_map))
  
  # Debug: Check after renaming
  message("Column names after renaming: ", paste(names(data), collapse = ", "))
  
  # Select only columns matching the template and fill missing columns with NA
  standardized_data <- template_df %>%
    bind_rows(data) %>%
    select(names(column_types)) %>%
    mutate(across(
      everything(),
      ~ {
        col_type <- column_types[[cur_column()]]
        # Handle factor levels explicitly if needed
        if (col_type == "factor" && is.character(.)) {
          factor(., levels = levels(template_df[[cur_column()]]))
        } else {
          as(., col_type)
        }
      }
    ))
  
  # Debug: Check the structure of standardized data
  message("Final standardized data structure:")
  glimpse(standardized_data)
  
  return(standardized_data)
}

# Step 4: Define a rename map for your dataset
rename_map <- c(
  spp = "species_name",
  abundance = "abundance",
  habitat = "forest_type",
  sample_year = "sampling_year",
  Lat = "latitude",
  Long = "longitude",
  site = "site",
  transect = "transect",
  trap = "point"
)

# Step 5: Standardize the dataset
standardized_df <- standardize_data(df_gc, rename_map, column_types, template_df)
view(dfSummary(standardized_df))

# Step 6: Identify sites where longitude and latitude are NA
sites_with_missing_coords <- standardized_df %>%
  filter(is.na(longitude) | is.na(latitude)) %>% # Find rows with missing coordinates
  distinct(site) # Extract unique sites

#7. Remove data we don't want 
standardized_df <- standardized_df  %>%  
  filter(!is.na(forest_type)) # i.e. plantations 

# View the standardized dataframe
print(standardized_df)
