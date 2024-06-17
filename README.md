

#1. CombineDungBeetleDatasets.R
This code combines multiple dung beetle datasets collected in Sabah by GC, ES, TL, CF, FE and 
also adds relevant information on habitats and ages, while applying a common taxonomic backbone.


#2. FormatDBForNegBinom.R
This takes the output csv from 1, and formats the data so that it is prepared for model fitting 

#3. FitModel.R 
This fits a zero-inflated negative binomial model to the output from 2, which allows us to predict 
species level abundances in different habitat types of different ages. Here we carry out model checks.  

#4. PredictAbundanceByHab.R
This calculates abundance for each species and habitat and age combo, using the model from 3. We do 
not interpolate or extrapolate beyond the chronosequence of age*habitat that we sampled in the field - 
here instead we apply certain rules to fill in age-gradient gaps. 

#5.ProcessScenarioOutcomes.R
Using the outputs from 4 we then calculate scenario-level abundances of species, and produce our summary
dung beetle performance outputs. 