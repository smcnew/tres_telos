## Adults co-variate compilation script ##

# A script to assemble and clean telomere dataset for adults
# Sabrina McNew December 2020

library(dplyr)



# load data ---------------------------------------------------------------


ts_results <- read.csv("results/telomere_results.csv") # telomere data by extract #
extractions <- read.csv("covariates/extraction_list.csv") %>% rename(sample_id = extract_label) # Study F = 0 parents / 1 nestlings
phys <- read.csv("covariates/Captures_Hormone_Bleeding_Blood_DNA 11.18.2020.csv", stringsAsFactors = FALSE) %>%
          mutate(Site = case_when(Site == "Unit 1" ~ "Unit_1",
                                  Site == "Unit 2" ~ "Unit_2",
                                  TRUE ~ Site)) %>%
          mutate(band_date = paste(Individual_Band, Capture_Date, sep = "_"),
          site_nest = paste(Site, Nest, Exp_Year, sep = "_"))

nests <- read.csv("covariates/Nest_Records 11.18.2020.csv", stringsAsFactors = FALSE) %>%
          mutate(Site = case_when(Site == "CU1" ~ "Unit_1",
                                  Site == "CU2" ~ "Unit_2",
                                  Site == "CU4" ~ "Unit_4",
                                  TRUE ~ Site)) %>%
          mutate(site_nest = paste(Site, Nest, Exp_Year, sep = "_"))


# data compilation --------------------------------------------------------

# keep entire extraction list, remove standards, golds, etc., filter to just adults
comp_data <- merge(extractions, ts_results, all.x = T, all.y = F) %>%
  filter(study_f == 0) %>%
  mutate(band_date = paste(band, Capture_Date, sep = "_"))

comp_data <- merge(comp_data, phys, all.x = T ) # add phys
comp_data %>% #select(Individual_Band, band) %>% # check that phys$Individual_band matches extractions$band
  mutate(check = Individual_Band == band) %>% #summary(check) #9 NAs
  filter(is.na(check))  %>% select(Capture_Date, sample_id, band, year, day)

dim(comp_data)
comp_data <- merge(comp_data,
            select(nests, -c("Location","Site","Nest","Species","Exp_Year","Fecal_Sample_nb","Notes")),
            all.x = T, all.y = F)



# Consolidate treatments
comp_data <- comp_data %>% mutate(treatment = case_when(Individual_Treatment == "" ~ "control",
                                           Individual_Treatment == "Color_Dull" ~ "control",
                                           Individual_Treatment == "Color_Sham" ~ "control",
                                           Individual_Treatment == "control" ~ "control",
                                           Individual_Treatment == "Control" ~ "control",
                                           Individual_Treatment == "Control_Control" ~ "control",
                                           Individual_Treatment == "Control_Predator" ~ "control",
                                           Individual_Treatment == "Control_Tape" ~ "control",
                                           Individual_Treatment == "CORT_3d" ~ "cort",
                                           Individual_Treatment == "CORT_6d" ~ "cort",
                                           Individual_Treatment == "CORT_XFoster" ~ "cort",
                                           Individual_Treatment == "DEX_Treatment" ~ "cort",
                                           Individual_Treatment == "DMSO" ~ "sham_control",
                                           Individual_Treatment == "DMSO_6d" ~ "sham_control",
                                           Individual_Treatment == "DMSO_XFoster" ~ "sham_control",
                                           Individual_Treatment == "Dull_Control" ~ "control",
                                           Individual_Treatment == "Dull_Predator" ~ "control",
                                           Individual_Treatment == "Dull_Tape" ~ "control",
                                           Individual_Treatment == "high" ~ "cort",
                                           Individual_Treatment == "High_Tape" ~ "control",
                                           Individual_Treatment == "long" ~ "cort",
                                           Individual_Treatment == "low" ~ "cort",
                                           Individual_Treatment == "Low_Tape" ~ "control",
                                           Individual_Treatment == "Medium_Tape" ~ "control",
                                           Individual_Treatment == "None" ~ "control",
                                           Individual_Treatment == "Predator" ~ "control"
                                           ))

table(comp_data$treatment)


# data exploration --------------------------------------------------------

table(comp_data$Site)

head(comp_data)
boxplot(ts_ratio4 ~ treatment, comp_data)
# Zuur-style checks for weird data (based on Zuur et al. 2010 Methods Ecol Evol)
# 1. Look for outliers in xs and ys using boxplots and cleveland dotplots
# Identify key variables
response <- c(
  "number_chicks_fledged_from_rear_nest",
  "day_14_tarsus_length",
  "day_14_weight",
  "chick_survival_to_first_breed_season"
)
par(mfrow = c(4,2))
for(i in 1:length(response)){
  tit2 <- tit
  tit2$chick_survival_to_first_breed_season <- as.numeric(as.character(tit$chick_survival_to_first_breed_season))
  dplyr::select(tit2, response[i]) %>% boxplot(., main = response[i]) #boxplot
  dplyr::select(tit2, response[i]) %>%
    pull(.) %>% plot (x = ., y = 1:length(.), ylab="order of data", xlab= response[i]) #Dotplot index ~ value
}

# Response variables look pretty good, decent spread, few outliers. One baby has a very small tarsus ( < 13)

# Now check predictors, including treatment variables and covariates
predictors <- c(
  "net_rearing_manipulation",
  "rear_nest_trt",
  "rear_Cs_at_start_of_rearing",
  "d14_rear_nest_brood_size",
  "hatch_nest_LD"
)
par(mfrow = c(5,2), mar=c(4,2,2,0))
for(i in 1:length(predictors)){
  tit2 <- tit
  tit2$rear_nest_trt <- as.numeric(as.character(tit2$rear_nest_trt))
  dplyr::select(tit2, predictors[i]) %>% boxplot(., main = predictors[i]) #boxplot
  dplyr::select(tit2, predictors[i]) %>%
    pull(.) %>% plot (x = ., y = 1:length(.), ylab="order of data", xlab= predictors[i]) #Dotplot index ~ value
}

dplyr::select(tit, rear_nest_trt) %>% pull %>% as.character %>% as.numeric %>% boxplot()
# weird patterns in LD because sheet was sorted in chronological order

# 2. Look for homogeneity of variance: sure, looks like trt 5 is smaller than others
dev.off()
boxplot(day_14_tarsus_length ~ rear_nest_trt, tit) #boxes are about the same size, more outliers in 5
boxplot(day_14_weight ~ rear_nest_trt, tit) #boxes are about the same size, more outliers in 5

# 3. Check for normality, why not
qqnorm(tit$day_14_tarsus_length[tit$rear_nest_trt==5]) # 5 and 7 not super normal but let's move on
qqline(tit$day_14_tarsus_length[tit$rear_nest_trt==5])


# Check for covariance among predictors
head(tit)
tit %>% select(predictors, response) %>% pairs() # clutch, brood, and fledglings correlated

#