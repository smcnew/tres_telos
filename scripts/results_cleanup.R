# QC script
# A script to remove bad triplicates from a run and calculate cq_mean and cq_sd
# Function asks for path to file, the date of the PCR run, and the primer (e.g. telo or gadph)
# There will be a warning message related to NAs created by removing "undetermined" from Cq results
library(dplyr)


# Set file and indicate where standards and gold are- double check with lab book
  file <- "2020_TRES_Telomere_Labwork/qpcr_results_raw/2020-12-04_telo.csv"
  run <- read.csv(file)
  date  <- strsplit(file, split = "/")[[1]][[3]] %>%
    strsplit(., "_") %>% unlist() %>% .[1]
  primer <-  strsplit(file, split = "/")[[1]][[3]] %>%
    strsplit(., "_") %>% unlist() %>% .[2] %>% sub(".csv", "", .)

  T_087 <- 21
  h20 <- 122
  std1 <- 123
  std2 <- 124
  std3 <- 125
  std4 <- 126
  std5 <- 127
  gold <- 128

# CLEAN UP
run$Cq <-run$Cq %>% as.character() %>% as.numeric #Change "undetermined" to NAs, this will create a warning message

# Filtering method
# 1. For all triplicates where sd < 0.25, calculate mean and sd based on all 3 replicates
# 2. For triplicates where sd > 0.25, calculate distance matrix between all three.
# Take the lowest sd, and if it's < 0.25, throw out the third point.

results <-
  data.frame(
    pcr_date = date,
    pcr_primer = primer,
    sample = sort(unique(run$Sample)),
    sample_id = "",
    Cq = NA,
    Cq_sd = NA,
    nsamps = NA
  ) # Data frame to save processed results in

for (i in 1:nrow(results)) {
  test <- dplyr::filter(run, Sample == i)
  standard_dev <- sd(test$Cq, na.rm = T) #what is the starting standard dev of the triplicate
  standard_dev <- ifelse(is.na(standard_dev), 1, standard_dev) #NAs cause errors; change to 1
  if(standard_dev <= 0.25) { #if sd < 0.25 calculate mean and sd on all 3
    results$Cq[i] <- mean(test$Cq, na.rm = T)
    results$Cq_sd[i] <- sd(test$Cq, na.rm = T)
    results$nsamps[i] <- sum(!is.na(test$Cq))
  }

  else { #find a bad triplicate
    sd_matrix  <- sqrt(((dist(test$Cq)/2)^2)*2) # manually calculate sd in distance matrix form
    if (all(is.na(sd_matrix)) == TRUE) {        # if all triplicates failed, return NA
      good_samples <- NA
    }
    else {
    best_pair <- which(sd_matrix == min (sd_matrix, na.rm = T)) * # find the pair with the lowest sd
      as.numeric(min (sd_matrix, na.rm = T) < 0.25) # make sure that the lowest sd is within 0.25, return 0 if not.

    good_samples <-
      case_when(best_pair == 1 ~ test$Cq[-3], # exclude c
                best_pair == 2 ~ test$Cq[-2], # exclude b
                best_pair == 3 ~ test$Cq[-1], # exclude a
                best_pair == 0 ~ as.numeric(NA)) #all are bad
    }
    results$Cq[i] <- mean(good_samples)
    results$Cq_sd[i] <- sd(good_samples)
    results$nsamps[i] <- case_when (all(is.na(good_samples)) == TRUE ~ as.numeric(NA),
                                    all(is.na(good_samples)) == FALSE ~ 2) # add number of samples used in calculation, NA if all were thrown out
  }
}

results[1:21,]

# If gold or 087 fails QC, inspect manually:
# gold_cq <- filter(run, Sample == gold) %>% select(Cq) %>% pull %>% .[2]

gold_cq <- results$Cq[results$sample == gold]
std3_cq <- results$Cq[results$sample == std3 ]
T_087_cq <- results$Cq[results$sample == T_087]

results$Cq_cv <- (results$Cq_sd / results$Cq) * 100
results$Cq_relative <- results$Cq / gold_cq
results$Cq_relative_sd3 <- results$Cq / std3_cq
results$Cq_relative_087 <- results$Cq / T_087_cq


# Calculate efficiency

# Slope of standards vs. concentration
# 1. Using just the averaged results
concentrations <- c(9.27, 4.78, 2.42, 0.993, .44) # concentrations of standards C
#concentrations <- c(10.3, 4.73, 1.79, 0.76, 0.3) # concentrations of standards B

cqs <- results$Cq[std1:std5]
slope <- lm(cqs ~ log10(concentrations))$coefficients[2] %>% as.numeric()
efficiency <-  (10 ^ (-1/slope))-1; efficiency

# 2. Using all standard raw values
all_cqs <- filter(run, Sample %in% std1:std5) %>% arrange(Sample) %>% select(Cq) %>% pull()
all_con <- rep(concentrations, 3) %>% sort(decreasing = T)
slope2 <- lm(all_cqs ~ log10(all_con))$coefficients[2] %>% as.numeric()
efficiency2 <-  (10 ^ (-1/slope2))-1; efficiency2

# 3. Using expected concentration (not qubited values)
concentrations <- c(10, 10/2, 10/4, 10/8, 10/16)
#concentrations <- c(10, 10/2, 10/4, 10/8) #exclude last standard

all_con <- rep(concentrations, 3) %>% sort(decreasing = T)
slope3 <- lm(all_cqs ~ log10(all_con))$coefficients[2] %>% as.numeric()
efficiency3 <- (10 ^ (-1/slope3))-1

print(noquote(c("Efficiency 1:", round(efficiency, 3)*100)))
print(noquote(c("Efficiency 2:", round(efficiency2, 3)*100)))
print(noquote(c("Efficiency 3:", round(efficiency3, 3)*100)))


sum(is.na(results$Cq[1:20])) # how many samples failed QC
dplyr::filter(run, Sample == h20) %>% select(Cq) #inspect water, should be NA
date #double check date and primer are right
primer
tail(results)
write.csv(results, paste("2020_TRES_Telomere_Labwork/processed_results/", date, "_", primer, "_processed.csv", sep = ""), row.names = FALSE)
