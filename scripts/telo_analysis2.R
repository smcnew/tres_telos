# Telo_analysis November

library(reshape2)
library(tidyr)
library(rptR)
library(dplyr)
# Read files

file_names <- dir("2020_TRES_Telomere_Labwork/processed_results/use_for_analysis/", full.names = T) #where you have your files
results <- do.call(rbind,lapply(file_names, function (x) read.csv(x, stringsAsFactors=F)))
results <- results[!is.na(results$Cq),] # drop rows that didn't pass QC
results$pcr_plate <- paste(results$pcr_date, results$pcr_primer, sep = "_")

results %>% select(pcr_date) %>% unique
plate_metadata <- read.csv("plate_metadata.csv")

#results %>% filter(sample_id %in% c("gold", "T_087", "Std3")) %>% select(pcr_plate, Cq, sample_id) %>%
#  pivot_wider(id_cols = pcr_plate, names_from = sample_id, values_from = Cq) %>% write.csv("plate_metadata_new.csv", row.names = FALSE)

results <- merge(results, plate_metadata)

# Inspect calibration samples to identify most consistent ones
results %>% filter(sample_id == "gold", pcr_primer == "gadph") %>% select(Cq) %>% pull %>% sd # 0.47
results %>% filter(sample_id == "Std3", pcr_primer == "gadph") %>% select(Cq) %>% pull %>% sd # 0.22
results %>% filter(sample_id == "T_087", pcr_primer == "gadph") %>% select(Cq) %>% pull %>% sd # 0.37
gadph_standards <- results %>% filter(pcr_primer == "gadph") %>% filter(sample_id %in% c("gold","Std3","T_087"))

results %>% filter(sample_id == "gold", pcr_primer == "telo") %>% select(Cq) %>% pull %>% sd # 0.50
results %>% filter(sample_id == "Std3", pcr_primer == "telo") %>% select(Cq) %>% pull %>% sd # 0.56
results %>% filter(sample_id == "T_087", pcr_primer == "telo") %>% select(Cq) %>% pull %>% sd # 0.35

telo_standards <- results %>% filter(pcr_primer == "telo") %>% filter(sample_id %in% c("gold","Std3","T_087"))

# Efficiency summaries: "efficiency_a" is using observed (qubited concentration); "efficiency_mean" is taking the mean of the efficiencies.
results %>% filter(pcr_primer == "telo") %>% select(pcr_date, efficiency_a, efficiency_mean) %>% unique %>%
  summarize(across(c(efficiency_a, efficiency_mean), mean))


results %>% filter(pcr_primer == "gadph") %>% select(pcr_date, efficiency_a, efficiency_mean) %>% unique %>%
  summarize(across(c(efficiency_a, efficiency_mean), mean))

# Repeatability, just standards
rpt(Cq ~ (1|sample_id), data = telo_standards, grname = "sample_id", datatype = "Gaussian", nboot= 100, npermut = 0) # Repeatability = 0.879
rpt(Cq ~ (1|sample_id), data = gadph_standards, grname = "sample_id", datatype = "Gaussian", nboot= 100, npermut = 0) # Repeatability = 0.77 (but they're all very similar)

# All samples
results %>% filter(pcr_primer == "telo") %>% group_by(sample_id) %>% filter(n() > 1) %>% # grab duplicated rows
  rpt(Cq ~ (1|sample_id), data = ., grname = "sample_id", datatype = "Gaussian", nboot= 100, npermut = 0) # run repeatability = 0.80

results %>% filter(pcr_primer == "gadph") %>% group_by(sample_id) %>% filter(n() > 1) %>% # grab duplicated rows
  rpt(Cq ~ (1|sample_id), data = ., grname = "sample_id", datatype = "Gaussian", nboot= 100, npermut = 0) # run repeatability = 0.72

dup_gadph <- results %>% filter(pcr_primer == "gadph") %>% group_by(sample_id) %>% filter(n() > 1) %>% droplevels()
table(dup_gadph$sample_id)
dup_gadph %>% as.data.frame %>% select(pcr_date) %>% unique


# T/S Ratio
# T/S ratio = 2^−(Cq_telo – Cq_telo_golden] – [Cq_gadph – Cq_gadph_golden)])
# T/S ratio = (Efficiency_telo)^(Cq_telo – Cq_telo_golden]/(Efficiency_gadph)^ [Cq_gadph – Cq_gadph_golden)])
#

ts_results <- results %>% select(pcr_primer, sample_id, Cq, Std3, gold, T_087, efficiency_a, efficiency_mean) %>%
  distinct(sample_id, pcr_primer, .keep_all = TRUE ) %>%
  pivot_wider(id_cols = sample_id, names_from = pcr_primer, values_from = c(Cq, Std3, gold, T_087, efficiency_a, efficiency_mean)) %>% as.data.frame

ts_results <- mutate(ts_results, ts_ratio1 = 2^-((Cq_telo - T_087_telo)-(Cq_gadph - gold_gadph))) # T/S ratio method Morinha #1
ts_results <- mutate(ts_results, ts_ratio2 = ((efficiency_a_telo/100)^(Cq_telo - T_087_telo))/((efficiency_a_gadph/100)^(Cq_gadph - gold_gadph))) # T/S ratio method 2, efficiency
ts_results <- mutate(ts_results, ts_ratio3 = 2 ^ ((T_087_telo-gold_gadph) - (Cq_telo - Cq_gadph))) #T/S ratio Britt
ts_results <- mutate(ts_results, ts_ratio4 = ((telo_eff/100)^(T_087_telo-Cq_telo))/((gadph_eff/100)^(gold_gadph-Cq_gadph))) # T/S ratio method 2


telo_eff <- results %>% filter(pcr_primer == "telo") %>% select(pcr_date, efficiency_a, efficiency_mean) %>% unique %>%
  select(efficiency_a) %>% pull %>% mean

gadph_eff <- results %>% filter(pcr_primer == "gadph") %>% select(pcr_date, efficiency_a, efficiency_mean) %>% unique %>%
  select(efficiency_a) %>% pull %>% mean

ts_results$gadph_eff_all <- gadph_eff
ts_results$telo_eff_all <- telo_eff


test <- ts_results$ts_ratio1 == ts_results$ts_ratio3
na.omit(test) %>% all()

plot(ts_ratio4 ~ ts_ratio3, ts_results)

head(ts_results)
plot(ts_ratio1 ~ ts_ratio2, ts_results)

write.csv(ts_results, "ts_results.csv")
ts_results %>% arrange(sample_id) %>% write.csv("telomere_results.csv")
ts_results$Cq_telo %>% fivenum

