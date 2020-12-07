
library(dplyr)
library(reshape2)
library(tidyr)
library(rptR)
#Data
telo <- read.csv("telo_results.csv")
head(telo)
filter(telo, sample == "Gold", primer == "gadph")


aggregate(Cq_sd ~ primer, median, data =  telo)
xduplicated(telo$sample)

replicated <- telo %>% filter(primer == "telo") %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter( !sample %in% c("Gold", "H20")) %>% data.frame() %>%
  pivot_wider(id_cols =sample, names_from = pcr_date, values_from = Cq_mu) %>% as.data.frame()
replicated$first <- apply(replicated[,2:5], 1, function(x) sum(x, na.rm= T))
replicated$first[replicated$first > 20 ] <- NA


plot(`31-Jul` ~ `28-Jul`, data = replicated_gadph )

replicated_gadph %>% select(`31-Jul`,first )%>% na.omit() %>% plot()
replicated_gadph %>% select(`31-Jul`,`28-Jul` )%>% na.omit() %>% plot(, main = "Gadph")


replicated_gadph <- telo %>% filter(primer == "gadph") %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter( !sample %in% c("Gold", "H20")) %>% data.frame() %>%
  pivot_wider(id_cols =sample, names_from = pcr_date, values_from = Cq_mu) %>% as.data.frame()
replicated_gadph$first <- apply(replicated_gadph[,2:5], 1, function(x) sum(x, na.rm= T))
replicated_gadph$first[replicated_gadph$first > 30 ] <- NA
plot(`31-Jul` ~ first, data = replicated_gadph[replicated_gadph$`31-Jul`!= 0,] )

std_conc <- c(12.7, 5.6, 2.76, 1.58, 1.36)*3
#amplification calculator use raw Cq not relative

lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),5] ~ log10(std_conc)) %>% summary()


par(mfrow=c(3,2))
plot(x = log(std_conc), y = replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),2],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),3],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),4],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),5],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),6],
     ylab = "cq_rel")

# AMPLIFICATION EFFICIENCIES: telo 77 - 115%
#lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4" ),2] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4" ),3] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
#lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4" ),4] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4" ),5] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated[replicated$sample %in% c("Std1", "Std2","Std3", "Std4" ),6] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]



par(mfrow=c(3,2))
plot(x = log(std_conc), y = replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),2],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),3],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),4],
     ylab = "cq_rel")
plot(x = log(std_conc), y = replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4", "Std5"),5],
     ylab = "cq_rel")
plot(x = log(std_conc)[1:4], y = replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),5],
     ylab = "cq_rel")


# AMPLIFICATION EFFICIENCIES: GADPH 83- 93%
lm(replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),2] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),3] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),4] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),5] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]
lm(replicated_gadph[replicated_gadph$sample %in% c("Std1", "Std2","Std3", "Std4" ),6] ~ log10(std_conc)[1:4]) %>% summary() %>% .$coefficients %>% .[2,1]




# Repeatabilities ---------------------------------------------------------
library(lme4)
## Repeatability GADPH of sample set
gadph_repeats <- telo[telo$primer == "gadph",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., !sample %in% c("Gold", "H20", "Std1","Std1","Std2","Std3", "Std4", "Std5", "golden_2")) %>%
  as.data.frame() %>% droplevels()  #group_by(sample) %>% tally() %>% as.data.frame() # 2 samples per group

rpt(Cq_mu ~ (1|sample), grname = "sample", data = gadph_repeats, datatype = "Gaussian", nboot = 100, npermut = 0) #low repeatability with all samples
rpt(Cq_mu ~ (1|sample), grname = "sample", data = gadph_repeats[gadph_repeats$pcr_date %in% c("28-Jul", "31-Jul"),],
    datatype = "Gaussian", nboot = 100, npermut = 0) #98% repeatability
rpt(Cq_mu ~ (1|sample), grname = "sample", data = gadph_repeats[gadph_repeats$pcr_date %in% c("21-Jul", "31-Jul"),],
    datatype = "Gaussian", nboot = 100, npermut = 0) #0 repeatability

## Repeatability telo of sample set (VERY LOW)
telo_repeats <- telo[telo$primer == "telo",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., !sample %in% c("Gold", "H20", "Std1","Std1","Std2","Std3", "Std4", "Std5", "golden_2")) %>%
  as.data.frame() %>% droplevels()  #group_by(sample) %>% tally() %>% as.data.frame() # 2 samples per group

rpt(Cq_mu ~ (1|sample), grname = "sample", data = telo_repeats, datatype = "Gaussian", nboot = 100, npermut = 0) # 0.09
rpt(Cq_mu ~ (1|sample), grname = "sample", data = telo_repeats[telo_repeats$pcr_date %in% c("28-Jul", "31-Jul"),],
    datatype = "Gaussian", nboot = 100, npermut = 0) #15% repeatability
rpt(Cq_mu ~ (1|sample), grname = "sample", data = telo_repeats[telo_repeats$pcr_date %in% c("21-Jul", "31-Jul"),],
    datatype = "Gaussian", nboot = 100, npermut = 0) #0 repeatability

# plots
telo_repeats %>% pivot_wider(id_cols = sample, names_from = pcr_date, values_from = Cq_mu) %>% as.data.frame() %>%
  select(`31-Jul`, `21-Jul`) %>%
  na.omit() %>% plot()

# Repeatability of golden samples

telo[telo$primer == "gadph",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., sample %in% c("Gold", "H20")) %>% as.data.frame() %>% droplevels() %>%
rpt(Cq_mu ~ (1|sample), grname = "sample", data = ., datatype = "Gaussian", nboot = 100, npermut = 0) ## 99% repeatability vs. water

telo[telo$primer == "gadph",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., sample %in% c("Gold", "golden_2")) %>% as.data.frame() %>% droplevels() %>%
  rpt(Cq_mu ~ (1|sample), grname = "sample", data = ., datatype = "Gaussian", nboot = 100, npermut = 0) ## 77% repeatability vs. golden 2


telo[telo$primer == "gadph",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., sample %in% c("Std1","Std1","Std2","Std3", "Std4", "Std5")) %>% as.data.frame() %>% droplevels() %>%
  rpt(Cq_mu ~ (1|sample), grname = "sample", data = ., datatype = "Gaussian", nboot = 100, npermut = 0) ## 94% repeatability vs. golden 2


telo[telo$primer == "telo",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., sample %in% c("Gold", "golden_2")) %>% as.data.frame() %>% droplevels() %>%
  rpt(Cq_mu ~ (1|sample), grname = "sample", data = ., datatype = "Gaussian", nboot = 100, npermut = 0) ## 12% repeatability vs. golden 2

telo[telo$primer == "telo",] %>%
  group_by(sample) %>% filter(n() > 1) %>%
  filter(., sample %in% c("Std1","Std1","Std2","Std3", "Std4", "Std5")) %>% as.data.frame() %>% droplevels() %>%
  rpt(Cq_mu ~ (1|sample), grname = "sample", data = ., datatype = "Gaussian", nboot = 100, npermut = 0) ##19% repeatability vs. golden 2

testrpt <- read.delim(pipe('pbpaste'))

rpt(Cq_mu ~ (1|sample), data = testrpt, grname = "sample", datatype = "Gaussian", nboot = 100, npermut= 0)


testdat <- read.table(pipe("pbpaste"))
colnames(testdat) <- c("sample", "cq", "cq_rel")
rpt(cq_rel ~ (1|sample), data = testdat, grname = "sample", datatype = "Gaussian", nboot= 100, npermut = 0)

dim(testdat)
colnames(testdat) <- c("cq", "round")
testdat$sample <- rep(1:118,2)
rpt(cq ~ (1|sample), data = testdat, grname = "sample", datatype = "Gaussian", nboot = 100, npermut = 0)
