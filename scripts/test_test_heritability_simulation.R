## This script is written to test an issue with traits that have low heritability in 
  # univariate models having moderate support for heritability when run in bivariate
  # models. 

  # Conor Taff. Last updated 28 January 2020

  # Run in RStudio with R version 1.1.463
  # R version 3.6.1

## Load packages
  library(MCMCglmm)         # running version 2.29
  library(MasterBayes)      # running version 2.55
  
## Make a simulated data set. This is similar in structure to a real dataset that I am
  # working with. Here there are nestlings from families of 3, 4, or 5. Extrapair
  # paternity is around 50%, so that nest-mates can be full or half siblings. Before
  # hatching about 30% of eggs were swapped between nests so that some full and
  # half siblings were raised in different nests. Parentage of all nestlings in 
  # this simulation is known (in the real dataset some nestlings have unknown fathers
  # and are not included in heritability estimation). For purposes of this example
  # no other covariates are included.
  
  # Set a random seed to get the same simulation results
    #set.seed(135)
    # Note that with a smallish sample size similar to what we have in the study,
    # the random seed chosen here makes a pretty big difference and you can get
    # apparently pretty high heritability even when the 'true' value is set very
    # low. Presumably this would be less of an issue with + sample, EPP, and swaps.
    
  ## Set parameters for simulation
    # Level of genetic correlation between two traits
      set_rho <- 0
    # Heritability of trait 1
      set_h2_1 <- 0.15
    # Heritability of trait 2
      set_h2_2 <- 0.1
    # Percent influence of nest environment on nestling phenotype trait 1
      set_Ve_1 <- 0.4
    # Percent influence of nest environment on nestling phenotype trait 2
      set_Ve_2 <- 0
    # Correlatoin between effect of nest on traits 1 and 2
      set_rho_nest <- .0
    # Number of nests (determines number of pairs of parents)
      set_num <- 150
    # Min and max nestlings per brood
      brood_min <- 3
      brood_max <- 6
    # Rate of extrapair paternity
      epp_rate <- 38
    # Percent of nestlings swapped into different nests
      swap_rate <- 30
      
  ## Set parameters for markov chain models
      nitt_bi_sim <- 1e5                # iterations for bivariate models
      thin_bi_sim <- 100                # thinning for biviarate models
      burn_bi_sim <- 5e4                # burnin for bivariate models
      nitt_un_sim <- 5e4                # iterations for univariate models
      thin_un_sim <- 100                # thinning for univariate models
      burn_un_sim <- 1e4                # burnin for univariate models     
    
  # Function to make random variables correlated at fixed amount
    complement <- function(y, rho, x) {
      if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
      y.perp <- residuals(lm(x ~ y))
      rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
  
  # Simulate parents
    d_parents <- as.data.frame(as.factor(paste0(seq(1, set_num), "_nest")))
    colnames(d_parents) <- "nest"
    d_parents$mother <- as.factor(paste0(seq(1, nrow(d_parents), 1), "_mom"))
    d_parents$father <- as.factor(paste0(seq(1, nrow(d_parents), 1), "_dad"))
    d_parents$clutch <- round(runif(nrow(d_parents), min = brood_min, max = brood_max))
    d_parents$mother_tr1 <- rnorm(nrow(d_parents), 0, 1)
    d_parents$mother_tr2 <- complement(d_parents$mother_tr1, set_rho)
    d_parents$father_tr1 <- rnorm(nrow(d_parents), 0, 1)
    d_parents$father_tr2 <- complement(d_parents$father_tr1, set_rho)
    d_parents$nest_tr1 <- rnorm(nrow(d_parents), 0, 1)
    d_parents$nest_tr2 <- complement(d_parents$nest_tr1, set_rho_nest)
  
  # Create data frame with each nestling in a row
    d_nestlings <- d_parents[rep(seq_len(nrow(d_parents)), d_parents$clutch), ]
    
  # Randomly assign xx% of nestlings as extra pair young. This isn't actually exactly
    # resulting in xx% epy because the randomization might pick the original father
    # in some cases. Its close enough though.
      for(i in 1:nrow(d_nestlings)){
        choose <- runif(1, 0, 100)
        if(choose > epp_rate){
          #d_nestlings$father[i] <- d_nestlings$father[round(runif(1, 0, nrow(d_nestlings)), 0)]
          d_nestlings[i, c("father", "father_tr1", "father_tr2")] <-
            d_nestlings[round(runif(1, 1, nrow(d_nestlings)), 0), c("father", "father_tr1", "father_tr2")]
        }
      }
 
  # Now randomly assign xx% to have been cross fostered and raised in a different nest.
    # Again this isn't exact, but should be close enough.
      for(i in 1:nrow(d_nestlings)){
        choose <- runif(1, 0, 100)
        if(choose > (100 - swap_rate)){
          #d_nestlings$nest[i] <- d_nestlings$nest[round(runif(1, 0, nrow(d_nestlings)), 0)]
          d_nestlings[i, c("nest", "nest_tr1", "nest_tr2")] <-
            d_nestlings[round(runif(1, 1, nrow(d_nestlings)), 0), c("nest", "nest_tr1", "nest_tr2")]
        }
      }
    
  # Calculate midparent values
      d_nestlings$mid_tr1 <- (d_nestlings$mother_tr1 + d_nestlings$father_tr1) / 2
      d_nestlings$mid_tr2 <- (d_nestlings$mother_tr2 + d_nestlings$father_tr2) / 2
    
  # Now create two traits for nestlings. These are randomized, but are based on the heritability
    # and environment effect values defined above. 
      #d_nestlings$trait1 <- rnorm(nrow(d_nestlings), 0, 1)   # early version, no heritability
      #d_nestlings$trait2 <- rnorm(nrow(d_nestlings), 0, 1)   # early version, no heritability
      d_nestlings$trait1 <- (set_h2_1 * d_nestlings$mid_tr1) +
                            (set_Ve_1 * d_nestlings$nest_tr1) +
                            ((1 - set_h2_1 - set_Ve_1) * rnorm(nrow(d_nestlings), 0, 1))
      d_nestlings$trait2 <- (set_h2_2 * d_nestlings$mid_tr2) +
                            (set_Ve_2 * d_nestlings$nest_tr2) +
                            ((1 - set_h2_2 - set_Ve_2) * rnorm(nrow(d_nestlings), 0, 1))
      
  # Give unique labels to each nestling
      d_nestlings$animal <- as.factor(paste0(seq(1, nrow(d_nestlings), 1), "_n"))
      
## Create a pedigree from the simulated data
    ped_sim <- d_nestlings[, c("animal", "mother", "father")]
    colnames(ped_sim) <- c("animal", "dam", "sire")
    ped_sim <- orderPed(ped_sim)
    ped_sim <- insertPed(ped_sim)
    
    cor(d_nestlings$mid_tr1, d_nestlings$trait1)
    cor(d_nestlings$mid_tr2, d_nestlings$trait2)
      
## Using the simulated data set created above, fit univariate heritability models for each trait.
    # Note that I've tried this with each of the priors listed below and with some 
    # changes to the strength of these priors. Modify the model below to choose which prior.
      
    # The prior in the Wilson et al tutorial
      uni_prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002),
                                  G2 = list(V = 1, nu = 0.002)), 
                                  R = list(V = 1, nu = 0.002))
      
    # the bivariate prior from Wilson et al tutorial (what Jocelyn used)
      pr_bi1 = list(G = list(G1 = list(V = diag(2), n = 1.002), G2 = list(V = diag(2), n = 1.002)),
                    R = list(V = diag(2), nu = 1.002))
    
    # A prior using parameter expansion
      uni_prior2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                  G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                                  R = list(V = 1, nu = 0.002))
      
    # A bivariate prior with parameter expansion
      pr_bi2 = list(G = list(G1 = list(V = diag(2), nu = 1.002, alpha.mu = c(0, 0), alpha.V = 1000*diag(2)), 
                             G2 = list(V = diag(2), nu = 1.002, alpha.mu = c(0, 0), alpha.V = 1000*diag(2))),
                    R = list(V = diag(2), nu = 1.002))
## Fit models      
    # The univariate model for trait 1 with parameter expansion
      m_un_tr1 <- MCMCglmm(trait1 ~ 1, random = ~animal + nest,
                           pedigree = ped_sim, prior = uni_prior2, verbose = FALSE,
                           data = d_nestlings, nitt = nitt_un_sim, thin = thin_un_sim, burnin = burn_un_sim)
    # The univariate model for trait 1 with original prior
      m_un_tr1o <- MCMCglmm(trait1 ~ 1, random = ~animal + nest,
                           pedigree = ped_sim, prior = uni_prior1,
                           data = d_nestlings, nitt = nitt_un_sim, thin = thin_un_sim, burnin = burn_un_sim)
      
      # Check on chains and model diagnostics. The mixing would surely be better with a longer
        # run and more thinning, but this seems OK for illustration.
          #plot(m_un_tr1$VCV)
          autocorr(m_un_tr1$VCV)
          effectiveSize(m_un_tr1$VCV)
          
      # Check posterior means and estimate heritability
          posterior.mode(m_un_tr1$VCV)
          tr1_herit <- m_un_tr1$VCV[, "animal"] / (m_un_tr1$VCV[, "animal"] + m_un_tr1$VCV[, "nest"] +
                                        m_un_tr1$VCV[, "units"])
          mean(tr1_herit)
          HPDinterval(tr1_herit, 0.95)
          
          posterior.mode(m_un_tr1o$VCV)
          tr1_herito <- m_un_tr1o$VCV[, "animal"] / (m_un_tr1o$VCV[, "animal"] + m_un_tr1o$VCV[, "nest"] +
                                                     m_un_tr1o$VCV[, "units"])
          mean(tr1_herito)
          HPDinterval(tr1_herito, 0.95)
          
    # The univariate model for trait 2 parameter expanded
          m_un_tr2 <- MCMCglmm(trait2 ~ 1, random = ~animal + nest,
                               pedigree = ped_sim, prior = uni_prior2,
                               data = d_nestlings, nitt = nitt_un_sim, thin = thin_un_sim, burnin = burn_un_sim)
    # The univariate model for trait 2 original prior
          m_un_tr2o <- MCMCglmm(trait2 ~ 1, random = ~animal + nest,
                               pedigree = ped_sim, prior = uni_prior1,
                               data = d_nestlings, nitt = nitt_un_sim, thin = thin_un_sim, burnin = burn_un_sim)
          
      # Check on chains and model diagnostics. The mixing would surely be better with a longer
        # run and more thinning, but this seems OK for illustration.
          #plot(m_un_tr2$VCV)
          autocorr(m_un_tr2$VCV)
          effectiveSize(m_un_tr2$VCV)
          
      # Check posterior means and estimate heritability
          posterior.mode(m_un_tr2$VCV)
          tr2_herit <- m_un_tr2$VCV[, "animal"] / (m_un_tr2$VCV[, "animal"] + m_un_tr2$VCV[, "nest"] +
                                                     m_un_tr2$VCV[, "units"])
          mean(tr2_herit)
          HPDinterval(tr2_herit, 0.95)
          
          posterior.mode(m_un_tr2o$VCV)
          tr2_herito <- m_un_tr2o$VCV[, "animal"] / (m_un_tr2o$VCV[, "animal"] + m_un_tr2o$VCV[, "nest"] +
                                                     m_un_tr2o$VCV[, "units"])
          mean(tr2_herito)
          HPDinterval(tr2_herito, 0.95)
        
## The heritability for both of these traits in univariate models is very close to 0, as it should
    # be. Now fit a bivariate model with the two traits together.
          
      # Fit the bivariate model with Jocelyn's prior
          m_bivo <- MCMCglmm(cbind(trait1, trait2) ~ trait - 1,
                             random = ~us(trait):animal + us(trait):nest, 
                             rcov = ~us(trait):units,
                             family = c("gaussian", "gaussian"), pedigree = ped_sim,
                             data = d_nestlings, prior = pr_bi1, 
                             nitt = nitt_bi_sim, thin = thin_bi_sim, burnin = burn_bi_sim)
      # Fit the bivariate model with parameter expansion
          m_biv <- MCMCglmm(cbind(trait1, trait2) ~ trait - 1,
                            random = ~us(trait):animal + us(trait):nest, 
                            rcov = ~us(trait):units,
                            family = c("gaussian", "gaussian"), pedigree = ped_sim,
                            data = d_nestlings, prior = pr_bi2, 
                            nitt = nitt_bi_sim, thin = thin_bi_sim, burnin = burn_bi_sim)
    
      # Look at chains and diagnostics
          #plot(m_biv$VCV)
          autocorr(m_biv$VCV)
          
      # calculate heritabilities from posterior means
          bi_tr1 <- m_biv$VCV[, "traittrait1:traittrait1.animal"] /
              (m_biv$VCV[, "traittrait1:traittrait1.animal"] + 
               m_biv$VCV[, "traittrait1:traittrait1.nest"] +
               m_biv$VCV[, "traittrait1:traittrait1.units"])
          mean(bi_tr1)
          HPDinterval(bi_tr1)
          
          bi_tr2 <- m_biv$VCV[, "traittrait2:traittrait2.animal"] /
              (m_biv$VCV[, "traittrait2:traittrait2.animal"] + 
               m_biv$VCV[, "traittrait2:traittrait2.nest"] +
               m_biv$VCV[, "traittrait2:traittrait2.units"])
          mean(bi_tr2)
          HPDinterval(bi_tr2)
          
          bi_cov <- m_biv$VCV[, "traittrait1:traittrait2.animal"] /
              (m_biv$VCV[, "traittrait1:traittrait2.animal"] +
               m_biv$VCV[, "traittrait1:traittrait2.nest"] +
               m_biv$VCV[, "traittrait1:traittrait2.units"])
          mean(bi_cov)
          HPDinterval(bi_cov)
          
          
          bi_tr1o <- m_bivo$VCV[, "traittrait1:traittrait1.animal"] /
            (m_bivo$VCV[, "traittrait1:traittrait1.animal"] + 
               m_bivo$VCV[, "traittrait1:traittrait1.nest"] +
               m_bivo$VCV[, "traittrait1:traittrait1.units"])
          mean(bi_tr1o)
          HPDinterval(bi_tr1o)
          
          bi_tr2o <- m_bivo$VCV[, "traittrait2:traittrait2.animal"] /
            (m_bivo$VCV[, "traittrait2:traittrait2.animal"] + 
               m_bivo$VCV[, "traittrait2:traittrait2.nest"] +
               m_bivo$VCV[, "traittrait2:traittrait2.units"])
          mean(bi_tr2o)
          HPDinterval(bi_tr2o)
          
          bi_covo <- m_bivo$VCV[, "traittrait1:traittrait2.animal"] /
            (m_bivo$VCV[, "traittrait1:traittrait2.animal"] +
               m_bivo$VCV[, "traittrait1:traittrait2.nest"] +
               m_bivo$VCV[, "traittrait1:traittrait2.units"])
          mean(bi_covo)
          HPDinterval(bi_covo)
          
      
      
      
      
      
      
      

    