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
    set.seed(135)
  
  # Simulate parents
    d_parents <- as.data.frame(as.factor(paste0(seq(1, 150), "_nest")))
    colnames(d_parents) <- "nest"
    d_parents$mother <- as.factor(paste0(seq(1, nrow(d_parents), 1), "_mom"))
    d_parents$father <- as.factor(paste0(seq(1, nrow(d_parents), 1), "_dad"))
    d_parents$clutch <- round(runif(nrow(d_parents), min = 3, max = 6))
  
  # Create data frame with each nestling in a row
    d_nestlings <- d_parents[rep(seq_len(nrow(d_parents)), d_parents$clutch), ]
    
  # Randomly assign 50% of nestlings as extra pair young. This isn't actually exactly
    # resulting in 50% epy because the randomization might pick the original father
    # in some cases. Its close enough though.
      for(i in 1:nrow(d_nestlings)){
        choose <- runif(1, 0, 100)
        if(choose > 50){
          d_nestlings$father[i] <- d_nestlings$father[round(runif(1, 0, nrow(d_nestlings)), 0)]
        }
      }
 
  # Now randomly assign 30% to have been cross fostered and raised in a different nest.
    # Again this isn't exact, but should be close enough.
      for(i in 1:nrow(d_nestlings)){
        choose <- runif(1, 0, 100)
        if(choose > 70){
          d_nestlings$nest[i] <- d_nestlings$nest[round(runif(1, 0, nrow(d_nestlings)), 0)]
        }
      }
    
  # Now create two traits with mean 0 and SD 1. These are entirely random, so they should 
    # not be heritable at all and should have no correlation with each other.
      d_nestlings$trait1 <- rnorm(nrow(d_nestlings), 0, 1)
      d_nestlings$trait2 <- rnorm(nrow(d_nestlings), 0, 1)
      
  # Give unique labels to each nestling
      d_nestlings$animal <- as.factor(paste0(seq(1, nrow(d_nestlings), 1), "_n"))
      
## Create a pedigree from the simulated data
    ped_sim <- d_nestlings[, c("animal", "mother", "father")]
    colnames(ped_sim) <- c("animal", "dam", "sire")
    ped_sim <- orderPed(ped_sim)
    ped_sim <- insertPed(ped_sim)
      
## Using the simulated data set created above, fit univariate heritability models for each trait.
    # Note that I've tried this with each of the priors listed below and with some 
    # changes to the strength of these priors. Modify the model below to choose which prior.
      
    # The prior in the Wilson et al tutorial
      uni_prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002),
                                  G2 = list(V = 1, nu = 0.002)), 
                                  R = list(V = 1, nu = 0.002))
    
    # A prior using parameter expansion
      uni_prior2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                                  G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                                  R = list(V = 1, nu = 0.002))
      
    # The univariate model for trait 1
      m_un_tr1 <- MCMCglmm(trait1 ~ 1, random = ~animal + nest,
                           pedigree = ped_sim, prior = uni_prior2,
                           data = d_nestlings, nitt = 25000, thin = 100, burnin = 5000,
                           verbose = FALSE)
      
      # Check on chains and model diagnostics. The mixing would surely be better with a longer
        # run and more thinning, but this seems OK for illustration.
          plot(m_un_tr1$VCV)
          autocorr(m_un_tr1$VCV)
          effectiveSize(m_un_tr1$VCV)
          
      # Check posterior means and estimate heritability
          posterior.mode(m_un_tr1$VCV)
          tr1_herit <- m_un_tr1$VCV[, "animal"] / (m_un_tr1$VCV[, "animal"] + m_un_tr1$VCV[, "nest"] +
                                        m_un_tr1$VCV[, "units"])
          mean(tr1_herit)
          HPDinterval(tr1_herit, 0.95)
          
    # The univariate model for trait 2
          m_un_tr2 <- MCMCglmm(trait2 ~ 1, random = ~animal + nest,
                               pedigree = ped_sim, prior = uni_prior2,
                               data = d_nestlings, nitt = 25000, thin = 100, burnin = 5000,
                               verbose = FALSE)
          
      # Check on chains and model diagnostics. The mixing would surely be better with a longer
        # run and more thinning, but this seems OK for illustration.
          plot(m_un_tr2$VCV)
          autocorr(m_un_tr2$VCV)
          effectiveSize(m_un_tr2$VCV)
          
      # Check posterior means and estimate heritability
          posterior.mode(m_un_tr2$VCV)
          tr2_herit <- m_un_tr2$VCV[, "animal"] / (m_un_tr2$VCV[, "animal"] + m_un_tr2$VCV[, "nest"] +
                                                     m_un_tr2$VCV[, "units"])
          mean(tr2_herit)
          HPDinterval(tr2_herit, 0.95)
        
## The heritability for both of these traits in univariate models is very close to 0, as it should
    # be. Now fit a bivariate model with the two traits together.
          
      ## Specify priors.
        # The original prior from the Wilson et al tutorial
          bi_prior1 = list(G = list(G1 = list(V = diag(2), nu = 1.002), 
                                    G2 = list(V = diag(2), nu = 1.002)),
                                    R = list(V = diag(2), nu = 1.002))
          
        # Another prior I found in a recent paper
          biv_mc_var = matrix(c(var(d_nestlings$trait1), 0, 0, var(d_nestlings$trait2)), 2, 2)
          bi_prior2 = list(G = list(G1 = list(V = biv_mc_var/2, n = 2),
                                        G2 = list(V = biv_mc_var/2, n = 2)),
                                        R = list(V = biv_mc_var/2, n = 2))
          
      # Fit the bivariate model
          m_biv <- MCMCglmm(cbind(trait1, trait2) ~ trait - 1,
                             random = ~us(trait):animal + us(trait):nest, 
                             rcov = ~us(trait):units,
                             family = c("gaussian", "gaussian"), pedigree = ped_sim,
                             data = d_nestlings, prior = bi_prior1, nitt = 50000, thin = 100, burnin = 5000)
    
      # Look at chains and diagnostics
          plot(m_biv$VCV)
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
          
      
      
      
      
      
      
      

    