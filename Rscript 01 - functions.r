################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# PARAMETER FUNCTION
# functions used in estimation of parameters for warmouth model
#-------------------------------------------------------------------------------
# Growth functions

# Function to est growth based on VB relationship
# can use form with t0 or form with L0
# use of either form based on specification of t0 or L0
# can only specify 1 and other must be NA
VB.growth <- function(t, Linf, K, t0 = NA, L0 = NA, names = FALSE) {
     
    if(is.na(L0[1])){
        L <- Linf * (1 - exp(-K * (t - t0)))  # t0 form
    } else if (is.na(t0[1])) {
        L <- Linf - (Linf - L0) * exp(-K * t) # L0 form
    }
    if(names == TRUE) names(L) <- paste0("L", t)
    L
}

# function to est age from VB growth function
# can use form with t0 or form with L0
# use of either form based on specification of t0 or L0
# can only specify 1 and other must be NA
agebyVB <- function(L, Linf, K, t0 = NA, L0 = NA){
    
    if(is.na(L0[1])){
        t <- -1 / K * log( 1 - L / Linf) + t0    # t0 form
    } else if (is.na(t0[1])) {
        t <- log(-((L - Linf) / (Linf - L0)))/-K # L0 form
    }
    t
}

#-------------------------------------------------------------------------------
# Estimate population growth rate and predict confidence/prediction interval
# Based on data and fit from Randall et al (1995) and adjustemnt to relationship
# from Randall and Minns (2000)

r_max_pred <- function(w,         # prediction weight
                       h,         # prediction habitat (1 = lake; 2 = river)
                       interval){ # confidence or prediction interval
    
    # Data from Randall et al. (1995)
    data <- as.data.frame(list(
        Production = c(30.3, 3.8, 2.0, 228.5, 15.5, 4.7, 21.5, 30.8, 5.7, 85.7, 8.6,
                       270.0, 596.0, 198.0, 140.0, 100.0,  74.0,  52.0, 30.0, 174.5,
                       280.3, 243.5, 226.4, 375.4, 379.1, 217.6, 500.6, 297.2, 110.3,
                       128.8, 161.4, 83.4, 807.2, 228.0, 329.0, 168.7, 57.2, 147.3,
                       54.8, 76.5, 30.7, 35.3, 261.5, 98.3, 39.3, 33.3, 47.5, 26.1,
                       31.2, 53.0, 456.2, 80.3, 29.4),
        biomass = c(73.0, 7.3, 9.2, 372.2, 32.2, 11.2, 28.9, 11.6, 24.4, 277.8, 29.5,
                    161.2, 198.0,  95.0,  75.0,  59.0,  51.0, 129.0, 62.0, 124.4,
                    274.8, 228.4, 149.9, 376.0, 278.7, 127.5, 217.5, 187.2, 111.2,
                    241.0, 307.5, 50.2, 310.5, 142.5, 86.6, 45.6, 10.8, 40.9, 49.8,
                    42.5, 38.4, 37.4, 173.0, 71.0, 30.5, 21.3, 38.5, 21.5, 27.1, 36.0,
                    233.4, 80.2, 15.4),
        weight = c(18.0246914, 3.5540409, 25.8426966, 6.6786291, 24.2105263, 37.8378378,
                   17.8947368, 5.7058534, 219.8198198, 13.6712598, 64.2701525, 1.24,
                   0.2544987, 1.4843750, 3.0, 2.4583333, 5.6666667, 18.4285714, 15.5,
                   2.4777915, 2.3869707, 1.4815487, 0.7197976, 3.7586469, 1.3582997,
                   2.0675564, 1.6603053, 0.6602313, 2.8333376, 13.9872316, 6.1428743,
                   2.1183222, 9.3047648, 4.2462529, 4.4194948, 0.4892704, 1.6410880,
                   0.8880686, 7.9935795,   0.9173916, 5.9222702, 10.0700054, 0.8397202,
                   1.3179144, 1.9223497, 2.1759117, 3.3804548, 2.6022755, 3.7215051,
                   2.7067669, 5.1344098, 12.0819524, 0.6273679),
        habitat = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1)
    ))
    
    # fit regression model
    mod = lm(log10(Production/biomass) ~ log10(weight) + habitat, data = data)
    
    # x prediction values
    pred.data = data.frame(weight = w, habitat = h)
    
    # make prediction with desired interval estiamte
    rmax = (predict(mod, pred.data, interval = interval))
    
    # retransform and multiply by 2 (Randall and Minns 2000)
    10^rmax*2
}


#-------------------------------------------------------------------------------
# OTHER

# function to calc geometric mean between sequenital values in a vector
geom_mean_seq <- function(vec) {
    n = length(vec)
    sqrt(vec[2:n]*vec[1:(n-1)])
}

# geometric mean of a vector
geom_mean <- function(vec, na.rm = FALSE) {
    exp(mean(log(vec), na.rm = na.rm))
}

#-------------------------------------------------------------------------------
# For parameterizing projection matrix from expression 

pmx_eval = function(mx, vars){
  matrix(sapply(mx, eval, vars),  sqrt(length(mx)), sqrt(length(mx)), byrow = T)
}

#===============================================================================
# Distribution functions

# Function to estimate shape parameters of the beta distribution
# from the defined mean and variance of data
beta_val <- function(mean, var){
  
  a <- mean * ((mean * (1 - mean)) / var - 1)
  b <- (1 - mean) * ((mean * (1 - mean)) / var - 1)
  
  list(a = a, b = b)
}

# Function to estimate shape parameters of the streched beta distribution
# from the defined mean and variance of data
# Streched beta is a beta dist'n scaled to outside (0,1) limits
# shape parameters give regulare beat dist'n which must be scaled with
# (max - min) + min when generating values with rbeta
beta_stretch_val <- function(mean, var = NA, sd = NA){
  
  if(is.na(var[1])) var = sd^2    # calc variance if not given
  if(is.na(sd[1])) sd = sqrt(var) # calc sd if not given
  
  # estimate mins and max of distribution
  min = qnorm(1e-6, mean, sd = sd)
  max = qnorm(1 - 1e-6, mean, sd = sd)
  
  # Scale means and variance to between 0 and 1
  mean.b <- (mean - min) / (max - min)
  var.b <- var * (1 / (max - min))^2
  
  # estimate beta dist'n shape parameters
  a <- mean.b * ((mean.b * (1 - mean.b)) / var.b - 1)
  b <- (1 - mean.b) * ((mean.b * (1 - mean.b)) / var.b - 1)
  
  # output
  list(a = a, b = b, min = min, max = max)
}

#===============================================================================
# Initial population size vector
# est. pop vector from desire adult population size

init_pop <- function(mx,   # Population Matrix 
                     Nt,   # desire adult density
                     p.rep # maturity vector - 1: mature; 0: Juvenle; length Tmax + 1
) 
{
  
  SS <- stable.stage(mx)                           # Stable stage distn
  Na <- SS * p.rep / sum(SS * p.rep) * Nt          # Calc adult pop
  Nj <- SS * (1-p.rep)/ sum(SS * p.rep) * Nt       # Calc juv pop
  N <- Nj + Na   
  N
}

#===============================================================================

# Function to project population size including density dependence effects
Projection <- function(mx,          # projection matrix expression
                       data,        # life history data
                       lh_mean,     # mean life history data
                       N0,          # initial age1+ population size pop size 
                       years,       # years to run simulation
                       p.cat,       # Probability of castastrophe
                       density_dependence = FALSE, # Boolean; if density-depedence included
                       demographic_stoch = FALSE,  # Boolean: if dempgraphic stochasticity included
                       allee = FALSE               # Boolean: if allee effects included
){
  #-----------------------------------------------------------------------------
  
  # needs popbio library
  require(popbio)
  
  # Castastrophes
  Catastrophe = sample(c(1,0), years, replace = TRUE,
                       prob = c(p.cat / data$gen.time, 1 - p.cat  / data$gen.time))
  
  # effect on castastrophe on pop size (percent recuction) - scaled Beta dist'n fit from Reed et al. (2003)
  e.cat <- sapply(Catastrophe, function(x) {
    ifelse(x == 0, NA, rbeta(1, shape1 = 0.762, shape2 = 1.5) * (1 - .5) + .5)
  })
  
  # Initialize parameters
  lh <- lh_mean 
  
  fs <- f_rand(with(lh, c(f1, f2, f3, f4)), sd = data$f.logsd)
  lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
  
  ss <- s_rand(with(lh, c(s0, s1, s2, s3)), cv = data$M.cv)
  lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
  
  # Population matrix
  A <- pmx_eval(mx, lh)
  
  # Initial population structure
  N <- as.vector(stable.stage(A) * N0 ) 
  #N <- as.vector(A %*% init_pop(A, Na, data$p.rep) ) 
  
  # initialize output vectors 
  Nvec <- rep(NA, years + 1); Nvec[1] <- sum(N) # * data$p.rep)        # Adult pop vector
  Evec <- rep(NA, years + 1);
  s0s <- NULL       # density dependent parameters
  Ns <- list(N)              # age-specific annual population size 
  lambdas <- rep(NA, years)  # population growth rate
  
  # stochastic parameters
  fts <- list(fs)
  sts <- list(ss$ss)
  
  #-----------------------------------------------------------------------------
  # loop through years
  for(t in 1:years){
    
    #---------------------------------------------------------------------------
    # Density Dependence
    
    # Allee Effects
    if(allee) {
    
      lh$a.f <- sum(N*data$p.rep) ^ 2 / (data$a ^ 2 + sum(N*data$p.rep) ^ 2)
      if(is.finite(lh$a.f) == F) lh$a.f = 0
      
    }
    
    # Survival Density-Dependence
    Et <- E_est(lh, N) # Egg count
    Evec[t] <- Et      # save egg count in vector
    if(density_dependence) {
      
      #s0.d <- data$s0.max * exp(-data$b * Et)
      s0.d <- data$s0.max / (1 + data$b * Et) 
      if(s0.d < 1e-300) s0.d <- 1e-300
      
      # Mean survival rate
      lh_mean$s0 <- s0.d
    } 
    
    # Demographic Stochasticity
    # simlulates random variation in vital rates when pop size is small.
    # make N random drawns for each age class for vital rates. as N increase the
    # draw will be closer to the mean, samller Ns will be more stochastic
    if(demographic_stoch) {
      # Survival rate
      # Binomial distribution - Sum N drawn and divide by N to give mean Survival
      lh$s0 <- lh_mean$s0
      lh$s1 <- lh$s2 <- lh$s3 <- ifelse(sum(N) > 1 & sum(N) < 500, 
                                        sum(rbinom(sum(N), 1, lh_mean$s1))/sum(N), 
                                        lh_mean$s1)
      
      # Fecundity
      # Poison dist. - take mean of N/2 (females) draws to give mean fecundity
      lh$f1 <- lh$f2 <- lh$f3 <- lh$f4 <- ifelse(sum(N*data$p.rep) > 2 & sum(N*data$p.rep) < 500, 
                                                 mean(rpois(sum(N*data$p.rep)/2, lh_mean$f2)), 
                                                 lh_mean$f2) 
    }else { 
      lh <- lh_mean 
    }
    
    #---------------------------------------------------------------------------
    # Stochasticity
    
    # Fecundity
    fs <- f_rand(with(lh, c(f1, f2, f3, f4)), sd = data$f.logsd)
    lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
    
    # Survival
    ss <- s_rand(with(lh, c(s0, s1, s2, s3)), cv = data$M.cv)
    lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
    
    #---------------------------------------------------------------------------
    # Population
    
    # Projection matrix
    A <- pmx_eval(mx, lh) 
    
    # project the population 1 year ahead.
    if(Catastrophe[t] == 1){
      N <- N * (1 - e.cat[t])
    } else {
      N <- as.vector(A %*% N )  
    }
    
    Nvec[t + 1] <- sum(N)         # Number of mature fish in pop
    s0s <- c(s0s, lh_mean$s0)
    
    Ns[t + 1] <- list(N)
    
    lambdas[t] <- Nvec[t + 1]/Nvec[t] # pop growth rate
    
    fts[t + 1] <- list(fs)
    sts[t + 1] <- list(ss$ss)
    
  }
  
  #-----------------------------------------------------------------------------
  # output
  list(pop = as.data.frame(list(year = 0:years,
                                N = Nvec,
                                E = Evec,
                                s0 = c(s0s, NA))),
       
       N = do.call(rbind, Ns),
       
       lambdas = lambdas,
       vars = list(
         ft = do.call(rbind, fts),
         st = do.call(rbind, sts)
       ),
       Cat = as.data.frame(list(
         Cat = Catastrophe,
         e.cat = e.cat
       ))
       
  )
}

#===============================================================================
# Repatriation

repatriation <- function(mx,                # projection matrix expression
                         data,              # life history data
                         lh_mean,           # mean life history data
                         Ninit = 0,         # initial adult population size pop size
                         Nt.vec,            # transfered population as vector of ages 0:tmax
                         s.t,               # survvial rate during translocation
                         t.time,            # when stocking takes place ="pre.spawn" OR "post.spawn" 
                         stock.month,       # num. months post spawned stocking takes place
                         stock.max,         # number of stocking events to take place
                         years,             # years to run simulation
                         p.cat              # Probability of castastrophe
){
  #-----------------------------------------------------------------------------
  
  # needs popbio library
  require(popbio)
  
  # Castastrophes
  Catastrophe = sample(c(1,0), years, replace = TRUE,
                       prob = c(p.cat / data$gen.time, 1 - p.cat  / data$gen.time))
  
  # effect on castastrophe on pop size (percent recuction) - scaled Beta dist'n fit from Reed et al. (2003)
  e.cat <- sapply(Catastrophe, function(x) {
    ifelse(x == 0, NA, rbeta(1, shape1 = 0.762, shape2 = 1.5) * (1 - .5) + .5)
  })
  
  # Initialize parameters
  lh <- lh_mean 
  
  fs <- f_rand(with(lh_mean , c(f1, f2, f3, f4)), sd = data$f.logsd)
  lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
  
  ss <- s_rand(with(lh_mean , c(s0, s1, s2, s3)), cv = data$M.cv)
  lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
  
  # Population matrix
  A <- pmx_eval(mx, lh)
  
  # Initial population structure
  N <- as.vector(stable.stage(A) * Ninit ) 
  
  # initialize output vectors 
  Nvec <- rep(NA, years);    #Adult pop vector
  Ns <- list(N)              # age-specific annual population size 
  
  # stock counter
  stock.count <- 0
  
  #-----------------------------------------------------------------------------
  # loop through years
  for(t in 1:years){
   
    # stocking pre spawning
    if(t.time == "pre.spawn" & stock.count < stock.max){
      N <- N + Nt.vec * s.t 
      stock.count <- stock.count + 1
    }
    
    #---------------------------------------------------------------------------
    # Density Dependence
    
    # Allee
    lh$a.f <- sum(N*data$p.rep) ^ 2 / (data$a ^ 2 + sum(N*data$p.rep) ^ 2)
    if(is.finite(lh$a.f) == F) lh$a.f = 0
    
    # SURVIVAL
    Et <- E_est(lh, N)
    
    #s0.d <- data$s0.max * exp(-data$b * Et)
    s0.d <- data$s0.max / (1 + data$b * Et) 
    if(s0.d < 1e-300) s0.d <- 1e-300
    
    # Mean survival rate
    lh$s0 <- s0.d

    #---------------------------------------------------------------------------
    # Demographic Stochasticity
    # simlulates random variation in vital rates when pop size is small.
    # make N random drawns for each age class for vital rates. as N increase the
    # draw will be closer to the mean, samller Ns will be more stochastic
    
    # Survival rate
    # Binomial distribution - Sum N drawn and divide by N to give mean Survival
    lh$s1 <- lh$s2 <- lh$s3 <- ifelse(sum(N) > 1 & sum(N) < 500, 
                                      sum(rbinom(sum(N), 1, lh_mean$s1))/sum(N), 
                                      lh_mean$s1)
    
    # Fecundity
    # Poison dist. - take mean of N/2 (females) draws to give mean fecundity
    lh$f1 <- lh$f2 <- lh$f3 <- lh$f4 <- ifelse(sum(N*data$p.rep) > 2 & sum(N*data$p.rep) < 500, 
                                               mean(rpois(sum(N*data$p.rep)/2, lh_mean$f2)), 
                                               lh_mean$f2) 
    
    #---------------------------------------------------------------------------
    # Stochasticity
    
    # Fecundity
    fs <- f_rand(with(lh, c(f1, f2, f3, f4)), sd = data$f.logsd)
    lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
    
    # Survival
    ss <- s_rand(means = with(lh, c(s0, s1, s2, s3)), cv = data$M.cv, X = NA)
    lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
    
    #---------------------------------------------------------------------------
    # Population
    
    # Projection matrix
    A <- pmx_eval(mx, lh) 
    
    # project the population 1 year ahead.
    if(Catastrophe[t] == 1){
      N <- N * (1 - e.cat[t])
    } else {
      N <- as.vector(A %*% N )  
    }
    
    # Stocking - post spawning
    # Stocked fish must survival the year before spawning - included stochasticity
    if(t.time == "post.spawn" & stock.count < stock.max){
      
      # Demographic stochasticity
      lh.s <- lh
      lh.s$s1 <- lh.s$s2 <- lh.s$s3 <- ifelse(sum(Nt.vec) > 1 & sum(Nt.vec) < 500, 
                                        sum(rbinom(sum(Nt.vec), 1, lh_mean$s1))/sum(Nt.vec), 
                                        lh_mean$s1)
      
      # Environmental stochasticity - same residuals as in matrix (X = ss$X)
      ss <- s_rand(with(lh.s, c(s0, s1, s2, s3)), cv = data$M.cv, X = ss$X)
      lh.s$s0 <- ss$ss[1]; lh.s$s1 <- ss$ss[2]; lh.s$s2 <- ss$ss[3]; lh.s$s3 <- ss$ss[4]; 
      
      # add stocked fish to pop vector - include mortality scaled to time remaiing in year
      N <- N + Nt.vec * s.t * with(lh.s, c(s1, s2, s3, 0)) ^ ((12 - stock.month)/12)
      
      # increment stocking counter
      stock.count <- stock.count + 1
    }
    
    # outputs
    Nvec[t] <- sum(N)         # Number of mature fish in pop

    Ns[t] <- list(N)          # Pop vector for output
    
    
  }
  #-----------------------------------------------------------------------------
  # output
  list(pop = as.data.frame(list(year = 1:years,
                                N = Nvec)),
       
       N = do.call(rbind, Ns)
  )
}

#===============================================================================
# Source Population Harm

# Function to project population size including density dependence effects
source_harm <- function(mx,          # projection matrix expression
                        data,        # life history data
                        lh_mean,     # mean life history data
                        Ninit ,      # initial adult population size pop size
                        Nt.vec,      # transfered population as vector of ages 0:tmax
                        t.time,      # when stocking takes place ="pre.spawn" OR "post.spawn" 
                        stock.max,   # number of stocking events to take place
                        years,       # years to run simulation
                        p.cat        # Probability of castastrophe "Length" and or "Survival"
){
  #-----------------------------------------------------------------------------
  
  # needs popbio library
  require(popbio)
  
  # Castastrophes
  Catastrophe = sample(c(1,0), years, replace = TRUE,
                       prob = c(p.cat / data$gen.time, 1 - p.cat  / data$gen.time))
  
  # effect on castastrophe on pop size (percent recuction) - scaled Beta dist'n fit from Reed et al. (2003)
  e.cat <- sapply(Catastrophe, function(x) {
    ifelse(x == 0, NA, rbeta(1, shape1 = 0.762, shape2 = 1.5) * (1 - .5) + .5)
  })
  
  # Initialize parameters
  lh <- lh_mean 
  
  fs <- f_rand(with(lh_mean , c(f1, f2, f3, f4)), sd = data$f.logsd)
  lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
  
  ss <- s_rand(with(lh_mean , c(s0, s1, s2, s3)), cv = data$M.cv)
  lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
  
  # Population matrix
  A <- pmx_eval(mx, lh)
  
  # Initial population structure
  N <- as.vector(stable.stage(A) * Ninit) 
  
  # initialize output vectors 
  Nvec <- rep(NA, years + 1); Nvec[1] <- sum(N)# * data$p.rep)        # Adult pop vector
  Ns <- list(N)              # age-specific annual population size 

  # stock counter
  stock.count <- 0
  
  #-----------------------------------------------------------------------------
  # loop through years
  for(t in 1:years){
    
    # stocking pre spawning
    if(t.time == "pre.spawn" & stock.count < stock.max){
      N <- N - Nt.vec 
      N <- ifelse(N > 1, N, 0)
      stock.count <- stock.count + 1
    }

    #---------------------------------------------------------------------------
    # Density Dependence
    
    # Allee Effects
    lh$a.f <- sum(N*data$p.rep) ^ 2 / (data$a ^ 2 + sum(N*data$p.rep) ^ 2)
    if(is.finite(lh$a.f) == F) lh$a.f = 0
    
    # Survival Density-Dependence
    Et <- E_est(lh, N) # Egg count
    s0.d <- data$s0.max / (1 + data$b * Et) 
    if(s0.d < 1e-300) s0.d <- 1e-300
      
    # Mean survival rate
    lh_mean$s0 <- s0.d
    
    # Demographic Stochasticity
    # simlulates random variation in vital rates when pop size is small.
    # make N random drawns for each age class for vital rates. as N increase the
    # draw will be closer to the mean, samller Ns will be more stochastic

    # Survival rate
    # Binomial distribution - Sum N drawn and divide by N to give mean Survival
    lh$s0 <- lh_mean$s0
    lh$s1 <- lh$s2 <- lh$s3 <- ifelse(sum(N) > 1 & sum(N) < 500, 
                                      sum(rbinom(sum(N), 1, lh_mean$s1))/sum(N), 
                                      lh_mean$s1)
    
    # Fecundity
    # Poison dist. - take mean of N/2 (females) draws to give mean fecundity
    lh$f1 <- lh$f2 <- lh$f3 <- lh$f4 <- ifelse(sum(N*data$p.rep) > 2 & sum(N*data$p.rep) < 500, 
                                               mean(rpois(sum(N*data$p.rep)/2, lh_mean$f2)), 
                                               lh_mean$f2) 
    
    #---------------------------------------------------------------------------
    # Stochasticity
    
    # Fecundity
    fs <- f_rand(with(lh_mean , c(f1, f2, f3, f4)), sd = data$f.logsd)
    lh$f1 <- fs[1]; lh$f2 <- fs[2]; lh$f3 <- fs[3]; lh$f4 <- fs[4]; 
    
    # Survival
    ss <- s_rand(with(lh_mean , c(s0, s1, s2, s3)), cv = data$M.cv)
    lh$s0 <- ss$ss[1]; lh$s1 <- ss$ss[2]; lh$s2 <- ss$ss[3]; lh$s3 <- ss$ss[4]; 
    
    #---------------------------------------------------------------------------
    # Population
    
    # Projection matrix
    A <- pmx_eval(mx, lh) 
    
    # project the population 1 year ahead.
    if(Catastrophe[t] == 1){
      N <- N * (1 - e.cat[t])
    } else {
      N <- as.vector(A %*% N )  
    }
    N <- ifelse(N > 1, N, 0)
    
    #output
    Nvec[t + 1] <- sum(N) # Number of age-1+ fish in pop
    Ns[t + 1] <- list(N)  # pop vector
    
    # Stocking - post spawning
    if(t.time == "post.spawn" & stock.count < stock.max){
      N <- N - Nt.vec 
      N <- ifelse(N > 1, N, 0)
      stock.count <- stock.count + 1
    }
    
  }
  
  #-----------------------------------------------------------------------------
  # output
  list(pop = as.data.frame(list(year = 0:years,
                                N = Nvec
                                )),
       
       N = do.call(rbind, Ns)
  )
}

#===============================================================================