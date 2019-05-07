################################################################################
################################################################################
# Eastern Sand Darter
# Repatriation model
#-------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())  ### clear the workspace

# Set working directory
setwd("")

# Load Libraries
library(pacman)
p_load(popbio)
p_load(ggplot2)
p_load(parallel)
p_load(RColorBrewer)

#-------------------------------------------------------------------------------
# Functions
# Files to source 
source("Rscript 01 - functions.r")  # sub routines for parameter fitting

#-------------------------------------------------------------------------------
# ggplot theme

theme_me <- theme_bw() +
      theme(axis.title = element_text(size = 11, family = "sans", face = "bold"),
            axis.text.x = element_text(size = 10, family = "sans", colour = "black"),
            axis.text.y = element_text(size = 10, family = "sans", hjust = 0.6,
                                       angle = 0, colour = "black"),
            legend.title = element_text(size = 10, family = "sans"),
            legend.text = element_text(size = 8, family = "sans"),
            strip.text = element_text(size = 11, family = "sans", face = "bold"))

# Colour blind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#-------------------------------------------------------------------------------
# Data
LWdata <- read.csv("ESD-Length.weight.csv", header = T)

#===============================================================================
# Parameters
# Initalize list to store parameter info
lh_data <- NULL   # life history data; parameters to est. growth, fecundity, and survival

#-------------------------------------------------------------------------------
# GROWTH
#-------------------------------------------------------------------------------
lh_data$H0 <- 5.6 # mm; Simon et al. 1992; COSEWIC 2009

# Von Bert growth function - Finch 2013
lh_data$Linf <- 55.52
lh_data$k <- 1.59
lh_data$t0 <- -0.47
#-------------------------------------------------------------------------------
# Length-Weight Relationship
#-------------------------------------------------------------------------------

mod <- nls(Weight ~ a * Length ^ b, data = LWdata, start = list(a = 1, b = 3))
L <- seq(20, 70, 0.25)
pred_data <- as.data.frame(list(
  Length = L,
  Weight = predict(mod, newdata = data.frame(Length = L))
))

# Plot
p <- ggplot(LWdata, aes(Length, Weight))
p <- p + geom_point(aes(), size = 1.5)
p <- p + geom_path(data = pred_data, aes(Length, Weight), size = 1)
p <- p + labs(x = "Total Length (mm)", y = "Weight (g)")           
p <- p + theme_me

#png("Figures/von Bert.png", width = 5, height = 3, unit = "in", res = 300)
p
#dev.off()

lh_data$a.LW <- coef(mod)[[1]]
lh_data$b.LW <- coef(mod)[[2]]

# clean up
rm(mod, L, pred_data, LWdata)

#-------------------------------------------------------------------------------
# Reproduction
#-------------------------------------------------------------------------------
# Fish of both sexes mature in the spring following their first growing season
# at age-1, but some females may not spawn until theri second year (COSEWIC 2009)
# Spawn in late spring and summer - water temp 14.4 - 25.5 C.
# Eggs laid several times during protracted spanwing season.
#-------------------------------------------------------------------------------

# Maturity ogive 
#lh_data$MO1 <- 0.911 # Finch et al. 2018

# Sex ratio 
#lh_data$R <- 0.5

#Spawning periodicity - 2 clutches per year - Finch 2018
#lh_data$P <- 2

# Fecdundity - mean clutch size - Finch 2018
lh_data$f <- 71.5

# Vector of proportion mature in each age class (0 to 4)
lh_data$p.rep = c(0, lh_data$MO1, lh_data$MO2, 1, 1)

#-------------------------------------------------------------------------------
# SURVIVAL
#-------------------------------------------------------------------------------
lh_data$S0 <- NA
lh_data$Sa <- 0.386 # Finch 2018

#===============================================================================
# Define the matrix structure 
#-------------------------------------------------------------------------------

# Projection matrix structure of a single population
# Matrix structure - Vital rates

# Tmat = 1; Tmax = 4
#A.esd.post <- expression(
#   MO1*P*R*f1*s0*a.f,  P*R*f2*s1*a.f,  P*R*f3*s2*a.f,  P*R*f4*s3*a.f,  0,
#   s0,                 0,              0,              0,              0,
#   0,                  s1,             0,              0,              0,
#   0,                  0,              s2,             0,              0,
#   0,                  0,              0,              s3,             0)
 
# PRE-SPAWN
# Tmat = 1; Tmax = 4
A.esd <- expression(
  MO1*P*R*f1*s0*a.f,  P*R*f2*s0*a.f,  P*R*f3*s0*a.f,  P*R*f4*s0*a.f,
  s1,                 0,              0,              0,          
  0,                  s2,             0,              0,      
  0,                  0,              s3,             0)

#===============================================================================
# Functions
# fucntion to parameterize matrix with different life history characteristics
# 8 matrices possible:
# Tmax - age 3 ro 4
# Tmat - age 1 or 2
# clutches - 2 or 3 per year
lh_func <- function(Tmax, Tmat, clutch, data) {
  # Tmax = 4, Tmat = 1, clutch = 2
  lh <- with(data, list(s0 = Syoy, s1 = Sa, s2 = Sa, s3 = Sa,
                        f1 = f, f2 = f, f3 = f, f4 = f, 
                        MO1= 0.911, P = 2, R = 0.5, a.f = 1))
  # if Tmax = 3
  if(Tmax == 3) {
    lh$s3 <- 0     # change survival rate at age-3 = 0%
    lh$f4 <- 0     # change fecundity at age-4 to 0
  } 
  # if tmat = 2
  if(Tmat == 2) {
    lh$MO1 <- 0     # change prop mature at age-1 tp 0
    lh$f1 <- 0      # change fecundity at age-1 to 0
  } 
  # if clutch = 3
  if(clutch == 3) {
    lh$P <- 3       # change number of annual spawns to 3
  } 
  lh
} 

# Egg count
E_est <- function(data, N) {
  with(data, MO1*P*R*f1*a.f * N[1] + P*R*f2*a.f * N[2] + P*R*f3*a.f * N[3] + P*R*f4*a.f * N[4])
}

#===============================================================================
# LIFE HISTORY PARAMETER UNCERTANTY
#===============================================================================

lh_params <- list(list(Tmax = 4, Tmat = 1, clutch = 3),
                  list(Tmax = 3, Tmat = 2, clutch = 2))
lh_names <- c(paste("Tmax -", 4,"Tmat -", 1,"clutch -", 3), 
              paste("Tmax -", 3,"Tmat -", 2,"clutch -", 2))
names(lh_params) <- lh_names

#===============================================================================
# MAXIMUM POPULATION GROWTH RATES
#===============================================================================

# Max lambda (Randal and Minns 2000) - rmax = 2.64*Wmat^-0.35 (take lower prediction interval)
L_init <- with(lh_data, VB.growth(0:4, Linf = Linf, K = k, t0 = t0, names = T)) # length at age
Wmat <- with(lh_data, a.LW * L_init ^ b.LW)[3]                                  # Weight at maturity (age-2)
lambda.max <- round(exp(r_max_pred(Wmat, 0, interval = "prediction"))[[2]],2)   # Max lambda est from randall and minns lower prediction interval

# clean up
rm(Wmat, L_init)

# Max lambda sequence - for diff. simulations
max.lambdas <- round(seq(1, lambda.max, length.out = 4)[-1], 2)
L.names <- c(paste("Lambda -", max.lambdas[1]), paste("Lambda -", max.lambdas[2]),
             paste("Lambda -", max.lambdas[3]))

#===============================================================================
# DETERMINISTIC OPTIMIXATIONS
#===============================================================================
# Syoy - age-0 survival rate
# At lambda = 1 and lambda = lambda max

# Optimization furntion - solve for syoy
s_optim = function (syoy, mx, data, Tmax, Tmat, clutch, target.lambda){
  
  data$Syoy <- syoy
  
  # Population matrix
  A <- pmx_eval(mx, lh_func(Tmax = Tmax, Tmat = Tmat, clutch = clutch, data = data))
  
  # lambda
  lambda <- lambda(A)
  
  # golden search algorithm
  (lambda - target.lambda)^2
}

#-------------------------------------------------------------------------------

# optimization - lambda = 1
lh_data$s0.1.det <- lapply(1:length(lh_params), function(x) {
  s_optim_min <- optimize(s_optim, c(0, 1), tol = 1e-16, mx = A.esd, 
                          data = lh_data, Tmax = lh_params[[x]]$Tmax, 
                          Tmat = lh_params[[x]]$Tmat, clutch = lh_params[[x]]$clutch,
                         target.lambda = 1)$minimum    # lambda = min
})
names(lh_data$s0.1.det)  <- lh_names

# optimization - lambda = max
lh_data$s0.max.det <- lapply(max.lambdas, function(lambda) {
  lhs <- lapply(1:length(lh_params), function(x) {
    s_optim_min <- optimize(s_optim, c(0, 1), tol = 1e-16, mx = A.esd, 
                            data = lh_data, Tmax = lh_params[[x]]$Tmax, 
                            Tmat = lh_params[[x]]$Tmat, clutch = lh_params[[x]]$clutch,
                            target.lambda = lambda)$minimum    # lambda = min
  })
  names(lhs)  <- lh_names
  lhs
})
names(lh_data$s0.max.det) <- L.names

#-------------------------------------------------------------------------------

# mean projection matrix - lambda = 1
pmx.1.det <- lapply(1:length(lh_params), function(x) {
  data <- lh_data
  data$Syoy <- lh_data$s0.1.det[[x]]
  pmx <- pmx_eval(A.esd, lh_func(Tmax = lh_params[[x]]$Tmax, 
                                 Tmat = lh_params[[x]]$Tmat,  
                                 clutch = lh_params[[x]]$clutch, data))
  
})
names(pmx.1.det) <- lh_names

# mean projection matrix - lambda = max lambdas
pmx.max.det <- lapply(1:length(max.lambdas), function(l) {
  lhs <- lapply(1:length(lh_params), function(x) {
    data <- lh_data
    data$Syoy <- lh_data$s0.max.det[[l]][[x]]
    pmx <- pmx_eval(A.esd, lh_func(Tmax = lh_params[[x]]$Tmax, 
                                   Tmat = lh_params[[x]]$Tmat,  
                                   clutch = lh_params[[x]]$clutch, data))
  })
  names(lhs)  <- lh_names
  lhs
  })
names(pmx.max.det) <- L.names

#-------------------------------------------------------------------------------
# generation time

lh_data$gen.time.lh <- median(sapply(1:length(lh_params), function(x) {
  pmx <- pmx.1.det[[x]]
  generation.time(pmx, c = c(1:4))
}))
#names(lh_data$gen.time.lh) <- lh_names

#-------------------------------------------------------------------------------
# Variability
#-------------------------------------------------------------------------------

# Fecundity variability
# Log standard deviation for log-normal distribution
lh_data$f.logsd <- 0.05

# Survival rate
# age-specific CV for Instantaneous mortality as a normal distribution
lh_data$M.cv <- c(0.1, 0.2, 0.2, 0.2)

# Function to est. stochatic fecundity 
f_rand <- function(means, sd) {
  
  X <- rnorm(1, mean = 0, sd = 1)  # Est standard normal residual
  pX <- pnorm(X, mean = 0, sd = 1) # Normal probability distribution
  
  # Assign residual to each age; with age-specigic distribution
  fs <- sapply(means, function(x){
    
    if(x == 0){                             # if 0 keep at 0
      0
    } else {
      qlnorm(pX, meanlog = log(x), sd = sd) # convert non-0s to log-normal distribution
    }
  })
  fs
}

# Function to est. stochatic survival
s_rand <- function(means, cv, X = NA) {
  if(is.finite(X)){
    X = X
  } else {
    X <- rnorm(1, mean = 0, sd = 1)  # Est standard normal residual
  }
  
  pX <- pnorm(X, mean = 0, sd = 1) # Normal probability distribution
  
  # Assign residual to each age; with age-specigic distribution
  ss <- sapply(1:length(means), function(x){
    
    if(means[x] == 0){ # if 0 keep at 0
      0
    } else if(means[x] >= 1){
      1
    } else {
      params = beta_stretch_val(mean = -log(means[x]), sd = -log(means[x]) * cv[x]) # identify stretched beta distn params
      M <- qbeta(pX, shape1 = params$a, shape2 = params$b) * (params$max - params$min) + params$min # convert non-0s to beta distribution for M
      exp(-M) # Convert M to survival rate
    }
  })
  list(ss = ss, X = X) # output survival rate
}

#===============================================================================
# STOCHASTIC OPTIMIZATION
#===============================================================================
# Optimization function

# solve for YOY survival rate to acheive a desires lambda under environmental stochasticity
optim_f <- function(v, target.lambda, type, reps, N0, data, 
                    Tmax, Tmat, clutch) {
  data$Syoy <- v
  lh_mean <- lh_func(Tmax = Tmax, Tmat = Tmat, clutch = clutch, data = data)
  data$p.rep = c(lh_mean$MO1, 1, 1, 1)
  seed <- 123
  set.seed(seed)
  pop_data <- lapply(1:reps, function(i) {
    res <- Projection(mx = A.esd,     # projection matrix expression
                      data = data,    # life history data
                      lh_mean = lh_mean,
                      N0 = N0,        # initial pop size as stage-structure vector
                      years = 150,    # years to run simulation
                      p.cat = 0       # Probability of castastrophe
    )
    res$pop[-(1:51), 1:2]
  })
  
  # est mean population growth rate
  if(type == "logmu") {
    pop_data <- do.call(rbind, pop_data)
    pop_data <- pop_data[pop_data$N > 0, ]
    pop_data <- pop_data[is.finite(pop_data$N), ]
    mod <- lm(log(N) ~ year, data = pop_data)
    lambda.est <- exp(coef(mod)[[2]])
  } else if(type == "loglambda") {
    log.lambdas <- sapply(pop_data, function(i) {
      i <- i[is.finite(i$N), ]
      (log(i$N[dim(i)[1]]) - log(i$N[1])) / dim(i)[1]
    })
    lambda.est <- exp(mean(log.lambdas))
  }
  (lambda.est - target.lambda)^2
}

#-------------------------------------------------------------------------------
# YOY SURVIVAL RATE
# s0 @ lambda = 1

if(FALSE) {
  # s0 @ lambda = 1 - in parallel over CVs, cor.rhos and maxlambdas
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, ls())          # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  s_optim.1 <- lapply(1:length(lh_params), function(x) {
    data <- lh_data
    data$gen.time <- lh_data$gen.time.lh
    s_optim_min <- optimize(optim_f, c(0, 1), tol = 1e-16,  target.lambda = 1, 
                            type = "loglambda", reps = 50, N0 = 10000,
                            data = data,  Tmax = lh_params[[x]]$Tmax, 
                            Tmat = lh_params[[x]]$Tmat, 
                            clutch = lh_params[[x]]$clutch)    
    s_optim_min
  })

  names(s_optim.1) <- lh_names
  stopCluster(cl) # close clusters
  s_optim.1
  
  # s0 @ lambda = max - in parallel over CVs, cor.rhos and maxlambdas
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, ls())          # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  s_optim.max <- parLapply(cl, max.lambdas, function(lambda) {
    lhs <- lapply(1:length(lh_params), function(x) {
      data <- lh_data
      data$gen.time <- lh_data$gen.time.lh[[x]]
      data$b <- 0
      s_optim_min <- optimize(optim_f, c(0, 1), tol = 1e-16,  target.lambda = lambda, 
                              type = "loglambda", reps = 5000, N0 = 10000,
                              data = data,  Tmax = lh_params[[x]]$Tmax, 
                              Tmat = lh_params[[x]]$Tmat, 
                              clutch = lh_params[[x]]$clutch)   
      s_optim_min
    })
    names(lhs) <- lh_names
    lhs
  })
  names(s_optim.max ) <- L.names
  stopCluster(cl) # close clusters
  s_optim.max
  
   # clean up
  rm(s_optim.1, s_optim.max)
}

# YOY Survival at lambda = max for each CV and rho
# maxLambda = 1.3, 1.6, 1.9
lh_data$s0.1.sto <- list(
 'Tmax - 4 Tmat - 1 clutch - 3' = 0.006092559,
 'Tmax - 3 Tmat - 2 clutch - 2' = 0.02593462
)

# YOY Survival at lambda = max for each CV and rho
# maxLambda = 1.3, 1.6, 1.9
lh_data$s0.max.sto <- list(
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.0117631,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.07130559),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.01766163,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.1414977),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.02349616,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.2340981)
)
names(lh_data$s0.max.sto) <- L.names  

#-------------------------------------------------------------------------------
# Density Dependence
#-------------------------------------------------------------------------------
# Use optimization to solve for the density dependence parameter
# such that the geometric mean population size is K

if(FALSE){
  b_optim_f <- function(b.d, N0, reps, data, lh_mean) {
    data$b <- b.d
    seed <- 123
    set.seed(seed)
    pop_data <- lapply(1:reps, function(i) {
      res <- Projection(mx = A.esd,     # projection matrix expression
                        data = data,    # life history data
                        lh_mean = lh_mean,
                        N0 = N0,        # initial pop size as stage-structure vector
                        years = 200,    # years to run simulation
                        p.cat = 0,      # Probability of castastrophe
                        density_dependence = TRUE
      )
      res$pop$N[-c(1:101)] # extract adult pop size removing 100 year burn in
    })
    
    # geommetric mean pop size
    N <- mean(sapply(pop_data, function(x) geom_mean(x)))
    
    # golden search algorithm - compare N to Na
    (N/N0 - 1)^2
  }
  
  #-----------------------------------------------------------------------------
  # Run Optimization
  
  N0 <- 10000 # target pop size - Carrying capacity
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, ls())          # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  
  b.dd.est <- parLapply(cl, 1:length(max.lambdas), function(l) { # loop over max almbdas
    lhs <- lapply(1:length(lh_params), function(x) { # loop over life history params  
      data <- lh_data                           # data
      data$gen.time <- lh_data$gen.time.lh      # assign generatioin time - based on LH data
      data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
      data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
      
      # assign LH mean values based on the lh parameters 
      lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                         Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                         clutch = lh_params[[x]]$clutch, # number of clutches
                         data = data)
      
      # Create proportion adult vector
      data$p.rep = c(lh_mean$MO1, 1, 1, 1)
      
      # run optimization
      b_optim <- optimize(b_optim_f, c(0, 1), tol = 1e-16,
                          N0 = N0, reps = 1000, data = data, 
                          lh_mean = lh_mean)
      b_optim$minimum * N0 # unscale b - K independent
    
    })
    names(lhs) <- lh_names
    lhs
  })
  names(b.dd.est) <- L.names
  stopCluster(cl) # close clusters
  b.dd.est
  
  # clean up
  rm(b.dd.est, b_optim_f)
}

# b - survival density dependence parameters per individual
# at lambda = max for each CV and rho
# maxLambda = 1.3, 1.6, 1.9
lh_data$b_dd <-list(
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.00882628,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.06895557),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.01865467,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.1789011),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 0.02842251,
       'Tmax - 3 Tmat - 2 clutch - 2' = 0.3243718)
)
names(lh_data$b_dd) <- L.names  

#-----------------------------------------------------------------------------
# Density-dependence fucntion plots

# Density dependence function
Et <- seq(100, 100000, 100)
DD_data <- as.data.frame(list(
  E = Et,
  s0 = lh_data$s0.max.sto[[3]][[1]] / (1 + lh_data$b_dd[[3]][[1]]/1000 * Et) 
))

p <- ggplot(DD_data, aes(x = E, y = s0)) 
p <- p + geom_path(aes())
p <- p + labs(x = "Egg Density", y = "YOY Survival")           
p <- p + theme_me 

#png("Figures/density-dependence function.png", width = 4, height = 3.5, unit = "in", res = 300)
p
#dev.off()

# Allee effect fucntion
Nt <- seq(0, 1000, 5)
AE_data <- as.data.frame(list(
  Na = rep(Nt, 2),
  af = c(Nt ^ 2 / (50 ^ 2 + Nt ^ 2), Nt ^ 2 / (100 ^ 2 + Nt ^ 2)),
  a = unlist(lapply(c(50, 100), rep, length(Nt)))
))

p <- ggplot(AE_data, aes(x = Na, y = af)) 
p <- p + geom_path(aes(linetype = as.factor(a)))
p <- p + labs(x = "Adult Density", y = "Allee Effect", linetype = "a")           
p <- p + theme_me 

#png("Figures/allee effect function.png", width = 5, height = 3.5, unit = "in", res = 300)
p
#dev.off()

# Clean up
rm(Et, DD_data, Nt, AE_data, p)

#===============================================================================
# DISTRIBUTION PLOTS
#===============================================================================
# plot variability of vital rates and lambda

if(FALSE) {
  reps <- 1000
  Na <- 10000
  years <- 100
  var_sim <- {
    x <- 1
    l <- 2 # max lambda 2
    data <- lh_data                           # data
    data$gen.time <- lh_data$gen.time.lh      # assign generatioin time - based on LH data
    data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
    data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
    
    # assign LH mean values based on the lh parameters 
    lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                       Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                       clutch = lh_params[[x]]$clutch, # number of clutches
                       data = data)
    
    # Create proportion adult vector
    data$p.rep = c(lh_mean$MO1, 1, 1, 1)
  
    # extract density dependenc paramter and scale to Na
    data$b <- data$b_dd[[l]][[x]]/Na
    
    out <- lapply(1:reps, function(i) {                # replicates
      
      res <- Projection(mx = A.esd,                    # projection matrix expression
                        data = data,                   # life history data
                        lh_mean = lh_mean,             # mean vital rate 
                        Na = Na,                       # initial pop size as stage-structure vector
                        years = years,                 # years to run simulation
                        p.cat = 0,                     # propability of catastrophe
                        density_dependence = c("Survival") # type of density dependence
      )
      
      list(N = res$pop$N,    # output adult pop size
           ft = res$vars$ft, # output fecundity
           st = res$vars$st) # output survival
    })
  }
  
  #-------------------------------------------------------------------------------
  # Lambda 
  
  lambdas <- {
      l <- lapply(1:reps, function(i){
        N <- var_sim[[i]]$N
        lambda1 = sapply(1:100, function(n) exp(log(N[n + 1]) - log(N[n])))                   # Annual lambda
        lambda10 = sapply(seq(1, 91, 10), function(n) exp((log(N[n + 10]) - log(N[n])) / 10)) # 10 years lambda 
        lambda100 = exp((log(N[101]) - log(N[1])) / 100)                                      # 100 year lambda
        as.data.frame(list(
          max_lambda = max.lambdas[2],
          LH = x,
          lambda = c(lambda1, lambda10, lambda100),
          time_frame = c(rep("Annual", length(lambda1)), rep("10 years", length(lambda10)), "100 years")
        ))
      })
      do.call(rbind, l)
  }

  lambdas$time_frame <- factor(lambdas$time_frame, levels = unique(lambdas$time_frame))
  
  p <- ggplot(lambdas, aes(x = lambda)) 
  p <- p + geom_density(aes(), fill = "grey", alpha = 0.5)
  p <- p + facet_wrap( ~ time_frame, scales = "free_y")
  p <- p + scale_x_continuous(limits = c(0.2, 3))
  p <- p + scale_fill_manual(values = cbPalette)
  p <- p + labs(x = ~lambda, y = "Density", fill = "Life History")           
  p <- p + theme_me + theme(axis.text.y = element_blank())
  
  png("Figures/Lambda dist - density-dependence.png", width = 7, height = 3, unit = "in", res = 300)
  p
  dev.off()
  
  # clean up
  rm(lambdas, p, no_cores, cl)
  
  #-------------------------------------------------------------------------------
  # Fecundity

  ft <- { # life histyory characteristics
    f <- lapply(1:reps, function(i){              # Reps
      as.data.frame(list(
        max_lambda = max.lambdas[2],
        LH = x,
        var = as.vector(var_sim[[i]]$ft[,2])
      ))
    })
    do.call(rbind, f)
  }
  
  p <- ggplot(ft, aes(x = var)) 
  p <- p + geom_density(fill = "grey", alpha = 0.5)
  p <- p + labs(x = "Length (mm)", y = "Density", fill = "Life History")           
  p <- p + theme_me +
    theme(axis.text.y = element_blank())
  
  #png("Figures/Fecundity.png", width = 3.5, height = 3, unit = "in", res = 300)
  p
  #dev.off()
  
  # clean up
  rm(ft, p)
  
  #-------------------------------------------------------------------------------
  # Survival
  
  st <- { # life histyory characteristics
    s <- lapply(1:reps, function(i){              # Reps
      as.data.frame(list(
        max_lambda = max.lambdas[2],
        LH = x,
        Age = sort(rep(0:1, 101)),
        var = as.vector(-log(var_sim[[i]]$st[,c(1,2)]))
      ))
    })
    do.call(rbind, s)
  }

  # assign names
  st$Age <- factor(st$Age, levels = c(0, 1), labels = c( "YOY", "Adult"))
  
  p <- ggplot(st, aes(x = var)) 
  p <- p + geom_density(aes(fill = as.factor(Age)), alpha = 0.5)
  p <- p + scale_x_continuous(limits = c(0, 10))
  p <- p + scale_fill_manual(values = cbPalette)
  p <- p + labs(x = "Instantaneous Mortality", y = "Density", fill = "Stage")           
  p <- p + theme_me +
    theme(axis.text.y = element_blank())
  png("Figures/Mortality dist.png", width = 3.5, height = 3, unit = "in", res = 300)
  p
  dev.off()
  
  # clean up
  rm(st, p)
  
  #-------------------------------------------------------------------------------
  # clean up
  rm(var_sim)
}

#===============================================================================
# MVP SIMULATIONS
#===============================================================================
# estimate minimum viabpe population size
# conduct PVA for ESD - sim population with density dependence over 100 years
# incorporate stochastic environmental effects and catastrophe (decline in total pop
# size of 50 to 100%)

years <- 100  # years to sim
reps <- 10000 # num. of replicates

# propability of catastrophe - per generation
p.cat <- 0.1

# Allee effect paraters
allee <- c(50, 100)

# carrying capacity - inital pop size
CC <- exp(seq(log(1000), log(100000), length.out = 20))

# Run simulation
if(FALSE) {
  
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, list("lh_data", "Projection", "reps", "years", "s_rand", 
                         "f_rand", "E_est", "lh_func", "lh_params", "max.lambdas",
                         "pmx_eval", 'A.esd', "init_pop", "beta_stretch_val", 
                         "lh_names", "CC", "p.cat", "allee")) # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  clusterSetRNGStream(cl, iseed = 123)
  
  MVPsim.sDD <- parLapply(cl, 1:length(max.lambdas), function(l) { # population growth rate
    lhs <- lapply(1:length(lh_params), function(x) {               # life history
      lapply(allee, function(a.f) { 
        lapply(CC, function(N0) {                                    # population carrying capacity
         
          data <- lh_data                           # data
          data$gen.time <- lh_data$gen.time.lh      # assign generatioin time - based on LH data
          
          data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
          data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
          
          # assign LH mean values based on the lh parameters 
          lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                             Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                             clutch = lh_params[[x]]$clutch, # number of clutches
                             data = data)
          
          data$p.rep = c(lh_mean$MO1, 1, 1, 1) # Create proportion adult vector 
          
          data$b <- data$b_dd[[l]][[x]]/N0     # extract density dependent param and scale to carrying capacity
          data$a <- a.f                        # Allee effect parameters - NA at 50% fecundity
          
          out <- lapply(1:reps, function(i) {            # replicates
            
            res <- Projection(mx = A.esd,                # projection matrix expression
                              data = data,               # life history data
                              lh_mean = lh_mean,         # mean vital rates for LH
                              N0 = N0,                   # initial pop size as stage-structure vector
                              years = years,             # years to run simulation
                              p.cat = p.cat,             # propability of catastrophe
                              density_dependence = TRUE, # Include density dependence
                              demographic_stoch = TRUE,  # Include dempgraphic stochasticity
                              allee = TRUE               # include Allee effects        
            )
            N = res$pop$N
          }) # reps
        }) # carrying capacity
      }) # Allee effects
    }) # Life history
    names(lhs) <- lh_names
    save(lhs, file = paste0("Results/Results 01 - MVPsim - lambda ", max.lambdas[l],".R")) # save file - one per pop growth rate
    NULL
  }) # pop growth rate
  stopCluster(cl)
  rm(MVPsim.sDD)
}

#-------------------------------------------------------------------------------
# MVP ESTIAMTEION
                         
# Load and extract data
if(FALSE){
  PVA_data <- do.call(rbind, lapply(max.lambdas, function(l){             # Loop through lambdas
    load(file = paste0("Results/Results 01 - MVPsim - lambda ", l, ".R")) # load R data
    do.call(rbind, lapply(1:length(lhs), function(x) {                    # Loop thorugh life history
      do.call(rbind, lapply(1:length(allee), function(a.f){               # loop through allee effects
        do.call(rbind, lapply(1:length(lhs[[x]][[a.f]]), function(cc){           # loop through carrying capacity
          as.data.frame(list(                    # extract data
            lambda = l,                          # LAMBDA
            LH = x,                              # life history
            allee = allee[a.f],                  # allee effect
            pop = CC[cc],                        # Carrying capacity
            rep = 1:reps,                        # replicate indicator
            P = sapply(lhs[[x]][[a.f]][[cc]], function(x) {as.numeric(min(x) < 2)}) # Determing If pop goes extinct - N < 2 adults
          ))
        })) # Carrying capacity
      })) # Allee effect
    })) # life history
  })) # lambda
  row.names(PVA_data) <- 1:nrow(PVA_data) # renames rows

  # Fit MVP model - Logistic regression of P[extinct] ~ log(carrying capacity)
  mod.PVA <-  lapply(max.lambdas, function(l) { # max lambdas
    lapply(1:2, function(x) {                   # LH
      lapply(allee, function(a.f){
        data <- PVA_data[PVA_data$lambda == l & PVA_data$LH == x & PVA_data$allee == a.f,]
        glm(P ~ log10(pop), data = data, family = binomial)
      })
    })
  })

  # MVP - est MVP
  MVP <- function(P, a, b) 10^(-((log(1/P - 1) + a) / b))
  
  MVP_sDD <- do.call(rbind, lapply(1:length(max.lambdas), function(l) { # max lambdas
    do.call(rbind, lapply(1:2, function(x) {                            # Life History
      do.call(rbind, lapply(1:length(allee), function(a.f) {            # allee effect
        do.call(rbind, lapply(c(0.05, 0.01), function(mvp){               # persistence prob
          as.data.frame(list(
            lambda = max.lambdas[l],
            LH = x,
            allee = allee[a.f],
            Pext = mvp,
            MVP = MVP(mvp, coef(mod.PVA[[l]][[x]][[a.f]])[[1]], coef(mod.PVA[[l]][[x]][[a.f]])[[2]])
          ))
        }))
      }))
    }))
  }))

  #-----------------------------------------------------------------------------
  # Plot
  
  # Data fro ploting
  # aggregate reps to est. prob of ext. per carrying capacity
  plot_data <- with(PVA_data, aggregate(P, by = list(pop, allee, LH, lambda), FUN = function(x) {sum(x)/reps}))
  names(plot_data) <- c("pop", "allee", "LH", "lambda", "P")
  
  plot_data$LH <- paste0("A", plot_data$LH)
  
  # Prediction data
  # predict prob of extiction lines for plot - for each lambda and llife history
  Na <- seq(min(CC), max(CC), length.out = 1000)           # Na to predict over
  pred_data <- do.call(rbind, lapply(1:length(max.lambdas), function(l) { # max lambdas
    do.call(rbind, lapply(1:2, function(x) {                              # life history
      do.call(rbind, lapply(1:length(allee), function(a.f){
      
        mod <- mod.PVA[[l]][[x]][[a.f]]
        pred <- predict(mod, list(pop = Na), type = "response") # predict prob of extinction
        
        as.data.frame(list(
          lambda = max.lambdas[l],
          LH = paste0("A", x), 
          allee = allee[a.f],
          pop = Na,
          P = pred
        ))
      }))
    }))
  }))
  
  # Plot - adult pop size vs. prob of extirpation
  options(scipen = 100)
  p <- ggplot(plot_data, aes(pop, P))
  p <- p + geom_hline(aes(yintercept = 0.01), colour = "black", linetype = "dashed")
  p <- p + geom_hline(aes(yintercept = 0.05), colour = "black", linetype = "dashed")
  p <- p + geom_point(aes(colour = as.factor(lambda), shape = as.factor(allee)), size = 1.5)
  p <- p + geom_path(data = pred_data, aes(pop, P, colour = as.factor(lambda), linetype = as.factor(allee)))
  p <- p + scale_x_continuous(trans = "log10")
  p <- p + scale_colour_manual(values = cbPalette)
  p <- p + facet_wrap(~LH, ncol = 2, scales = "free_x")
  p <- p + labs(x ="Adult population size ", y = "Probability of Extirpation", colour = ~lambda, 
                linetype = "Allee", shape = "Allee") 
  p <- p + theme_me + theme(legend.position = "top")
  
  png("Figures/MVP.png", width = 7.5, height = 4.5, unit = "in", res = 300)
  p
  dev.off()
}

MVP_sDD <- list(
  list('Tmax - 4 Tmat - 1 clutch - 3' = 
         list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                              MVP = c(7239, 27019))), 
              'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                               MVP = c(13203, 43186)))), 
       'Tmax - 3 Tmat - 2 clutch - 2' =
       list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                            MVP = c(8347, 28711))), 
            'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                             MVP = c(14927, 44018))))),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 
         list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                              MVP = c(3335, 11465))), 
              'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                               MVP = c(6523, 20866)))), 
       'Tmax - 3 Tmat - 2 clutch - 2' =
         list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                              MVP = c(4234, 14322))), 
              'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                               MVP = c(8116, 24430))))),
  list('Tmax - 4 Tmat - 1 clutch - 3' = 
         list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                              MVP = c(2451, 8403))), 
              'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                               MVP = c(4798, 15663)))), 
       'Tmax - 3 Tmat - 2 clutch - 2' =
         list('allee 50' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                              MVP = c(2884, 9354))), 
              'allee 100' = as.data.frame(list(P_ext = c(0.05, 0.01),
                                               MVP = c(5750, 17596)))))
)
names(MVP_sDD) <- L.names  

# clean up
rm(reps, years, p.cat, CC)

#===============================================================================
# REPATRIATION MODEL
#===============================================================================
# Assess the potential success of translocation of ESD 
# Run PVA model for ESD pop with inital density of 0. Add individuals age-1 to 
# Tmax for syears (1 to 10). Vary number of indiviudals translocated.
# a portion will die during and soon after stocking (translocation mortality).
# model included 3 levels of allee effects (weak to strong)
# model includes density dependence on YOY survival - 3 level of pop growth  
# 2 suites of life-history characters included
# run for 500 reps for each combinatio of variables 
#-------------------------------------------------------------------------------

reps <- 5000 # num. of replicates

# propability of catastrophe - per generation
p.cat <- 0.1

# Allee effect paraters
allee <- c(50, 100)

# Stocking density - age 1 to Tmax ESD
N.t <- round(exp(seq(log(10), log(2000), length.out = 8)))

# Stocking duration
stock.years <- c(1, 5, 10)

# Translocation mortality
M.trans <- seq(0.1, 0.9, length.out = 5)

#-------------------------------------------------------------------------------
# Run simulation

if(FALSE) {
for(stock_time in c("pre.spawn", "post.spawn")){ # Loop through stock times
  
  print(stock_time)
  print(Sys.time())
  
  # Set file path for loading results for pre/post-spawn
  file.path <- ifelse(stock_time == "pre.spawn",
                      paste0("Results/Results 02 - repatriation - ", stock_time),
                      paste0("Results/Results 03 - repatriation - ", stock_time))

  # Run simulation
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, list("lh_data", "repatriation",  "s_rand", "f_rand", "E_est", "lh_func",
                         "max.lambdas", "lh_params", "reps", "allee", "p.cat", "N.t", 
                         "stock.years", "M.trans", "pmx_eval", 'A.esd', "init_pop", 
                         "beta_stretch_val", "MVP_sDD", "pmx.1.det", "stock_time", "file.path")) # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  clusterSetRNGStream(cl, iseed = 123)
  
  rPatsim <- parLapply(cl, 1:length(max.lambdas), function(l) { # population growth rate
    lapply(1:length(lh_params), function(x) {                   # life history
      A <- lapply(1:length(allee), function(a.f) {                        # Allee effect
        lapply(N.t, function(Nt) {                              # Stocking density
          lapply(stock.years, function(syears) {                # years of stocking
            lapply(M.trans, function(M.t) {                     # translocation mortality
              
              data <- lh_data                           # data
              data$gen.time <- lh_data$gen.time.lh      # assign generatioin time - based on LH data
              
              data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
              data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
              
              # assign LH mean values based on the lh parameters 
              lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                                 Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                                 clutch = lh_params[[x]]$clutch, # number of clutches
                                 data = data)
              
              data$p.rep = c(lh_mean$MO1, 1, 1, 1) # Create proportion adult vector 
              
              K <- MVP_sDD[[l]][[x]][[a.f]]$MVP[2] # Carrying capacity set to MVP with 1% P[ext.]
              data$b <- data$b_dd[[l]][[x]]/K      # extract density dependent param and scale to carrying capacity
              data$a <- allee[a.f]                 # Allee effect parameters - NA at 50% fecundity
              
              # Translostion pop vectory - Age structure of translocated fish - Juv and Adults
              Nt.vec <- stable.stage(pmx.1.det[[x]]) * Nt
              
              out <- lapply(1:reps, function(i) {  # replicates
                
                res <- repatriation(mx = A.esd,          # projection matrix expression
                                    data = data,         # life history data
                                    lh_mean = lh_mean,   # Mean vital rates
                                    Ninit = 0,           # initial adult population size pop size
                                    Nt.vec = Nt.vec,     # transfered population as vector of ages 0:tmax
                                    s.t = 1 - M.t,       # translocation survival rate
                                    t.time = stock_time, # when stocking takes place ="pre.spawn" OR "post.spawn" 
                                    stock.month = 0,     # num. months p0st spawned stocking takes place
                                    stock.max = syears,  # number of stocking events to take place
                                    years = 50,          # years to run simulation
                                    p.cat = p.cat        # propability of catastrophe

                )
                N = res$pop$N 
                rm(res); gc()
                N
              }) # Replicates
            }) # translocation mortality
          }) # years of stocking
        }) # Stocking density
      }) # Allee
      save(A, file = paste0(file.path, " - lambda ", max.lambdas[l]," - LH ", x, ".R")) # save file - one per pop growth rate
      NULL
    }) # Life history
  }) # pop growth rate
  stopCluster(cl)
  
  rm(rPatsim); gc()
}
}

#-------------------------------------------------------------------------------
# Repatriation Results

# Load and extract data
for(stock_time in c("pre.spawn", "post.spawn")){
    
  print(stock_time)

  # Function to determine if repatriation was successful
  # Unseccessful if: Pop drops below 50 after stocking stops
  #                  If average pop size is < MVP.05 over last 15 years of sim
  FUN = function(x, s.y, MVP) {
    crit1 <- as.numeric(min(x[-(1:s.y)]) > 2)
    crit2 <- as.numeric(geom_mean(x[36:50]) > MVP)
    crit1 * crit2
  }

  # Set file path for loading results for pre/post-spawn
  file.path <- ifelse(stock_time == "pre.spawn",
                  paste0("Results/Results 02 - repatriation - ", stock_time),
                  paste0("Results/Results 03 - repatriation - ", stock_time))
  
  # Extract data
  # Determine if successful introduction
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, list("MVP_sDD", "max.lambdas", "allee", "N.t", "stock.years",
                         "M.trans", "reps", "FUN", "geom_mean", "stock_time", "file.path")) # send data to clusters
  clusterEvalQ(cl, library(ggplot2))# load library in clusters
  
  rpat_data <- do.call(rbind, lapply(1:length(max.lambdas), function(l){    # Maximum lambda 
    do.call(rbind, lapply(1:2, function(x){                             # Life history
      
      # Load simulatin data
      load(file = paste0(file.path, " - lambda ", 
                         max.lambdas[l]," - LH ", x, ".R"))
      
      do.call(rbind, lapply(1:length(allee), function(a.f) {                # Allee effect
        
        # extract MVP value at 05% extinction risk 
        MVP.05 <- MVP_sDD[[l]][[x]][[a.f]]$MVP[1]
      
        do.call(rbind, lapply(1:length(N.t), function(Nt) {               # Stocking density
          do.call(rbind, lapply(1:length(stock.years), function(syears) { # years of stocking
            do.call(rbind, lapply(1:length(M.trans), function(M.t) {      # Translocation mortality
              as.data.frame(list(                    # extract data
                lambda = max.lambdas[l],             # LAMBDA
                LH = x,                              # life history
                allee = allee[a.f],                  # allee effect
                Nt = N.t[Nt],                        # Stocking density
                s.years = stock.years[syears],       # years of stocking
                Mt = M.trans[M.t],                   # Translocation mortality
                rep = 1:reps,                        # replicate indicator
                P = sapply(A[[a.f]][[Nt]][[syears]][[M.t]], FUN,  s.y = stock.years[syears], MVP = MVP.05) # Success or not - Boolean
              ))
            })) # Translocation mortality
          })) # Years stocked
        })) # Stocking density
      })) # Allee effect
    })) # Life history
  })) # Max growth rate
  stopCluster(cl)
  
  #-----------------------------------------------------------------------------
  # Fit logistic regression to binary repatriation success
  
  rpat_mod <- lapply(1:2, function(x){  # Life history
    lapply(allee, function(a.f) {       # Allee effect             
            
      # subset data 
      subData <- rpat_data[rpat_data$LH == x & rpat_data$allee == a.f, ]
      
      # it logistic regression
      # success a funciton of stocking density + lambda + stocking yeras and Translocation mortality
      mod <- glm(P ~ log(Nt) + log(lambda) + (s.years) + (Mt), data = subData, family = binomial)
      mod
    }) # Allee effect
  }) # life history
  
  save(rpat_mod, file = paste0(file.path, " - logistic model.R"), row.names = F)
  
  # Output logistic regression parameters
  mod.out <- lapply(1:2, function(x){  # Life history
    a <- lapply(1:length(allee), function(a.f) {       # Allee effect        
      
      coefs <- as.data.frame(summary(rpat_mod[[x]][[a.f]])$coefficients)
      coefs$Vars <- as.vector(rownames(coefs))
      coefs$allee <- allee[a.f]
      coefs$LH <- x
      coefs
    })
    do.call(rbind, a)
  })
  mod.out <- do.call(rbind, mod.out)
  
  write.csv(mod.out, file = paste0(file.path, " - model.csv"), row.names = F)
  
  #-----------------------------------------------------------------------------
  # Plot Results
  # Logistic regression of success prob again stocking density with other variables
  
  # Aggregate Boolean rep data to prob of success
  plot_data <- with(rpat_data, aggregate(P, by = list(Mt, s.years, Nt, allee, LH, lambda), FUN = function(x) {sum(x)/reps}))
  names(plot_data) <- c("Mt", "s.years", "Nt", "allee", "LH", "lambda", "P")
  
  # Prediction data for plotting
  Nx <- seq(10,2000, 2) # x valse for prediction
  pred_data <- lapply(1:length(max.lambdas), function(l){     # max lambda
    X <- lapply(1:2, function(x){                             # life history
      A <- lapply(1:length(allee), function(a.f) {            # Allee effect                    
        Y <- lapply(1:length(stock.years), function(syears) { # Year sof stocking                   
          M <- lapply(1:length(M.trans), function(M.t) {      # Translocation mortality
            as.data.frame(list(              # Pridiction data frame
              lambda = max.lambdas[l],       # LAMBDA
              LH = x,                        # life history
              allee = paste0("Allee Effect: ", allee[a.f]),            # Carrying capacity
              s.years = paste0("Years Stocked: ",stock.years[syears]), # years of stocking
              Mt = M.trans[M.t],             # Translocation mortality
              Nt = Nx,                       # Sotkcing density
              P = predict(rpat_mod[[x]][[a.f]], # Predicted prob of success
                          data.frame(Nt = Nx, lambda = max.lambdas[l], 
                                     s.years = stock.years[syears], Mt = M.trans[M.t]), 
                          type = "response")
            ))
          }) # Translocationmortality
          do.call(rbind, M)
        }) # yeras of stocking
        do.call(rbind, Y)
      }) # allee effect
      do.call(rbind, A)
    }) # life history
    do.call(rbind, X)
  }) # max lambda
  pred_data <- do.call(rbind, pred_data)
  pred_data$allee <- factor(pred_data$allee, levels = unique(pred_data$allee), labels = unique(pred_data$allee))
  pred_data$s.years <- factor(pred_data$s.years, levels = unique(pred_data$s.years), labels = unique(pred_data$s.years))
  
  
  # Plot output - one fore each lambda LH combo; faceted by Mt and sYears
  lapply(max.lambdas, function(l){ # max lambda
    X <- lapply(1:2, function(x){  # Life history
      
      subPlot_data <- plot_data[plot_data$lambda == l & plot_data$LH == x,] # subset plot data
      subPred_data <- pred_data[pred_data$lambda == l & pred_data$LH == x,] # subset pred data
      
      # plot
      p <- ggplot(subPred_data, aes(Nt, P))
      p <- p + geom_path(data = subPred_data, aes(Nt, P, colour = as.factor(Mt)))
      p <- p + facet_grid(allee ~ s.years)
      p <- p + scale_x_continuous(trans = "log10", breaks = c(10, 100, 1000),
                                  minor_breaks = c(seq(20, 90, 10), seq(200, 900, 100), 2000))
      p <- p + scale_colour_manual(values = cbPalette)
      p <- p + labs(x ="Fish Stocked Annually", y = "Probability of Success", colour = "Translocation Mortality") 
      p <- p + theme_me + theme(legend.position = "top"
      )
      p
      
      # Save plots
      file = paste0("Figures/Repatriation - ", stock_time, " - lambda ", l," - LH ", x, ".R") # save file - one per pop growth rate
      png(paste0(file, ".png"), width = 11, height = 7.5, unit = "in", res = 300)
      print(p)
      dev.off()
      
      NULL
    }) # 
  })
  
  #-----------------------------------------------------------------------------
  # Find pop size when success is likely
  
  # Optimization functi to find Nt at 90% success rate
  optim_f <- function(Nx, mod, lambda, syears, Mt, target){
  
    P <- predict(mod, data.frame(Nt = Nx, lambda = lambda,  s.years = syears, Mt = Mt), 
                type = "response")
    
    (P - target)^2
  }
  
  # find success prob
  rpat_Nt <- lapply(1:length(max.lambdas), function(l){       # Max lambda
    X <- lapply(1:2, function(x){                             # Life history
      A <- lapply(1:length(allee), function(a.f) {            # Allee effect                    
        Y <- lapply(1:length(stock.years), function(syears) { # Stocking years    
          M <- lapply(1:length(M.trans), function(M.t) {      # Translocation mortality
            as.data.frame(list(
              lambda = max.lambdas[l],          # Lambda
              LH = x,                           # life history
              allee = paste0("Allee Effect: ", allee[a.f]), # Carrying capacity
              s.years = paste0("Years Stocked: ",stock.years[syears]),# Stocking years
              Mt = M.trans[M.t],                # translocation mortality
              N = optimize(optim_f, c(10, 1e6), tol = 1e-16, mod = rpat_mod[[x]][[a.f]],
                           lambda = max.lambdas[l], syears = stock.years[syears], M.trans[M.t],
                           target = 0.9)$minimum  
            ))
          }) # Translocation mortality
          do.call(rbind, M)
        }) # Yera stocking
        do.call(rbind, Y)
      }) # Allee effect
      do.call(rbind, A)
    }) # Life history
    do.call(rbind, X)
  }) # Max lambda
  rpat_Nt <- do.call(rbind, rpat_Nt)
  
  rpat_Nt$allee <- factor(rpat_Nt$allee, levels = unique(rpat_Nt$allee), labels = unique(rpat_Nt$allee))
  rpat_Nt$s.years <- factor(rpat_Nt$s.years, levels = unique(rpat_Nt$s.years), labels = unique(rpat_Nt$s.years))
  
  
  # set very small or large values to constants
  #rpat_Nt$N[rpat_Nt$N <= 1] = 1    # Small
  #rpat_Nt$N[rpat_Nt$N >= 1e4] = NA # large
  
  # Save output as CSV
  write.csv(rpat_Nt, file = paste0(file.path, ".csv"), row.names = FALSE)
  
  # Plot - stocking density at 75% success again variables
  lapply(1:2, function(x) { # Life history
    # subset data
    subPlot_data <- rpat_Nt[rpat_Nt$LH == x,]
    
    # plot
    p <- ggplot(subPlot_data, aes(lambda, N))
    p <- p + geom_point(aes(colour = as.factor(Mt)), size = 1.5)
    p <- p + geom_path(aes(colour = as.factor(Mt)), size = 1)
    p <- p + scale_x_continuous(breaks = seq(1.5, 2.75, 0.25))
    p <- p + scale_colour_manual(values = cbPalette)
    p <- p + facet_grid(allee~s.years)
    p <- p + labs(x ="Population Growth Rate", y = "Stocking Density", colour = "Translocation mortality") 
    #p <- p + scale_y_continuous(limits= c(0, 10000)) 
    p <- p + theme_me + theme(legend.position = "top")
    
    # save
    file = paste0("Figures/Repatriation Density - ", stock_time, " - LH ", x, ".R") # save file - one per pop growth rate
    png(paste0(file, ".png"), width = 11, height = 7.5, unit = "in", res = 300)
    print(p)
    dev.off()
    
    NULL
    
  })
  
} # clost stock time loop

#===============================================================================
# SOURCE POPULATION HARM 
#===============================================================================
# Assess harm to source population 
# Run PVA model of source population. Inital pop size ranges form MVP to 10*MVP
# Harm scenarios same and repatriation model: various stocking denisties &
# stocking durations.
# Nt individuals remove from pop for syears with viability tested later
#-------------------------------------------------------------------------------

reps <- 5000 # num. of replicates

# propability of catastrophe - per generation
p.cat <- 0.1

# Allee effect paraters
allee <- c(50, 100)

# Stocking density - age 1 to Tmax ESD
N.t <- round(exp(seq(log(10), log(2000), length.out = 8)))

# range of carrying capacities to loop trough 
K <- round(exp(seq(log(5000), log(50000), length.out = 5)))     

# Stocking duration
stock.years <- c(1, 5, 10)

#-------------------------------------------------------------------------------
# Run simulation
if(FALSE) {
for(stock_time in c("pre.spawn", "post.spawn")){ # Loop through stock times
  
  print(stock_time)  
  
  # Set file path for loading results for pre/post-spawn
  file.path <- ifelse(stock_time == "pre.spawn",
                      paste0("Results/Results 04 - Harm - ", stock_time),
                      paste0("Results/Results 05 - Harm - ", stock_time))

  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, list("lh_data", "source_harm",  "s_rand", "f_rand", "E_est", "lh_func",
                         "max.lambdas", "lh_params", "reps",  "p.cat", "N.t", "stock_time",
                         "stock.years", "allee", "pmx_eval", 'A.esd', "init_pop", "file.path",
                         "beta_stretch_val", "MVP_sDD", "pmx.1.det", "K")) # send data to clusters
  clusterEvalQ(cl, library(popbio))# load library in clusters
  clusterSetRNGStream(cl, 123)
  
  Harmsim <- parLapply(cl, 1:length(max.lambdas), function(l) { # population growth rate
    lapply(1:length(lh_params), function(x) {                   # life history
      D <- lapply(1:length(allee), function(a.f){               # Allee effect
        lapply(N.t, function(Nt) {                              # Stocking density
          lapply(stock.years, function(syears) {                # years of stocking
            lapply(K, function(k){                              # Carrying Capacity
              data <- lh_data                           # data
              data$gen.time <- lh_data$gen.time.lh      # assign generatioin time - based on LH data
              
              data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
              data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
              
              # assign LH mean values based on the lh parameters 
              lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                                 Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                                 clutch = lh_params[[x]]$clutch, # number of clutches
                                 data = data)
              
              data$p.rep = c(lh_mean$MO1, 1, 1, 1) # Create proportion adult vector 
              
              data$b <- data$b_dd[[l]][[x]]/k      # extract density dependent param and scale to carrying capacity
              data$a <- allee[a.f]                 # Allee effect parameters - NA at 50% fecundity
              
              # Translostion pop vectory - Age structure of translocated fish - Juv and Adults
              Nt.vec <- stable.stage(pmx.1.det[[x]]) * Nt
              
              out <- lapply(1:reps, function(i) {  # replicates
                res <- source_harm(mx = A.esd,                       # projection matrix expression
                                   data = data,                      # life history data
                                   lh_mean = lh_mean,                # Mean vital rates
                                   Ninit = k,                        # initial adult population size pop size
                                   Nt.vec = Nt.vec,                  # transfered population as vector of ages 0:tmax
                                   t.time = stock_time,              # when stocking takes place ="pre.spawn" OR "post.spawn" 
                                   stock.max = syears,               # number of stocking events to take place
                                   years = 50,                       # years to run simulation
                                   p.cat = 0.1                       # propability of catastrophe
                )
                N = res$pop$N 
                rm(res); gc()
                N
              }) 
            })
          }) # years of stocking
        }) # Stocking density
      }) # allee effects
      save(D, file = paste0(file.path, " - lambda ", max.lambdas[l]," - LH ", x, ".R")) # save file - one per pop growth rate
      NULL
    }) # Life history
  }) # pop growth rate
  stopCluster(cl)
  
  rm(Harmsim)
} # Stock time
} # if

#-------------------------------------------------------------------------------
# Results
# Load and extract data
for(stock_time in c("pre.spawn", "post.spawn")){
  
  print(stock_time)
  
  # Set file path for loading results for pre/post-spawn
  file.path <- ifelse(stock_time == "pre.spawn",
                      paste0("Results/Results 04 - Harm - ", stock_time),
                      paste0("Results/Results 05 - Harm - ", stock_time))
  
  no_cores <- detectCores() - 1    # number of cores
  cl <- makeCluster(no_cores)      # create clusters
  clusterExport(cl, list("MVP_sDD", "max.lambdas", "N.t", "stock.years", "K",
                         "allee", "reps", "file.path")) # send data to clusters
  clusterEvalQ(cl, library(ggplot2))# load library in clusters
  
  Harm.data <- do.call(rbind,  parLapply(cl, 1:length(max.lambdas), function(l){
    do.call(rbind,  lapply(1:2, function(x){
      
      # Load simulatin data
      load(file = paste0(file.path, " - lambda ", max.lambdas[l]," - LH ", x, ".R"))
      
      do.call(rbind, lapply(1:length(allee), function(a.f){
        do.call(rbind,  lapply(1:length(N.t), function(Nt) {                                 # Stocking density
          do.call(rbind,  lapply(1:length(stock.years), function(syears) {                   # years of stocking
            do.call(rbind,  lapply(1:length(K), function(k) {
              as.data.frame(list(                    
                lambda = max.lambdas[l],             # LAMBDA
                LH = x,                              # life history
                allee = allee[a.f],
                Nt = N.t[Nt],  
                s.years = stock.years[syears],
                K = K[k],
                P = as.numeric(sapply(D[[a.f]][[Nt]][[syears]][[k]], function(i) any(i < 2))) # Population size
              )) 
            })) # Carry capacity
          })) # Years of stocking
        })) # Stocking density - harm
      })) # allee effect
    })) # Life History
  })) # max lambda
  stopCluster(cl)
  
  # Add simulation for 5000 transplants - frome extra simulation
  Harm.data.5000 <- do.call(rbind,  lapply(1:length(max.lambdas), function(l){
    do.call(rbind,  lapply(1:2, function(x){
      
      # Load simulatin data
      load(file = paste0(file.path, " - lambda ", max.lambdas[l]," - LH ", x, "- N.t 5000.R"))
      
      do.call(rbind, lapply(1:length(allee), function(a.f){
        do.call(rbind,  lapply(5000, function(Nt) {                                 # Stocking density
          do.call(rbind,  lapply(1:length(stock.years), function(syears) {                   # years of stocking
            do.call(rbind,  lapply(1:length(K), function(k) {
              as.data.frame(list(                    
                lambda = max.lambdas[l],             # LAMBDA
                LH = x,                              # life history
                allee = allee[a.f],
                Nt = Nt,  
                s.years = stock.years[syears],
                K = K[k],
                P = as.numeric(sapply(D[[a.f]][[1]][[syears]][[k]], function(i) any(i < 2))) # Population size
              )) 
            })) # Carry capacity
          })) # Years of stocking
        })) # Stocking density - harm
      })) # allee effect
    })) # Life History
  })) # max lambda
  
  # combine data frame
  Harm.data <- as.data.frame(cbind(Harm.data, Harm.data.5000))
  
  # clean up
  rm(Harm.data.5000)
  
  #-----------------------------------------------------------------------------
  # Fit logistic regression to binary repatriation success

  harm_mod <- lapply(1:length(max.lambdas), function(l){ # max lambda
    lapply(1:2, function(x){                          # Life history
      lapply(1:length(allee), function(a.f){
        lapply(stock.years, function(syears) {        # Years of stocking             
              
          # subset data 
          subData <- Harm.data[Harm.data$lambda == max.lambdas[l] & Harm.data$LH == x &            
                               Harm.data$allee == allee[a.f] & Harm.data$s.years == syears, ]
              
          # it logistic regression
          mod <- glm(P ~ log(Nt) + log(K), data = subData, family = binomial)
          mod
        }) # Yeras of stocking
      }) # allee efects
    }) # life history
  }) # max lambda
  
  save(harm_mod, file = paste0(file.path, " - logistic model.R"), row.names = F)
  
  print(lapply(1:length(max.lambdas), function(l){    # max lambda
    lapply(1:2, function(x){                          # life history
      lapply(1:length(allee), function(a.f){          # allee effects
          P = predict(harm_mod[[l]][[x]][[a.f]][[3]], # Predicted prob of success
                      data.frame(Nt = 1000, K = 25000), type = "response")*100
      })
    })
  }))
  
  # Optimization furntion - solve for syoy
  s_optim = function (syoy, mx, data, Tmax, Tmat, clutch, target.lambda){
    
    data$Syoy <- syoy
    
    # Population matrix
    A <- pmx_eval(mx, lh_func(Tmax = Tmax, Tmat = Tmat, clutch = clutch, data = data))
    
    # lambda
    lambda <- lambda(A)
    
    # golden search algorithm
    (lambda - target.lambda)^2
  }
  
  #-------------------------------------------------------------------------------
  
  # optimization - lambda = 1
  lh_data$s0.1.det <- lapply(1:length(lh_params), function(x) {
    s_optim_min <- optimize(s_optim, c(0, 1), tol = 1e-16, mx = A.esd, 
                            data = lh_data, Tmax = lh_params[[x]]$Tmax, 
                            Tmat = lh_params[[x]]$Tmat, clutch = lh_params[[x]]$clutch,
                            target.lambda = 1)$minimum    # lambda = min
  })
  names(lh_data$s0.1.det)  <- lh_names
  
  # Output logistic regression parameters
  mod.out <- do.call(rbind, lapply(1:length(max.lambdas), function(l){ # max lambda
    do.call(rbind, lapply(1:2, function(x){                            # Life history
      do.call(rbind, lapply(1:length(allee), function(a.f){
        do.call(rbind, lapply(1:length(stock.years), function(syears) {  # Years of stocking             
  
          coefs <- as.data.frame(summary(harm_mod[[l]][[x]][[a.f]][[syears]])$coefficients)
          coefs$Vars <- as.vector(rownames(coefs))
          coefs$s.years <- stock.years[syears]
          coefs$allee <- allee[a.f]
          coefs$LH <- x
          coefs$lambda <- max.lambdas[l]
          coefs
        }))
      }))
    }))
  }))

  write.csv(mod.out, file = paste0(file.path, " - model.csv"), row.names = F)
  
  # Prediction data for plotting
  Nx <- seq(0, 2000, 10) # x valse for prediction
  Kx <- 10^(seq(log10(5000), log10(60000), length.out = 100))
  xvals <- expand.grid(Nx, round(Kx))
  
  pred_data <- do.call(rbind, lapply(1:length(max.lambdas), function(l){ # max lambda
    do.call(rbind, lapply(1:2, function(x){                              # life history
      do.call(rbind, lapply(1:length(allee), function(a.f){                             # allee effects
        do.call(rbind, lapply(1:length(stock.years), function(syears) {  # Year sof stocking                   
          as.data.frame(list(              # Pridiction data frame
            lambda = max.lambdas[l],       # LAMBDA
            LH = x,                        # life history
            allee = allee[a.f],            # allee effect
            s.years = stock.years[syears], # years of stocking
            Nt = xvals[, 1],               # Sotkcing density
            K = xvals[, 2],                # Carrying Capacity
            P = predict(harm_mod[[l]][[x]][[a.f]][[syears]], # Predicted prob of success
                        data.frame(Nt = xvals[,1], K = xvals[,2]), type = "response")
          ))
        })) # yeras of stocking
      })) # allee effects
    })) # life history
  })) # max lambda

  # Plot - stocking density at 75% success again variables
  lapply(1:2, function(x) { # Life history
    lapply(1:length(allee), function(a.f){
      # subset data
      subPlot_data <- pred_data[pred_data$LH == x & pred_data$allee == allee[a.f],]
      subPlot_data$s.years <- paste0("Years of removals: ", subPlot_data$s.years)
      subPlot_data$s.years <- factor(subPlot_data$s.years, levels = unique(subPlot_data$s.years),
                                     labels = unique(subPlot_data$s.years), ordered = T)
      
      subPlot_data$K <- subPlot_data$K/1000
      
      MVP.val <- do.call(rbind, lapply(1:length(MVP_sDD), function(l) {
        as.data.frame(list(lambda = max.lambdas[l], s.years = c(1, 5, 10), MVP = MVP_sDD[[l]][[x]][[a.f]][2,2]))
      }))
      MVP.val$s.years <- paste0("Years of removals: ", MVP.val$s.years)
      MVP.val$s.years <- factor(MVP.val$s.years, levels = unique(MVP.val$s.years),
                                     labels = unique(MVP.val$s.years), ordered = T)
      
      # plot
      p <- ggplot(subPlot_data, aes(Nt, K))
      p <- p + geom_raster(aes(fill = P), interpolate = T)
      p <- p + geom_hline(data = MVP.val, aes(yintercept = MVP/1000))
      p <- p + scale_fill_distiller(palette= "Spectral", direction=-1, limits = c(0, 1))
      p <- p + scale_y_continuous(trans = "log10", breaks = c(5, seq(10, 50, 10)), 
                                  minor_breaks = c(seq(60, 90, 10), seq(20, 70, 10)),
                                  expand = c(0, 0))
      p <- p + scale_x_continuous(expand = c(0, 0))
      p <- p + facet_grid(lambda ~ s.years)
      p <- p + labs(x ="Removal Density", y = "Carrying Capacity (x 1000)", fill = "Extirpation Probability") 
      p <- p + theme_me + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
      
      # save
      file = paste0("Figures/Harm Ext prob - ", stock_time, " - LH ", x, "- allee ", allee[a.f], ".R") # save file - one per pop growth rate
      png(paste0(file, ".png"), width = 11, height = 7.5, unit = "in", res = 300)
      print(p)
      dev.off()
      
      NULL
    })
  })
  
}

#-------------------------------------------------------------------------------
# Success Vs. Harm
#-------------------------------------------------------------------------------

# Allee effect paraters
allee <- c(50, 100)

# Stocking density - age 1 to Tmax ESD
N.t <- round(exp(seq(log(10), log(2000), length.out = 8)))

# Stocking duration
stock.years <- c(1, 5, 10)

# Translocation mortality
M.trans <- seq(0.1, 0.9, length.out = 5)

for(stock_time.r in c("pre.spawn", "post.spawn")){   # Stock time repatiration
  for(stock_time.h in c("pre.spawn", "post.spawn")){ # Stock tim harm
    print(paste0("r - ", stock_time.r))
    print(paste0("h - ", stock_time.h))
    
    # Set file path for loading results for pre/post-spawn
    file.path.rpat <- ifelse(stock_time.r == "pre.spawn",
                        paste0("Results/Results 02 - repatriation - ", stock_time.r),
                        paste0("Results/Results 03 - repatriation - ", stock_time.r))
    
    file.path.harm <- ifelse(stock_time.h == "pre.spawn",
                        paste0("Results/Results 04 - Harm - ", stock_time.h),
                        paste0("Results/Results 05 - Harm - ", stock_time.h))
    
    load(paste0(file.path.rpat, " - logistic model.R"))
    load(paste0(file.path.harm, " - logistic model.R"))
    
    #-----------------------------------------------------------------------------
    # predict K for source pop that give 1% extinction prob when repatriation is 90% successfull
    
    # Optimization function to find Nt needed for 90% repat success
    optim_rpat <- function(N, mod, lambda, s.years, Mt, target){
      
      pred = predict(mod, data.frame(Nt = N, lambda = lambda, s.years = s.years, Mt = Mt), 
                     type = "response")
      
      # golden search algorithm
      (pred - target)^2
    }
    
    # Optimization function to find K needed to 99% presistence prob when Nt90 removed
    optim_harm <- function(C, mod, Nt, target){
      
      pred = predict(mod, data.frame(Nt = Nt, K = C), type = "response")
      
      # golden search algorithm
      (pred - target)^2
    }
    
    pred_data <- do.call(rbind, lapply(1:length(max.lambdas), function(l){ 
      do.call(rbind,lapply(1:length(stock.years), function(syears) {
        
        # Repariation model
        mod.rpat <- rpat_mod[[1]][[2]]  #lambda = l, LH1, allee = 100
        
        # Harm model
        mod.harm <- harm_mod[[l]][[1]][[2]][[syears]]  # lambda = l, LH1, allee = 100, stock.years = syears
        
        # Est stocking density requred to acheive 90% success
        N <- optimize(optim_rpat, c(1, 1e4), tol = 1e-16, mod = mod.rpat, 
                 lambda = max.lambdas[l], s.years = stock.years[syears], Mt = 0.5, target = 0.9)$minimum  
        
        # Est. K need to maintain pop with N removals
        source.K <- optimize(optim_harm, c(1, 1e6), tol = 1e-16, mod = mod.harm, Nt = N,
                 target = 0.01)$minimum  
        
        # Output data
        as.data.frame(list(
          harm_time = stock_time.h,     # When harm occurs
          rpat_time = stock_time.r,     # When stocking occurs
          lambda = max.lambdas[l],      # pop growth rate
          LH = 1,                       # Life history type
          allee = 100,                  # Allee effect
          s.years = stock.years[syears],# years of stocking/removal
          Mt = 0.5,                     # tanslocation mortality
          stock.density = N,            # stocking/removal density giving 90% success
          source.K = source.K           # Carrying capacity giving 99% persistence with N removals
        ))
      }))
    }))
    
    # save data
    write.csv(pred_data, file = paste0("Results/Results 06 - Harm K Repat N ", stock_time.r, " - ", stock_time.h, ".csv"),
              row.names = F)
    
    # Carrying capacity required for source pop to have 1% extinction prob when
    # 500 fish removed for 10 years and lambda = 2.13
    print(paste("Source K - 500 removals over 10 year with 2.13 lambda:",
      optimize(optim_harm, c(1, 1e6), tol = 1e-16, mod = harm_mod[[2]][[1]][[2]][[3]], 
               Nt = 500, target = 0.01)$minimum)
    )
    
    #-------------------------------------------------------------------------------
        # Prediction data for plotting
    # Lambda = 2.3
    # Life History 1
    # stock years = 10
    # Translocation mortality = 0.5
    # Allee = 50 OR 100
    # K = MVP.01
  
    Nx <- exp(seq(log(50), log(5000), length.out = 100)) # Stocking density (harm)
    Kx <- c(10000, 20000, 30000, 50000)                 # Carrying capacity
    xvals <- expand.grid(Nx, round(Kx))
    
    pred_data <- do.call(rbind, lapply(1:length(max.lambdas), function(l) {
      do.call(rbind, lapply(1:length(stock.years), function(y) {
        do.call(rbind, lapply(1:length(allee), function(a.f) {            # Allee effect                    
          as.data.frame(list(              # Pridiction data frame
            lambda = max.lambdas[l],       # LAMBDA
            LH = 1,                        # life history
            allee = paste0("Allee Effect: ", allee[a.f]),            # Carrying capacity
            s.years = paste0("Years Stocked: ", stock.years[y]), # years of stocking
            Mt = M.trans[3],             # Translocation mortality
            K = xvals[,2],
            Nt = xvals[,1],                       # Sotkcing density
            P.succ = predict(rpat_mod[[1]][[a.f]], # Predicted prob of success
                        data.frame(Nt = xvals[,1], lambda = max.lambdas[l], 
                                   s.years = stock.years[y], Mt = M.trans[3]), 
                        type = "response"),
            P.harm = predict(harm_mod[[l]][[1]][[a.f]][[y]], # Predicted prob of success
                        data.frame(Nt = xvals[,1], K = xvals[,2]), type = "response")
          ))
        }))
      }))
    })) # max lambda
    pred_data$allee <- factor(pred_data$allee, levels = unique(pred_data$allee), labels = unique(pred_data$allee))
    pred_data$s.years <- factor(pred_data$s.years, levels = unique(pred_data$s.years), labels = unique(pred_data$s.years))
    
    lab_data <- with(pred_data, aggregate(P.harm, by = list(K, allee, s.years, lambda), FUN = min))
    colnames(lab_data) <- c("K", "allee", "s.years", "lambda", "y")
   
    lab_data$x <- with(pred_data, aggregate(P.succ, by = list(K, allee, s.years, lambda), FUN = min))[,5]
    
    lines <- list(
      data.frame(list(x = c(0.9, 0.9), y = c(0.001, 0.01), K = NA)),
      data.frame(list(x = c(0.9, 1), y = c(0.01, 0.01)), K = NA)
    )
    
    # Plot
    p <- ggplot(pred_data[pred_data$allee == "Allee Effect: 100", ], aes(P.succ, P.harm, group = K))
    p <- p + geom_line(data = lines[[1]], aes(x = x, y = y), size = 1)
    p <- p + geom_line(data = lines[[2]], aes(x = x, y = y), size = 1)
    p <- p + geom_path(aes(colour = (Nt)), size = 1)
    p <- p + geom_text(data = lab_data[lab_data$allee == "Allee Effect: 100", ], 
                        aes(x = x, y = y, label = K/1000), nudge_x = -0.05, vjust = 0.5, check_overlap = T)
    p <- p + scale_colour_gradientn(colours = brewer.pal(3,"Set2"))
    p <- p + labs(x ="Probability of Success", y = "Probability of Extirpation", colour = "Stocking \nDensity") 
    p <- p + facet_grid(lambda ~ s.years, scales = "free_x")
    p <- p + scale_y_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1)) 
    p <- p + scale_x_continuous(breaks = seq(0, 1, 0.2)) 
    p <- p + theme_me + theme(legend.position = "right")
   
    file = paste0("Figures/Success ", stock_time.r, " vs Harm ", stock_time.h, ".R") # save file - one per pop growth rate
    png(paste0(file, ".png"), width = 8.5, height = 8, unit = "in", res = 300)
    print(p)
    dev.off()
  }
}

#################################################################################
#################################################################################

reps = 1000
data <- lh_data
syears <- 10
l <- 1
x <- 2
Ninit =  MVP_sDD[[l]][[x]]$MVP[2] 
Nt <- 000
data$gen.time <- lh_data$gen.time.lh[[x]] # assign generatioin time - based on LH data
data$Syoy = data$s0.1.sto[[x]]            # Syoy at lambda = 1
data$s0.max = data$s0.max.sto[[l]][[x]]   # Syoy at lambda = max
# assign LH mean values based on the lh parameters 
lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax,     # Longevity
                   Tmat = lh_params[[x]]$Tmat,     # age-at maturity
                   clutch = lh_params[[x]]$clutch, # number of clutches
                   data = data)

data$p.rep = c(lh_mean$MO1, 1, 1, 1) # Create proportion adult vector 
data$b <- data$b_dd[[l]][[x]]/Ninit      # extract density dependent param and scale to carrying capacity

# Translostion pop vectory - Age structure of translocated fish
Nt.vec <- stable.stage(pmx.1.det[[x]]) * Nt
out <- lapply(1:reps, function(i) {  # replicates
  res <- source_harm(mx = A.esd,                                  # projection matrix expression
                      data = data,                                # life history data
                      lh_mean = lh_mean,                          # Mean vital rates
                      Ninit = Ninit,                              # initial adult population size pop size
                      Nt.vec = Nt.vec,                            # transfered population as vector of ages 0:tmax
                      t.time = "post.spawn",                      # when stocking takes place ="pre.spawn" OR "post.spawn" 
                      stock.max = syears,                         # number of stocking events to take place
                      years = 20,                                 # years to run simulation
                      p.cat = 0.1,                                # propability of catastrophe
                     density_dependence = "Survival"
  )
  N = res$pop$N
  rm(res); N
})

sum(sapply(out, function(i){
  any(i < 50)
}))/reps

mean(sapply(out, function(i){
  i[20]
}))

log.lambdas <- sapply(out, function(i) {
  i <- i[is.finite(i)]
  (log(i[length(i)]) - log(i[1])) / length(i)
})
length(which(exp(log.lambdas) < 1))/reps
length(which(exp(log.lambdas) ==0))/reps

hist(exp(log.lambdas), breaks = 50)
N0 <- 1000
reps <- 1000

no_cores <- detectCores() - 1    # number of cores
cl <- makeCluster(no_cores)      # create clusters
clusterExport(cl, list("lh_data", "Projection", "N0", "reps", "s_rand", 
                       "f_rand", "E_est", "lh_func", "lh_params", "max.lambdas",
                       "pmx_eval", 'A.esd', "init_pop", "beta_stretch_val", 
                       "lh_names", "geom_mean")) # send data to clusters
clusterEvalQ(cl, library(popbio))# load library in clusters

pop<- parLapply(cl, 1:length(max.lambdas), function(l) { # population growth rate
  lhs <- lapply(1:length(lh_params), function(x) {              # life history
    data <- lh_data
    data$gen.time <- data$gen.time.lh[[x]]
    data$Syoy = data$s0.1.sto[[x]]
    data$s0.max = data$s0.max.sto[[l]][[x]]
    lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax, 
                       Tmat = lh_params[[x]]$Tmat, 
                       clutch = lh_params[[x]]$clutch, data = data)
    data$p.rep <- c(lh_mean$MO1, 1, 1, 1)
    data$b <- data$b_dd[[l]][[x]]/N0
        
    out <- lapply(1:reps, function(i) {                # replicates
      
      res <- Projection(mx = A.esd,                        # projection matrix expression
                        data = data,                       # life history data
                        lh_mean = lh_mean,
                        N0 = N0,                  # initial pop size as stage-structure vector
                        years = 200,              # years to run simulation
                        p.cat = 0,                # propability of catastrophe
                        density_dependence = TRUE # type of density dependence
                        
      )
      res$pop$N[-(1:100)]
    })
    as.data.frame(list(
      pop = mean(sapply(out, geom_mean)),
      lambda = exp(mean(sapply(out, function(i) (log(i[length(i)]) - log(i[1]))/length(i))))
    ))
    
  })
})

stopCluster(cl)
pop

N0 <- 10000
reps <- 1000

no_cores <- detectCores() - 1    # number of cores
cl <- makeCluster(no_cores)      # create clusters
clusterExport(cl, list("lh_data", "Projection", "N0", "reps", "s_rand", 
                       "f_rand", "E_est", "lh_func", "lh_params", "max.lambdas",
                       "pmx_eval", 'A.esd', "init_pop", "beta_stretch_val", 
                       "lh_names", "geom_mean")) # send data to clusters
clusterEvalQ(cl, library(popbio))# load library in clusters

pop<- parLapply(cl, 1:length(max.lambdas), function(l) { # population growth rate
  lhs <- lapply(1:length(lh_params), function(x) {              # life history
    
    data <- lh_data
    data$gen.time <- data$gen.time.lh
    data$Syoy = data$s0.max.sto[[l]][[x]]
    data$s0.max = data$s0.max.sto[[l]][[x]]
    lh_mean <- lh_func(Tmax = lh_params[[x]]$Tmax, 
                       Tmat = lh_params[[x]]$Tmat, 
                       clutch = lh_params[[x]]$clutch, data = data)
    data$p.rep <- c(lh_mean$MO1, 1, 1, 1)
    data$b <- data$b_dd[[l]][[x]]/N0
    
    out <- lapply(1:reps, function(i) {                # replicates
      
      res <- Projection(mx = A.esd,                    # projection matrix expression
                        data = data,                   # life history data
                        lh_mean = lh_mean,
                        N0 = N0,                       # initial pop size as stage-structure vector
                        years = 200,                   # years to run simulation
                        p.cat = 0,                     # propability of catastrophe
                        density_dependence = FALSE # type of density dependence
                        
      )
      res$pop$N[-(1:100)]
    })
    as.data.frame(list(
      pop = mean(sapply(out, geom_mean)),
      lambda = exp(mean(sapply(out, function(i) (log(i[length(i)]) - log(i[1]))/length(i))))
    ))
    
  })
})

stopCluster(cl)
pop
