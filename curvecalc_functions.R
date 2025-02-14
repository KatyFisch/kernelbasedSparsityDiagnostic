
# dose-response curve
get_drc <- function(doses, Dset, Odat){
  # simulate truth
  true_drc <- foreach(i = doses, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- node("A.", distr = "rconst", const = i)
    Dset <- Dset + action(paste("A+", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("A+", i, sep="")), n = 1000000, rndseed = 1, 
                              verbose=F)
    mean(dat_int[[paste("A+", i, sep="")]]$Y) 
  }
  df_drc <- data.frame(dose=doses, truth=true_drc)
  # estimate with g-formula
  estimates_drc <- foreach(i = doses, .combine = 'c', .packages = c("CICI", "simcausal", "lmtp")) %dopar% { 
    est <- gformula(X=Odat[c("L", "A", "Y")], Ynodes = c("Y"), Anodes = c("A"), 
                    abar=i, verbose=F, Yform="GAM", seed = 1)
    est$results$psi
    #f <- function(data, trt){rep(i, length(data[[trt]]))}
    #est <- lmtp_sub(Odat, "A", "Y", baseline="L", shift=f, outcome_type="continuous", folds=1,
    #                learners = "SL.gam")
    #est$theta
  }
  df_drc["est"] <- estimates_drc
  return(df_drc)
}


# deterministic lmtp
get_src_det <- function(shifts, Dset, sl_lib, Odat){
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "rconst", const = A)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- node("A.", distr = "rconst", const = A + i)
    Dset <- Dset + action(paste("Ap", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ap", i, sep="")), n = 1000000, rndseed = 1, 
                              verbose=F)
    mean(dat_int[[paste("Ap", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
    f <- function(data, trt){ data[[trt]] + i }
    lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f, 
                        mtp = TRUE, folds = 2, outcome_type = "continuous",
                        learners_trt = sl_lib, learners_outcome = sl_lib)
    result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
    list(theta = result$theta, low = result$conf.low, high = result$conf.high)
  }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# stochastic lmtp
get_src_sto <- function(shifts, Dset, sl_lib, Odat){
  # function to estimate conditional density
  dens <- function (n, a, var, shift=0) { 
    dat <- data.frame(var)
    dat["a"] <- a
    dens <- lm(a~., data = dat)
    means <- predict(dens, dat)
    draws <- rnorm(n = n, mean = means + shift, sd = summary(dens)$sigma)
    return(draws)
  }
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "dens", a=A, var=L, shift=0)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal"), .export="dens") %do% {
    Aplus1 <- node("A.", distr = "dens", a=A, var=L, shift=i)
    Dset <- Dset + action(paste("As", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("As", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("As", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f_si <- function(data, trt){
       dat <- data[,2:which(colnames(data)==trt)-1]
       dat <- dat[ , grepl( "L" , names( dat ) ) ]
       dat <- data.frame(dat)
       dat["a"] <- data[[trt]]
       dens <- lm(a~., data = dat)
       means <- predict(dens, dat)
       set.seed(1)
       draws <- rnorm(n = length(dat$a), mean = means + i, sd = summary(dens)$sigma)
       draws
     }
     lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f_si, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# dynamic deterministic lmtp
get_src_detdy <- function(shifts, Dset, sl_lib, Odat){
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "rconst", const = A)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shifted intervention
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- node("A.", distr = "rconst", const = ifelse(L==1, A + i, A))
    Dset <- Dset + action(paste("Ap", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ap", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("Ap", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shifted intervention
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ ifelse(data[["L"]] == 1, data[[trt]] + i, data[[trt]]) }
     lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# dynamic stochastic lmtp
get_src_stody <- function(shifts, Dset, sl_lib, Odat){
  # function to estimate conditional density
  dens <- function (n, a, var, shift=0) { 
    dat <- data.frame(var)
    dat["a"] <- a
    dens <- lm(a~., data = dat)
    means <- predict(dens, dat)
    draws <- rnorm(n = n, mean = ifelse(var==1, means + shift, means), sd = summary(dens)$sigma)
    return(draws)
  }
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "dens", a=A, var=L, shift=0)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, verbose = F)
  # simulate truth
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal"), .export="dens") %do% {
    Aplus1 <- node("A.", distr = "dens", a=A, var=L, shift=i)
    Dset <- Dset + action(paste("As", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("As", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("As", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f_si <- function(data, trt){
       dat <- data[,2:which(colnames(data)==trt)-1]
       dat <- dat[ , grepl( "L" , names( dat ) ) ]
       dat <- data.frame(dat)
       dat["a"] <- data[[trt]]
       dens <- lm(a~., data = dat)
       means <- predict(dens, dat)
       set.seed(1)
       draws <- rnorm(n = length(dat$a), mean = ifelse(data[['L']]==1, means + i, means), sd = summary(dens)$sigma)
       draws
     }
     lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f_si, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
  }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# deterministic threshold lmtp
get_src_detth <- function(thresholds, Dset, sl_lib, Odat, left=TRUE){
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "rconst", const = A)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shifted intervention
  responses_thresh <- foreach(i = thresholds, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    if (left){Athresh <- node("A.", distr = "rconst", const = ifelse(A<i, i, A))}
    else {Athresh <- node("A.", distr = "rconst", const = ifelse(A>i, i, A))}
    Dset <- Dset + action(paste("At", i, sep=""), nodes=Athresh)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("At", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("At", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(threshold=thresholds, response=responses_thresh)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shifted intervention
  estimates_shift <- foreach(i = thresholds, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ 
       if(left){ ifelse(data[[trt]] < i, i, data[[trt]]) }
       else { ifelse(data[[trt]] > i, i, data[[trt]]) }
    }
     lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# deterministic shift lmtp
get_src_detsh <- function(shifts, Dset, sl_lib, Odat){
  maxAL0 <- max(Odat[Odat$L==0,"A"])
  maxAL1 <- max(Odat[Odat$L==1,"A"])
  minAL0 <- min(Odat[Odat$L==0,"A"])
  minAL1 <- min(Odat[Odat$L==1,"A"])
  u_func_pos <- function(L){    ifelse(L==0, maxAL0, maxAL1)}
  u_func_neg <- function(L){    ifelse(L==0, minAL0, minAL1)}
  # simulate truth
  # natural intervention
  Anat <- node("A.", distr = "rconst", const = A)
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shifted intervention
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- node("A.", distr = "rconst", const = ifelse(A <= ifelse(L==0, maxAL0, maxAL1) - i & A >= ifelse(L==0, minAL0, minAL1) - i, A + i, A))
    Dset <- Dset + action(paste("Ap", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ap", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("Ap", i, sep="")]]$Y) - mean(dat_nat[["A_"]]$Y)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shifted intervention
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ ifelse(data[[trt]] <= u_func_pos(data[["L"]]) - i & data[[trt]] >= u_func_neg(data[["L"]]) - i, data[[trt]] + i, data[[trt]]) }
     lmtp_si <- lmtp_sdr(Odat, "A", "Y", baseline = "L", shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


## longitudinal

# dose-response curve
get_drc_long <- function(doses0, doses1, Dset, Odat){
  df_drc <-data.frame()
  # simulate truth
  for (j in doses0){
    true_drc <- foreach(i = doses1, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
      Aplus1 <- c(node("A.", t=0, distr = "rconst", const = j), node("A.", t=1, distr = "rconst", const = i))
      Dset <- Dset + action(paste("A+", j, i, sep=""), nodes=Aplus1)
      dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("A+", j, i, sep="")), n = 1000000, rndseed = 1, 
                                verbose=F)
      mean(dat_int[[paste("A+", j, i, sep="")]]$Y_1)
    }
    df_drc <- rbind(df_drc, data.frame(dose0=rep(j, length(doses1)), dose1=doses1, truth=true_drc))
  }
  # estimate with g-formula
  for (j in doses0){
    estimates_drc <- foreach(i = doses1, .combine = 'c', .packages = c("CICI", "simcausal", "lmtp")) %dopar% { 
      est <- gformula(X=Odat[c("L_0", "A_0", "Y_0", "L_1", "A_1", "Y_1")], 
                      Ynodes = c("Y_0", "Y_1"), Anodes = c("A_0", "A_1"), Lnodes = c("L_1"),
                      abar=matrix(c(j,i), nrow=1), verbose=F, Yform="GAM", seed = 1)
      est$results$psi[2]
    }
    mask <- df_drc["dose0"]==j 
    df_drc[mask, "est"] <- estimates_drc
  }
  
  return(df_drc)
}


# deterministic lmtp
get_src_det_long <- function(shifts, Dset, sl_lib, Odat){
  # simulate truth
  # natural intervention
  Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- c(node("A.", t=0, distr = "rconst", const = A_0 + i), node("A.", t=1, distr = "rconst", const = A_1 + i))
    Dset <- Dset + action(paste("Ap", i, i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ap", i, i, sep="")), n = 1000000, rndseed = 1, 
                              verbose=F)
    mean(dat_int[[paste("Ap", i, i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)

  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")),  shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ print(trt);data[[trt]] + i }
     lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), c("Y_1"), time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     result
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
    
    dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
    dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
    dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])

  return(dat_shift)
}


# stochastic lmtp
get_src_sto_long <- function(shifts, Dset, sl_lib, Odat){
  # function to estimate conditional density
  dens <- function (n, a, var, shift=0) { 
    dat <- data.frame(var)
    dat["a"] <- a
    dens <- lm(a~., data = dat)
    means <- predict(dens, dat)
    draws <- rnorm(n = n, mean = means + shift, sd = summary(dens)$sigma)
    return(draws)
  }
  # simulate truth
  # natural intervention
  #Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Anat <- c(node("A.", t=0, distr = "dens", a=A_0, var=L_0, shift=0), node("A.", t=1, distr = "dens", a=A_1, var=c(L_0, A._0, Y_0, L_1) , shift=0))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal"), .export="dens") %do% {
    Aplus1 <- c(node("A.", t=0, distr = "dens", a=A_0, var=L_0, shift=i), node("A.", t=1, distr = "dens", a=A_1, var=c(L_0, A._0, Y_0, L_1) , shift=i))
    Dset <- Dset + action(paste("As", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("As", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("As", i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  Sys.time()
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
   f_si <- function(data, trt){
     dat <- data[,2:which(colnames(data)==trt)-1]
     #dat <- dat[ , grepl( "L" , names( dat ) ) ]
     dat <- data.frame(dat)
     dat["a"] <- data[[trt]]
     dens <- lm(a~., data = dat)
     means <- predict(dens, dat)
     set.seed(1)
     draws <- rnorm(n = length(dat$a), mean = means + i, sd = summary(dens)$sigma)
     draws
   }
   lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f_si, 
                       mtp = TRUE, folds = 2, outcome_type = "continuous",
                       learners_trt = sl_lib, learners_outcome = sl_lib)
   result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
   list(theta = result$theta, low = result$conf.low, high = result$conf.high)
 }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}

# dynamic deterministic lmtp
get_src_detdy_long <- function(shifts, Dset, sl_lib, Odat){
  # simulate truth
  # natural intervention
  Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- c(node("A.", t=0, distr = "rconst", const = ifelse(L_0==1, A_0 + i, A_0)), node("A.", t=1, distr = "rconst", const = ifelse(L_1==1, A_1 + i, A_1)))
    Dset <- Dset + action(paste("Adp", i, i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Adp", i, i, sep="")), n = 1000000, rndseed = 1, 
                              verbose=F)
    mean(dat_int[[paste("Adp", i, i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")),  shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ ifelse(data[[paste("L", substring(trt,2,3), sep="")]] == 1, data[[trt]] + i, data[[trt]]) }
     lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  
  return(dat_shift)
}


# dynamic stochastic lmtp
get_src_stody_long <- function(shifts, Dset, sl_lib, Odat){
  # function to estimate conditional density
  dens <- function (n, a, var, shift=0) {
    var <- data.frame(var)
    dat <- var
    dat["a"] <- a
    dens <- lm(a~., data = dat)
    means <- predict(dens, dat)
    if (is.null(ncol(var))){
      cond_l <- var
    } else {
      cond_l <- var[length(var)]
    }
    draws <- rnorm(n = n, mean = ifelse(cond_l==1, means + shift, means), sd = summary(dens)$sigma) 
    dat['a_new'] <- draws
    return(draws)
  }
  # simulate truth
  # natural intervention
  #Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Anat <- c(node("A.", t=0, distr = "dens", a=A_0, var=L_0, shift=0), node("A.", t=1, distr = "dens", a=A_1, var=c(L_0, A._0, Y_0, L_1) , shift=0))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shift interventions
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal"), .export="dens") %do% {
    Aplus1 <- c(node("A.", t=0, distr = "dens", a=A_0, var=L_0, shift=i), node("A.", t=1, distr = "dens", a=A_1, var=c(L_0, A._0, Y_0, L_1) , shift=i))
    Dset <- Dset + action(paste("Ads", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ads", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("Ads", i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shift interventions
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f_si <- function(data, trt){
       dat <- data[,2:which(colnames(data)==trt)-1]
       #dat <- dat[ , grepl( "L" , names( dat ) ) ]
       dat <- data.frame(dat)
       dat["a"] <- data[[trt]]
       dens <- lm(a~., data = dat)
       means <- predict(dens, dat)
       set.seed(1)
       draws <- rnorm(n = length(dat$a), mean = ifelse(data[[paste("L", substring(trt,2,3), sep="")]]==1, means + i, means), sd = summary(dens)$sigma)
       draws
     }
     lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f_si, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# deterministic threshold lmtp
get_src_detth_long <- function(thresholds, Dset, sl_lib, Odat, left=TRUE){
  # simulate truth
  # natural intervention
  Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shifted intervention
  responses_thresh <- foreach(i = thresholds, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    if (left){Athresh <- c(node("A.", t=0, distr = "rconst", const = ifelse(A_0<i, i, A_0)), 
                           node("A.", t=1, distr = "rconst", const = ifelse(A_1<i, i, A_1)))}
    else {Athresh <- c(node("A.", t=0, distr = "rconst", const = ifelse(A_0>i, i, A_0)),
                       node("A.", t=1, distr = "rconst", const = ifelse(A_1>i, i, A_1)))}
    Dset <- Dset + action(paste("At", i, sep=""), nodes=Athresh)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("At", i, sep="")), n = 1000000, rndseed = 1,
                              verbose = F)
    mean(dat_int[[paste("At", i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(threshold=thresholds, response=responses_thresh)
  
  # estimate shift-response curve
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shifted intervention
  estimates_shift <- foreach(i = thresholds, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ 
       if(left){ ifelse(data[[trt]] < i, i, data[[trt]]) }
       else { ifelse(data[[trt]] > i, i, data[[trt]]) }
     }
     lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}


# deterministic shift lmtp
get_src_detsh_long <- function(shifts, Dset, sl_lib, Odat){
  # simulate truth
  maxA0L0 <- max(Odat[Odat$L_0==0,"A_0"])
  minA0L0 <- min(Odat[Odat$L_0==0,"A_0"])
  maxA0L1 <- max(Odat[Odat$L_0==1,"A_0"])
  minA0L1 <- min(Odat[Odat$L_0==1,"A_0"])
  maxA1L0 <- max(Odat[Odat$L_0==0,"A_1"])
  minA1L0 <- min(Odat[Odat$L_0==0,"A_1"])
  maxA1L1 <- max(Odat[Odat$L_0==1,"A_1"])
  minA1L1 <- min(Odat[Odat$L_0==1,"A_1"])
  # natural intervention
  Anat <- c(node("A.", t=0, distr = "rconst", const = A_0), node("A.", t=1, distr = "rconst", const = A_1))
  Dset <- Dset  + action("A_", nodes=Anat)
  dat_nat <- simcausal::sim(DAG = Dset, actions = "A_", n = 1000000, rndseed = 1, 
                            verbose=F)
  # shifted intervention
  responses_shift <- foreach(i = shifts, .combine = 'c', .packages = c("lmtp", "simcausal")) %do% {
    Aplus1 <- c(node("A.", t=0, distr = "rconst", const = ifelse(A_0 <= ifelse(L_0==0,maxA0L0,maxA0L1) - i & A_0 >= ifelse(L_0==0,minA0L0,minA0L1) - i, A_0 + i, A_0)),
                node("A.", t=1, distr = "rconst", const = ifelse(A_1 <= ifelse(L_1==0,maxA1L0,maxA1L1) - i & A_1 >= ifelse(L_1==0,minA1L0,minA1L1) - i, A_1 + i, A_1)) )
    Dset <- Dset + action(paste("Ap", i, sep=""), nodes=Aplus1)
    dat_int <- simcausal::sim(DAG = Dset, actions = c(paste("Ap", i, sep="")), n = 1000000, rndseed = 1, 
                              verbose = F)
    mean(dat_int[[paste("Ap", i, sep="")]]$Y_1) - mean(dat_nat[["A_"]]$Y_1)
  }
  dat_shift <- data.frame(shift=shifts, response=responses_shift)
  # estimate shift-response curve
  u_func_pos <- function(L, A){    ifelse(L==0,max(Odat[Odat$L_1==0,A]),max(Odat[Odat$L_1==1,A]))}
  u_func_neg <- function(L, A){    ifelse(L==0,min(Odat[Odat$L_1==0,A]),min(Odat[Odat$L_1==1,A]))}
  # natural intervention
  lmtp_natural <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = NULL, 
                           mtp = TRUE, folds = 2, outcome_type = "continuous",
                           learners_trt = sl_lib, learners_outcome = sl_lib)
  # shifted intervention
  estimates_shift <- foreach(i = shifts, .combine = 'comb', .multicombine=T, 
                             .init=list(list(),list(),list()), .packages = c("lmtp")) %dopar% {
     f <- function(data, trt){ 
       ifelse(data[[trt]] <= u_func_pos(data[[paste("L", substring(trt,2,3), sep="")]], paste("A", substring(trt,2,3), sep="")) - i & 
              data[[trt]] >= u_func_neg(data[[paste("L", substring(trt,2,3), sep="")]], paste("A", substring(trt,2,3), sep="")) - i, 
              data[[trt]] + i, data[[trt]]) }
     lmtp_si <- lmtp_sdr(Odat, c("A_0", "A_1"), "Y_1", time_vary = list(c("L_0"), c("L_1","Y_0")), shift = f, 
                         mtp = TRUE, folds = 2, outcome_type = "continuous",
                         learners_trt = sl_lib, learners_outcome = sl_lib)
     result <- lmtp_contrast(lmtp_si, ref=lmtp_natural)$vals
     list(theta = result$theta, low = result$conf.low, high = result$conf.high)
   }
  dat_shift["lmtp"] <- unlist(estimates_shift[[1]])
  dat_shift["lmtp_low"] <- unlist(estimates_shift[[2]])
  dat_shift["lmtp_high"] <- unlist(estimates_shift[[3]])
  return(dat_shift)
}



























## other helpful functions

# truncate A to make sure all values lie on the graph 
# (where the real dose-response curve is not extreme)
rnorm_trunc <- function(n, mean, sd, minval = 0, maxval=1) {
  out <- rnorm(n = n, mean = mean, sd = sd) 
  minval <- minval 
  out[out < minval] <- minval 
  out[out > maxval] <- maxval 
  return(out) 
}

# combine outputs in dopar
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
