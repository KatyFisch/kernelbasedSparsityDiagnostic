library(simcausal)
library(lmtp)
library(SuperLearner)
library(ggplot2)
library(CICI)
library(patchwork)
library(cowplot)
#install_github("ehkennedy/npcausal")
library(npcausal)

source("curvecalc_functions.R")
source("plot_functions.R")
source("pdf_functions_wdiagnostic.R")
source("diagnostic_function.R")

#parallelization
library(foreach)
parallel::detectCores() #12
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# general settings
sl_lib <- c("SL.glm", "SL.gam", "SL.glmnet", "SL.earth")
length_out <- 50 #50
length_out_dia <- 9 #9
n_datapoints <- 350 #350

# polynomial linear
y_func <- function(A.){ 
  result <- 1.1 - 23.29167*A. + 154.0625*A.**2 - 392.1875*A.**3 + 429.6875*A.**4 - 169.2708*A.**5
  if (A. <= 0.11){
    point <- 0.11
    a <-  1.1 - 23.29167*point + 154.0625*point**2 - 392.1875*point**3 + 429.6875*point**4 - 169.2708*point**5
    deriv <- -23.29167 + 2*154.0625*point - 3*392.1875*point**2 + 4*429.6875*point**3 - 5*169.2708*point**4
    result <- deriv*A. + a - deriv*point
  }
  if (A. >=0.9){
    point <- 0.9
    a <-  1.1 - 23.29167*point + 154.0625*point**2 - 392.1875*point**3 + 429.6875*point**4 - 169.2708*point**5
    deriv <- -23.29167 + 2*154.0625*point - 3*392.1875*point**2 + 4*429.6875*point**3 - 5*169.2708*point**4
    result <- deriv*A. + a - deriv*point
  }
  return(result)
}
D.polynomial_linear <- DAG.empty() + 
  node("L", distr = "rbern", prob = 0.3) +
  node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
  node("A.",  distr = "rconst", const = A) +
  node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)

make_pdfs_stp(sl_lib = sl_lib,
              D.base = D.polynomial_linear,
              doses = seq(0, 1, length.out = length_out),
              kennedy_search = seq(.015, 0.1, length.out = length_out),
              shifts_det_lpv = seq(-0.3, 0.3, length.out = length_out), 
              shifts_sto_lpv = seq(-0.3, 0.3, length.out = length_out), 
              shifts_det_d_lpv = seq(-0.3, 0.3, length.out = length_out),  
              shifts_sto_d_lpv = seq(-0.3, 0.3, length.out = length_out),  
              thresholds_lpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_lpv = seq(-0.8, 0.8, length.out = length_out),
              shifts_det_mpv = seq(-0.3, 0.3, length.out = length_out),
              shifts_sto_mpv = seq(-0.3, 0.3, length.out = length_out), 
              shifts_det_d_mpv = seq(-0.3, 0.3, length.out = length_out),  
              shifts_sto_d_mpv = seq(-0.3, 0.3, length.out = length_out),  
              thresholds_mpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_mpv = seq(-0.5, 0.5, length.out = length_out), 
              shifts_det_spv = seq(-0.3, 0.3, length.out = length_out), 
              shifts_sto_spv = seq(-0.3, 0.3, length.out = length_out), 
              shifts_det_d_spv = seq(-0.3, 0.3, length.out = length_out),  
              shifts_sto_d_spv = seq(-0.3, 0.3, length.out = length_out),  
              thresholds_spv = seq(0, 1, length.out = length_out),  
              shifts_detsh_spv = seq(-0.4, 0.4, length.out = length_out),  
              descr_str = "polynomial-linear", 
              length_out_dia = length_out_dia, 
              n = n_datapoints
)

# polynomial
D.polynomial <- DAG.empty() + 
  node("L", distr = "rbern", prob = 0.3) +
  node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) +
  node("A.",  distr = "rconst", const = A) +
  node("Y",  distr = "rnorm", sd=0.05,
       mean = 0.2*L + 1.1 - 23.29167*A. + 154.0625*A.**2 - 392.1875*A.**3 + 429.6875*A.**4 - 169.2708*A.**5)
make_pdfs_stp(sl_lib = sl_lib,
              D.base = D.polynomial,
              doses = seq(0, 1, length.out = length_out),
              kennedy_search = seq(.015, 0.1, length.out = length_out), 
              shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_lpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
              shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
              shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_mpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
              shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_spv = seq(0, 1, length.out = length_out),  
              shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
              descr_str = "polynomial", 
              length_out_dia = length_out_dia, 
              n = n_datapoints)


# linear piece 2
y_func <- function(A.){ 
  result <- 0.4 + 3* A.
  if (A. >= 0.1){
    result <- result - 3.5*(A.-0.1)
  }
  if (A. >= 0.2){
    result <- result + 2*(A.-0.2)
  }
  if (A. >= 0.3){
    result <- result - 3*(A.-0.3)
  }
  if (A. >= 0.4){
    result <- result - 2.5*(A.-0.4)
  }
  if (A. >= 0.5){
    result <- result + 6*(A.-0.5)
  }
  if (A. >= 0.6){
    result <- result - 1*(A.-0.6)
  }
  if (A. >= 0.7){
    result <- result - 1.8*(A.-0.7)
  }
  if (A. >= 0.8){
    result <- result + 2.5*(A.-0.8)
  }
  if (A. >= 0.9){
    result <- result - 2.5*(A.-0.9)
  }
  return(result)
}
D.linear_piece2 <- DAG.empty() + 
  node("L", distr = "rbern", prob = 0.3) +
  node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
  node("A.",  distr = "rconst", const = A) +
  node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
make_pdfs_stp(sl_lib = sl_lib,
              D.base = D.linear_piece2,
              doses = seq(0, 1, length.out = length_out),
              kennedy_search = seq(.015, 0.1, length.out = length_out), 
              shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_lpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
              shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
              shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_mpv = seq(0, 1, length.out = length_out), 
              shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
              shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
              shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
              shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
              thresholds_spv = seq(0, 1, length.out = length_out),  
              shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
              descr_str = "linear-piece2", 
              length_out_dia = length_out_dia, 
              n = n_datapoints)




if(TRUE){
  # linear
  D.linear <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) +
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", sd=0.05, mean = 0.2*L - 0.3 + 1.4*A.)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.linear,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out),
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "linear", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)
  
  # linear piece 1
  y_func <- function(A.){ 
    result <- 0
    if (A. <= 0.2){
      result <- result + 0.5*A.
    }
    if (A. <= 0.4){
      result <- result - 2*A.
    }
    if (A. <= 0.6){
      result <- result + 3*A.
    }
    if (A. <= 0.8){
      result <- result - 5*A.
    }
    {
      result <- result + 2*A.
    }
    return(result)
  }
  D.linear_piece1 <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.linear_piece1,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out), 
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "linear-piece1", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)
  

  

  
  # polynomial flat
  y_func <- function(A.){ 
    A. <- A.*1.1
    result <- 1.1 - 23.29167*A. + 154.0625*A.**2 - 392.1875*A.**3 + 429.6875*A.**4 - 169.2708*A.**5
    if (A. <= 0.2){
      if (A. <=0){result <- 0}
      else {
        result <- (10*A.)**3*A.*A./0.7/2/2
      }
    }
    #if (A. >=0.91 && A. <0.94){
    #  result <- 0.997 - 0.02*(10*A.)**2 + 8.1*0.2
    #}
    if (A. >=0.94){
      if (A. >=1){result <- 0.4}
      else{
        result <- 0.4 + (10*(1-A.))**3/0.85
      }
    }
    return(result)
  }
  D.polynomial_flat <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.polynomial_flat,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out), 
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "polynomial-flat", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)
  
  # sine
  y_func <- function(A.){ 
    A. <- A.*20
    result <- sin(A.)*0.5+0.5
    return(result)
  }
  D.sine <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.sine,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out), 
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "sine", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)
  
  # step 1
  y_func <- function(A.){ 
    if (A. <= 0.2){
      result <- 0.5
    }
    else if (A. <= 0.4){
      result <- 0.8
    }
    else if (A. <= 0.6){
      result <- 0.2
    }
    else if (A. <= 0.8){
      result <- 0.6
    }
    else {
      result <- 0
    }
    return(result)
  }
  D.step1 <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.step1,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out), 
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "step1", 
                length_out_dia = length_out_dia, 
                n = n_datapoints) 
  
  # step2
  y_func <- function(A.){ 
    if (A. <= 0.2){
      result <- 0.5
    }
    else if (A. <= 0.22){
      result <- 1.1
    }
    else if (A. <= 0.4){
      result <- 0.6
    }
    else if (A. <= 0.45){
      result <- 0.2
    }
    else if (A. <= 0.6){
      result <- -0.1
    }
    else if (A. <= 0.7){
      result <- 0.6
    }
    else {
      result <- 0
    }
    return(result)
  }
  D.step2 <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.step2,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out), 
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "step2", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)
  

  
  # step3
  y_func <- function(A.){ 
    if (A. <= 0.2){
      result <- 0
    }
    else if (A. <= 0.4){
      result <- 0.25
    }
    else if (A. <= 0.6){
      result <- 0.5
    }
    else if (A. <= 0.8){
      result <- 0.75
    }
    else {
      result <- 1
    }
    return(result)
  }
  D.step3 <- DAG.empty() + 
    node("L", distr = "rbern", prob = 0.3) +
    node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14) + 
    node("A.",  distr = "rconst", const = A) +
    node("Y",  distr = "rnorm", mean = 0.2*L+y_func(A.), sd=0.05)
  make_pdfs_stp(sl_lib = sl_lib,
                D.base = D.step3,
                doses = seq(0, 1, length.out = length_out),
                kennedy_search = seq(.015, 0.1, length.out = length_out),
                shifts_det_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_lpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_lpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_lpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_lpv = seq(-1, 1, length.out = length_out),
                shifts_det_mpv = seq(-0.2, 0.2, length.out = length_out),
                shifts_sto_mpv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_mpv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_mpv = seq(0, 1, length.out = length_out), 
                shifts_detsh_mpv = seq(-1, 1, length.out = length_out), 
                shifts_det_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_sto_spv = seq(-0.2, 0.2, length.out = length_out), 
                shifts_det_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                shifts_sto_d_spv = seq(-0.3, 0.15, length.out = length_out),  
                thresholds_spv = seq(0, 1, length.out = length_out),  
                shifts_detsh_spv = seq(-1, 1, length.out = length_out),  
                descr_str = "step3", 
                length_out_dia = length_out_dia, 
                n = n_datapoints)

}
