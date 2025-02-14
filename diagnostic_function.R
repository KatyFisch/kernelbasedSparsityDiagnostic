library(Rfast)

library(simcausal)
library(lmtp)
library(SuperLearner)
library(randomForest)
library(ggplot2)
library(ggtext)
library(ltmle)
library(CICI)
library(rJava)
library(bartMachine)
library(patchwork)
#source("curvecalc_functions.R")
#source("plot_functions.R")
#source("pdf_functions.R")
#install_github("ehkennedy/npcausal")
library(npcausal)
library(cowplot)
library(ks)

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






## CREATE DAG
# D <- DAG.empty()
# # D.base <- D + node("L", t=0, distr = "rbern", prob = 0.3)
# # D.base <- D.base + node("A", t=0,  distr = "rnorm", mean = (0.21)*L[0] + 0.4, sd=0.13)
# # D.base <- D.base +
# #   node("Y", t=0,  distr = "rnorm", sd=0.05,
# #        mean = 0.2*L[0] + 1.1 - 23.29167*A[0] + 154.0625*A[0]**2 - 392.1875*A[0]**3 + 429.6875*A[0]**4 - 169.2708*A[0]**5)
# # #t=1
# # D.base <- D.base + node("L", t=1, distr = "rbern", prob = ifelse(A[0]>0.5, 0.5, 0.3))
# # D.base <- D.base + node("A", t=1,  distr = "rnorm", mean = (0.1)*L[1] + 0.3*A[0] + 0.3, sd=0.13)
# # D.base <- D.base +
# #   node("Y", t=1,  distr = "rnorm", sd=0.05,
# #        mean = 0.2*L[1] + Y[0] + 1.1 - 23.29167*A[1] + 154.0625*A[1]**2 - 392.1875*A[1]**3 + 429.6875*A[1]**4 - 169.2708*A[1]**5)
# D.base <- D + node("L", distr = "rbern", prob = 0.3)
# D.base <- D.base + node("A", distr = "rnorm", mean = (0.21)*L + 0.4, sd=0.13)
# D.base <- D.base +
#   node("Y", distr = "rnorm", sd=0.05,
#        mean = 0.2*L + 1.1 - 23.29167*A + 154.0625*A**2 - 392.1875*A**3 + 429.6875*A**4 - 169.2708*A**5)
# Dset <- simcausal::set.DAG(D.base)
# # get data for plotting
# Odat <- simcausal::sim(DAG = Dset, n = 100, rndseed = 1, verbose=F)
# data <- subset(Odat, select = -c(Y, ID) )




# # example intervention function
# int_func <- function(data, param){
#   data_int <- data
#   data_int$A_0 <- data_int$A_0 + param
#   data_int$A_1 <- data_int$A_1 + param
#   return(data_int)
# }

# kernel
rbf.normal <- function(x, point, disthalf) {
  gamma <- -log(0.5)/(disthalf)^2
  exp(-gamma * (x-point)^2 )  
}

# calc diagnostic
singlediagnostic <- function(data_observed, obs_intervened, disthalf_vec, kernel=rbf.normal){
  foreach(i=1:nrow(data_observed), .combine='+') %:%
    foreach(coln=colnames(data_observed), .combine='*') %do% {
      if (coln %in% names(disthalf_vec)){
        res <- kernel(data_observed[i, coln], obs_intervened[[coln]], disthalf=unname(disthalf_vec[coln]))
      }
      else{1}
    }
}



# calc diagnostic harmonic mean
singlediagnostic_hm <- function(data_observed, obs_intervened, disthalf_vec, kernel=rbf.normal){
  hm_fin_func <- function(val){
    n <- length(disthalf_vec)
    return(n/val)
  }
  foreach(i=1:nrow(data_observed), .combine='+') %:%
    foreach(coln=colnames(data_observed), .combine='+', .final=hm_fin_func) %do% {
      if (coln %in% names(disthalf_vec)){
        res <- 1/kernel(data_observed[i, coln], obs_intervened[[coln]], disthalf=unname(disthalf_vec[coln]))
      }
      else{1}
    }
}

# calc diagnostic with minimum values
singlediagnostic_mv <- function(data_observed, obs_intervened, disthalf_vec, minval_vec, kernel=rbf.normal){
  foreach(i=1:nrow(data_observed), .combine='+') %:%
    foreach(coln=colnames(data_observed), .combine='*') %do% {
      if (coln %in% names(disthalf_vec)){
        res <- kernel(data_observed[i, coln], obs_intervened[[coln]], disthalf=unname(disthalf_vec[coln]))
        minval <- minval_vec[coln]
        #scale
        res <- res*(1-minval)+minval
      }
      else{1}
    }
}





#########

f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
o <- function(x) {
  subset(x, x < quantile(x, 0.05) | quantile(x, 0.95) < x)
}

# the diagnostic
my_diagnostic <- function(int_func, int_func_param, data, disthalf_vec, col="darkgreen", type="multiplicative",
                          minval_vec=NA){
  
  original_data <- data
  #data <- subset(data, select = -c(A.,ID) )
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]
  data[is.na(data)] <- Inf #Na handling might not be optimal?

  # actual function
  result_df <- data.frame(matrix(NA, nrow=0, ncol=3))
  names(result_df) <- c("shift", "diagnostic", "observation")
  for (param in int_func_param){
    # data
    data_intervened <- int_func(original_data, param)
    data_intervened <- data_intervened[vartoconsider]
    

    if (type=="minval"){
      dia <- apply(data_intervened, 1, function(x){singlediagnostic_mv(data, x, disthalf_vec, 
                                                                       minval_vec, kernel=rbf.normal)})
    }
    else if (type=="harmonicmean"){
      dia <- apply(data_intervened, 1, function(x){singlediagnostic_hm(data, x, disthalf_vec, kernel=rbf.normal)})
    }
    else{
      dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
      dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
      dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    }


    #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec, kernel=rbf.normal)})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )
    result_df <- rbind(result_df, new_data)
  }
  
  result_plot <- ggplot(result_df, color=col, aes(x=shift, y=diagnostic, group=shift)) + 
    stat_summary(fun.data = f, geom="boxplot") + xlab("") + ylab("EDP") +
    stat_summary(fun = o, geom="point")
  

  return(result_plot)
}

my_plot <- function(diagnostic_object){

}


# the diagnostic
my_diagnostic_dat <- function(data_int, data, disthalf_vec, col="darkgreen", type="multiplicative",
                          minval_vec=NA){
  
  
  #data <- subset(data, select = -c(A.,ID) )
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]
  #data[is.na(data)] <- Inf #Na handling might not be optimal?
  data_int <- data_int[vartoconsider]
  
  # actual function
  result_df <- data.frame(matrix(NA, nrow=0, ncol=3))
  names(result_df) <- c("shift", "diagnostic", "observation")
  # data
  data_intervened <- data_int
  
  
  if (type=="minval"){
    dia <- apply(data_intervened, 1, function(x){singlediagnostic_mv(data, x, disthalf_vec, 
                                                                     minval_vec, kernel=rbf.normal)})
  }
  else if (type=="harmonicmean"){
    dia <- apply(data_intervened, 1, function(x){singlediagnostic_hm(data, x, disthalf_vec, kernel=rbf.normal)})
  }
  else{
    dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
    dat_int_scaled[sapply(dat_int_scaled, is.infinite)] <- 1 #what's happening here?
    dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
    dat_scaled[sapply(dat_scaled, is.infinite)] <- 1
    dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
  } 
    
  #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec, kernel=rbf.normal)})
  result_df <- data.frame(shift= rep(1, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )

  result_plot <- ggplot(result_df, color=col, aes(x=shift, y=diagnostic)) + 
    stat_summary(fun.data = f, geom="boxplot") + xlab("") + 
    stat_summary(fun = o, geom="point")
  return(result_plot)
}
  

my_diagnostic_values <- function(int_func, int_func_param, data, disthalf_vec, col="darkgreen", type="multiplicative",
                          minval_vec=NA){
  
  
  #data <- subset(data, select = -c(A.,ID) )
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]
  data[is.na(data)] <- Inf #Na handling might not be optimal?
  
  # actual function
  result_df <- data.frame(matrix(NA, nrow=0, ncol=3))
  names(result_df) <- c("shift", "diagnostic", "observation")
  for (param in int_func_param){
    # data
    data_intervened <- int_func(data, param)
    
    
    if (type=="minval"){
      dia <- apply(data_intervened, 1, function(x){singlediagnostic_mv(data, x, disthalf_vec, 
                                                                       minval_vec, kernel=rbf.normal)})
    }
    else if (type=="harmonicmean"){
      dia <- apply(data_intervened, 1, function(x){singlediagnostic_hm(data, x, disthalf_vec, kernel=rbf.normal)})
    }
    else{
      dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
      dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
      dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    }
    
    
    #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec, kernel=rbf.normal)})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )
    result_df <- rbind(result_df, new_data)
  }
  return(result_df)
}


# the extended diagnostic
my_extended_diagnostic <- function(int_func, int_func_param, data, disthalf_vec_Q, disthalf_vec_g, disthalf_vec_w, col="darkgreen"){
  
  #data <- subset(data, select = -c(A.,ID) )
  
  
  # actual function
  result_df <- data.frame(matrix(NA, nrow=0, ncol=4))
  names(result_df) <- c("shift", "diagnostic", "observation", "type")
  data[is.na(data)] <- Inf #Na handling might not be optimal?
  # go over different
  for (param in int_func_param){
    # data
    data_intervened <- int_func(data, param)
    
    # Q
    dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec_Q, FUN='*')
    dat_scaled <- sweep(data, 2, 1/disthalf_vec_Q, FUN='*')
    dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec_Q, kernel=rbf.normal)})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia), 
                           type= rep("Q", length(dia)))
    result_df <- rbind(result_df, new_data)
    # g
    dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec_g, FUN='*')
    dat_scaled <- sweep(data, 2, 1/disthalf_vec_g, FUN='*')
    dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec_g, kernel=rbf.normal)})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia), 
                           type= rep("g", length(dia)))
    result_df <- rbind(result_df, new_data)
    # weights
    dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec_w, FUN='*')
    dat_scaled <- sweep(data, 2, 1/disthalf_vec_w, FUN='*')
    dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    #dia <- apply(data_intervened, 1, function(x){min(200,1/singlediagnostic(data, x, disthalf_vec_w, kernel=rbf.normal))})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia), 
                           type= rep("<i>w</i>", length(dia)))
    result_df <- rbind(result_df, new_data)
  }
  
  result_df$type <- factor(result_df$type, levels=c("Q", "g", paste0("<i>w</i>")))
  
  scaling_factor <- max(dat[dat$type!="<i>w</i>", "diagnostic"])/max(dat[dat$type=="<i>w</i>", "diagnostic"])
  
  
  result_plot <- ggplot(result_df, 
                   aes(x=type,
                       y = ifelse(type == "<i>w</i>", diagnostic*scaling_factor, diagnostic), 
                       color = ifelse(type == "<i>w</i>", "violet", "black"),
                       fill=shift, group=type)) + 
    stat_summary(fun.data = f, geom="boxplot", fill="white") +
    stat_summary(fun = o, geom="point") + facet_wrap(~shift, nrow=1) +  
    scale_y_continuous(sec.axis = sec_axis(transform = ~./scaling_factor,
                                           breaks = (seq(0,200,50)))) +
    scale_color_identity() +
    guides(fill="none") +
    labs(y = "y") +
    theme(axis.text.y.right = element_text(color = "violet")) + theme(axis.text.x = element_markdown())
  
  #result_plot <- ggplot(result_df, aes(x=type, y=diagnostic, fill=shift, group=type)) + 
  #  stat_summary(fun.data = f, geom="boxplot", fill="white") + xlab("") + #scale_y_log10() +
  #  stat_summary(fun = o, geom="point") + facet_wrap(~shift, nrow=1) + theme(panel.spacing = unit(0.2, "lines")) + 
  #  guides(fill="none") + theme(axis.text.x = element_markdown()) #to do make third boxplots different color, different axis
  return(result_plot)
}




# # Install and load the ks package
# 
# 
# 
# 
# 
# # example intervention function
# int_func <- function(data, param){
#   data_int <- data
#   data_int$A <- data_int$A + param
#   return(data_int)
# }
# 
# 
# disthalf_vec_Q <- c(L=0.5, A=0.1)
# 
# disthalf_vec_g <- c(L=0.5)
# 
# disthalf_vec_w <- c(A=0.01)
#   
# int_func_param <- 0.1#seq(0,1,0.1)
# 
# #int_func(data, int_func_param)
# 
# my_extended_diagnostic(int_func, int_func_param, data, disthalf_vec_Q, disthalf_vec_g, disthalf_vec_w)
# 
# 
# # library(ggplot2)
# # ggplot(result_df) + geom_boxplot(aes(x=shift, y=diagnostic, group=shift))
# # #plot(data_intervened$A_1, dia)
# # data[1,]
# 
# a1 <- rnorm(n=100, mean=5, sd=20)
# a2 <- rnorm(n=100, mean=6, sd=20)
# a3 <- rnorm(n=100, mean=8, sd=20)
# 
# b <- rnorm(n=100, mean=11, sd=10)
# 
# c1 <- rnorm(n=100, mean=500, sd=80)
# c2 <- rnorm(n=100, mean=600, sd=80)
# c3 <- rnorm(n=100, mean=800, sd=80)
# 
# letters <- c(rep("allgemein", 300))
# type <- rep(c(rep(1,100), rep(2,100), rep(3,100)),3)
# 
# dat <- data.frame(y=c(a1,a2,a3,b,b,b,c1,c2,c3), letters, type)
# 
# 
# f <- function(x) {
#   r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   r
# }
# o <- function(x) {
#   subset(x, x < quantile(x, 0.05) | quantile(x, 0.95) < x)
# }
# ggplot(dat, aes(x=letters, y=y)) + 
#   stat_summary(fun.data = f, geom="boxplot", fill="white") +
#   stat_summary(fun = o, geom="point") +  
#   guides(fill="none")
# 
# 
# scaling_factor <- max(dat[dat$type!="<i>w</i>", "diagnostic"])/max(dat[dat$type=="<i>w</i>", "diagnostic"])
# 
# 
# ggplot(dat, 
#        aes(x=type,
#            y = ifelse(type == "<i>w</i>", diagnostic*scaling_factor, diagnostic), 
#            color = ifelse(type == "<i>w</i>", "violet", "black"),
#            fill=shift, group=type)) + 
#   stat_summary(fun.data = f, geom="boxplot", fill="white") +
#   stat_summary(fun = o, geom="point") + facet_wrap(~shift, nrow=1) +  
#   scale_y_continuous(sec.axis = sec_axis(transform = ~./scaling_factor,
#                                          breaks = (seq(0,200,50)))) +
#   scale_color_identity() +
#   guides(fill="none") +
#   labs(y = "y") +
#   theme(axis.text.y.right = element_text(color = "violet")) + theme(axis.text.x = element_markdown())
# 
# 
# 

# df1 <- data.frame(a=c(1,1,1,1,1,1,1,1), b=c(1,2,3,1,2,3,1,2), c=c(0,0,1,0,0,1,0,0))
# df2 <- data.frame(a=c(1,2,1,2,1,2,1,2), b=c(1,2,3,1,2,3,1,2), c=c(0,0,1,0,0,1,0,0))
# 
# n <- nrow(df1)
# p <- ncol(df1)
# 
# rbf.normal <- function(point, x) {
#   exp(log(0.5) * (x-point)^2 )  
# }
# 
# result <-c()
# 
# for (i in 1:n){
#   row1 <- df1[i,]
#   # sum over all df2
#   result_i <- 0
#   for (j in 1:n){
#     row2 <- df2[j,]
#     # calculate the kernel_function for each column and multiply the results
#     result_j <- 1
#     for (d in 1:p){
#       result_kernel <- rbf.normal(row1[[d]], row2[[d]])
#       result_j <- result_j*result_kernel
#     }
#     result_i <- result_i + result_j
#   }
#   result <- c(result, result_i)
# }
# 
# result
# 
# 
# rbf.normal = \(x, y) exp(-.6931472 *(x-y)^2L); rowSums(combn(cbind(df1, df2), 2L, FUN=Reduce, f=rbf.normal))
# 
# 
# 
# 



























