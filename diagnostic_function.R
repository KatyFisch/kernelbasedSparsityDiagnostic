library(Rfast)
library(ggtext)
library(ks)



# kernel
rbf.normal <- function(x, point, disthalf) {
  gamma <- -log(0.5)/(disthalf)^2
  if (gamma == Inf){
    result <- ifelse(x==point, 1, 0)
  } else {
    result <- exp(-gamma * (x-point)^2 ) 
  }
  result
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

# helper functions for adjusted boxplot
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
o <- function(x) {
  subset(x, x < quantile(x, 0.05) | quantile(x, 0.95) < x)
}

# the diagnostic
my_diagnostic <- function(int_func, int_func_param, data, disthalf_vec, col="darkgreen", type="Rfast",
                          minval_vec=NA, kernel="rbf.normal", plot.out=TRUE){
  
  if (kernel=="rbf.normal"){
    kernel_used <- rbf.normal
  } else {
    kernel_used <- kernel
  }
  
  original_data <- data
  #data <- subset(data, select = -c(A.,ID) )
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]

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
    else if (type=="Rfast") {
      dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
      dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
      dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)))
    } else {
      dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec,  
                                                                    kernel=kernel_used)})
    }


    #dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec, kernel=rbf.normal)})
    new_data <- data.frame(shift= rep(param, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )
    result_df <- rbind(result_df, new_data)
  }
  
  if (!plot.out){
    return(result_df)
  }
  
  result_plot <- ggplot(result_df, color=col, aes(x=shift, y=diagnostic, group=shift)) + 
    stat_summary(fun.data = f, geom="boxplot") + xlab("") + ylab("EDP") +
    stat_summary(fun = o, geom="point")
  

  return(result_plot)
}



