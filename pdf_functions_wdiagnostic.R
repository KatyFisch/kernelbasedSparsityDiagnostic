make_pdfs_stp <- function(sl_lib, D.base, doses, kennedy_search, shifts_det_lpv, shifts_sto_lpv, 
                          shifts_det_d_lpv, shifts_sto_d_lpv, thresholds_lpv, 
                          shifts_detsh_lpv, shifts_det_mpv, shifts_sto_mpv, 
                          shifts_det_d_mpv, shifts_sto_d_mpv, thresholds_mpv, 
                          shifts_detsh_mpv, shifts_det_spv, shifts_sto_spv, 
                          shifts_det_d_spv, shifts_sto_d_spv, thresholds_spv, 
                          shifts_detsh_spv, descr_str, length_out_dia, n, 
                          disthalf_vec=c(A = 0.15)
                          ){
  ####################################################
  # little positivity violation
  ####################################################
  
  # function that changes the ranges for the diagnostic
  dia_range <- function(sequence, n=length_out_dia){ 
    seq(min(sequence), max(sequence), length.out=n)
  }
  length_out_dia_cdrc <- length_out_dia*2
  
  
  ## CREATE DAG
  Dset <- simcausal::set.DAG(D.base)
  # get data for plotting
  Odat <- simcausal::sim(DAG = Dset, n = n, rndseed = 1, verbose=F) 
  
  ## dose-response curve with CICI
  dat_drc <- get_drc(doses, Dset, Odat)
  ## dose-response curve with Kennedy
  l <- matrix(c(Odat$L, rep(0,n)), nrow=n, ncol=2)
  ce.res <- ctseff(y=Odat$Y, a=Odat$A, x=l, bw.seq=kennedy_search, sl.lib=sl_lib)
  # check that bandwidth choice is minimizer
  plot(ce.res$bw.risk$bw, ce.res$bw.risk$risk)
  # plot drc
  gg_drc_lpv <- plot_drc(dat_drc, ce.res, Odat)
  gg_drc_lpv_l <- plot_drc_legend(dat_drc, ce.res, Odat)
  # add diagnostic
  drc_intervention <- function(data, param){
    data["A"] <- param
    return(data)
  }
  disthalf_vec <- c(A = 0.5*sd(Odat$A))
  gg_drc_lpv_dia <- my_diagnostic(int_func=drc_intervention, int_func_param=dia_range(doses, length_out_dia_cdrc),
                                  data=Odat, disthalf_vec=disthalf_vec, col="black")
  gg_drc_lpv_noest <- plot_drc_no_estimation(dat_drc, Odat)
  gg_drc_lpv_noest_l <- plot_drc_no_estimation_legend(dat_drc, Odat)



  # deterministic shift shifts, Dset, sl_lib, Odat
  dat_src_det <- get_src_det(shifts=shifts_det_lpv, Dset=Dset, sl_lib=sl_lib, Odat=Odat)
  # stochastic shift
  dat_src_sto <- get_src_sto(shifts_sto_lpv, Dset, sl_lib, Odat)
  # plot src
  gg_src_lpv <- plot_src(dat_src_det, dat_src_sto)
  # add diagnostic
  src_det_intervention <- function(data, param){
    data["A"] <- data["A"] + param
    return(data)
  }
  gg_src_lpv_dia <- my_diagnostic(int_func=src_det_intervention, int_func_param=dia_range(shifts_sto_lpv),
                                  data=Odat, disthalf_vec=disthalf_vec)

  
  # deterministic dynamic shift
  dat_src_detdy <- get_src_detdy(shifts_det_d_lpv, Dset, sl_lib, Odat)
  # stochastic dynamic shift
  dat_src_stody <- get_src_stody(shifts_sto_d_lpv, Dset, sl_lib, Odat)
  # plot src dynamic
  gg_src_d_lpv <- plot_src_d(dat_src_stody, dat_src_detdy)
  # add diagnostic
  src_detdy_intervention <- function(data, param){
    data["A"] <- data["A"] + ifelse(data["L"]==1, param, 0)
    return(data)
  }
  gg_src_d_lpv_dia <- my_diagnostic(int_func=src_detdy_intervention, int_func_param=dia_range(shifts_sto_d_lpv),
                                    data=Odat, disthalf_vec=disthalf_vec)


  ## deterministic "threshold" lmtp
  dat_src_detth <- get_src_detth(thresholds_lpv, Dset, sl_lib, Odat, left=TRUE)
  # plot "threshold"
  gg_src_th_lpv <- plot_src_th(dat_src_detth)
  # add diagnostic
  src_detth_intervention <- function(data, param){
    data["A"] <- ifelse(data[["A"]] < param, param, data[["A"]])
    return(data)
  }
  gg_src_th_lpv_dia <- my_diagnostic(int_func=src_detth_intervention, int_func_param=dia_range(thresholds_lpv),
                                     data=Odat, disthalf_vec=disthalf_vec)


  # deterministic "shift" lmtp
  dat_src_detsh <- get_src_detsh(shifts_detsh_lpv, Dset, sl_lib, Odat)
  # plot "shift"
  gg_src_sh_lpv <- plot_src_sh(dat_src_detsh)
  # add diagnostic
  src_detsh_intervention <- function(data, param){
    maxAL0 <- max(data[data$L==0,"A"])
    maxAL1 <- max(data[data$L==1,"A"])
    minAL0 <- min(data[data$L==0,"A"])
    minAL1 <- min(data[data$L==1,"A"])
    u_func_pos <- function(L){    ifelse(L==0, maxAL0, maxAL1)}
    u_func_neg <- function(L){    ifelse(L==0, minAL0, minAL1)}
    data["A"] <- ifelse(data[["A"]] <= u_func_pos(data[["L"]]) - param & data[["A"]] >= u_func_neg(data[["L"]]) - param, data[["A"]] + param, data[["A"]])
    return(data)
  }
  gg_src_sh_lpv_dia <- my_diagnostic(int_func=src_detsh_intervention, int_func_param=dia_range(shifts_detsh_lpv),
                                     data=Odat, disthalf_vec=disthalf_vec)

  ### save final plot with diagnostic
  layout <- "
  AAC
  BBD
  EGI
  FHJ
  "
  gg_full <- wrap_elements(full = gg_drc_lpv) + wrap_elements(full = gg_drc_lpv_dia) +
    add_double_legend(gg_src_lpv, dat_src_sto) + gg_src_lpv_dia +
    gg_src_d_lpv + gg_src_d_lpv_dia +
    gg_src_th_lpv + gg_src_th_lpv_dia +
    gg_src_sh_lpv  + gg_src_sh_lpv_dia +
    plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_dia_little.pdf', sep=''), gg_full, device="pdf", width = 12, height = 10)
  dev.off()

  ### save final plot without diagnostic
  layout <- "
  AAB
  CDE
  "
  gg_full <- wrap_elements(full = gg_drc_lpv) + add_double_legend(gg_src_lpv, dat_src_sto) + gg_src_d_lpv + gg_src_th_lpv + gg_src_sh_lpv  + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_little.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()


  ####################################################
  # medium positivity violation
  ####################################################

  ## modify DAG
  D.base <- D.base + node("A",  distr = "rnorm", mean = (0.55)*L + 0.26, sd=0.08)
  Dset <- simcausal::set.DAG(D.base)
  # get data for plotting
  Odat <- simcausal::sim(DAG = Dset, n = n, rndseed = 1, verbose=F)
  disthalf_vec$A <- 0.5*sd(Odat$A)

  ## dose-response curve with CICI
  dat_drc <- get_drc(doses, Dset, Odat)
  ## dose-response curve with Kennedy
  l <- matrix(c(Odat$L, rep(0,n)), nrow=n, ncol=2)
  ce.res <- ctseff(y=Odat$Y, a=Odat$A, x=l, bw.seq=kennedy_search, sl.lib=sl_lib)
  # check that bandwidth choice is minimizer
  plot(ce.res$bw.risk$bw, ce.res$bw.risk$risk)
  # plot drc
  gg_drc_mpv <- plot_drc(dat_drc, ce.res, Odat)
  gg_drc_mpv_l <- plot_drc_legend(dat_drc, ce.res, Odat)
  # add diagnostic
  gg_drc_mpv_dia <- my_diagnostic(int_func=drc_intervention, int_func_param=dia_range(doses, n=length_out_dia_cdrc),
                                  data=Odat, disthalf_vec=disthalf_vec, col="black")
  gg_drc_mpv_noest <- plot_drc_no_estimation(dat_drc, Odat)
  gg_drc_mpv_noest_l <- plot_drc_no_estimation_legend(dat_drc, Odat)

  # deterministic shift
  dat_src_det <- get_src_det(shifts_det_mpv, Dset, sl_lib, Odat)
  # stochastic shift
  dat_src_sto <- get_src_sto(shifts_sto_mpv, Dset, sl_lib, Odat)
  # plot src
  gg_src_mpv <- plot_src(dat_src_det, dat_src_sto)
  # add diagnostic
  gg_src_mpv_dia <- my_diagnostic(int_func=src_det_intervention, int_func_param=dia_range(shifts_sto_mpv),
                                  data=Odat, disthalf_vec=disthalf_vec)


  # deterministic dynamic shift
  dat_src_detdy <- get_src_detdy(shifts_det_d_mpv, Dset, sl_lib, Odat)
  # stochastic dynamic shift
  dat_src_stody <- get_src_stody(shifts_sto_d_mpv, Dset, sl_lib, Odat)
  # plot src dynamic
  gg_src_d_mpv <- plot_src_d(dat_src_stody, dat_src_detdy)
  # add diagnostic
  gg_src_d_mpv_dia <- my_diagnostic(int_func=src_detdy_intervention, int_func_param=dia_range(shifts_sto_d_mpv),
                                    data=Odat, disthalf_vec=disthalf_vec)


  ## deterministic "threshold" lmtp
  dat_src_detth <- get_src_detth(thresholds_mpv, Dset, sl_lib, Odat, left=TRUE)
  # plot "threshold"
  gg_src_th_mpv <- plot_src_th(dat_src_detth)
  # add diagnostic
  gg_src_th_mpv_dia <- my_diagnostic(int_func=src_detth_intervention, int_func_param=dia_range(thresholds_mpv),
                                     data=Odat, disthalf_vec=disthalf_vec)


  # deterministic "shift" lmtp
  dat_src_detsh <- get_src_detsh(shifts_detsh_mpv, Dset, sl_lib, Odat)
  # plot "shift"
  gg_src_sh_mpv <- plot_src_sh(dat_src_detsh)
  # add diagnostic
  gg_src_sh_mpv_dia <- my_diagnostic(int_func=src_detsh_intervention, int_func_param=dia_range(shifts_detsh_mpv),
                                     data=Odat, disthalf_vec=disthalf_vec)


  ### save final plot with diagnostic
  layout <- "
  AAC
  BBD
  EGI
  FHJ
  "
  gg_full <- wrap_elements(full = gg_drc_mpv) + wrap_elements(full = gg_drc_mpv_dia) +
    add_double_legend(gg_src_mpv, dat_src_sto) + gg_src_mpv_dia +
    gg_src_d_mpv + gg_src_d_mpv_dia +
    gg_src_th_mpv + gg_src_th_mpv_dia +
    gg_src_sh_mpv  + gg_src_sh_mpv_dia +
    plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_dia_medium.pdf', sep=''), gg_full, device="pdf", width = 12, height = 10)
  dev.off()

  ### save final plot without diagnostic
  layout <- "
  AAB
  CDE
  "
  gg_full <- wrap_elements(full = gg_drc_mpv) + add_double_legend(gg_src_mpv, dat_src_sto) + gg_src_d_mpv + gg_src_th_mpv + gg_src_sh_mpv  + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_medium.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()
  
  # ### poster plot
  # layout <- "
  # AACE
  # BBDF
  # "
  # gg_full <- wrap_elements(full = gg_drc_mpv) + wrap_elements(full = gg_drc_mpv_dia) +
  #   add_double_legend(gg_src_mpv, dat_src_sto) + gg_src_mpv_dia +
  #   gg_src_th_mpv + gg_src_th_mpv_dia +
  #   plot_layout(design = layout)
  # gg_full <- gg_full + plot_annotation(title=" ",
  #                                      theme = theme(plot.title = element_text(size = 18)))
  # pdf()
  # ggsave(paste('figures/', descr_str, '_posterplot_medium.pdf', sep=''), gg_full, device="pdf", width = 16, height = 7)
  # dev.off()

  ####################################################
  # strong positivity violation
  ####################################################
  
  ## modify DAG
  D.base <- D.base + node("A",  distr = "rnorm", mean = (0.64)*L + 0.18, sd=0.05)
  Dset <- simcausal::set.DAG(D.base)
  # get data for plotting
  Odat <- simcausal::sim(DAG = Dset, n = n, rndseed = 1, verbose=F)
  disthalf_vec$A <- 0.5*sd(Odat$A)
  
  
  ## dose-response curve with CICI
  dat_drc <- get_drc(doses, Dset, Odat)
  ## dose-response curve with Kennedy
  l <- matrix(c(Odat$L, rep(0,n)), nrow=n, ncol=2)
  ce.res <- ctseff(y=Odat$Y, a=Odat$A, x=l, bw.seq=kennedy_search, sl.lib=sl_lib)
  # check that bandwidth choice is minimizer
  plot(ce.res$bw.risk$bw, ce.res$bw.risk$risk)
  # plot drc
  gg_drc_spv <- plot_drc(dat_drc, ce.res, Odat)
  gg_drc_spv_l <- plot_drc_legend(dat_drc, ce.res, Odat)
  gg_drc_spv_lh <- plot_drc_legend_high(dat_drc, ce.res, Odat)
  # add diagnostic
  gg_drc_spv_dia <- my_diagnostic(int_func=drc_intervention, int_func_param=dia_range(doses, n=length_out_dia_cdrc), 
                                  data=Odat, disthalf_vec=disthalf_vec, col="black")
  gg_drc_spv_noest <- plot_drc_no_estimation(dat_drc, Odat) 
  gg_drc_spv_noest_l <- plot_drc_no_estimation_legend(dat_drc, Odat) 
  
  # deterministic shift
  dat_src_det <- get_src_det(shifts_det_spv, Dset, sl_lib, Odat)
  # stochastic shift
  dat_src_sto <- get_src_sto(shifts_sto_spv, Dset, sl_lib, Odat)
  # plot src
  gg_src_spv <- plot_src(dat_src_det, dat_src_sto)
  # add diagnostic
  gg_src_spv_dia <- my_diagnostic(int_func=src_det_intervention, int_func_param=dia_range(shifts_sto_spv), 
                                  data=Odat, disthalf_vec=disthalf_vec)
  
  
  # deterministic dynamic shift
  dat_src_detdy <- get_src_detdy(shifts_det_d_spv, Dset, sl_lib, Odat)
  # stochastic dynamic shift
  dat_src_stody <- get_src_stody(shifts_sto_d_spv, Dset, sl_lib, Odat)
  # plot src dynamic 
  gg_src_d_spv <- plot_src_d(dat_src_stody, dat_src_detdy)
  # add diagnostic
  gg_src_d_spv_dia <- my_diagnostic(int_func=src_detdy_intervention, int_func_param=dia_range(shifts_sto_d_spv), 
                                    data=Odat, disthalf_vec=disthalf_vec)
  
  
  ## deterministic "threshold" lmtp
  dat_src_detth <- get_src_detth(thresholds_spv, Dset, sl_lib, Odat, left=TRUE)
  # plot "threshold"
  gg_src_th_spv <- plot_src_th(dat_src_detth)
  # add diagnostic
  gg_src_th_spv_dia <- my_diagnostic(int_func=src_detth_intervention, int_func_param=dia_range(thresholds_spv), 
                                     data=Odat, disthalf_vec=disthalf_vec)
  
  
  # deterministic "shift" lmtp
  dat_src_detsh <- get_src_detsh(shifts_detsh_spv, Dset, sl_lib, Odat)
  # plot "shift"
  gg_src_sh_spv <- plot_src_sh(dat_src_detsh)
  # add diagnostic
  gg_src_sh_spv_dia <- my_diagnostic(int_func=src_detsh_intervention, int_func_param=dia_range(shifts_detsh_spv), 
                                     data=Odat, disthalf_vec=disthalf_vec)
  
  ### save final plot with diagnostic
  layout <- "
  AAC
  BBD
  EGI
  FHJ
  "
  gg_full <- wrap_elements(full = gg_drc_spv) + wrap_elements(full = gg_drc_spv_dia) +
    add_double_legend(gg_src_spv, dat_src_sto) + gg_src_spv_dia +
    gg_src_d_spv + gg_src_d_spv_dia +
    gg_src_th_spv + gg_src_th_spv_dia +
    gg_src_sh_spv  + gg_src_sh_spv_dia +
    plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_dia_strong.pdf', sep=''), gg_full, device="pdf", width = 12, height = 10)
  dev.off()
  
  ### save final plot without diagnostic
  layout <- "
  AAB
  CDE
  "
  gg_full <- wrap_elements(full = gg_drc_spv) + add_double_legend(gg_src_spv, dat_src_sto) + gg_src_d_spv + gg_src_th_spv + gg_src_sh_spv  + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/', descr_str, '_strong.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()
  
  
  
  
  #######################################
  # additional plots (without diagnostic) # well let's make it with diagnostic
  #######################################
  
  ## save plots by type not strength of positivity violation

  # dose response curves full with dia
  layout <- "
  A
  B
  C
  D
  E
  F
  "

  gg_drc_lpv <- gg_drc_lpv + ggtitle("overlapping data")
  gg_drc_mpv <- gg_drc_mpv + ggtitle("adjacent data")
  gg_drc_spv <- gg_drc_spv + ggtitle("divided data")

  gg_full <- (gg_drc_lpv_l + ggtitle("overlapping data")) + gg_drc_lpv_dia + gg_drc_mpv + gg_drc_mpv_dia +
             gg_drc_spv + gg_drc_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_drc.pdf', sep=''), gg_full, device="pdf", width = 8, height = 11.5)
  dev.off()
  
  # dose response curves without dia and actually also no curves
  layout <- "
  A
  B
  C
  "
  

  gg_drc_lpv_noest <- gg_drc_lpv_noest + ggtitle("overlapping data")
  gg_drc_mpv_noest <- gg_drc_mpv_noest + ggtitle("adjacent data")
  gg_drc_spv_noest <- gg_drc_spv_noest + ggtitle("divided data")

  
  gg_full <- gg_drc_lpv_noest_l + gg_drc_mpv_noest + gg_drc_spv_noest + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  gg_full
  
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_drc_wellnodia.pdf', sep=''), gg_full, device="pdf", width = 8, height = 9)
  dev.off()
  
  # dose response curves small with dia
  layout <- "
  ABC
  DEF
  "
  
  gg_drc_lpv <- gg_drc_lpv + ggtitle("overlapping data") + theme(legend.position="none")
  gg_drc_mpv <- gg_drc_mpv + ggtitle("adjacent data") + theme(legend.position="none")
  gg_drc_spv <- gg_drc_spv + ggtitle("divided data") + theme(legend.position="none")
  
  
  gg_full <- gg_drc_lpv + gg_drc_mpv + (gg_drc_spv_lh + ggtitle("divided data")) + 
             gg_drc_lpv_dia + gg_drc_mpv_dia + gg_drc_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_drc_small.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()



  # shift response curves
  layout <- "
  ABC
  DEF
  "

  gg_src_lpv <- gg_src_lpv + ggtitle("overlapping data")
  gg_src_mpv <- gg_src_mpv + ggtitle("adjacent data")
  gg_src_spv <- add_double_legend(gg_src_spv + ggtitle("divided data"), dat_src_sto)

  gg_full <- gg_src_lpv + gg_src_mpv + gg_src_spv + 
             gg_src_lpv_dia + gg_src_mpv_dia + gg_src_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_src.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()


  # dynamic shift response curves
  layout <- "
  ABC
  DEF
  "

  gg_src_d_lpv <- gg_src_d_lpv + ggtitle("overlapping data")
  gg_src_d_mpv <- gg_src_d_mpv + ggtitle("adjacent data")
  gg_src_d_spv <- add_double_legend(gg_src_d_spv + ggtitle("divided data"), dat_src_sto)

  gg_full <- gg_src_d_lpv + gg_src_d_mpv + gg_src_d_spv + 
             gg_src_d_lpv_dia + gg_src_d_mpv_dia + gg_src_d_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_src_d.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()


  # "threshold" response curves
  layout <- "
  ABC
  DEF
  "

  gg_src_th_lpv <- gg_src_th_lpv + ggtitle("overlapping data")
  gg_src_th_mpv <- gg_src_th_mpv + ggtitle("adjacent data")
  gg_src_th_spv <- add_single_legend(gg_src_th_spv + ggtitle("divided data"), dat_src_detth)

  gg_full <- gg_src_th_lpv + gg_src_th_mpv + gg_src_th_spv + 
             gg_src_th_lpv_dia + gg_src_th_mpv_dia + gg_src_th_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_src_th.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()


  # "shift" response curves
  layout <- "
  ABC
  DEF
  "

  gg_src_sh_lpv <- gg_src_sh_lpv + ggtitle("overlapping data")
  gg_src_sh_mpv <- gg_src_sh_mpv + ggtitle("adjacent data")
  gg_src_sh_spv <- add_single_legend(gg_src_sh_spv + ggtitle("divided data"), dat_src_detth)

  gg_full <- gg_src_sh_lpv + gg_src_sh_mpv + gg_src_sh_spv + 
             gg_src_sh_lpv_dia + gg_src_sh_mpv_dia + gg_src_sh_spv_dia + plot_layout(design = layout)
  gg_full <- gg_full + plot_annotation(title=" ",
                                       theme = theme(plot.title = element_text(size = 18)))
  pdf()
  ggsave(paste('figures/by_type/', descr_str, '_dia_src_sh.pdf', sep=''), gg_full, device="pdf", width = 12, height = 6)
  dev.off()
}

