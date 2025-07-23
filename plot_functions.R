
## SINGLE TIME POINT

plot_drc <- function(dat_drc, ce.res, Odat){
  ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth), linetype="solid") + 
    geom_line(data=dat_drc, aes(x=dose, y=est), linetype="dashed") +
    geom_line(data=ce.res$res, aes(x=a.vals, y=est), linetype="dotted") + # Kennedy
    guides(fill="none", color = "none", linetype = "none", shape = "none") +
    scale_color_manual(values = c("#E69F00", "#0072B2"))
}

plot_drc_legend <- function(dat_drc, ce.res, Odat){
  ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    geom_line(data=dat_drc, aes(x=dose, y=est, linetype="g-formula estimate")) +
    geom_line(data=ce.res$res, aes(x=a.vals, y=est, linetype="doubly robust estimate")) + # Kennedy
    ggtitle("Dose RC") + 
    theme(legend.position = c(0.65, 1.155), legend.direction = "horizontal", legend.box = "horizontal") +
    labs(color='L:') +
    scale_linetype_manual("",
                         breaks = c("truth", "g-formula estimate", "doubly robust estimate"),
                         values = c("solid", "dashed", "dotted")) +
    scale_color_manual(values = c("#E69F00", "#0072B2"))
}

plot_drc_legend_high <- function(dat_drc, ce.res, Odat){
  ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    geom_line(data=dat_drc, aes(x=dose, y=est, linetype="g-formula estimate")) +
    geom_line(data=ce.res$res, aes(x=a.vals, y=est, linetype="doubly robust estimate")) + # Kennedy
    ggtitle("Dose RC") + #c(0.65, 1.105)
    theme(legend.position = c(0.2, 1.205), legend.direction = "horizontal", legend.box = "horizontal") +
    labs(color='L:') +
    scale_linetype_manual("",
                          breaks = c("truth", "g-formula estimate", "doubly robust estimate"),
                          values = c("solid", "dashed", "dotted")) +
    scale_color_manual(values = c("#E69F00", "#0072B2"))
}

plot_drc_no_estimation <- function(dat_drc, Odat){
  ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth), linetype="solid") + 
    ggtitle("Dose RC") + guides(fill="none", color = "none", linetype = "none", shape = "none") +
    scale_color_manual(values = c("#E69F00", "#0072B2"))
}

plot_drc_no_estimation_legend <- function(dat_drc, Odat){
  ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    ggtitle("Dose RC") + 
    theme(legend.position = c(0.85, 1.105), legend.direction = "horizontal", legend.box = "horizontal") +
    labs(color='L:') +
    scale_linetype_manual("",
                         breaks = c("truth"),
                         values = c("solid")) +
    scale_color_manual(values = c("#E69F00", "#0072B2"))
}


plot_src <- function(dat_src_det, dat_src_sto){
  ggplot(data=dat_src_sto, aes(x=shift, y=response, group=1)) + 
  geom_line(col="#009E73", linetype="solid") + 
  geom_line(aes(x=shift, y=lmtp, group=2), linetype="dashed", col="#009E73") +
  geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#009E73") +
  geom_line(data=dat_src_det, aes(x=shift, y=response, group=4), col="#56B4E9") + 
  geom_line(data=dat_src_det, aes(x=shift, y=lmtp, group=5), col="#56B4E9", linetype="dashed") +
  geom_ribbon(data=dat_src_det, aes(ymin=lmtp_low, ymax=lmtp_high, group=6), alpha=0.3, fill="#56B4E9") +
  ggtitle("relative Shift RC")
}


plot_src_d <- function(dat_src_stody, dat_src_detdy){
  ggplot(data=dat_src_stody, aes(x=shift, y=response, group=1)) + geom_line(col="#009E73") + 
    geom_line(aes(x=shift, y=lmtp, group=2), col="#009E73", linetype="dashed") +
    geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#009E73") +
    geom_line(data=dat_src_detdy, aes(x=shift, y=response, group=4), linetype=1, alpha=1, col="#56B4E9") + 
    geom_line(data=dat_src_detdy, aes(x=shift, y=lmtp, group=5), linetype=2, alpha=1, col="#56B4E9") +
    geom_ribbon(data=dat_src_detdy, aes(ymin=lmtp_low, ymax=lmtp_high, group=6), alpha=0.3, fill="#56B4E9") +
    ggtitle("relative SRC: dynamic (only L=1)") 
}


plot_src_th <- function(dat_src_detth){
  ggplot(data=dat_src_detth, aes(x=threshold, y=response, group=1)) + geom_line(col="#56B4E9") + 
  geom_line(aes(x=threshold, y=lmtp, group=2), linetype="dashed", col="#56B4E9") +
  geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#56B4E9") + 
  #ggtitle("relative threshold RC")
  ggtitle("relative Threshold RC")
}


plot_src_sh <- function(dat_src_detsh){
  ggplot(data=dat_src_detsh, aes(x=shift, y=response, group=1)) + geom_line(col="#56B4E9") + 
  geom_line(aes(x=shift, y=lmtp, group=2), linetype="dashed", col="#56B4E9") +
  geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#56B4E9") +
  ggtitle("relative cond. shift RC")
}

# legend(ary) hacks

add_double_legend <- function(plot, dat_src_sto){
  # src plot again to steal the legend from it
  src_plot <- ggplot(data=dat_src_sto, aes(x=shift, y=response, group=1)) + 
    geom_line(aes(col="stochastic", linetype="truth")) + 
    geom_line(aes(x=shift, y=lmtp, group=2, linetype="estimate"), col="#009E73") +
    geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#009E73") +
    geom_line(data=dat_src_sto, aes(x=shift, y=response, group=4, col="deterministic")) + 
    geom_line(data=dat_src_sto, aes(x=shift, y=lmtp, group=5), col="#56B4E9", linetype="dashed") +
    geom_ribbon(data=dat_src_sto, aes(ymin=lmtp_low, ymax=lmtp_high, group=6), alpha=0.3, fill="#56B4E9") +
    ggtitle("relative SRC") +
    scale_color_manual("", breaks=c('deterministic', 'stochastic'),
                       values=c('#56B4E9', '#009E73')) +
    scale_linetype_manual("", breaks = c("truth", "estimate"),
                          values = c("solid", "dashed")) +
    guides(linetype = guide_legend(order = 1, override.aes = list(color = c('azure4', 'azure4'))),
           color = guide_legend(order = 2)) + 
    theme(legend.position = "top", legend.direction = "vertical", legend.box = "horizontal",
          legend.margin = ggplot2::margin(-1.9,-0.1,0.2,-0.1, unit="cm"))

  p2 <- my_get_legend(src_plot) 
  
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1, l = 4.5, b = 1, r = 10) # l used to be 2.5
  )
  return(plot + p2 + plot_layout(design = layout))
  
}

add_single_legend <- function(plot, dat_src_detth){
  #"threshold" plot again to steal from
  thresh_plot <- ggplot(data=dat_src_detth, aes(x=threshold, y=response, group=1)) + 
    geom_line(aes(linetype="truth"), col="#56B4E9") + 
    geom_line(aes(x=threshold, y=lmtp, group=2, linetype="estimate"), col="#56B4E9") +
    geom_ribbon(aes(ymin=lmtp_low, ymax=lmtp_high, group=3), alpha=0.3, fill="#56B4E9") + 
    ggtitle("relative threshold RC") +
    scale_linetype_manual("", breaks = c("truth", "estimate"), values = c("solid", "dashed")) + 
    theme(legend.position = c(0.65,0.5), legend.direction = "vertical", legend.box = "horizontal",
          legend.margin = ggplot2::margin(-1.9,-0.1,0.2,-0.1, unit="cm"))
  
  p2 <- my_get_legend(thresh_plot) 
  
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1, l = 6, b = 1, r = 10)
  )
  return(plot + p2 + plot_layout(design = layout))
}

add_drc_2_legend <- function(plot, dat_drc, Odat){
  drc_plot <- ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    theme(legend.position = c(0.6, 1.105), legend.direction = "horizontal", legend.box = "horizontal") + 
    labs(color='L:')
  p2 <- my_get_legend(drc_plot) 
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1, l = 4, b = 1, r = 10)
  )
  return(plot + p2 + plot_layout(design = layout))
}

add_drc_3_legend <- function(plot, dat_drc, Odat){
  drc_plot <- ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    theme(legend.position = c(0.6, 1.105), legend.direction = "horizontal", legend.box = "horizontal") + 
    labs(color='L:') +
    scale_linetype_manual("",
                          breaks = c("truth"),
                          values = c("solid"))
  
  p2 <- my_get_legend(drc_plot) 
  
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1, l = 4, b = 1, r = 10)
  )
  return(plot + p2 + plot_layout(design = layout))
}


add_drc_5_legend <- function(plot, dat_drc, Odat){
  drc_plot <- ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    geom_line(data=dat_drc, aes(x=dose, y=est, linetype="g-formula estimate")) +
    geom_line(data=ce.res$res, aes(x=a.vals, y=est, linetype="doubly robust estimate")) + # Kennedy
    ggtitle("Dose RC") + 
    theme(legend.position = c(0.6, 1.105), legend.direction = "horizontal", legend.box = "horizontal") + 
    labs(color='L:') +
    scale_linetype_manual("",
                          breaks = c("truth", "g-formula estimate", "doubly robust estimate"),
                          values = c("solid", "dashed", "dotted"))
  
  p2 <- my_get_legend(drc_plot) 
  
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1.5, l = 2.5, b = 2, r = 10)
  )
  return(plot + p2 + plot_layout(design = layout))
}


add_drc_5_legend_small <- function(plot, dat_drc, Odat){
  drc_plot <- ggplot() + geom_point(data=Odat, aes(x=A, y=Y, color=factor(L), group=1), alpha=0.5) +
    geom_line(data=dat_drc, aes(x=dose, y=truth, linetype="truth")) + 
    geom_line(data=dat_drc, aes(x=dose, y=est, linetype="g-formula estimate")) +
    geom_line(data=ce.res$res, aes(x=a.vals, y=est, linetype="doubly robust estimate")) + # Kennedy
    ggtitle("Dose RC") + 
    theme(legend.position = c(0.6, 1.105), legend.direction = "horizontal", legend.box = "horizontal") + 
    labs(color='L:') +
    scale_linetype_manual("",
                          breaks = c("truth", "g-formula estimate", "doubly robust estimate"),
                          values = c("solid", "dashed", "dotted"))
  
  p2 <- my_get_legend(drc_plot) 
  
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 10),
    area(t = 1, l = 4, b = 1, r = 10)
  )
  return(plot + p2 + plot_layout(design = layout))
}

## helper

my_get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  return(legend)
}


























