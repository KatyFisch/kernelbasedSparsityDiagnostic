library(simcausal)
source("diagnostic_function.R")
library(tidyverse)



D.unif <- DAG.empty() +   node("A",  distr = "rnorm", mean=0, sd=1, t=1) +
  node("L", distr = "rnorm", mean=0, sd=1, t=1:1001)

Dset <- simcausal::set.DAG(D.unif)
# get data for plotting
data <- simcausal::sim(DAG = Dset, n = 1000, rndseed = 1, verbose=F)
data["ID"] <- NULL

drc_intervention <- function(data, param){
  data["A_1"] <- param
  return(data)
}

dia_list <- list()
j <- 1
# use diagnostic
for (i in seq(2,10,1)){
  data_sub <- data[,1:i]
  disthalf_vec=c("A_1"=0.5, setNames(rep(0.5, i-1), paste0("L_", 1:(i-1))))
  dia <- my_diagnostic(drc_intervention, 0, data_sub, disthalf_vec, plot.out=FALSE)
  dia_list[[j]] <- dia$diagnostic
  j <- j + 1
}

df_long <- enframe(dia_list, name = "Group", value = "Value") %>%
  unnest(cols = c(Value))





sim_plot <- ggplot(df_long, aes(x = factor(Group), y = Value)) +
  geom_boxplot(fill = "#F0E442", color = "#0072B2", width = 0.6, outlier.shape = 16, outlier.size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    x = "number of variables in adjustment set",
    y = "Effective Data Points (EDP)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0)
  )

pdf()
ggsave('figures/add_sim.pdf', sim_plot, device="pdf", width = 5, height = 4)
dev.off()
