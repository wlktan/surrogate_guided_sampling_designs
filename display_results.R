####################################################################################
# Displays results for data example illustration for Biostatistics manuscript
# Surrogate-guided sampling designs for classification of rare outcomes from electronic medical records data
# Tan and Heagerty, 2020
####################################################################################

library(dplyr)
library(magrittr)
library(xtable)

sim.res1 <- read.csv(paste0("p1_lasso_n100.csv")) %>%
  mutate(n_train = 100)
sim.res2 <- read.csv(paste0("p1_lasso_n250.csv")) %>%
  mutate(n_train = 250)
sim.res3 <- read.csv(paste0("p1_lasso_n500.csv")) %>%
  mutate(n_train = 500)

bind_rows(sim.res1, sim.res2, sim.res3) %>%
  group_by(n_train) %>%
  summarise(srs_est = mean(srs),
            srs_ci_lower = quantile(srs, probs = 0.05),
            srs_ci_upper = quantile(srs, probs = 0.95),
            sgs_est = mean(sgs),
            sgs_ci_lower = quantile(sgs, probs = 0.05),
            sgs_ci_upper = quantile(sgs, probs = 0.95)) %>%
  transmute(n_train,
            AUC_srs = paste0(round(srs_est, 2), " (", round(srs_ci_lower, 2), ", ", round(srs_ci_upper, 2), ")"),
            AUC_sgs = paste0(round(sgs_est, 2), " (", round(sgs_ci_lower, 2), ", ", round(sgs_ci_upper, 2), ")")) %>%
  as.data.frame() %>%
  set_colnames(c("Training sample size",
                 "$\\hat{AUC}(\\mathbf{D}^{SRS}(n))$",
                 "$\\hat{AUC}(\\mathbf{D}^{SGS}(n))$")) %>%
  xtable(.,
         caption = "Average validation AUC (95\\% C.I.) for various training sample sizes, 
         based on B=100 bootstrap resamples, for illustration of surrogate-guided sampling (SGS) designs on radiology reports drawn from the LIRE data set.",
         label = "sgs_tb:data_results",
         align = c("c", "c", "c", "c")) %>%
  print(., table.placement = "h!",
        caption.placement = "top",
        include.rownames = FALSE,
        sanitize.text.function=function(x){x})
