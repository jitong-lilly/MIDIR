#################################################################################
# FILENAME    : MIDIR_analysis_simgbda_pw.R
# AUTHOR      : Jitong Lou <lou_jitong@lilly.com>
# DATE        : 06/10/2024
# DESCRIPTION : Proposed direct estimation method for PW imputation, tested on simulated data based on the real data application in the article
#################################################################################

rm(list = ls())

# Load packages ####
library(dplyr)
library(Matrix) # bdiag
library(purrr) # map
library(ggplot2)
options(digits = 6)
options(dplyr.summarise.inform = FALSE) 

# Load functions ####
source(file.path("./MIDIR_functions_simgbda.R"))

# Generate dataset ####
set.seed(54523)
gbdalong = simgbda(min_rd = 3)

gbdalong1 = gbdalong %>% mutate(time_order=factor(time)) 
gbdadata_complete=tidyr::complete(gbdalong1, id, time_order, fill=list(miss=1)) %>% group_by(id) %>% 
  mutate(base = rep(mean(base, na.rm = TRUE), length(base)), TRTSORT2 = rep(mean(TRTSORT2, na.rm = TRUE), length(TRTSORT2))) %>%
  mutate(time=time_order) %>% # time for formula, time_order for waves argument
  mutate(subject=id) # subject for id argument

df_fit_tf_list=list(
  control = gbdadata_complete %>% filter(TRTSORT2==4) %>%
    mutate(value_pw_chg = value_raw_chg), 
  drug1 = gbdadata_complete %>% filter(TRTSORT2==1) %>%
    mutate(value_pw_chg = value_raw_chg),
  drug2 = gbdadata_complete %>% filter(TRTSORT2==2) %>%
    mutate(value_pw_chg = value_raw_chg)
) # create a list of data for each treatment group

df_fit_tf_completer_list=list(
  control = gbdadata_complete %>% filter(TRTSORT2==4) %>%
    mutate(value_pw_chg = value_raw_chg), 
  drug1 = gbdadata_complete %>% filter(TRTSORT2==1) %>%
    group_by(subject) %>%
    mutate(last_miss = as.numeric(is.na(last(value_raw)))) %>%
    mutate(
      value_raw_pw = value_raw
    ) %>%
    mutate(
      value_pw = value_raw_pw * ifelse(
        miss == 1 | last_miss == 1, NA, 1)
    ) %>%
    mutate(
      value_raw_pw_chg = value_raw_pw - base,
      value_pw_chg = value_pw - base
    ) %>%
    ungroup(subject) %>%
    filter(time_order != 0),
  drug2 = gbdadata_complete %>% filter(TRTSORT2==2) %>%
    group_by(subject) %>%
    mutate(last_miss = as.numeric(is.na(last(value_raw)))) %>%
    mutate(
      value_raw_pw = value_raw
    ) %>%
    mutate(
      value_pw = value_raw_pw * ifelse(
        miss == 1 | last_miss == 1, NA, 1)
    ) %>%
    mutate(
      value_raw_pw_chg = value_raw_pw - base,
      value_pw_chg = value_pw - base
    ) %>%
    ungroup(subject) %>%
    filter(time_order != 0)
) # create a list of data of completers for each treatment group

N_list=lapply(df_fit_tf_list,function(x){n_distinct(x$id)}) %>% as.list # extract number of subjects for each treatment group

K=length(df_fit_tf_list)-1 # set number of treatment groups, not including the placebo arm

# Parameter estimation ####
result_betaREF = fit_data_PW_beta(
  data = df_fit_tf_list$control, # all groups use the same reference
  N = N_list$control,
  K = K, beta_family = gaussian(), beta_corstr = "unstructured"
) # MMRM for reference group
range(rowSums(result_betaREF$U)) # value of estimating equations, should be close to 0

result_beta = mapply(
  FUN = fit_data_PW_beta,
  data = df_fit_tf_completer_list,
  N = N_list, 
  MoreArgs = list(K = K, beta_family = gaussian(), beta_corstr = "unstructured"),
  SIMPLIFY = FALSE
) # MMRM for MAR for each group
sapply(result_beta, function(x){
  range(rowSums(x$U))
}) # value of estimating equations, should be close to 0

result_pi = mapply(
  FUN = fit_data_PW_pi,
  data = df_fit_tf_list,
  MoreArgs = list(K = K),
  SIMPLIFY = FALSE
) # proportion of missingness for each group
sapply(result_pi, function(x){
  range(sum(x$U))
}) # value of estimating equations, should be close to 0

result_nuMISS = mapply(
  FUN = fit_data_PW_nuMISS,
  data = df_fit_tf_list,
  MoreArgs = list(K = K),
  SIMPLIFY = FALSE
) # baseline for subjects with missing outcomes for each group
sapply(result_nuMISS, function(x){
  range(sum(x$U))
}) # value of estimating equations, should be close to 0

result_nu = mapply( 
  FUN = fit_data_PW_nu,
  data = df_fit_tf_list,
  MoreArgs = list(K = K),
  SIMPLIFY = FALSE
) # baseline for each group
sapply(result_nu, function(x){
  range(sum(x$U))
}) # value of estimating equations, should be close to 0


nu_hat = result_nu %>% sapply(function(x){x$coef})
mean_base = ifelse(
  (class(nu_hat))[1] == "matrix", 
  rowSums(nu_hat %*% (bdiag(N_list)/J)), 
  weighted.mean(nu_hat, unlist(N_list))
)
mean_base # mean for baseline for all subjects

var_base = df_fit_tf_list %>%
  do.call(rbind, .) %>%
  filter(time == K) %>%
  pull(base) %>%
  var
var_base  # variance of baseline for all subjects

result_delta = mapply( 
  FUN = calc_est_PW_delta_v1,
  beta_result = result_beta, pi_result = result_pi, 
  nuMISS_result = result_nuMISS, nu_result = result_nu, 
  reference = as.list(names(df_fit_tf_list) == "control"),
  MoreArgs = list(betaREF_result = result_betaREF),
  SIMPLIFY = FALSE
) # treatment effect for each group


l_delta = cbind(-1, diag(2)); colnames(l_delta) = names(result_delta); 
l_delta # comparison matrix
names_comp = apply(l_delta, 1, function(x){
  paste(x[x!=0], names(x[x!=0]), sep = "*", collapse = " + ")
})

result_theta = calc_est_PW_theta_v1(
  result_delta, l_delta = l_delta, mean_base, var_base,
  tip = FALSE, complete = FALSE
)  # treatment difference for each comparison

# Tipping point analysis ####
tip_delta = expand.grid(
  control = seq(-10,0,0.2),
  drug1 = seq(10,0,-0.2)
) %>% mutate(drug2 = drug1) # set lower and upper limits for potential tipping points

result_tip = cbind(
  tip_delta,
  calc_est_PW_TIP_v1(
    beta_result = result_beta, betaREF_result = result_betaREF, pi_result = result_pi,
    nuMISS_result = result_nuMISS, nu_result = result_nu,
    l_delta, mean_base, var_base, tip = TRUE, tip_delta, complete = FALSE)
) # obtain tipping points and statistics

# Plot tipping point analysis result ####
data=as.data.frame(cbind(
  z_v=result_tip$z_unadj2, # z stat of drug2 vs control
  delta_t_v=result_tip$drug2, # delta adding to drug 2
  delta_c_v=result_tip$control
)) # z_v is the vector for z stats; delta_t_v and delta_c_v are tipping point vectors for treatment and comparator

# define the number of breaks
mybreaks = c(
  min(data$z_v),qnorm(0.0001/2),qnorm(0.001/2),qnorm(0.01/2),
  qnorm(0.025/2),qnorm(0.05/2), max(data$z_v)
)

# Output contour plot 
gbda_pw_tip=ggplot(data, aes(x = delta_t_v, y = delta_c_v, z = z_v)) +
  geom_contour_filled(breaks= mybreaks, show.legend = TRUE) +
  scale_fill_manual(palette=mycolors, values=breaklabel(6), name="P value", drop=FALSE) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-10, 0)) +
  theme(legend.position = "right")+
  geom_hline(yintercept=0, color='black',linetype='dashed')+
  geom_vline(xintercept=0.0, color='black',linetype='dashed')+
  annotate("point", x = 0.0, y = 0.0, colour = "black")+
  xlab(expression(paste(italic(delta)," of Drug2",sep=''))) +
  ylab(expression(paste(italic(delta)," of Placebo",sep=''))) +
  labs(title = "2-way tipping point analysis", 
       subtitle = "PW imputation")
gbda_pw_tip

# Output of estimates ####
result_temp = rbind(
  result_theta$delta %>% do.call(cbind,.),
  result_theta$theta %>% do.call(cbind,.) 
)
rownames(result_temp)[-(1:length(result_delta))] = names_comp
result = data.frame(
  treatment = rownames(result_temp),
  N = c(unlist(N_list), rep("NA", nrow(l_delta))),
  pi = c(unlist(map(result_pi, "coef")), rep("NA", nrow(l_delta))),
  result_temp
) %>%
  dplyr::select(-est,-var,-se,-est_adj,-var_adj,-se_adj) %>%
  arrange(treatment)
rownames(result) = NULL
result

## END OF FILE

