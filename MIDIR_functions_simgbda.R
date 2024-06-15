#################################################################################
# FILENAME    : MIDIR_functions_simgbda.R
# AUTHOR      : Jitong Lou <lou_jitong@lilly.com>
# DATE        : 06/10/2024
# DESCRIPTION : Used functions for the proposed direct estimation method for all imputations, tested on simulated data based on the real data application in the article
#################################################################################

# functions for data generation ####
simgbda = function(min_rd = 3){
  # minrd: minimum number of retrieved dropouts in any treatment groupu. This argument must be >= 3.
  if (min_rd < 3){
    print("Please make sure the minimum number of retrieved dropouts is not less than 3.")
  } else {
    ## generate data from MVN
    ## 3 treatment groups
    ## dula 0.75mg
    meanvec1=c(8.05,6.77,6.74)
    covmat1=matrix(c(1.546,0.854,0.701,0.85361,0.89603,0.77214,0.70148,0.77214,0.93348),nrow=3)
    
    ## dula 1.5
    meanvec2=c(8.10,6.53,6.53)
    covmat2=matrix(c(1.80,0.82,0.76,0.82,0.70,0.61,0.76,0.61,0.77),nrow=3)
    
    ## Placebo
    meanvec4=c(8.06,7.50,7.38)
    covmat4=matrix(c(1.71,1.01,0.528,1.01,1.136,0.655,0.528,0.655,0.989),nrow=3)
    
    ## simulate data until numbers of RDs in all treatments are sufficient
    minrd=-Inf
    while (minrd<min_rd){
      #sample size
      n1=280
      n2=279
      n4=141
      
      #hba1c for each of three groups
      data1=MASS::mvrnorm(n = n1, 
                    mu = meanvec1,  
                    Sigma = covmat1)
      data2=MASS::mvrnorm(n=n2,
                    mu = meanvec2,  
                    Sigma = covmat2)
      data4=MASS::mvrnorm(n=n4,
                    mu = meanvec4,  
                    Sigma = covmat4)
      
      #treatment indicators
      trt1=rep(1,n1)
      trt2=rep(2,n2)
      trt4=rep(4,n4)
      
      #adherence indicator for each of three groups at time 2 only, proportions from actual GBDA data
      a1=rbinom(n1,1,0.9143)
      a2=rbinom(n2,1,0.9104)
      a4=rbinom(n4,1,0.8794)
      
      #function to generate a bernoulli conditional on 1 vs. 0. 
      #Probability of missing=0.65 if patient not on treatment, .005 if patient is on treatment at time 2 only
      
      genmiss=function(x){
        if (x==0){return(rbinom(1,1,0.65))}
        else if (x==1){return(rbinom(1,1,0.005))}
      }
      
      #generate missing indicators for last visit (time 2 only)
      miss1=sapply(a1,genmiss)
      miss2=sapply(a2,genmiss)
      miss4=sapply(a4,genmiss)
      
      #generate id
      id1=seq(1,n1)
      id2=seq(n1+1,(n1+n2))
      id4=seq(n1+n2+1,n1+n2+n4)
      
      #start piecing them together
      gbdawide=data.frame(rbind(cbind(id1,trt1,data1,a1,miss1),cbind(id2,trt2,data2,a2,miss2),cbind(id4,trt4,data4,a4,miss4)))
      names(gbdawide)=c("id","trtsort2","a1c0","a1c1","a1c2","ONTRTFLAG","miss")
      
      #Convert wide to long and calculate change from baseline
      
      gbdalong <- gbdawide %>% 
        tidyr::pivot_longer(
          cols = `a1c0`:`a1c2`, 
          names_to = "timechr",
          values_to = "a1c"
        ) %>% mutate(time=as.numeric(substring(timechr,4,4)))
      
      #join with baseline data to calculate change
      gbdalong1=left_join(gbdalong,gbdawide%>%dplyr::select(id,a1c0),by="id") %>%
        mutate(value_raw_chg=ifelse( (miss==1) & (time==2),NA,a1c-a1c0),value_raw=ifelse((miss==1)&(time==2),NA,a1c)) %>%
        dplyr::filter(time>0) %>% rename(base=a1c0,TRTSORT2=trtsort2) %>% 
        mutate(miss=ifelse(time==1,0,miss),ONTRTFLAG=ifelse(time==1,1,ONTRTFLAG)) # Patient always on treatment at the first time point, always no missing at first time point
      #if miss=1 and time=2 set a1c to be missing
      
      ## get number of RD's by each treatment group. RD's are those who do not adhere to treatment and have non-missing primary
      minrd=gbdalong1 %>% dplyr::filter((time==2)&(ONTRTFLAG==0)&(miss==0)) %>%
        group_by(TRTSORT2) %>%
        summarise(n_rd=n()) %>% 
        summarise(minrd=min(n_rd)) %>% as.numeric()
    } # end of while (minrd<min_rd)
    
    print("Numbers of retrieved dropouts for each treatment group are:")
    print(gbdalong1 %>% dplyr::filter((time==2)&(ONTRTFLAG==0)&(miss==0)) %>%
            group_by(TRTSORT2) %>%
            summarise(n_rd=n()))
    return(gbdalong1)
  } # end of the else condition of if (min_rd < 3)
} # end of the function

# function for all models ####
working_correlation_v2 = function(K, corstr, corr){
  ## Function for working correlation matrix in MMRM estimation
  
  # https://statisticaloddsandends.wordpress.com/2020/02/07/generating-correlation-matrix-for-ar1-model/
  if (corstr == "independence"){
    return(diag(K))
  } else if (corstr == "ar1"){
    return(outer(seq(K), seq(K), function(i, j){corr^abs(i-j)}))
  } else if (corstr == "exchangeable"){
    return(stats::toeplitz(corr^(c(0,rep(1,K-1)))))
  } else if (corstr == "unstructured"){
    output = diag(K)
    output[lower.tri(output)] = corr
    output[upper.tri(output)] = t(output)[upper.tri(t(output))]
    return(output)
  } # consider incorporate "userdefined" later
}

mycolors <- function(x) {
  ## Function to return the desired number of colors in tipping point analysis
  colors<-c(colorRampPalette(c("dark green","green"))(x-1), "red")
  return(colors)
}

breaklabel <- function(x){
  ## Function to create labels for legends in tipping point analysis
  labels<- c(
    "<.0001 superiority of LY","<.001 superiority of LY", "<.01 superiority of LY", 
    "<.025 superiority of LY", "<.05 superiority of LY","superiority of LY not achieved"
  )
  return(labels[1:x])
}

# functions for RTB ####
fit_data_RTB_beta = function(data, N, K, beta_family, beta_corstr, eps = 0){
  ## Function for estimation of MMRM
  
  ## fit GEE, data must be sorted by subject, can check attributes(beta_fit$data)$groups to verify
  beta_fit = geepack::geeglm(
    formula = value_rtb_chg ~ factor(time)*base, # change from baseline
    id = subject,
    data = data, 
    waves = time_order, # waves argument handles missing value
    scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
    family = beta_family, 
    # corstr = "unstructured"
    corstr = ifelse(K >= 3, beta_corstr, "ar1")
  )
  # https://faculty.washington.edu/heagerty/Courses/b571/homework/geepack-paper.pdf
  # https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm
  
  beta_coef = coef(beta_fit) # extract coefficients
  p_beta = length(beta_coef) # dimension of coefficients
  beta_cov = vcov(beta_fit) # extract covariance matrix of coefficients, https://rdrr.io/cran/merDeriv/man/vcov.lmerMod.html
  beta_R = working_correlation_v2(K, beta_corstr, beta_fit$geese$alpha)
  beta_residual = beta_fit$y-as.vector(beta_fit$fitted.values) # beta_fit$fitted.values is a matrix, this value is different from beta_fit$residuals
  
  ## derive score equations for main effects
  # block diagnol matrix, https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/bdiag.html
  # convert matrix to list, https://r-lang.com/how-to-convert-r-matrix-to-list/#:~:text=How%20to%20Convert%20R%20Matrix%20to%20List%201,into%20the%20list%20of%20vectors.%203%20See%20also
  
  beta_nonmissing_factor = beta_fit$data %>%
    group_by(subject) %>%
    summarise(n=sum(!is.na(value_rtb_chg))) %>% #  change from baselina
    # summarise(n=sum(!is.na(value_rtb))) %>% # outcome
    pull(n) %>%
    mapply(rep, x=1:N, times=., SIMPLIFY = FALSE) %>%
    unlist
  
  beta_V_idx = split(
    beta_fit$data$time_order[!is.na(beta_fit$data$value_rtb_chg)], # change from baseline
    # beta_fit$data$time_order[!is.na(beta_fit$data$value_rtb)], # raw outcome
    beta_nonmissing_factor
  ) # a list, each sublist records the time orders of nonmissing measurements for each subject with at least 1 post-baseline measurement
  beta_S_idx = split(
    1:length(beta_nonmissing_factor), 
    beta_nonmissing_factor
  ) # a list, each sublist records the row indices of nonmissing measurements for each subject with at least 1 post-baseline measurement
  
  beta_D = as.matrix(modelr::model_matrix(beta_fit$data, beta_fit)) * 
    as.vector(beta_fit$family$mu.eta(beta_fit$linear.predictors)) # [K*(N_ctrl+N_drug)]*p_beta if no missing
  # mu.eta: if the inverse-link function is mu = ginv(eta) where eta is the value of the linear predictor, then this function returns d(ginv(eta))/d(eta) = d(mu)/d(eta).
  beta_A_sqrt = beta_fit$fitted.values %>%
    beta_fit$family$variance() %>%
    as.vector() %>%
    split(x = ., f = 1:length(.)) %>%
    bdiag %>%
    sqrt
  beta_V_inv = solve(
    beta_fit$geese$gamma * # dispersion parameter
      beta_A_sqrt %*% 
      bdiag(lapply(beta_V_idx, function(x){beta_R[x,x]})) %*%
      beta_A_sqrt + diag(nrow(beta_A_sqrt))*eps # add a diagonal matrix (optional) to avoid singular matrix
  ) # (\phi*A^(1/2)RA^(1/2))^(-1), [K*(N_ctrl+N_drug)]*[K*(N_ctrl+N_drug)] if no missing
  beta_S = bdiag(lapply(beta_S_idx, function(x){beta_residual[x]})) # [K*(N_ctrl+N_drug)]*[N_ctrl+N_drug] if non missing
  
  ## estimating equations
  beta_U = matrix(0, nrow = p_beta, ncol = N)
  beta_U[,unique(beta_nonmissing_factor)] = as.matrix(t(beta_D) %*% beta_V_inv %*% beta_S) # dgCMatrix can't be directly set to a matrix
  
  ## information matrix
  eps = 0
  beta_I = - t(beta_D) %*% beta_V_inv %*% beta_D + eps*diag(p_beta) # in case this matrix is singular
  
  return(list(coef=beta_coef, U=beta_U, I=beta_I, fit=beta_fit))
}

fit_data_RTB_pi = function(data, K){
  ## Function for estimation of missingness
  
  pi_df = data %>% 
    filter(time_order == K)
  
  pi_R = 1-pi_df$miss # missing(=0) or not(=1)
  pi_T = rep(1, length(pi_R)) # treatment indicator
  
  ## coefficients
  pi_coef = mean(pi_R)
  
  ## estimating equations
  pi_U = pi_T * (pi_R - pi_coef) # change order will change the result
  
  ## information matrix
  pi_I = sum(-pi_T)
  
  return(list(coef=pi_coef, U=pi_U, I=pi_I))
}

fit_data_RTB_nu = function(data, K){
  ## Function for estimation of baseline covariates
  
  nu_df = data %>% 
    filter(time_order == K)
  
  nu_X = nu_df$base # observed baseline covariate
  nu_T = rep(1, length(nu_X)) # treatment indicator
  
  # coefficients
  nu_coef = mean(nu_X)
  
  # estimating equations
  nu_U = nu_T * (nu_X - nu_coef) # change order will change the result
  
  # information matrix
  nu_I = sum(-nu_T)
  
  return(list(coef=nu_coef, U=nu_U, I=nu_I))
}

calc_est_RTB_delta_v2 = function(beta_result, pi_result, nu_result){
  ## Function for estimation of treatment effect
  
  ## Extract coefficients
  beta_coef = beta_result$coef
  pi_coef = pi_result$coef
  nu_coef = nu_result$coef
  
  ## Sandwich variance estimator for joint distribution of all parameters
  joint_U = rbind(beta_result$U, pi_result$U, nu_result$U)
  joint_I = bdiag(beta_result$I, pi_result$I, nu_result$I)
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  
  ## mu, observed response (linear predictors of beta, nu)
  l_beta = c(1, rep(0,K-2), 1, nu_coef, rep(0,K-2), nu_coef) # value_rtb_chg ~ factor(time)*base
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## delta, mu*pi
  delta_hat = mu_hat*pi_hat # change from baseline
  
  ## se(delta)
  deriv_delta_params = c(
    l_beta*pi_hat, # beta
    mu_hat, # pi
    sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*pi_hat # nu, value_rtb_chg ~ factor(time)*base
    # sum(beta_coef[c(paste0("factor(time)",K,":base"))])*pi_hat # nu, value_rtb_chg ~ factor(time):base
  )
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(pi_coef)), # pi
    rep(1, length(nu_coef)) # nu
  )
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    mu = mu_hat
  ))
}

calc_est_RTB_theta_v3 = function(delta_result_list, l_delta, mean_base, var_base){
  ## Function for estimation of treatment difference
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  var_delta_hat_joint = map(delta_result_list, "var") %>% bdiag
  cov_nu_delta_joint = map(delta_result_list, "cov_nu_delta") %>% bdiag
  var_nu_inv_joint = map(delta_result_list, "var_nu") %>% bdiag %>% solve
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  return(list(
    theta = list(
      est = theta_hat, var = var_theta_hat, se = se_theta_hat,
      est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
      est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
    ),
    delta = list(
      est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
      est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
      est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
    )
  ))
}

calc_est_RTB_delta_v3.1 = function(beta_result, pi_result, nu_result, tip_delta = 0){
  ## Function for estimation of treatment effect in tipping point analysis 

  ## Extract coefficients
  beta_coef = beta_result$coef
  pi_coef = pi_result$coef
  nu_coef = nu_result$coef
  
  ## Sandwich variance estimator for joint distribution of all parameters
  joint_U = rbind(beta_result$U, pi_result$U, nu_result$U)
  joint_I = bdiag(beta_result$I, pi_result$I, nu_result$I)
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  
  ## mu, observed response (linear predictors of beta, nu)
  l_beta = c(1, rep(0,K-2), 1, nu_coef, rep(0,K-2), nu_coef) # value_rtb_chg ~ factor(time)*base
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## delta, mu*pi
  delta_hat = mu_hat*pi_hat # change from baseline
  vec_coef_pi = c( # for the coefficient of phi_hat in deriv_delta_params
    l_beta, # beta
    0, # pi
    sum(beta_coef[c("base", paste0("factor(time)",K,":base"))]) # nu
  )
  vec_coef_constant = c( # for the coefficient not related to phi_hat in deriv_delta_params
    rep(0, length(l_beta)), # beta
    mu_hat,
    0 # nu
  )
  deriv_delta_params = vec_coef_pi*pi_hat + vec_coef_constant
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(pi_coef)), # pi
    rep(1, length(nu_coef)) # nu
  )
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*(1-pi_hat)
  deriv_delta_params_tip = deriv_delta_params + c(
    rep(0, length(beta_coef)), # beta
    -rep(tip_delta, length(pi_coef)), # pi
    rep(0, length(nu_coef)) # nu
  )
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip
  ))
}

calc_est_RTB_theta_v4 = function(delta_result_list, l_delta, mean_base, var_base, tip = FALSE, complete = TRUE){
  ## Function for estimation of treatment difference in tipping point analysis
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  var_delta_hat_joint = map(delta_result_list, "var") %>% bdiag
  cov_nu_delta_joint = map(delta_result_list, "cov_nu_delta") %>% bdiag
  var_nu_inv_joint = map(delta_result_list, "var_nu") %>% bdiag %>% solve
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  if (tip){ # tipping point analysis is required
    delta_hat_joint_tip = map(delta_result_list, "est_tip") %>% unlist
    var_delta_hat_joint_tip = map(delta_result_list, "var_tip") %>% bdiag
    
    theta_hat_tip = l_delta %*% delta_hat_joint_tip %>% as.vector
    var_theta_hat_tip = (l_delta %*% var_delta_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_hat_tip = sqrt(var_theta_hat_tip)
    
    delta_adj_hat_joint_tip = (delta_hat_joint_tip + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
      as.vector
    var_delta_adj_hat_joint_tip = var_delta_hat_joint_tip - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_hat_tip = l_delta %*% delta_adj_hat_joint_tip %>% as.vector
    var_theta_adj_hat_tip = (l_delta %*% var_delta_adj_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_hat_tip = sqrt(var_theta_adj_hat_tip)
    
    delta_adj_unc_hat_joint_tip = delta_adj_hat_joint_tip
    var_delta_adj_unc_hat_joint_tip = var_delta_adj_hat_joint_tip + 
      1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_unc_hat_tip = theta_adj_hat_tip
    var_theta_adj_unc_hat_tip = (l_delta %*% var_delta_adj_unc_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_unc_hat_tip = sqrt(var_theta_adj_unc_hat_tip)
    
    if (complete){ # complete result is required
      return(list(
        theta = list(
          est = theta_hat, var = var_theta_hat, se = se_theta_hat,
          est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
          est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
        ),
        delta = list(
          est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
          est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
          est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
        ),
        theta_tip = list(
          est = theta_hat_tip, var = var_theta_hat_tip, se = se_theta_hat_tip,
          est_adj = theta_adj_hat_tip, var_adj = var_theta_adj_hat_tip, se_adj = se_theta_adj_hat_tip,
          est_adj_unc = theta_adj_unc_hat_tip, var_adj_unc = var_theta_adj_unc_hat_tip, se_adj_unc = se_theta_adj_unc_hat_tip
        )
      ))
    } else { # pvalue only
      return(list(
        z_unadj = theta_hat_tip/se_theta_hat_tip,
        z_adj = theta_adj_hat_tip/se_theta_adj_hat_tip,
        z_adj_unc = theta_adj_unc_hat_tip/se_theta_adj_unc_hat_tip
      ))
    }
  } else { # tipping point analysis is not required
    return(list(
      theta = list(
        est = theta_hat, var = var_theta_hat, se = se_theta_hat,
        est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
        est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
      ),
      delta = list(
        est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
        est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
        est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
      )
    ))
  }
}

calc_est_RTB_TIP_v1 = function(
    beta_result, pi_result, nu_result, 
    l_delta, mean_base, var_base, tip = FALSE, tip_delta = NULL, complete = TRUE){
  ## general function for tipping point analysis
  
  if (tip & !complete){
    apply(tip_delta, 1, function(z){
      calc_est_RTB_theta_v4(
        delta_result_list = mapply( 
          FUN = calc_est_RTB_delta_v3.1,
          beta_result = beta_result, pi_result = pi_result, nu_result = nu_result, 
          tip_delta = as.list(z),
          SIMPLIFY = FALSE
        ),
        l_delta, mean_base, var_base, tip, complete
      ) %>% unlist
      # only consider complete = FALSE in the current version
      # a vector, each C elements are the p values of C treatment effects, where C is the number of rows in l_delta
    }) %>% t # matrix, each row corresponds to each set of tipping points
  }
}

# functions for J2R ####
fit_data_J2R_beta = function(data, N, K, beta_family, beta_corstr, eps = 0){
  ## Function for estimation of MMRM
  
  ## fit GEE, data must be sorted by subject, can check attributes(beta_fit$data)$groups to verify
  if (K >= 2){
    beta_fit = geepack::geeglm(
      formula = value_j2r_chg ~ factor(time)*base, # change from baseline
      id = subject,
      data = data, 
      waves = time_order, # waves argument handles missing value
      scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
      family = beta_family, 
      corstr = ifelse(K >= 3, beta_corstr, "ar1")
    )
  } else {
    beta_fit = geepack::geeglm(
      formula = value_j2r_chg ~ base, # change from baseline
      id = subject,
      data = data, 
      waves = time_order, # waves argument handles missing value
      scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
      family = beta_family, 
      # corstr = "unstructured"
      corstr = ifelse(K >= 3, beta_corstr, "ar1")
    )
  }
  # https://faculty.washington.edu/heagerty/Courses/b571/homework/geepack-paper.pdf
  # https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm
  
  beta_coef = coef(beta_fit) # extract coefficients
  p_beta = length(beta_coef) # dimension of coefficients
  beta_cov = vcov(beta_fit) # extract covariance matrix of coefficients, https://rdrr.io/cran/merDeriv/man/vcov.lmerMod.html
  beta_R = working_correlation_v2(K, beta_corstr, beta_fit$geese$alpha)
  beta_residual = beta_fit$y-as.vector(beta_fit$fitted.values) # beta_fit$fitted.values is a matrix, this value is different from beta_fit$residuals
  
  ## derive score equations for main effects
  # block diagnol matrix, https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/bdiag.html
  # convert matrix to list, https://r-lang.com/how-to-convert-r-matrix-to-list/#:~:text=How%20to%20Convert%20R%20Matrix%20to%20List%201,into%20the%20list%20of%20vectors.%203%20See%20also
  
  beta_nonmissing_factor = beta_fit$data %>%
    group_by(subject) %>%
    summarise(n=sum(!is.na(value_j2r_chg))) %>% #  change from baselina
    # summarise(n=sum(!is.na(value_j2r))) %>% # outcome
    pull(n) %>%
    mapply(rep, x=1:N, times=., SIMPLIFY = FALSE) %>%
    unlist
  
  beta_V_idx = split(
    beta_fit$data$time_order[!is.na(beta_fit$data$value_j2r_chg)], # change from baseline
    # beta_fit$data$time_order[!is.na(beta_fit$data$value_j2r)], # raw outcome
    beta_nonmissing_factor
  ) # a list, each sublist records the time orders of nonmissing measurements for each subject with at least 1 post-baseline measurement
  beta_S_idx = split(
    1:length(beta_nonmissing_factor), 
    beta_nonmissing_factor
  ) # a list, each sublist records the row indices of nonmissing measurements for each subject with at least 1 post-baseline measurement
  
  beta_D = as.matrix(modelr::model_matrix(beta_fit$data, beta_fit)) * 
    as.vector(beta_fit$family$mu.eta(beta_fit$linear.predictors)) # [K*(N_ctrl+N_drug)]*p_beta if no missing
  # mu.eta: if the inverse-link function is mu = ginv(eta) where eta is the value of the linear predictor, then this function returns d(ginv(eta))/d(eta) = d(mu)/d(eta).
  beta_A_sqrt = beta_fit$fitted.values %>%
    beta_fit$family$variance() %>%
    as.vector() %>%
    split(x = ., f = 1:length(.)) %>%
    bdiag %>%
    sqrt
  beta_V_inv = solve(
    beta_fit$geese$gamma * # dispersion parameter
      beta_A_sqrt %*% 
      bdiag(lapply(beta_V_idx, function(x){beta_R[x,x]})) %*%
      beta_A_sqrt + diag(nrow(beta_A_sqrt))*eps # add a diagonal matrix (optional) to avoid singular matrix
  ) # (\phi*A^(1/2)RA^(1/2))^(-1), [K*(N_ctrl+N_drug)]*[K*(N_ctrl+N_drug)] if no missing
  beta_S = bdiag(lapply(beta_S_idx, function(x){beta_residual[x]})) # [K*(N_ctrl+N_drug)]*[N_ctrl+N_drug] if non missing
  
  ## estimating equations
  beta_U = matrix(0, nrow = p_beta, ncol = N)
  beta_U[,unique(beta_nonmissing_factor)] = as.matrix(t(beta_D) %*% beta_V_inv %*% beta_S) # dgCMatrix can't be directly set to a matrix
  
  ## information matrix
  eps = 0
  beta_I = - t(beta_D) %*% beta_V_inv %*% beta_D + eps*diag(p_beta) # in case this matrix is singular
  
  return(list(coef=beta_coef, U=beta_U, I=beta_I, fit=beta_fit))
}

fit_data_J2R_pi = function(data, K){
  ## Function for estimation of missingness
  
  pi_df = data %>% 
    filter(time_order == K)
  
  pi_R = 1-pi_df$miss # missing(=0) or not(=1)
  pi_T = rep(1, length(pi_R)) # treatment indicator
  
  ## coefficients
  pi_coef = mean(pi_R)
  
  ## estimating equations
  pi_U = pi_T * (pi_R - pi_coef) # change order will change the result
  
  ## information matrix
  pi_I = sum(-pi_T)
  
  return(list(coef=pi_coef, U=pi_U, I=pi_I))
}

fit_data_J2R_nu = function(data, K){
  ## Function for estimation of baseline covariates
  
  nu_df = data %>% 
    filter(time_order == K)
  
  nu_X = nu_df$base # observed baseline covariate
  nu_T = rep(1, length(nu_X)) # treatment indicator
  
  # coefficients
  nu_coef = mean(nu_X)
  
  # estimating equations
  nu_U = nu_T * (nu_X - nu_coef) # change order will change the result
  
  # information matrix
  nu_I = sum(-nu_T)
  
  return(list(coef=nu_coef, U=nu_U, I=nu_I))
}

calc_est_J2R_delta_v8 = function(
    beta_result, betaREF_result, pi_result, 
    nuREF_result, nu_result, reference = FALSE, tip_delta = 0
){
  ## Function for estimation of treatment effect

  ## Extract coefficients
  beta_coef = beta_result$coef
  betaREF_coef = betaREF_result$coef; ncoef_betaREF = length(betaREF_coef)
  pi_coef = pi_result$coef
  nuREF_coef = nuREF_result$coef; ncoef_nuREF = length(nuREF_coef)
  nu_coef = nu_result$coef
  
  
  ## Sandwich variance estimator for joint distribution of all parameters
  if (reference){
    joint_U = rbind(
      betaREF_result$U, nuREF_result$U,
      beta_result$U, pi_result$U, nu_result$U
    )
  } else {
    nREF = ncol(betaREF_result$U)
    n = ncol(beta_result$U)
    joint_U = rbind(
      cbind(
        betaREF_result$U,
        matrix(0, nrow = ncoef_betaREF, ncol = nREF + n)
      ),
      c(rep(0, nREF), nuREF_result$U, rep(0, n)), 
      cbind(
        matrix(0, nrow = ncoef_betaREF, ncol = nREF*2),
        beta_result$U
      ), 
      c(rep(0, nREF), rep(0, nREF), pi_result$U), 
      # c(rep(0, nREF), nuMISS_result$U), 
      c(rep(0, nREF), rep(0, nREF), nu_result$U)
    )
  }
  joint_I = bdiag(
    betaREF_result$I, nuREF_result$I,
    beta_result$I, pi_result$I, nu_result$I
  )
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  joint_cov_shared = joint_cov[
    1:(ncoef_betaREF + ncoef_nuREF), 1:(ncoef_betaREF + ncoef_nuREF)
  ]
  joint_cov_distinct = joint_cov[
    -(1:(ncoef_betaREF + ncoef_nuREF)), -(1:(ncoef_betaREF + ncoef_nuREF))
  ]
  joint_cov_cov = joint_cov[
    -(1:(ncoef_betaREF + ncoef_nuREF)), 1:(ncoef_betaREF + ncoef_nuREF)
  ]
  
  ## mu, completer (linear predictors of muOBS)
  if (K >= 2){
    l_beta = c(
      1, rep(0,K-2), 1, nu_coef, # int, factor(time), base
      rep(0,K-2), nu_coef # factor(time)*base
    ) # value ~ factor(time)*base
  } else {
    l_beta = c(1, nu_coef) # ???, value_ad_rd_chg ~ factor(time)*base
  }
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## muMISS, subjects with missing outcome (linear predictors of betaREF, nuMISS)
  if (K >= 2){
    l_betaREF = c(1, rep(0,K-2), 1, nuREF_coef, 
                  rep(0,K-2), nuREF_coef)
  } else {
    l_betaREF = c(1, nuREF_coef)
  }
  muREF_hat = beta_result$fit$family$linkinv(sum(l_betaREF*betaREF_coef)) # can handle only continuous variable for now
  
  
  ## delta
  delta_hat = mu_hat*pi_hat + muREF_hat*(1-pi_hat) # change from baseline
  
  ## se(delta)
  if (K >= 2){
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base", paste0("factor(time)",K,":base"))]) # nu  
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muREF_hat, # pi
      0 # nu
    )
  } else {
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base")]) # nu
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muREF_hat, # pi
      0 # nu
    )
  }
  vec_coef_pi_shared = c(
    -l_betaREF, # betaREF
    -sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))]) # nuREF  
  )
  vec_coef_constant_shared = c(
    l_betaREF, # betaREF
    sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))]) # nuREF
  )
  # deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  # deriv_delta_params = c(vec_coef_pi_shared, vec_coef_pi_distinct) * pi_hat + c(vec_coef_constant_shaerd, vec_coef_pi_shared)
  deriv_delta_params_shared = vec_coef_pi_shared * pi_hat + vec_coef_constant_shared
  deriv_delta_params_distinct = vec_coef_pi_distinct * pi_hat + vec_coef_constant_distinct
  deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params_distinct = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(pi_coef)), # pi
    rep(1, length(nu_coef)) # nu
  )
  # deriv_nu_params_shared = c(
  #   rep(0, length(betaREF_coef)), # betaREF
  #   ifelse(reference, rep(1, length(nuREF_coef)), rep(0, length(nuREF_coef))) # nuREF
  # )
  deriv_nu_params_shared = c(
    rep(0, length(betaREF_coef)), # betaREF
    rep(0, length(nuREF_coef)) # nuREF
  )
  deriv_nu_params = c(deriv_nu_params_shared, deriv_nu_params_distinct)
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*pi_hat
  deriv_delta_params_distinct_tip = deriv_delta_params_distinct + c(
    rep(0, length(beta_coef)), # beta
    rep(tip_delta, length(pi_coef)), # pi
    rep(0, length(nu_coef)) # nu
  )
  deriv_delta_params_shared_tip = deriv_delta_params_shared
  deriv_delta_params_tip = c(deriv_delta_params_shared_tip, deriv_delta_params_distinct_tip)
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    cov_shared = joint_cov_shared, cov_distinct = joint_cov_distinct,
    cov = joint_cov, cov_cov = joint_cov_cov,
    deriv_delta_distinct = deriv_delta_params_distinct, deriv_delta_shared = deriv_delta_params_shared,
    deriv_nu_distinct = deriv_nu_params_distinct, deriv_nu_shared = deriv_nu_params_shared,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip, # tipping
    deriv_delta_distinct_tip = deriv_delta_params_distinct_tip, # tipping 
    deriv_delta_shared_tip = deriv_delta_params_shared_tip # tipping
  ))
}

calc_est_J2R_theta_v4 = function(delta_result_list, l_delta, mean_base, var_base){
  ## Function for estimation of treatment difference
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  deriv_delta_joint = rbind(
    map(delta_result_list, "deriv_delta_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_delta_distinct") %>% bdiag
  )
  
  cov_joint = rbind(
    cbind(
      delta_result_list[[1]]$cov_shared,
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.) %>% t
    ),
    cbind(
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.),
      map(delta_result_list, "cov_distinct") %>% bdiag
    )
  )
  
  deriv_nu_joint = rbind(
    map(delta_result_list, "deriv_nu_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_nu_distinct") %>% bdiag
  )
  var_nu_inv_joint = t(deriv_nu_joint) %*% cov_joint %*% deriv_nu_joint %>% solve
  cov_nu_delta_joint =  t(deriv_nu_joint) %*% cov_joint %*% deriv_delta_joint
  
  var_delta_hat_joint = t(deriv_delta_joint) %*% cov_joint %*% deriv_delta_joint
  
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  return(list(
    theta = list(
      est = theta_hat, var = var_theta_hat, se = se_theta_hat,
      est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
      est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
    ),
    delta = list(
      est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
      est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
      est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
    )
  ))
}

calc_est_J2R_delta_v8.2 = function(
    beta_result, betaREF_result, pi_result, 
    nuREF_result, nu_result, reference = FALSE, tip_delta = 0,
    base_var
){
  ## Function for estimation of treatment effect in tipping point analysis 

  ## Extract coefficients
  beta_coef = beta_result$coef
  betaREF_coef = betaREF_result$coef; ncoef_betaREF = length(betaREF_coef)
  pi_coef = pi_result$coef
  nuREF_coef = nuREF_result$coef; ncoef_nuREF = length(nuREF_coef)
  nu_coef = nu_result$coef
  
  
  ## Sandwich variance estimator for joint distribution of all parameters
  if (reference){
    joint_U = rbind(
      betaREF_result$U, nuREF_result$U,
      beta_result$U, pi_result$U, nu_result$U
    )
  } else {
    nREF = ncol(betaREF_result$U)
    n = ncol(beta_result$U)
    joint_U = rbind(
      cbind(
        betaREF_result$U,
        matrix(0, nrow = ncoef_betaREF, ncol = nREF + n)
      ),
      c(rep(0, nREF), nuREF_result$U, rep(0, n)), 
      cbind(
        matrix(0, nrow = ncoef_betaREF, ncol = nREF*2),
        beta_result$U
      ), 
      c(rep(0, nREF), rep(0, nREF), pi_result$U), 
      # c(rep(0, nREF), nuMISS_result$U), 
      c(rep(0, nREF), rep(0, nREF), nu_result$U)
    )
    # joint_U = rbind(
    #   cbind(
    #     betaREF_result$U,
    #     matrix(0, nrow = ncoef_betaREF, ncol = nREF + n)
    #   ),
    #   cbind(
    #     matrix(0, nrow = ncoef_nuREF, ncol = nREF),
    #     nuREF_result$U,
    #     matrix(0, nrow = ncoef_nuREF, ncol = n)
    #   ),
    #   cbind(
    #     matrix(0, nrow = ncoef_betaREF, ncol = nREF*2),
    #     beta_result$U
    #   ),
    #   c(rep(0, nREF), rep(0, nREF), pi_result$U),
    #   cbind(
    #     matrix(0, nrow = ncoef_nuREF, ncol = nREF*2),
    #     nu_result$U
    #   )
    # )
  }
  joint_I = bdiag(
    betaREF_result$I, nuREF_result$I,
    beta_result$I, pi_result$I, nu_result$I
  )
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  joint_cov_shared = joint_cov[
    1:(ncoef_betaREF + ncoef_nuREF), 1:(ncoef_betaREF + ncoef_nuREF)
  ]
  joint_cov_distinct = joint_cov[
    -(1:(ncoef_betaREF + ncoef_nuREF)), -(1:(ncoef_betaREF + ncoef_nuREF))
  ]
  joint_cov_cov = joint_cov[
    -(1:(ncoef_betaREF + ncoef_nuREF)), 1:(ncoef_betaREF + ncoef_nuREF)
  ]
  
  ## mu, completer (linear predictors of muOBS)
  if (K >= 2){
    l_beta = c(
      1, rep(0,K-2), 1, nu_coef, # int, factor(time), base, base1
      rbind(matrix(0, nrow = K-2, ncol = length(base_var)), nu_coef) %>% c # factor(time)*base, factor(time)*base1
    ) # value_chg ~ factor(time)*base + factor*base1
  } else {
    l_beta = c(1, nu_coef) # value_chg ~ factor(time)*base
  }
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## muMISS, subjects with missing outcome (linear predictors of betaREF, nuMISS)
  if (K >= 2){
    l_betaREF = c(
      1, rep(0,K-2), 1, nuREF_coef, # int, factor(time), base, base1
      rbind(matrix(0, nrow = K-2, ncol = length(base_var)), nuREF_coef) %>% c # factor(time)*base, factor(time)*base1
    ) # value ~ factor(time)*base + factor*base1
  } else {
    l_betaREF = c(1, nuREF_coef)
  }
  muREF_hat = beta_result$fit$family$linkinv(sum(l_betaREF*betaREF_coef)) # can handle only continuous variable for now
  
  
  ## delta
  delta_hat = mu_hat*pi_hat + muREF_hat*(1-pi_hat) # change from baseline
  
  ## se(delta)
  if (K >= 2){
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sapply(base_var, function(z){
        sum(beta_coef[c(z, paste0("factor(time)",K,":",z))]) # nu, length p_nu
      })
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muREF_hat, # pi
      rep(0, length(base_var)) # nu, length p_nu
    )
  } else {
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      beta_coef[base_var] # nu, length p_nu
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muREF_hat, # pi
      rep(0, length(base_var)) # nu, length p_nu
    )
  }
  vec_coef_pi_shared = c(
    -l_betaREF, # betaREF
    sapply(base_var, function(z){
      -sum(beta_coef[c(z, paste0("factor(time)",K,":",z))]) # nu, length p_nu
    })
  )
  vec_coef_constant_shared = c(
    l_betaREF, # betaREF
    sapply(base_var, function(z){
      sum(beta_coef[c(z, paste0("factor(time)",K,":",z))]) # nu, length p_nu
    })
  )
  # deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  # deriv_delta_params = c(vec_coef_pi_shared, vec_coef_pi_distinct) * pi_hat + c(vec_coef_constant_shaerd, vec_coef_pi_shared)
  deriv_delta_params_shared = vec_coef_pi_shared * pi_hat + vec_coef_constant_shared
  deriv_delta_params_distinct = vec_coef_pi_distinct * pi_hat + vec_coef_constant_distinct
  deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate, WRONG??? 10/4/2023
  deriv_nu_params_shared = matrix(0,
                                  nrow = length(betaREF_coef) + length(nuREF_coef), ncol = length(nuREF_coef)
  )
  deriv_nu_params_distinct = matrix(0,
                                    nrow = length(beta_coef) + length(pi_coef), ncol = length(nuREF_coef)
  ) %>%
    rbind(diag(length(nu_coef)))
  
  deriv_nu_params = rbind(deriv_nu_params_shared, deriv_nu_params_distinct)
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*(1-pi_hat)
  deriv_delta_params_distinct_tip = deriv_delta_params_distinct + c(
    rep(0, length(beta_coef)), # beta
    -rep(tip_delta, length(pi_coef)), # pi
    rep(0, length(nu_coef)) # nu
  )
  deriv_delta_params_shared_tip = deriv_delta_params_shared
  deriv_delta_params_tip = c(deriv_delta_params_shared_tip, deriv_delta_params_distinct_tip)
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    cov_shared = joint_cov_shared, cov_distinct = joint_cov_distinct,
    cov = joint_cov, cov_cov = joint_cov_cov,
    deriv_delta_distinct = deriv_delta_params_distinct, deriv_delta_shared = deriv_delta_params_shared,
    deriv_nu_distinct = deriv_nu_params_distinct, deriv_nu_shared = deriv_nu_params_shared,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip, # tipping
    deriv_delta_distinct_tip = deriv_delta_params_distinct_tip, # tipping 
    deriv_delta_shared_tip = deriv_delta_params_shared_tip # tipping
  ))
}

calc_est_J2R_theta_v5 = function(delta_result_list, l_delta, mean_base, var_base, tip = FALSE, complete = TRUE){
  ## Function for estimation of treatment difference in tipping point analysis
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  deriv_delta_joint = rbind(
    map(delta_result_list, "deriv_delta_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_delta_distinct") %>% bdiag
  )
  
  cov_joint = rbind(
    cbind(
      delta_result_list[[1]]$cov_shared,
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.) %>% t
    ),
    cbind(
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.),
      map(delta_result_list, "cov_distinct") %>% bdiag
    )
  )
  
  deriv_nu_joint = rbind(
    map(delta_result_list, "deriv_nu_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_nu_distinct") %>% bdiag
  )
  var_nu_inv_joint = t(deriv_nu_joint) %*% cov_joint %*% deriv_nu_joint %>% solve
  cov_nu_delta_joint =  t(deriv_nu_joint) %*% cov_joint %*% deriv_delta_joint
  
  var_delta_hat_joint = t(deriv_delta_joint) %*% cov_joint %*% deriv_delta_joint
  
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  if (tip){ # tipping point analysis is required
    delta_hat_joint_tip = map(delta_result_list, "est_tip") %>% unlist
    deriv_delta_joint_tip = rbind(
      map(delta_result_list, "deriv_delta_shared_tip") %>% do.call(cbind,.),
      map(delta_result_list, "deriv_delta_distinct_tip") %>% bdiag
    )
    var_delta_hat_joint_tip = t(deriv_delta_joint_tip) %*% cov_joint %*% deriv_delta_joint_tip
    
    theta_hat_tip = l_delta %*% delta_hat_joint_tip %>% as.vector
    var_theta_hat_tip = (l_delta %*% var_delta_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_hat_tip = sqrt(var_theta_hat_tip)
    
    delta_adj_hat_joint_tip = (delta_hat_joint_tip + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
      as.vector
    var_delta_adj_hat_joint_tip = var_delta_hat_joint_tip - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_hat_tip = l_delta %*% delta_adj_hat_joint_tip %>% as.vector
    var_theta_adj_hat_tip = (l_delta %*% var_delta_adj_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_hat_tip = sqrt(var_theta_adj_hat_tip)
    
    delta_adj_unc_hat_joint_tip = delta_adj_hat_joint_tip
    var_delta_adj_unc_hat_joint_tip = var_delta_adj_hat_joint_tip + 
      1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_unc_hat_tip = theta_adj_hat_tip
    var_theta_adj_unc_hat_tip = (l_delta %*% var_delta_adj_unc_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_unc_hat_tip = sqrt(var_theta_adj_unc_hat_tip)
    
    if (complete){ # complete result is required
      return(list(
        theta = list(
          est = theta_hat, var = var_theta_hat, se = se_theta_hat,
          est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
          est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
        ),
        delta = list(
          est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
          est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
          est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
        ),
        theta_tip = list(
          est = theta_hat_tip, var = var_theta_hat_tip, se = se_theta_hat_tip,
          est_adj = theta_adj_hat_tip, var_adj = var_theta_adj_hat_tip, se_adj = se_theta_adj_hat_tip,
          est_adj_unc = theta_adj_unc_hat_tip, var_adj_unc = var_theta_adj_unc_hat_tip, se_adj_unc = se_theta_adj_unc_hat_tip
        )
      ))
    } else { # pvalue only
      return(list(
        z_unadj = theta_hat_tip/se_theta_hat_tip,
        z_adj = theta_adj_hat_tip/se_theta_adj_hat_tip,
        z_adj_unc = theta_adj_unc_hat_tip/se_theta_adj_unc_hat_tip
      ))
    }
  } else { # tipping point analysis is not required
    return(list(
      theta = list(
        est = theta_hat, var = var_theta_hat, se = se_theta_hat,
        est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
        est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
      ),
      delta = list(
        est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
        est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
        est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
      )
    ))
  }
  
}


calc_est_J2R_TIP_v1 = function(
    beta_result, betaREF_result, pi_result, nuREF_result, nu_result, 
    l_delta, mean_base, var_base, tip = FALSE, tip_delta = NULL, complete = TRUE){
  ## general function for tipping point analysis
  
  if (tip & !complete){
    apply(tip_delta, 1, function(z){
      calc_est_J2R_theta_v5(
        delta_result_list = mapply( 
          FUN = calc_est_J2R_delta_v8.2,
          beta_result = beta_result, pi_result = pi_result, 
          nu_result = nu_result, tip_delta = as.list(z),
          MoreArgs = list(
            nuREF_result = nuREF_result, betaREF_result = betaREF_result, 
            reference = FALSE, base_var = "base"
          ),
          SIMPLIFY = FALSE
        ),
        l_delta, mean_base, var_base, tip, complete
      ) %>% unlist
      # only consider complete = FALSE in the current version
      # a vector, each C elements are the p values of C treatment effects, where C is the number of rows in l_delta
    }) %>% t # matrix, each row corresponds to each set of tipping points
  }
}

# functions for PW ####
fit_data_PW_beta = function(data, N, K, beta_family, beta_corstr, eps = 0){
  ## Function for estimation of MMRM
  
  ## fit GEE, data must be sorted by subject, can check attributes(beta_fit$data)$groups to verify
  if (K >= 2){
    beta_fit = geepack::geeglm(
      formula = value_pw_chg ~ factor(time)*base, # change from baseline
      id = subject,
      data = data, 
      waves = time_order, # waves argument handles missing value
      scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
      family = beta_family, 
      corstr = ifelse(K >= 3, beta_corstr, "ar1")
    )
  } else {
    beta_fit = geepack::geeglm(
      formula = value_pw_chg ~ base, # change from baseline
      id = subject,
      data = data, 
      waves = time_order, # waves argument handles missing value
      scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
      family = beta_family, 
      # corstr = "unstructured"
      corstr = ifelse(K >= 3, beta_corstr, "ar1")
    )
  }
  # https://faculty.washington.edu/heagerty/Courses/b571/homework/geepack-paper.pdf
  # https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm
  
  beta_coef = coef(beta_fit) # extract coefficients
  p_beta = length(beta_coef) # dimension of coefficients
  beta_cov = vcov(beta_fit) # extract covariance matrix of coefficients, https://rdrr.io/cran/merDeriv/man/vcov.lmerMod.html
  beta_R = working_correlation_v2(K, beta_corstr, beta_fit$geese$alpha)
  beta_residual = beta_fit$y-as.vector(beta_fit$fitted.values) # beta_fit$fitted.values is a matrix, this value is different from beta_fit$residuals
  
  ## derive score equations for main effects
  # block diagnol matrix, https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/bdiag.html
  # convert matrix to list, https://r-lang.com/how-to-convert-r-matrix-to-list/#:~:text=How%20to%20Convert%20R%20Matrix%20to%20List%201,into%20the%20list%20of%20vectors.%203%20See%20also
  
  beta_nonmissing_factor = beta_fit$data %>%
    group_by(subject) %>%
    summarise(n=sum(!is.na(value_pw_chg))) %>% #  change from baselina
    # summarise(n=sum(!is.na(value_pw))) %>% # outcome
    pull(n) %>%
    mapply(rep, x=1:N, times=., SIMPLIFY = FALSE) %>%
    unlist
  
  beta_V_idx = split(
    beta_fit$data$time_order[!is.na(beta_fit$data$value_pw_chg)], # change from baseline
    # beta_fit$data$time_order[!is.na(beta_fit$data$value_pw)], # raw outcome
    beta_nonmissing_factor
  ) # a list, each sublist records the time orders of nonmissing measurements for each subject with at least 1 post-baseline measurement
  beta_S_idx = split(
    1:length(beta_nonmissing_factor), 
    beta_nonmissing_factor
  ) # a list, each sublist records the row indices of nonmissing measurements for each subject with at least 1 post-baseline measurement
  
  beta_D = as.matrix(modelr::model_matrix(beta_fit$data, beta_fit)) * 
    as.vector(beta_fit$family$mu.eta(beta_fit$linear.predictors)) # [K*(N_ctrl+N_drug)]*p_beta if no missing
  # mu.eta: if the inverse-link function is mu = ginv(eta) where eta is the value of the linear predictor, then this function returns d(ginv(eta))/d(eta) = d(mu)/d(eta).
  beta_A_sqrt = beta_fit$fitted.values %>%
    beta_fit$family$variance() %>%
    as.vector() %>%
    split(x = ., f = 1:length(.)) %>%
    bdiag %>%
    sqrt
  beta_V_inv = solve(
    beta_fit$geese$gamma * # dispersion parameter
      beta_A_sqrt %*% 
      bdiag(lapply(beta_V_idx, function(x){beta_R[x,x]})) %*%
      beta_A_sqrt + diag(nrow(beta_A_sqrt))*eps # add a diagonal matrix (optional) to avoid singular matrix
  ) # (\phi*A^(1/2)RA^(1/2))^(-1), [K*(N_ctrl+N_drug)]*[K*(N_ctrl+N_drug)] if no missing
  beta_S = bdiag(lapply(beta_S_idx, function(x){beta_residual[x]})) # [K*(N_ctrl+N_drug)]*[N_ctrl+N_drug] if non missing
  
  ## estimating equations
  beta_U = matrix(0, nrow = p_beta, ncol = N)
  beta_U[,unique(beta_nonmissing_factor)] = as.matrix(t(beta_D) %*% beta_V_inv %*% beta_S) # dgCMatrix can't be directly set to a matrix
  
  ## information matrix
  eps = 0
  beta_I = - t(beta_D) %*% beta_V_inv %*% beta_D + eps*diag(p_beta) # in case this matrix is singular
  
  return(list(coef=beta_coef, U=beta_U, I=beta_I, fit=beta_fit))
}

fit_data_PW_pi = function(data, K){
  ## Function for estimation of missingness
  
  pi_df = data %>% 
    filter(time_order == K)
  
  pi_R = 1-pi_df$miss # missing(=0) or not(=1)
  pi_T = rep(1, length(pi_R)) # treatment indicator
  
  ## coefficients
  pi_coef = mean(pi_R)
  
  ## estimating equations
  pi_U = pi_T * (pi_R - pi_coef) # change order will change the result
  
  ## information matrix
  pi_I = sum(-pi_T)
  
  return(list(coef=pi_coef, U=pi_U, I=pi_I))
}

fit_data_PW_nuMISS = function(data, K){
  ## Function for estimation of baseline covariates for subjects with missing outcomes
  
  nuMISS_df = data %>% 
    filter(time_order == K)
  
  nuMISS_X = nuMISS_df$base # observed baseline covariates
  nuMISS_TR = nuMISS_df$miss   # missing indicator
  
  # coefficients
  nuMISS_coef = sum(nuMISS_X * nuMISS_TR)/sum(nuMISS_TR)
  
  # estimating equations
  nuMISS_U = nuMISS_TR * (nuMISS_X - nuMISS_coef) # change order will change the result
  
  # information matrix
  nuMISS_I = sum(-nuMISS_TR)
  
  return(list(coef=nuMISS_coef, U=nuMISS_U, I=nuMISS_I))
}

fit_data_PW_nu = function(data, K){
  ## Function for estimation of baseline covariates
  
  nu_df = data %>% 
    filter(time_order == K)
  
  nu_X = nu_df$base # observed baseline covariate
  nu_T = rep(1, length(nu_X)) # treatment indicator
  
  # coefficients
  nu_coef = mean(nu_X)
  
  # estimating equations
  nu_U = nu_T * (nu_X - nu_coef) # change order will change the result
  
  # information matrix
  nu_I = sum(-nu_T)
  
  return(list(coef=nu_coef, U=nu_U, I=nu_I))
}

calc_est_PW_delta_v1 = function(beta_result, betaREF_result, pi_result, nuMISS_result, nu_result, reference = FALSE, tip_delta = 0){
  ## Function for estimation of treatment effect
  
  ## Extract coefficients
  beta_coef = beta_result$coef
  betaREF_coef = betaREF_result$coef
  pi_coef = pi_result$coef
  nuMISS_coef = nuMISS_result$coef
  nu_coef = nu_result$coef
  
  
  ## Sandwich variance estimator for joint distribution of all parameters
  if (reference){
    joint_U = rbind(betaREF_result$U, beta_result$U, pi_result$U, nuMISS_result$U, nu_result$U)
  } else {
    nREF = ncol(betaREF_result$U)
    joint_U = rbind(
      cbind(betaREF_result$U, matrix(0, nrow = length(betaREF_coef), ncol = ncol(beta_result$U))), 
      cbind(matrix(0, nrow = length(beta_coef), ncol = nREF), beta_result$U), 
      c(rep(0, nREF), pi_result$U), 
      c(rep(0, nREF), nuMISS_result$U), 
      c(rep(0, nREF), nu_result$U)
    )
  }
  joint_I = bdiag(betaREF_result$I, beta_result$I, pi_result$I, nuMISS_result$I, nu_result$I)
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  joint_cov_shared = joint_cov[1:length(betaREF_coef), 1:length(betaREF_coef)]
  joint_cov_distinct = joint_cov[-(1:length(betaREF_coef)), -(1:length(betaREF_coef))]
  joint_cov_cov = joint_cov[-(1:length(betaREF_coef)), 1:length(betaREF_coef)]
  
  ## mu, completer (linear predictors of muOBS)
  l_nuOBS = c(-nu_result$I, nuMISS_result$I)/sum(c(-nu_result$I, nuMISS_result$I))
  if (K >= 2){
    l_beta = c(1, rep(0,K-2), 1, sum(l_nuOBS*c(nu_coef, nuMISS_coef)), 
               rep(0,K-2), sum(l_nuOBS*c(nu_coef, nuMISS_coef)))
  } else {
    l_beta = c(1, sum(l_nuOBS*c(nu_coef, nuMISS_coef)))
  }
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## muMISS, subjects with missing outcome (linear predictors of betaREF, nuMISS)
  l_nuMISS = c(0, 1)
  if (K >= 2){
    l_betaREF = c(1, rep(0,K-2), 1, sum(l_nuMISS*c(nu_coef, nuMISS_coef)), 
                  rep(0,K-2), sum(l_nuMISS*c(nu_coef, nuMISS_coef)))
  } else {
    l_betaREF = c(1, sum(l_nuMISS*c(nu_coef, nuMISS_coef)))
  }
  muMISS_hat = beta_result$fit$family$linkinv(sum(l_betaREF*betaREF_coef)) # can handle only continuous variable for now
  
  ## delta
  delta_hat = mu_hat*pi_hat + muMISS_hat*(1-pi_hat) # change from baseline
  
  ## se(delta)
  if (K >= 2){
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuOBS[2] -  
        sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[2],
      sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuOBS[1] +
        sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[1]
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muMISS_hat, # pi
      sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[2],
      0
    )
  } else {
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base")])*l_nuOBS[2] - 
        sum(betaREF_coef[c("base")])*l_nuMISS[2],
      sum(beta_coef[c("base")])*l_nuOBS[1] +
        sum(betaREF_coef[c("base")])*l_nuMISS[1]
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muMISS_hat, # pi
      sum(betaREF_coef[c("base")])*l_nuMISS[2],
      0
    )
    
  }
  vec_coef_pi_shared = c(
    -l_betaREF # betaREF
  )
  vec_coef_constant_shared = c(
    l_betaREF # betaREF
  )
  deriv_delta_params_shared = vec_coef_pi_shared * pi_hat + vec_coef_constant_shared
  deriv_delta_params_distinct = vec_coef_pi_distinct * pi_hat + vec_coef_constant_distinct
  deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params_distinct = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(pi_coef)), # pi
    rep(0, length(nuMISS_coef)), # nuMISS
    rep(1, length(nu_coef)) # nu
  )
  deriv_nu_params_shared = c(
    rep(0, length(betaREF_coef)) # betaREF
  )
  deriv_nu_params = c(deriv_nu_params_shared, deriv_nu_params_distinct)
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*pi_hat
  # deriv_delta_params_tip = deriv_delta_params + c(
  #   rep(0, length(beta_coef)), # beta
  #   rep(tip_delta, length(pi_coef)), # pi
  #   rep(0, length(nuMISS_coef)), # nuMISS
  #   rep(0, length(nu_coef)) # nu
  # )
  deriv_delta_params_distinct_tip = deriv_delta_params_distinct + c(
    rep(0, length(beta_coef)), # beta
    rep(tip_delta, length(pi_coef)), # pi
    rep(0, length(nuMISS_coef)), # nuMISS
    rep(0, length(nu_coef)) # nu
  )
  deriv_delta_params_shared_tip = deriv_delta_params_shared
  deriv_delta_params_tip = c(deriv_delta_params_shared_tip, deriv_delta_params_distinct_tip)
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    cov_shared = joint_cov_shared, cov_distinct = joint_cov_distinct,
    cov = joint_cov, cov_cov = joint_cov_cov,
    deriv_delta_distinct = deriv_delta_params_distinct, deriv_delta_shared = deriv_delta_params_shared,
    deriv_nu_distinct = deriv_nu_params_distinct, deriv_nu_shared = deriv_nu_params_shared,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip, # tipping
    deriv_delta_distinct_tip = deriv_delta_params_distinct_tip, # tipping 
    deriv_delta_shared_tip = deriv_delta_params_shared_tip # tipping
  ))
}

calc_est_PW_theta_v1 = function(delta_result_list, l_delta, mean_base, var_base, tip = FALSE, complete = TRUE){
  ## Function for estimation of treatment difference
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  deriv_delta_joint = rbind(
    map(delta_result_list, "deriv_delta_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_delta_distinct") %>% bdiag
  )
  
  cov_joint = rbind(
    cbind(
      delta_result_list[[1]]$cov_shared,
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.) %>% t
    ),
    cbind(
      map(delta_result_list, "cov_cov") %>% do.call(rbind,.),
      map(delta_result_list, "cov_distinct") %>% bdiag
    )
  )
  
  deriv_nu_joint = rbind(
    map(delta_result_list, "deriv_nu_shared") %>% do.call(cbind,.),
    map(delta_result_list, "deriv_nu_distinct") %>% bdiag
  )
  var_nu_inv_joint = t(deriv_nu_joint) %*% cov_joint %*% deriv_nu_joint %>% solve
  cov_nu_delta_joint =  t(deriv_nu_joint) %*% cov_joint %*% deriv_delta_joint
  
  var_delta_hat_joint = t(deriv_delta_joint) %*% cov_joint %*% deriv_delta_joint
  
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  if (tip){ # tipping point analysis is required
    delta_hat_joint_tip = map(delta_result_list, "est_tip") %>% unlist
    deriv_delta_joint_tip = rbind(
      map(delta_result_list, "deriv_delta_shared_tip") %>% do.call(cbind,.),
      map(delta_result_list, "deriv_delta_distinct_tip") %>% bdiag
    )
    var_delta_hat_joint_tip = t(deriv_delta_joint_tip) %*% cov_joint %*% deriv_delta_joint_tip
    
    theta_hat_tip = l_delta %*% delta_hat_joint_tip %>% as.vector
    var_theta_hat_tip = (l_delta %*% var_delta_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_hat_tip = sqrt(var_theta_hat_tip)
    
    delta_adj_hat_joint_tip = (delta_hat_joint_tip + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
      as.vector
    var_delta_adj_hat_joint_tip = var_delta_hat_joint_tip - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_hat_tip = l_delta %*% delta_adj_hat_joint_tip %>% as.vector
    var_theta_adj_hat_tip = (l_delta %*% var_delta_adj_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_hat_tip = sqrt(var_theta_adj_hat_tip)
    
    delta_adj_unc_hat_joint_tip = delta_adj_hat_joint_tip
    var_delta_adj_unc_hat_joint_tip = var_delta_adj_hat_joint_tip + 
      1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_unc_hat_tip = theta_adj_hat_tip
    var_theta_adj_unc_hat_tip = (l_delta %*% var_delta_adj_unc_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_unc_hat_tip = sqrt(var_theta_adj_unc_hat_tip)
    
    if (complete){ # complete result is required
      return(list(
        theta = list(
          est = theta_hat, var = var_theta_hat, se = se_theta_hat,
          est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
          est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
        ),
        delta = list(
          est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
          est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
          est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
        ),
        theta_tip = list(
          est = theta_hat_tip, var = var_theta_hat_tip, se = se_theta_hat_tip,
          est_adj = theta_adj_hat_tip, var_adj = var_theta_adj_hat_tip, se_adj = se_theta_adj_hat_tip,
          est_adj_unc = theta_adj_unc_hat_tip, var_adj_unc = var_theta_adj_unc_hat_tip, se_adj_unc = se_theta_adj_unc_hat_tip
        )
      ))
    } else { # pvalue only
      return(list(
        z_unadj = theta_hat_tip/se_theta_hat_tip,
        z_adj = theta_adj_hat_tip/se_theta_adj_hat_tip,
        z_adj_unc = theta_adj_unc_hat_tip/se_theta_adj_unc_hat_tip
      ))
    }
  } else { # tipping point analysis is not required
    return(list(
      theta = list(
        est = theta_hat, var = var_theta_hat, se = se_theta_hat,
        est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
        est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
      ),
      delta = list(
        est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
        est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
        est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
      )
    ))
  }
  
}

calc_est_PW_delta_v1.1 = function(beta_result, betaREF_result, pi_result, nuMISS_result, nu_result, reference = FALSE, tip_delta = 0){
  ## Function for estimation of treatment effect in tipping point analysis 
  
  ## Extract coefficients
  beta_coef = beta_result$coef
  betaREF_coef = betaREF_result$coef
  pi_coef = pi_result$coef
  nuMISS_coef = nuMISS_result$coef
  nu_coef = nu_result$coef
  
  
  ## Sandwich variance estimator for joint distribution of all parameters
  if (reference){
    joint_U = rbind(betaREF_result$U, beta_result$U, pi_result$U, nuMISS_result$U, nu_result$U)
  } else {
    nREF = ncol(betaREF_result$U)
    joint_U = rbind(
      cbind(betaREF_result$U, matrix(0, nrow = length(betaREF_coef), ncol = ncol(beta_result$U))), 
      cbind(matrix(0, nrow = length(beta_coef), ncol = nREF), beta_result$U), 
      c(rep(0, nREF), pi_result$U), 
      c(rep(0, nREF), nuMISS_result$U), 
      c(rep(0, nREF), nu_result$U)
    )
  }
  joint_I = bdiag(betaREF_result$I, beta_result$I, pi_result$I, nuMISS_result$I, nu_result$I)
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  joint_cov_shared = joint_cov[1:length(betaREF_coef), 1:length(betaREF_coef)]
  joint_cov_distinct = joint_cov[-(1:length(betaREF_coef)), -(1:length(betaREF_coef))]
  joint_cov_cov = joint_cov[-(1:length(betaREF_coef)), 1:length(betaREF_coef)]
  
  ## mu, completer (linear predictors of muOBS)
  l_nuOBS = c(-nu_result$I, nuMISS_result$I)/sum(c(-nu_result$I, nuMISS_result$I))
  if (K >= 2){
    l_beta = c(1, rep(0,K-2), 1, sum(l_nuOBS*c(nu_coef, nuMISS_coef)), 
               rep(0,K-2), sum(l_nuOBS*c(nu_coef, nuMISS_coef)))
  } else {
    l_beta = c(1, sum(l_nuOBS*c(nu_coef, nuMISS_coef)))
  }
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## pi, missingness (linear predictors of pi)
  pi_hat = pi_coef
  
  ## muMISS, subjects with missing outcome (linear predictors of betaREF, nuMISS)
  l_nuMISS = c(0, 1)
  if (K >= 2){
    l_betaREF = c(1, rep(0,K-2), 1, sum(l_nuMISS*c(nu_coef, nuMISS_coef)), 
                  rep(0,K-2), sum(l_nuMISS*c(nu_coef, nuMISS_coef)))
  } else {
    l_betaREF = c(1, sum(l_nuMISS*c(nu_coef, nuMISS_coef)))
  }
  muMISS_hat = beta_result$fit$family$linkinv(sum(l_betaREF*betaREF_coef)) # can handle only continuous variable for now
  
  ## delta
  delta_hat = mu_hat*pi_hat + muMISS_hat*(1-pi_hat) # change from baseline
  
  ## se(delta)
  if (K >= 2){
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuOBS[2] -  
        sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[2], # nuMISS, use mu^{R=0}
      sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuOBS[1] +
        sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[1] # nu, value_rtb_chg ~ factor(time)*base
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muMISS_hat, # pi
      sum(betaREF_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuMISS[2], # nuMISS, use mu^{R=0}
      0 # nu, value_rtb_chg ~ factor(time)*base
    )
  } else {
    vec_coef_pi_distinct = c(
      l_beta, # beta
      0, # pi
      sum(beta_coef[c("base")])*l_nuOBS[2] - 
        sum(betaREF_coef[c("base")])*l_nuMISS[2],
      sum(beta_coef[c("base")])*l_nuOBS[1] +
        sum(betaREF_coef[c("base")])*l_nuMISS[1]
    )
    vec_coef_constant_distinct = c(
      rep(0, length(l_beta)), # beta
      mu_hat - muMISS_hat, # pi
      sum(betaREF_coef[c("base")])*l_nuMISS[2],
      0
    )
    
  }
  vec_coef_pi_shared = c(
    -l_betaREF # betaREF
  )
  vec_coef_constant_shared = c(
    l_betaREF # betaREF
  )
  deriv_delta_params_shared = vec_coef_pi_shared * pi_hat + vec_coef_constant_shared
  deriv_delta_params_distinct = vec_coef_pi_distinct * pi_hat + vec_coef_constant_distinct
  deriv_delta_params = c(deriv_delta_params_shared, deriv_delta_params_distinct)
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params_distinct = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(pi_coef)), # pi
    rep(0, length(nuMISS_coef)), # nuMISS
    rep(1, length(nu_coef)) # nu
  )
  deriv_nu_params_shared = c(
    rep(0, length(betaREF_coef)) # betaREF
  )
  deriv_nu_params = c(deriv_nu_params_shared, deriv_nu_params_distinct)
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*(1-pi_hat)
  deriv_delta_params_distinct_tip = deriv_delta_params_distinct + c(
    rep(0, length(beta_coef)), # beta
    rep(tip_delta, length(pi_coef)), # pi
    rep(0, length(nuMISS_coef)), # nuMISS
    rep(0, length(nu_coef)) # nu
  )
  deriv_delta_params_shared_tip = deriv_delta_params_shared
  deriv_delta_params_tip = c(deriv_delta_params_shared_tip, deriv_delta_params_distinct_tip)
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    cov_shared = joint_cov_shared, cov_distinct = joint_cov_distinct,
    cov = joint_cov, cov_cov = joint_cov_cov,
    deriv_delta_distinct = deriv_delta_params_distinct, deriv_delta_shared = deriv_delta_params_shared,
    deriv_nu_distinct = deriv_nu_params_distinct, deriv_nu_shared = deriv_nu_params_shared,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip, # tipping
    deriv_delta_distinct_tip = deriv_delta_params_distinct_tip, # tipping 
    deriv_delta_shared_tip = deriv_delta_params_shared_tip # tipping
  ))
}

calc_est_PW_TIP_v1 = function(
    beta_result, betaREF_result, 
    pi_result, nuMISS_result, nu_result, 
    l_delta, mean_base, var_base, tip = FALSE, tip_delta = NULL, complete = TRUE
){
  ## general function for tipping point analysis
  
  if (tip & !complete){
    apply(tip_delta, 1, function(z){
      calc_est_PW_theta_v1(
        delta_result_list = mapply( 
          FUN = calc_est_PW_delta_v1.1,
          beta_result = beta_result, pi_result = pi_result, 
          nu_result = nu_result, nuMISS_result = nuMISS_result, tip_delta = as.list(z),
          MoreArgs = list(
            betaREF_result = betaREF_result, 
            reference = FALSE
          ),
          SIMPLIFY = FALSE
        ),
        l_delta, mean_base, var_base, tip, complete
      ) %>% unlist
      # only consider complete = FALSE in the current version
      # a vector, each C elements are the p values of C treatment effects, where C is the number of rows in l_delta
    }) %>% t # matrix, each row corresponds to each set of tipping points
  }
}

# functions for RD ####
fit_data_RD_beta = function(data, N, K, beta_family, beta_corstr, eps = 0){## for outcome of interest
  ## Function for estimation of MMRM
  
  ## fit GEE, data must be sorted by subject, can check attributes(beta_fit$data)$groups to verify
  beta_fit = geepack::geeglm(
    formula = value_ad_rd_chg ~ factor(time)*base, # change from baseline
    id = subject,
    data = data, 
    waves = time_order, # waves argument handles missing value
    scale.fix = T, # T or F depends on the type of response, https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture22.htm
    family = beta_family, 
    # corstr = "unstructured"
    corstr = ifelse(K >= 3, beta_corstr, "ar1")
  )
  # https://faculty.washington.edu/heagerty/Courses/b571/homework/geepack-paper.pdf
  # https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm
  
  beta_coef = coef(beta_fit) # extract coefficients
  p_beta = length(beta_coef) # dimension of coefficients
  beta_cov = vcov(beta_fit) # extract covariance matrix of coefficients, https://rdrr.io/cran/merDeriv/man/vcov.lmerMod.html
  beta_R = working_correlation_v2(K, beta_corstr, beta_fit$geese$alpha)
  beta_residual = beta_fit$y-as.vector(beta_fit$fitted.values) # beta_fit$fitted.values is a matrix, this value is different from beta_fit$residuals
  
  ## derive score equations for main effects
  # block diagnol matrix, https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/bdiag.html
  # convert matrix to list, https://r-lang.com/how-to-convert-r-matrix-to-list/#:~:text=How%20to%20Convert%20R%20Matrix%20to%20List%201,into%20the%20list%20of%20vectors.%203%20See%20also
  
  beta_nonmissing_factor = beta_fit$data %>%
    group_by(subject) %>%
    summarise(n=sum(!is.na(value_ad_rd_chg))) %>% #  change from baselina
    pull(n) %>%
    mapply(rep, x=1:N, times=., SIMPLIFY = FALSE) %>%
    unlist
  
  beta_V_idx = split(
    beta_fit$data$time_order[!is.na(beta_fit$data$value_ad_rd_chg)], # change from baseline
    beta_nonmissing_factor
  ) # a list, each sublist records the time orders of nonmissing measurements for each subject with at least 1 post-baseline measurement
  beta_S_idx = split(
    1:length(beta_nonmissing_factor), 
    beta_nonmissing_factor
  ) # a list, each sublist records the row indices of nonmissing measurements for each subject with at least 1 post-baseline measurement
  
  beta_D = as.matrix(modelr::model_matrix(beta_fit$data, beta_fit)) * 
    as.vector(beta_fit$family$mu.eta(beta_fit$linear.predictors)) # [K*(N_ctrl+N_drug)]*p_beta if no missing
  # mu.eta: if the inverse-link function is mu = ginv(eta) where eta is the value of the linear predictor, then this function returns d(ginv(eta))/d(eta) = d(mu)/d(eta).
  beta_A_sqrt = beta_fit$fitted.values %>%
    beta_fit$family$variance() %>%
    as.vector() %>%
    split(x = ., f = 1:length(.)) %>%
    bdiag %>%
    sqrt
  beta_V_inv = solve(
    beta_fit$geese$gamma * # dispersion parameter
      beta_A_sqrt %*% 
      bdiag(lapply(beta_V_idx, function(x){beta_R[x,x]})) %*%
      beta_A_sqrt + diag(nrow(beta_A_sqrt))*eps # add a diagonal matrix (optional) to avoid singular matrix
  ) # (\phi*A^(1/2)RA^(1/2))^(-1), [K*(N_ctrl+N_drug)]*[K*(N_ctrl+N_drug)] if no missing
  beta_S = bdiag(lapply(beta_S_idx, function(x){beta_residual[x]})) # [K*(N_ctrl+N_drug)]*[N_ctrl+N_drug] if non missing
  
  ## estimating equations
  beta_U = matrix(0, nrow = p_beta, ncol = N)
  beta_U[,unique(beta_nonmissing_factor)] = as.matrix(t(beta_D) %*% beta_V_inv %*% beta_S) # dgCMatrix can't be directly set to a matrix
  
  ## information matrix
  eps = 0
  beta_I = - t(beta_D) %*% beta_V_inv %*% beta_D + eps*diag(p_beta) # in case this matrix is singular
  
  return(list(coef=beta_coef, U=beta_U, I=beta_I, fit=beta_fit))
}

fit_data_RD_betaRD = function(data, K){
  ## Function for estimation of beta^{-}, for retrieved dropout

  betaRD_df = data %>% 
    filter(time_order == K)
  
  betaRD_X = betaRD_df$base %>%
    cbind(1, .) %>%
    t() # a matrix of dimension 2*N, each column for each subject
  
  betaRD_TAR = as.numeric(betaRD_df$adherence == 0 & betaRD_df$miss == 0) # a vector of retrieved dropout indicator, each element for each subject
  
  # coefficients
  betaRD_coef = solve(
    a = betaRD_X %*% diag(betaRD_TAR) %*% t(betaRD_X), # dimension of [I*(p_base+1)]*[I*(p_base+1)]
    b = betaRD_X %*% (betaRD_TAR*ifelse(is.na(betaRD_df$value_rd_chg), 0, betaRD_df$value_rd_chg)) # change from baseline
  )
  
  # estimating equations
  betaRD_U = t(
    (t(betaRD_X) * ifelse(is.na(betaRD_df$value_rd_chg), 0, betaRD_df$value_rd_chg) - # change from baseline
       t(apply(betaRD_X, 2, function(x){x %*% t(x) %*% betaRD_coef}))
    ) * betaRD_TAR
  ) # current method can only handle one baseline covariate
  
  # information matrix
  betaRD_I = - betaRD_X %*% diag(betaRD_TAR) %*% t(betaRD_X)
  
  return(list(coef=betaRD_coef, U=betaRD_U, I=betaRD_I))
}

fit_data_RD_pi = function(data, K){
  ## Function for estimation of proportion of A=0 and R=0
  
  param_df = data %>% 
    filter(time_order == K)
  
  param_A = (1 - param_df$adherence) * param_df$miss # if A=0 and R=0(miss=1), param_A=1
  param_T = rep(1, length(param_A)) # treatment indicator
  
  # coefficients
  param_coef = mean(param_A)
  
  # estimating equations
  param_U = param_T * (param_A - param_coef) # change order will change the result
  
  # information matrix
  param_I = sum(-param_T)
  
  return(list(coef=param_coef, U=param_U, I=param_I))
}

fit_data_RD_tau = function(data, K){
  ## Function for estimation of proportion of A=1 and R=0
  
  param_df = data %>% 
    filter(time_order == K)
  
  param_A = param_df$adherence * param_df$miss # if A=1 and R=0(miss=1), param_A=1
  param_T = rep(1, length(param_A)) # treatment indicator
  
  # coefficients
  param_coef = mean(param_A)
  
  # estimating equations
  param_U = param_T * (param_A - param_coef) # change order will change the result
  
  # information matrix
  param_I = sum(-param_T)
  
  return(list(coef=param_coef, U=param_U, I=param_I))
}

fit_data_RD_tauAD = function(data, K){
  ## Function for estimation of proportion of A=1 and R=1
  
  param_df = data %>% 
    filter(time_order == K)
  
  param_A = param_df$adherence * (1 - param_df$miss) # if A=1 and R=1(miss=0), param_A=1
  param_T = rep(1, length(param_A)) # treatment indicator
  
  # coefficients
  param_coef = mean(param_A)
  
  # estimating equations
  param_U = param_T * (param_A - param_coef) # change order will change the result
  
  # information matrix
  param_I = sum(-param_T)
  
  return(list(coef=param_coef, U=param_U, I=param_I))
}

fit_data_RD_nuNAD = function(data, K){
  ## Function for estimation of baseline covariates of non-adherers

  nuNAD_df = data %>% 
    filter(time_order == K)
  
  nuNAD_X = nuNAD_df$base # observed baseline covariates
  nuNAD_TAR = 1 - nuNAD_df$adherence # non-adherers indicator
  
  # coefficients
  nuNAD_coef = sum(nuNAD_X * nuNAD_TAR)/sum(nuNAD_TAR)
  
  # estimating equations
  nuNAD_U = nuNAD_TAR * (nuNAD_X - nuNAD_coef) # change order will change the result
  
  # information matrix
  nuNAD_I = sum(-nuNAD_TAR)
  
  return(list(coef=nuNAD_coef, U=nuNAD_U, I=nuNAD_I))
}

fit_data_RD_nu = function(data, K){
  ## Function for estimation of baseline covariates

  nu_df = data %>% 
    filter(time_order == K)
  
  nu_X = nu_df$base # observed baseline covariate
  nu_T = rep(1, length(nu_X)) # treatment indicator
  
  # coefficients
  nu_coef = mean(nu_X)
  
  # estimating equations
  nu_U = nu_T * (nu_X - nu_coef) # change order will change the result
  
  # information matrix
  nu_I = sum(-nu_T)
  
  return(list(coef=nu_coef, U=nu_U, I=nu_I))
}

calc_est_RD_delta_v3 = function(
    beta_result, betaRD_result, 
    pi_result, tau_result, tauAD_result, 
    nuNAD_result, nu_result,
    tip_delta = 0){
  ## Function for estimation of treatment effect

  ## Extract coefficients
  beta_coef = beta_result$coef
  betaRD_coef = betaRD_result$coef
  pi_coef = pi_result$coef
  tau_coef = tau_result$coef
  tauAD_coef = tauAD_result$coef
  nuNAD_coef = nuNAD_result$coef
  nu_coef = nu_result$coef
  
  ## Sandwich variance estimator for joint distribution of all parameters
  joint_U = rbind(beta_result$U, betaRD_result$U, 
                  pi_result$U, tau_result$U, tauAD_result$U, 
                  nuNAD_result$U, nu_result$U)
  joint_I = bdiag(beta_result$I, betaRD_result$I, 
                  pi_result$I, tau_result$I, tauAD_result$I, 
                  nuNAD_result$I, nu_result$I)
  joint_cov = solve(joint_I) %*% joint_U %*% t(joint_U) %*% solve(joint_I) # covariance matrix of parameters
  
  ## mu, adherers (linear predictors of beta, nu, nuNAD)
  l_nuAD = c(-nu_result$I, nuNAD_result$I)/sum(c(-nu_result$I, nuNAD_result$I))
  l_beta = c(1, rep(0,K-2), 1, sum(l_nuAD*c(nu_coef, nuNAD_coef)), 
             rep(0,K-2), sum(l_nuAD*c(nu_coef, nuNAD_coef))) # ???, value_ad_rd_chg ~ factor(time)*base
  mu_hat = beta_result$fit$family$linkinv(sum(l_beta*beta_coef)) # can handle only continuous variable for now
  
  ## phi, adherence (linear predictors of phi)
  phi_hat = tau_coef + tauAD_coef
  
  ## muNAD, non-adherers (linear predictors of betaRD, nuNAD)
  l_betaRD = c(1, nuNAD_coef)
  muNAD_hat = sum(l_betaRD * betaRD_coef)
  
  ## delta, mu*phi + muNAD*(1-phi)
  delta_hat = mu_hat*phi_hat + muNAD_hat*(1-phi_hat) # change from baseline
  
  ## se(delta)
  deriv_delta_params = c(
    l_beta*phi_hat, # beta
    l_betaRD*(1-phi_hat), # betaRD
    0, # pi_coef
    (mu_hat - muNAD_hat), # tau_coef
    (mu_hat - muNAD_hat), # tauAD_coef
    sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuAD[2]*phi_hat + betaRD_coef[2]*(1-phi_hat), # nuNAD
    sum(beta_coef[c("base", paste0("factor(time)",K,":base"))])*l_nuAD[1]*phi_hat # nu
  )
  var_delta_hat = (t(deriv_delta_params) %*% joint_cov %*% deriv_delta_params) %>% as.vector
  se_delta_hat = sqrt(var_delta_hat)
  
  ## adjusted for baseline covariate
  deriv_nu_params = c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(betaRD_coef)), # betaRD
    rep(0, length(pi_coef)), # pi
    rep(0, length(tau_coef)), # tau
    rep(0, length(tauAD_coef)), # tauAD
    rep(0, length(nuNAD_coef)), # nuNAD
    rep(1, length(nu_coef)) # nu
  )
  var_nu = t(deriv_nu_params) %*% joint_cov %*% deriv_nu_params # p_nu * p_nu
  cov_nu_delta =  t(deriv_nu_params) %*% joint_cov %*% deriv_delta_params # p_nu * 1
  
  ## values for tipping point analysis
  delta_hat_tip = delta_hat + tip_delta*(pi_coef + tau_coef)
  deriv_delta_params_tip = deriv_delta_params + c(
    rep(0, length(beta_coef)), # beta
    rep(0, length(betaRD_coef)), # betaRD
    rep(tip_delta, length(pi_coef)), # pi
    rep(tip_delta, length(tau_coef)), # pi
    rep(0, length(tauAD_coef)), # nuNAD
    rep(0, length(nuNAD_coef)), # nuNAD
    rep(0, length(nu_coef)) # nu
  )
  var_delta_hat_tip = (t(deriv_delta_params_tip) %*% joint_cov %*% deriv_delta_params_tip) %>% as.vector
  se_delta_hat_tip = sqrt(var_delta_hat_tip)
  
  
  return(list(
    est = delta_hat, var = var_delta_hat, se = se_delta_hat,
    var_nu = var_nu, cov_nu_delta = cov_nu_delta, nu_coef = nu_coef,
    est_tip = delta_hat_tip, var_tip = var_delta_hat_tip, se_tip = se_delta_hat_tip
  ))
}

calc_est_RD_theta_v2 = function(delta_result_list, l_delta, mean_base, var_base){
  ## Function for estimation of treatment difference

  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  var_delta_hat_joint = map(delta_result_list, "var") %>% bdiag
  cov_nu_delta_joint = map(delta_result_list, "cov_nu_delta") %>% bdiag
  var_nu_inv_joint = map(delta_result_list, "var_nu") %>% bdiag %>% solve
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  return(list(
    theta = list(
      est = theta_hat, var = var_theta_hat, se = se_theta_hat,
      est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
      est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
    ),
    delta = list(
      est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
      est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
      est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
    )
  ))
}

calc_est_RD_theta_v3 = function(delta_result_list, l_delta, mean_base, var_base, tip = FALSE, complete = TRUE){
  ## Function for estimation of treatment difference in tipping point analysis
  
  J = sum(unlist(N_list))
  I = length(N_list) - 1
  A = matrix(1, nrow = J, ncol = I+1) # J*(I+1)
  
  delta_hat_joint = map(delta_result_list, "est") %>% unlist
  var_delta_hat_joint = map(delta_result_list, "var") %>% bdiag
  cov_nu_delta_joint = map(delta_result_list, "cov_nu_delta") %>% bdiag
  var_nu_inv_joint = map(delta_result_list, "var_nu") %>% bdiag %>% solve
  nu_coef_joint = map(delta_result_list, "nu_coef") %>% unlist
  
  theta_hat = l_delta %*% delta_hat_joint %>% as.vector
  var_theta_hat = (l_delta %*% var_delta_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_hat = sqrt(var_theta_hat)
  
  delta_adj_hat_joint = (delta_hat_joint + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
    as.vector
  var_delta_adj_hat_joint = var_delta_hat_joint - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_hat = l_delta %*% delta_adj_hat_joint %>% as.vector
  var_theta_adj_hat = (l_delta %*% var_delta_adj_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_hat = sqrt(var_theta_adj_hat)
  
  delta_adj_unc_hat_joint = delta_adj_hat_joint
  var_delta_adj_unc_hat_joint = var_delta_adj_hat_joint + 
    1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
  theta_adj_unc_hat = theta_adj_hat
  var_theta_adj_unc_hat = (l_delta %*% var_delta_adj_unc_hat_joint %*% t(l_delta)) %>%
    diag
  se_theta_adj_unc_hat = sqrt(var_theta_adj_unc_hat)
  
  if (tip){ # tipping point analysis is required
    delta_hat_joint_tip = map(delta_result_list, "est_tip") %>% unlist
    var_delta_hat_joint_tip = map(delta_result_list, "var_tip") %>% bdiag
    
    theta_hat_tip = l_delta %*% delta_hat_joint_tip %>% as.vector
    var_theta_hat_tip = (l_delta %*% var_delta_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_hat_tip = sqrt(var_theta_hat_tip)
    
    delta_adj_hat_joint_tip = (delta_hat_joint_tip + t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% (mean_base - nu_coef_joint)) %>%
      as.vector
    var_delta_adj_hat_joint_tip = var_delta_hat_joint_tip - t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_hat_tip = l_delta %*% delta_adj_hat_joint_tip %>% as.vector
    var_theta_adj_hat_tip = (l_delta %*% var_delta_adj_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_hat_tip = sqrt(var_theta_adj_hat_tip)
    
    delta_adj_unc_hat_joint_tip = delta_adj_hat_joint_tip
    var_delta_adj_unc_hat_joint_tip = var_delta_adj_hat_joint_tip + 
      1/J^2 * var_base * t(cov_nu_delta_joint) %*% var_nu_inv_joint %*% t(A) %*% A %*% var_nu_inv_joint %*% cov_nu_delta_joint
    theta_adj_unc_hat_tip = theta_adj_hat_tip
    var_theta_adj_unc_hat_tip = (l_delta %*% var_delta_adj_unc_hat_joint_tip %*% t(l_delta)) %>%
      diag
    se_theta_adj_unc_hat_tip = sqrt(var_theta_adj_unc_hat_tip)
    
    if (complete){ # complete result is required
      return(list(
        theta = list(
          est = theta_hat, var = var_theta_hat, se = se_theta_hat,
          est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
          est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
        ),
        delta = list(
          est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
          est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
          est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
        ),
        theta_tip = list(
          est = theta_hat_tip, var = var_theta_hat_tip, se = se_theta_hat_tip,
          est_adj = theta_adj_hat_tip, var_adj = var_theta_adj_hat_tip, se_adj = se_theta_adj_hat_tip,
          est_adj_unc = theta_adj_unc_hat_tip, var_adj_unc = var_theta_adj_unc_hat_tip, se_adj_unc = se_theta_adj_unc_hat_tip
        )
      ))
    } else { # pvalue only
      return(list(
        z_unadj = theta_hat_tip/se_theta_hat_tip,
        z_adj = theta_adj_hat_tip/se_theta_adj_hat_tip,
        z_adj_unc = theta_adj_unc_hat_tip/se_theta_adj_unc_hat_tip
      ))
    }
  } else { # tipping point analysis is not required
    return(list(
      theta = list(
        est = theta_hat, var = var_theta_hat, se = se_theta_hat,
        est_adj = theta_adj_hat, var_adj = var_theta_adj_hat, se_adj = se_theta_adj_hat,
        est_adj_unc = theta_adj_unc_hat, var_adj_unc = var_theta_adj_unc_hat, se_adj_unc = se_theta_adj_unc_hat
      ),
      delta = list(
        est = delta_hat_joint, var = diag(var_delta_hat_joint), se = sqrt(diag(var_delta_hat_joint)),
        est_adj = delta_adj_hat_joint, var_adj = diag(var_delta_adj_hat_joint), se_adj = sqrt(diag(var_delta_adj_hat_joint)),
        est_adj_unc = delta_adj_unc_hat_joint, var_adj_unc = diag(var_delta_adj_unc_hat_joint), se_adj_unc = sqrt(diag(var_delta_adj_unc_hat_joint))
      )
    ))
  }
}

calc_est_RD_TIP_v2 = function(
    beta_result, betaRD_result,
    pi_result, tau_result, tauAD_result,
    nuNAD_result, nu_result, 
    l_delta, mean_base, var_base, 
    tip = FALSE, tip_delta = NULL, complete = TRUE){
  ## general function for tipping point analysis

  if (tip & !complete){
    apply(tip_delta, 1, function(z){
      calc_est_RD_theta_v3(
        delta_result_list = mapply( 
          FUN = calc_est_RD_delta_v3,
          beta_result = beta_result, betaRD_result = betaRD_result, 
          pi_result = pi_result, tau_result = tau_result, tauAD_result = tauAD_result,
          nuNAD_result = nuNAD_result, nu_result = nu_result, 
          tip_delta = as.list(z),
          SIMPLIFY = FALSE
        ),
        l_delta, mean_base, var_base, tip, complete
      ) %>% unlist
      # only consider complete = FALSE in the current version
      # a vector, each C elements are the p values of C treatment effects, where C is the number of rows in l_delta
    }) %>% t # matrix, each row corresponds to each set of tipping points
  }
}
