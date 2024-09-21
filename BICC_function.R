# BICC function ----------------------------------------------------------

# Arguments
# data = data frame of data
# niter = number of iterations
# nchain = number of Gibbs chains
# nburn = number of burn-in

BICC <- function(data, niter, nchain, nburn){
  
  #################### clmm  
  
  # convert data to run clmm 
  data$video = as.factor(data$video)
  data$rater = as.factor(data$rater)
  data$score_ord = as.factor(data$score_ord)
  
  
  # run clmm
  clmm_fm1 <- clmm(score_ord ~ predictor + (1 | video) + (1 | rater), data = data, Hess = TRUE, link = "probit")
  
  #print(summary(clmm_fm1))
  
  # convert data to run gibbs
  video = as.numeric(data$video)
  rater = as.numeric(data$rater)
  predictor = as.numeric(data$predictor)
  score = as.numeric(data$score_ord)
  
  # number of categories of the responses
  K <- length(unique(score))
  # number of subjects
  Nvideo <- length(unique(video))
  # number of raters
  Nrater <- length(unique(rater))
  # total observations
  Nobs <- length(score)
  
  # number of scores of each video
  numberoftimes <- c()
  for(v in 1:Nvideo){numberoftimes[v] <- length(subset(score, video==v))}
  scores_per_subject <- ifelse(var(numberoftimes)==0, numberoftimes[1], 'Unequal')
  
  # workload of each rater 
  raterwatched <- c()
  for(r in 1:Nrater){raterwatched[r] <- length(subset(score, rater==r))}
  rater_workload <- ifelse(var(raterwatched)==0, raterwatched[1], 'Unequal')
  
  # Data type
  if(rater_workload==Nvideo){
    type = 'Complete'
  } else if(var(raterwatched)==0 & rater_workload<Nvideo){
    type = 'Balanced'
  } else {
    type = 'Imbalanced'
  }
  
  
  gibbs_fn <- function(niter, video, rater, predictor, score
                       , beta.inital = unname(clmm_fm1$beta)
                       , delta.inital = unname(clmm_fm1$alpha)
                       , sigma_s_sq.init = unname((clmm_fm1$optRes$par[5])^2)
                       , sigma_r_sq.init = unname((clmm_fm1$optRes$par[6])^2)
  ) 
  {
    
    # sampling progress
    samp_prog <- seq(0, niter, 1000)
    
    # initialize parameters
    
    # beta 
    # hyperparameter of prior distribution of beta
    mu0 <- 0
    Sigma0 <- 1
    
    beta <- beta.inital
    
    # set sigma_sq = 1
    sigma_sq <- 1
    
    # sigma_s_sq
    sigma_s_sq <- sigma_s_sq.init
    
    # Sigma_r_sq
    sigma_r_sq <- sigma_r_sq.init
    
    # video effect
    video_gibbs <- rnorm(Nvideo, mean = 0, sd = sqrt(sigma_s_sq))
    
    # rater effect
    rater_gibbs <- rnorm(Nrater, mean = 0, sd = sqrt(sigma_r_sq))
    
    # cut-point - delta
    delta <- rep(0, K+1)
    delta[1] <- -Inf
    delta[K+1] <- Inf
    rand <- rnorm(K-1, mean = unname(quantile(sapply(1:Nobs, function(i) beta*predictor[i] + video_gibbs[video[i]] + rater_gibbs[rater[i]]), probs = c(0.25,0.50,0.75))), sd = 0.5)
    delta[2:K] <- rand[order(rand)] 
    
    # create empty matrices to store draws
    beta_sample <- matrix(0, nrow = niter, ncol = 1)
    delta_sample <- matrix(0, nrow = niter, ncol = max(score)-1)
    video_sample <- matrix(0, nrow = niter, ncol = Nvideo)
    rater_sample <- matrix(0, nrow = niter, ncol = Nrater)
    sigma_s_sq_sample <- matrix(0, nrow = niter, ncol = 1)
    sigma_r_sq_sample <- matrix(0, nrow = niter, ncol = 1)
    rho_sample <- matrix(0, nrow = niter, ncol = 1)
    
    # Metropolis-Hastings algorithm for selection of hyperparameters
    # create empty matrices to store draws
    # prior distributions of sigma_s_sq: IG(a_s, b_s) 
    # hyperparameter: a_s, b_s
    a_s <- matrix(0, nrow = niter, ncol = 1)
    b_s <- matrix(0, nrow = niter, ncol = 1)
    # prior distributions of sigma_r_sq: IG(a_r, b_r) 
    # hyperparameter: a_r, b_r
    a_r <- matrix(0, nrow = niter, ncol = 1)
    b_r <- matrix(0, nrow = niter, ncol = 1)
    
    # Initialize the hyperparameter values
    a_s[1,] <- 5
    b_s[1,] <- 5
    a_r[1,] <- 5
    b_r[1,] <- 5
    
    iter <- 2
    
    while (iter <= niter) {
      
      # store draws
      beta_sample[iter, ] <- beta
      delta_sample[iter, ] <- delta[2:K]
      video_sample[iter, ] <- video_gibbs
      rater_sample[iter, ] <- rater_gibbs
      sigma_s_sq_sample[iter, ] <- sigma_s_sq
      sigma_r_sq_sample[iter, ] <- sigma_r_sq
      
      # latent variable draw
      lower_b <- delta[score]
      upper_b <- delta[score+1]
      
      # sample latent y from the truncated normal distribution
      y_g <- sapply(1:Nobs,
                    function(i) truncnorm::rtruncnorm(n = 1,
                                                      mean = beta*predictor[i] + video_gibbs[video[i]] + rater_gibbs[rater[i]],
                                                      sd = sqrt(sigma_sq),
                                                      a = lower_b[i],
                                                      b = upper_b[i]))  
      
      # beta draw
      Sigma1 <- solve(solve(Sigma0) + (1/sigma_sq)*t(predictor)%*%predictor)
      y_g_prime <- sapply(1:Nobs,
                          function(i)  y_g[i] - video_gibbs[video[i]] - rater_gibbs[rater[i]])
      mu1 <- Sigma1%*%(solve(Sigma0)%*%(mu0) + (1/sigma_sq)*t(predictor)%*%y_g_prime)
      beta <- mvrnorm(n = 1, mu = mu1, Sigma = Sigma1)
      
      # cutoff points draw
      for (i in 2:K) {
        delta[i] <- runif(1, min = max(max(y_g[score==i-1]), delta[i-1], delta[1]), 
                          max = min(min(y_g[score==i]), delta[K+1])) 
      }
      
      # video effect draw
      sig_s_sq <- c()
      sig_s_sq <- sapply(1:Nvideo, function(v) ((numberoftimes[v]/sigma_sq) + (1/sigma_s_sq))^(-1))
      yy <- sapply(1:Nobs, function(i)  (y_g[i]  - beta*predictor[i] - rater_gibbs[rater[i]]))
      new_mat_s <- cbind(yy, video)
      sum_s <- rep(0, Nvideo)
      for(i in 1:Nvideo){sum_s[i] <- sum(subset(new_mat_s,video==i)[,1])}
      mean_s <- sapply(1:Nvideo, function(v)  (sig_s_sq[v]/sigma_sq)*sum_s[v])
      video_gibbs <- sapply(1:Nvideo, function(i) rnorm(1, mean = mean_s[i], sd = sqrt(sig_s_sq[i])) )
      
      # rater effect draw
      sig_r_sq <- c()
      sig_r_sq <- sapply(1:Nrater, function(r) ((raterwatched[r]/sigma_sq) + (1/sigma_r_sq))^(-1))
      yyy <- sapply(1:Nobs, function(i)  (y_g[i]  - beta*predictor[i] - video_gibbs[video[i]] ))
      new_mat_r <- cbind(yyy, rater)
      sum_r <- rep(0, Nrater)
      for(i in 1:Nrater){sum_r[i] <- sum(subset(new_mat_r,rater==i)[,1])}
      mean_r <- sapply(1:Nrater, function(r) (sig_r_sq[r]/sigma_sq)*sum_r[r])
      rater_gibbs <- sapply(1:Nrater, function(i) rnorm(1, mean = mean_r[i], sd = sqrt(sig_r_sq[i])))
      
      # MH algorithm for hyperparameter - a_s, b_s
      # pdf of IG distribution
      a_s_func <- function(x){(b_s[iter-1,]^x)*(sigma_s_sq^(-x-1))/gamma(x)}
      b_s_func <- function(x){(x^(a_s[iter-1,]))*exp((-x/sigma_s_sq))}
      
      # generate a candidate for a_s 
      z1=rexp(n = 1, rate = 0.5)
      # generate a candidate for b_s
      z2=rexp(n = 1, rate = 0.5)
      # calculate probabilities of acceptance
      # for component 1 a_s
      rat1 = exp(log(a_s_func(z1))+log(dexp(x = a_s[iter-1,], rate = 0.5, log = FALSE))-log(a_s_func(a_s[iter-1,]))-log(dexp(x = z1, rate = 0.5, log = FALSE)))
      rho1 = min(rat1,1)
      # for component 2 beta_v
      rat2 = exp(log(b_s_func(z2))+log(dexp(x = b_s[iter-1,], rate = 0.5, log = FALSE))-log(b_s_func(b_s[iter-1,]))-log(dexp(x = z2, rate = 0.5, log = FALSE)))
      rho2 = min(rat2,1)
      ##decide whether accept the candidates or not
      coin1 = rbinom(1,1,rho1)
      coin2 = rbinom(1,1,rho2)
      if(coin1==1){a_s[iter,]=z1}else{a_s[iter,]=a_s[iter-1,]}
      if(coin2==1){b_s[iter,]=z2}else{b_s[iter,]=b_s[iter-1,]}
      
      
      # sigma_s_sq draw
      new_a_s <- a_s[iter,] + (Nvideo/2)
      new_b_s <- b_s[iter,] + (0.5*(sum(video_gibbs^2)))
      sigma_s_sq <- rinvgamma(1, shape = new_a_s, rate = new_b_s)
      
      # MH algorithm for hyperparameter - a_r, b_r
      # pdf of IG distribution
      a_r_func <- function(x){(b_r[iter-1,]^x)*(sigma_r_sq^(-x-1))/gamma(x)}
      b_r_func <- function(x){(x^(a_r[iter-1,]))*exp((-x/sigma_r_sq))}
      
      # generate a candidate for alpha_r 
      z3=rexp(n = 1, rate = 0.5)
      # generate a candidate for beta_r
      z4=rexp(n = 1, rate = 0.5)
      # calculate probabilities of acceptance
      # for component 1 alpha_r
      rat3 = exp(log(a_r_func(z3))+log(dexp(x = a_r[iter-1,], rate = 0.5, log = FALSE))-log(a_r_func(a_r[iter-1,]))-log(dexp(x = z3, rate = 0.5, log = FALSE)))
      rho3 = min(rat3,1)
      # for component 2 beta_r
      rat4 = exp(log(b_r_func(z4))+log(dexp(x = b_r[iter-1,], rate = 0.5, log = FALSE))-log(b_r_func(b_r[iter-1,]))-log(dexp(x = z4, rate = 0.5, log = FALSE)))
      rho4 = min(rat4,1)
      # decide whether accept the candidates or not
      coin3 = rbinom(1,1,rho3)
      coin4 = rbinom(1,1,rho4)
      if(coin3==1){a_r[iter,]=z3}else{a_r[iter,]=a_r[iter-1,]}
      if(coin4==1){b_r[iter,]=z4}else{b_r[iter,]=b_r[iter-1,]}
      
      # sigma_r_sq draw
      new_a_r <- a_r[iter,] + (Nrater/2)
      new_b_r <- b_r[iter,] + (0.5*(sum(rater_gibbs^2)))
      sigma_r_sq <- rinvgamma(1, shape = new_a_r, rate = new_b_r)
      
      # calculate rho
      rho_g <- sigma_s_sq/(sigma_sq+sigma_s_sq+sigma_r_sq)
      
      # store draws
      rho_sample[iter, ] <- rho_g
      
      iter <- iter + 1
      
      # print the progress
      if (iter %in% samp_prog) {
        cat('Sampling iteration', iter, 'complete','\n')
      }
      
    }
    
    out <- list(beta = beta_sample[,1],
                sigma_s_sq = sigma_s_sq_sample[,1],
                sigma_r_sq = sigma_r_sq_sample[,1],
                delta_1 = delta_sample[,1],
                delta_2 = delta_sample[,2],
                delta_3 = delta_sample[,3],
                rho = rho_sample[,1]
    )
    
    return(out)
    
  }
  
  
  #################### Gibbs sampling
  
  cat('-----------------------------------------------------------------','\n',sep = '')
  
  cat('Gibbs sampling progress', '\n','\n')
  
  # run gibbs function
  for(k in 1:nchain){
    cat('Chain', k, '\n')
    assign(paste0('BICC_chain',k),gibbs_fn(niter = niter, video = video, rater = rater, predictor = predictor, score = score))
  }
  
  warmup = nburn
  
  start = warmup+1
  end = niter
  
  
  ######## Gibbs parameters names #########
  parameternames <- c("beta", 
                      "sigma_s_sq", 
                      "sigma_r_sq", 
                      "delta_1",
                      "delta_2",
                      "delta_3",
                      "rho")
  
  performance_metrics <- c()
  
  for(l in 1:length(parameternames)){
    
    
    # create empty vector to store all samples of each parameter
    assign(paste0('BICC_',parameternames[l],'_samples'), c())
    
    for(k in 1:nchain){
      # store all samples of each parameter
      assign(paste0('BICC_',parameternames[l],'_samples'), cbind(get(paste0('BICC_',parameternames[l],'_samples')),
                                                                 get(paste0('BICC_chain',k))[[l]][start:end]))
    }
    
    # draw samples histogram
    #hist(get(paste0('BICC_',parameternames[l],'_samples')), main = paste0('BICC_',parameternames[l],'_samples'))
    
    # store mean, SD, 95 lower bound, upper bound 
    assign(paste0('BICC_',parameternames[l],'_samp'),c(mean(get(paste0('BICC_',parameternames[l],'_samples'))[,1]),
                                                       sd(get(paste0('BICC_',parameternames[l],'_samples'))[,1]),
                                                       unname(quantile(get(paste0('BICC_',parameternames[l],'_samples'))[,1], na.rm=TRUE, probs = c(0.025, 0.975)))))
    
    performance_metrics <- rbind(performance_metrics,get(paste0('BICC_',parameternames[l],'_samp')))
    colnames(performance_metrics) <- c("Mean", "SD", "95% CI Lower bound", "Upper bound")
    
  }
  
  rownames(performance_metrics) <- c(parameternames)
  
  
  #######  Check convergence  ############################################################
  
  cat('-----------------------------------------------------------------','\n',sep = '')
  
  for(l in 1:length(parameternames)){
    
    # create empty vector to store all samples of each parameter
    assign(paste0('BICC_',parameternames[l],'_monitor'), c())
    
    for(k in 1:nchain){
      # store all samples of each parameter
      assign(paste0('BICC_',parameternames[l],'_monitor'), cbind(get(paste0('BICC_',parameternames[l],'_monitor')),
                                                                 get(paste0('BICC_chain',k))[[l]]))
    }
    cat('Check convergence of', parameternames[l], '\n','\n')
    monitor(get(paste0('BICC_',parameternames[l],'_monitor')), warmup = nburn)
    cat('-----------------------------------------------------------------------------','\n',sep = '')
    
  }
  
  
  # Report data information
  
  cat('Number of subjects:', Nvideo, '\n')
  cat('Number of raters:', Nrater, '\n')
  cat('Number of scores per subject:', scores_per_subject, '\n')
  cat('Workload of each rater:', rater_workload, '\n')
  cat('Total number of observations:', Nobs, '\n')
  cat('Data type:', type)
  
  cat('\n','-----------------------------------------------------------------','\n',sep = '')
  
  # Report parameter estimates and 95% Credible interval
  cat('Gibbs Sampling Results', '\n',sep = '')
  
  print(performance_metrics)
  
  cat('-----------------------------------------------------------------','\n',sep = '')
  
  
  
  #################     plot the samples    #################################################
  
  
  for(l in 1:length(parameternames)){
    
    # Get min, max to plot samples
    ylim_min <- min(get(paste0('BICC_',parameternames[l],'_samples')))
    ylim_max <- max(get(paste0('BICC_',parameternames[l],'_samples')))
    
    
    # plot samples
    plot(start:end, get(paste0('BICC_',parameternames[l],'_samples'))[,1],
         ylim=c(ylim_min,ylim_max),type="l", ylab = "", 
         main = paste('BICC', parameternames[l]), xlab = "Iteration")
    if(nchain>1){for(k in 1:nchain){points(start:end, get(paste0('BICC_',parameternames[l],'_samples'))[,k], type="l", col=k)}}
  }
}
