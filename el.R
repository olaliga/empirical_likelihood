
package <- c('numDeriv', 'matlib', 'ggplot2')
lapply(package, require, character.only = T)

import <- c('profile_likelihood.R', 'plot.R')
lapply(import, source)


# objective function : pi from lagrange, profile in likelihood and maximize likelihood 
## 2 version : 1 for nonparametric 2 for semiparametric

# nonparametric
#inner loop: apply function, fixed beta, minimized profile likelihood
inner <- function(beta)
{
  #step 0
  iteration = 0; epsilon = 10^-4; lambda = c(0, 0); lagrange = 0
  lambda_norm <- 1
  
  
  
  while(lambda_norm >= epsilon)
  {
    iteration <- iteration + 1
    print(paste('inner iteration:', iteration))
    #step 1
    
    #print(profile_loglikelihood_d2(lambda = lambda, beta = beta))
    #print(profile_loglikelihood_d1(lambda = lambda, beta = beta))
    print(paste('profile_ll_d2:', solve(profile_loglikelihood_d2(lambda = lambda, beta = beta) ) ) )
    print(paste('profile_ll_d1:',profile_loglikelihood_d1(lambda = lambda, beta = beta) ) )
    
    delta = solve(profile_loglikelihood_d2(lambda = lambda, beta = beta)) %*% profile_loglikelihood_d1(lambda = lambda, beta = beta)
    
    
    print(paste('delta', delta))
    
    #step 2
    eta = 1
    lambda_temp <- lambda - eta*delta ; lagrange_temp = profile_loglikelihood(lambda = lambda_temp, beta = beta)
    while(check_condition(lambda = lambda_temp, beta = beta) != 0 | lagrange_temp >= lagrange){
      eta = eta/2
      lambda_temp <- lambda - eta*delta; lagrange_temp = profile_loglikelihood(lambda = lambda_temp, beta = beta)
    }
    
    #step 3
    lambda_norm <- sqrt(t(lambda - lambda_temp)%*%(lambda - lambda_temp))
    print(paste('lagrange_temp:', lagrange_temp))
    print(paste('lagrange:', lagrange))
    print(paste('eta:', eta))
    lambda <- lambda_temp; lagrange = lagrange_temp
    
    print(paste('lambda_norm:', lambda_norm))
    print('--------------------------------------------------------')
  }
  
  
  need <- list('lambda' = lambda, 'lagrange' = lagrange, 'iteration' = iteration)
  return(need)
}

#outer loop:
outer <- function()
{
  #step 0
  iteration_outer = 0; sigma = 5; epsilon = 10^-4
  beta = mean(data) # mle of data
  
  #
  beta_norm = 1
  
  inner_beta <- inner(beta)
  lagrange = profile_loglikelihood(lambda = inner_beta$lambda, beta = beta)
  
  #save
  need_beta <- c()
  need_lagrange <- c()
  need_inner <- c()
  need_lambda <- c()
  need_iteration <- c()
  need_lambda_1 <- c()
  need_lambda_2 <- c()
  
  
  while(beta_norm >= epsilon)
  {
    iteration_outer = iteration_outer + 1
    print('*******************************')
    print(paste('iteration_outer:', iteration_outer))
    
    #step 1
    m_1 <- T_fcn(lambda = inner_beta$lambda, beta = beta) %*% inner_beta$lambda
    print(paste('m_1', m_1))
    
    m_2 <- -T_fcn(lambda = inner_beta$lambda, beta = beta) %*% 
      solve(profile_loglikelihood_d2(lambda = inner_beta$lambda, beta = beta)) %*% t( T_fcn(lambda = inner_beta$lambda, beta = beta))
    print(paste('m_2', m_2))
    
    delta = solve(m_2)%*%m_1
    print(paste('delta', delta))
    
    eta = 1; t = 0
    
    #step 2
    beta_temp = beta - eta*delta
    
    inner_beta_temp <- inner(beta_temp)
    lagrange_temp <- profile_loglikelihood(lambda = inner_beta_temp$lambda, beta = beta_temp)
    
    while(lagrange_temp <= lagrange & t != sigma)
    {
      t = t + 1; eta = eta/2
      beta_temp = beta - eta*delta
      inner_beta_temp <- inner(beta_temp)
      lagrange_temp <- profile_loglikelihood(lambda = inner_beta_temp$lambda, beta = beta_temp)
      
    }
    
    #step 3
    print(paste('beta_temp:', beta_temp))
    print(paste('beta:', beta))
    print(paste('lagrange_temp:', lagrange_temp))
    print(paste('lagrange:', lagrange))
    
    
    beta_norm = sqrt(t(beta_temp - beta)%*%(beta_temp-beta))
    
    beta = beta_temp; lagrange = lagrange_temp; inner_beta = inner_beta_temp
    
    need_beta[iteration_outer] <- beta; need_lagrange[iteration_outer] <- lagrange
    need_lambda_1[iteration_outer] <- inner_beta$lambda[1]
    need_lambda_2[iteration_outer] <- inner_beta$lambda[2]
    need_iteration[iteration_outer] <- inner_beta$iteration
    
  }
  
  need <- data.frame(beta = unlist(need_beta), lagrange = unlist(need_lagrange), 
                     lambda_1 = unlist(need_lambda_1), lambda_2 = unlist(need_lambda_2),
                     iteration = unlist(need_iteration))
  
  
  
  return(need)
  
}

beta_ls <- c()
apply_fcn <- function(x)
{
  source('profile_likelihood.R')
  result <- outer()
  return(result$beta[length(result$beta)])
}
beta_ls <- sapply(1:1000, apply_fcn)


beta_ls

