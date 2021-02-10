package <- c('numDeriv')
lapply(package, require, character.only = T)


###### setting for likelihood

# data
mu = 10
n = 1000
data <- rnorm(n, mu, sqrt(mu^2+1))
#data <-  c(-0.3835412,  0.3875394, -0.1183843, -0.4563267, -0.6981712,  0.8591746, -0.6643783, 
           #1.1393084, -0.2235991, -0.4846631,  1.0518586,  1.0743049,  1.0640797, -0.8561697, -0.1756521)
# constraint _inform
constraint_inform_1 <- function(x, beta){
  x - beta
}

constraint_inform_2 <- function(x, beta){
  x^2 - 2*beta^2 - 1
}

constraint_inform_1_d1 <- function(x, beta)
{
  # x[2] = beta
  constraint_inform_1_ins <- function(x){
    x[1]-x[2]
  }
  
  return(grad(constraint_inform_1_ins, c(x, beta))[2])
}
  
constraint_inform_2_d1 <- function(x, beta)
{
  #x[2] = beta
  constraint_inform_2_ins <- function(x){
    x[1]^2-2*x[2]^2 - 1
  }
  
  return(grad(constraint_inform_2_ins, c(x, beta))[2])
}  

# loglikelihood

# 0 means condition is ok 
check_condition <-  function(lambda, beta){
    
    g1 <- matrix(constraint_inform_1(data, beta), ncol = 1)
    g2 <- matrix(constraint_inform_2(data, beta), ncol = 1)
    g <- cbind(g1, g2)
    check_ls <- 1 - g%*%lambda
    need <- sum(check_ls <= 1/length(data) )
    
    return(need)
  }

profile_loglikelihood <- function(lambda, beta){
  
  g1 <- matrix(constraint_inform_1(data, beta), ncol = 1)
  g2 <- matrix(constraint_inform_2(data, beta), ncol = 1)
  g <- cbind(g1, g2)
  need <- sum(-log(1 - g%*%lambda))
  
  return(need)
}
  
profile_loglikelihood_d1 <- function(lambda, beta){
  
  g1 <- matrix(constraint_inform_1(data, beta), ncol = 1)
  g2 <- matrix(constraint_inform_2(data, beta), ncol = 1)
  g <- cbind(g1, g2)
  denom <- 1 - g%*%lambda 
  need <- t(cbind(sum(g1/denom), sum(g2/denom)))

  return(need)
}


profile_loglikelihood_d2 <- function(lambda, beta)
{
  g1 <- matrix(constraint_inform_1(data, beta), ncol = 1)
  g2 <- matrix(constraint_inform_2(data, beta), ncol = 1)
  g <- cbind(g1, g2)
  denom <- 1 - g%*%lambda
  
  
  gi_vec <- c(g1[1], g2[1])
  need <- gi_vec%*%t(gi_vec)/denom[1]^2
  
  for(i in c(2:length(data)))
  {
    gi_vec <- c(g1[i], g2[i])
    need <- need + gi_vec%*%t(gi_vec)/denom[i]^2
  }
  
  
  return(need)
}


T_fcn <- function(lambda, beta){
  
  g1 <- matrix(constraint_inform_1(data, beta) , ncol = 1)
  g2 <- matrix(constraint_inform_2(data, beta) , ncol = 1)
  g <- cbind(g1, g2)
  denom <- 1 - g%*%lambda
  
  dg1 <- c();dg2 <- c()
  
  for(i in c(1:length(data))){
    dg1[i] <- constraint_inform_1_d1(data[i], beta)
    dg2[i] <- constraint_inform_2_d1(data[i], beta)
  }

  need <- cbind(sum(dg1/denom), sum(dg2/denom))
  
  return(need)
}

  