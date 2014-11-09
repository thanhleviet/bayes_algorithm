#learn Gibbs sampler:
#sampling from bivarite normal distribution.
#To do: rewrite the code to implement multinormal cases.
library(mvtnorm)
library(mcmc)
gen_norm_gibbs <- function(rou, start_point, size = 1000, dim = 2){
  if(length(start_point) != dim) 
    return('the dimension of start_point is wrong')
  corr_matrix <- matrix(rou, nrow = dim, ncol = dim)
  diag(corr_matrix) <- 1
  samples <- rmvnorm(100, c(10, 20), corr_matrix)
  mean_sp <- colMeans(samples)
  theta <- matrix(NA, nrow = size, ncol = dim)
  theta[1, ] <- start_point
  for (i in 2:nrow){
    for(j in 1:dim)
    theta[i, j] <- rnorm(1, mean_sp[1] + rou * (theta[i - 1, 2] - mean_sp[2]) , sqrt(1 - rou ^ 2))
    theta[i, j] <- rnorm(1, mean_sp[2] + rou * (theta[i, 1] - mean_sp[1]), sqrt(1 - rou ^ 2))
  }
  return(theta)
}

#line plots for 3 sequences of simulated data starting from different point, with rou = 0.8
theta <- gen_norm_gibbs(rou = 0.8, start_point = c(0, 0), size = 200)
theta1 <- gen_norm_gibbs(rou = 0.8, start_point = c(10, 10), size = 200)
theta2 <- gen_norm_gibbs(rou = 0.8, start_point = c(-10, -10), size = 200)
plot(theta[1:nrow, 1], theta[1:nrow, 2], asp = 1, type = 'o')
lines(theta1[1:nrow, 1], theta1[1:nrow, 2], asp = 1, type = 'o', col = 'red')
lines(theta2[1:nrow, 1], theta2[1:nrow, 2], asp = 1, type = 'o', col = 'blue')


#line plots for 3 sequences of simulated data starting from different point, with rou = 0
theta <- gen_norm_gibbs(rou = 0, start_point = c(0, 0), size = 200)
theta1 <- gen_norm_gibbs(rou = 0, start_point = c(10, 10), size = 200)
theta2 <- gen_norm_gibbs(rou = 0, start_point = c(-10, -10), size = 200)
plot(theta[1:nrow, 1], theta[1:nrow, 2], asp = 1, type = 'o')
lines(theta1[1:nrow, 1], theta1[1:nrow, 2], asp = 1, type = 'o', col = 'red')
lines(theta2[1:nrow, 1], theta2[1:nrow, 2], asp = 1, type = 'o', col = 'blue')

#Maybe for two dimensional parameter, the effects of correlation on convergence rate isn't very significant, so I try to do experiment on higher dimensions, lets say 10

#I want to test two things. 1. parameters related to each other converges very
#slowly compared to independent parameters. 2. After proper transformation, 
#the convergence rate will become faster.



#Metropolis hasting:How to choose jumnping distribution
#The variance of jumping distribution. This controls the jumping distance
metropolis_norm <- function(size, scale, start_point = c(0, 0)){
  post_mean <- c(0, 0)
  post_corr <- diag(2)
  theta <- matrix(NA, nrow = size, ncol = 2)
  theta[1, ] <- start_point
  jump_corr <- post_corr / scale ^ 2
  for (i in 2:size){
    theta_star <- rmvnorm(1, theta[i-1, ], jump_corr)
    r <-  dmvnorm(theta_star, post_mean, post_corr) / dmvnorm(theta[i - 1, ], post_mean, post_corr)
    flag <- rbinom(1, size = 1, prob = min(r, 1))
    if(flag == 1)
      theta[i, ] <- theta_star
    else
      theta[i, ] <- theta[i - 1,]
    }
  return(theta)
}

#set scale parameter to 5, the sequence seems to be a random walk
theta1 <- metropolis_norm(200, 5)
theta2 <- metropolis_norm(200, 5, c(-1, 1))
theta3 <- metropolis_norm(200, 5, c(-1, -1))
theta4 <- metropolis_norm(200, 5, c(1, 1))
theta5 <- metropolis_norm(200, 5, c(1, -1))
plot(theta1[, 1], theta1[, 2], asp = 1, type = 'o')
lines(theta2[, 1], theta2[, 2], asp = 1, type = 'o', col = 'red')
lines(theta3[, 1], theta3[, 2], asp = 1, type = 'o', col = 'blue')
lines(theta4[, 1], theta4[, 2], asp = 1, type = 'o', col = 'yellow')
lines(theta5[, 1], theta5[, 2], asp = 1, type = 'o', col = 'purple')

#set scale parameter to 1, this seems to be appropriate
theta1 <- metropolis_norm(200, 1)
theta2 <- metropolis_norm(200, 1, c(-1, 1))
theta3 <- metropolis_norm(200, 1, c(-1, -1))
theta4 <- metropolis_norm(200, 1, c(1, 1))
theta5 <- metropolis_norm(200, 1, c(1, -1))
plot(theta1[, 1], theta1[, 2], asp = 1, type = 'o')
lines(theta2[, 1], theta2[, 2], asp = 1, type = 'o', col = 'red')
lines(theta3[, 1], theta3[, 2], asp = 1, type = 'o', col = 'blue')
lines(theta4[, 1], theta4[, 2], asp = 1, type = 'o', col = 'yellow')
lines(theta5[, 1], theta5[, 2], asp = 1, type = 'o', col = 'purple')

#set scale parameter to 1/1, the problem is rejection rate is too high
theta1 <- metropolis_norm(200, 1/5)
theta2 <- metropolis_norm(200, 1/5, c(-1, 1))
theta3 <- metropolis_norm(200, 1/5, c(-1, -1))
theta4 <- metropolis_norm(200, 1/5, c(1, 1))
theta5 <- metropolis_norm(200, 1/5, c(1, -1))
plot(theta1[, 1], theta1[, 2], type = 'o')
lines(theta2[, 1], theta2[, 2], asp = 1, type = 'o', col = 'red')
lines(theta3[, 1], theta3[, 2], asp = 1, type = 'o', col = 'blue')
lines(theta4[, 1], theta4[, 2], asp = 1, type = 'o', col = 'yellow')
lines(theta5[, 1], theta5[, 2], asp = 1, type = 'o', col = 'purple')


###test for simulated tempering
#First I need to create a multimodal distribution. I decdied to use a mixed guassian
post_dist <- function(theta){
  sigma1 <- diag(2)
  sigma2 <- diag(2) 
  mean1 <- c(0, 0)
  mean2 <- c(20, 20)
  dens <- dmvnorm(theta, mean1, sigma1) + dmvnorm(theta, mean2, sigma2)
  dens
}

#The contour of posterior density
x <- y <- seq(-3, 23, 0.1)
theta <- expand.grid(x, y)
dens_val <- post_dist(theta)
mat_dens <- matrix(dens_val, nrow = 261)
contour(x, y, mat_dens, asp = 1)

#Then I use traditional Metropolis to sample data from this distribution
metro <- function(start = c(0, 0), size = 1000){
  theta <- matrix(NA, size, 2)
  sigma <- diag(2) * 4
  theta[1, ] <- start
  for (i in 2:size){
    theta_star <- rmvnorm(1, theta[i - 1, ], sigma)
    r <- min(1, post_dist(theta_star) / post_dist(theta[i - 1, ]))
    flag <- rbinom(1, 1, prob = r)
    if (flag == 1)
      theta[i, ] <- theta_star
    else
      theta[i, ] <- theta[i - 1, ]
  }
  theta
}
theta <- metro()
plot(theta[, 1], theta[, 2])

#In the above case, the random walk is stuck around one mode, it can't jump to 
#another mode. What a pity. But fortunately, we can use simulated tempering.

#How many k we need? It's difficult to choose, but a vague idea is that we need 
#the ramdom walk jump from one mode to another mdoe more easily.

#In this case, I choose k = 50 because the flat area between this two 
mat_dens2 <- mat_dens ^ 0.1
contour(x, y, mat_dens2, asp = 1)

mat_dens1 <- mat_dens ^ 0.02
contour(x, y, mat_dens1, asp = 1)

post_distk <- function(temp, theta) return post_dist(theta) ^ temp
#let k = 20
sim_temp <- function(size = 10000){
  temp <- seq(1, 0.05, -0.05)
  alpha0
  theta <- matrix(NA, size)
  theta[1, ] <- c(0, 0)
  k <- rep(NA, size)
  k[1] <- 1
  for (i in 2:size){
    k_star <- k[1] + 1
    u <- runif(1)
    k[i] <- ifelse(u > alpha0, k_star, k[i - 1])
    
    
  }
}


