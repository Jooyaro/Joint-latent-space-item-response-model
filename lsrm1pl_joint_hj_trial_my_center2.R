library("Rcpp")
sourceCpp("joint_lsirm.cpp")

# male_young_center2
data1_center2_my <- as.matrix(read.table("center2_my_exp.txt"))
data2_center2_my <- as.matrix(read.table("center2_my_psy.txt"))
nrow(data1_center2_my)
nrow(data2_center2_my)

nsample_center2_my <- nrow(data1_center2_my)
nitem1 <- ncol(data1_center2_my)
nitem2 <- ncol(data2_center2_my)

ndim <- 2
niter    <- 30000
nburn    <- 5000
nthin    <- 5
nprint   <- 200

jump_beta1 <- 0.3
jump_beta2 <- 0.3
jump_theta <-2.5
jump_z <- 1
jump_w1 <- 0.25
jump_w2 <- 0.2

pr_mean_beta  <- 0
pr_sd_beta    <- 1
pr_mean_theta <- 0
prior_a1 <- 0.001
prior_a2 <- 0.001
prior_b1 <- 0.001
prior_b2 <- 0.001

output_center2_my <- lsrm_joint_cpp(data1=data1_center2_my, data2=data2_center2_my, nsample=nsample_center2_my, nitem1=nitem1, nitem2=nitem2, ndim=ndim, niter=niter,
                             nburn=nburn, nthin=nthin, nprint=nprint, jump_beta1=jump_beta1, jump_beta2=jump_beta2,
                             jump_theta=jump_theta, jump_z=jump_z, jump_w1=jump_w1, jump_w2=jump_w2,
                             pr_mean_beta=pr_mean_beta, pr_sd_beta=pr_sd_beta, prior_a1=prior_a1, prior_a2=prior_a2,
                             prior_b1=prior_b1, prior_b2=prior_b2, pr_mean_theta=pr_mean_theta)


nmcmc = as.integer((niter - nburn) / nthin)
max.address = min(which.max(output_center2_my$map))
w1.star = output_center2_my$w1[max.address,,]
w2.star = output_center2_my$w2[max.address,,]
z.star = output_center2_my$z[max.address,,]
w1.proc = array(0,dim=c(nmcmc,nitem1,ndim))
w2.proc = array(0,dim=c(nmcmc,nitem2,ndim))
z.proc = array(0,dim=c(nmcmc,nsample_center2_my,ndim))
library(MCMCpack)
for(iter in 1:nmcmc){
  z.iter = output_center2_my$z[iter,,]
  if(iter != max.address) z.proc[iter,,] = procrustes(z.iter,z.star)$X.new
  else z.proc[iter,,] = z.iter
  
  w1.iter = output_center2_my$w1[iter,,]
  if(iter != max.address) w1.proc[iter,,] = procrustes(w1.iter,w1.star)$X.new
  else w1.proc[iter,,] = w1.iter
  
  w2.iter = output_center2_my$w2[iter,,]
  if(iter != max.address) w2.proc[iter,,] = procrustes(w2.iter,w2.star)$X.new
  else w2.proc[iter,,] = w2.iter
}

w1.est = matrix(NA,nitem1,ndim)
for(i in 1:nitem1){
  for(j in 1:ndim){
    w1.est[i,j] = mean(w1.proc[,i,j])
  }
}

w2.est = matrix(NA,nitem2,ndim)
for(i in 1:nitem2){
  for(j in 1:ndim){
    w2.est[i,j] = mean(w2.proc[,i,j])
  }
}

z.est = matrix(NA,nsample_center2_my,ndim)
for(k in 1:nsample_center2_my){
  for(j in 1:ndim){
    z.est[k,j] = mean(z.proc[,k,j])
  }
}

beta1.estimate = apply(output_center2_my$beta1, 2, mean)
beta2.estimate = apply(output_center2_my$beta2, 2, mean)
theta1.estimate = apply(output_center2_my$theta1, 2, mean)
theta2.estimate = apply(output_center2_my$theta2, 2, mean)
sigma_theta1.estimate = mean(output_center2_my$sigma_theta1)
sigma_theta2.estimate = mean(output_center2_my$sigma_theta2)

output_center2_my <- list(beta1_estimate  = beta1.estimate,
            beta2_estimate  = beta2.estimate,
            theta1_estimate = theta1.estimate,
            theta2_estimate = theta2.estimate,
            sigma_theta1_estimate    = sigma_theta1.estimate,
            sigma_theta2_estimate    = sigma_theta2.estimate,
            z_estimate     = z.est,
            w1_estimate     = w1.est,
            w2_estimate     = w2.est,
            beta1           = output_center2_my$beta1,
            beta2           = output_center2_my$beta2,
            theta1          = output_center2_my$theta1,
            theta2          = output_center2_my$theta2,
            theta_sd1       = output_center2_my$sigma_theta1,
            theta_sd2       = output_center2_my$sigma_theta2,
            z              = z.proc,
            w1              = w1.proc,
            w2              = w2.proc,
            accept_beta1    = output_center2_my$accept_beta1,
            accept_beta2    = output_center2_my$accept_beta2,
            accept_theta1   = output_center2_my$accept_theta1,
            accept_theta2   = output_center2_my$accept_theta2,
            accept_w1       = output_center2_my$accept_w1,
            accept_w2       = output_center2_my$accept_w2,
            accept_z       = output_center2_my$accept_z)

t(output_center2_my$accept_beta1)
t(output_center2_my$accept_beta2)
output_center2_my$accept_theta1
output_center2_my$accept_theta2
output_center2_my$accept_z
t(output_center2_my$accept_w1)
t(output_center2_my$accept_w2)

for(i in 1:44){
  for(j in 1:2){
    ts.plot(output_center2_my$w1[,i,j],main=i)
  }
}

for(i in 1:21){
  for(j in 1:2){
    ts.plot(output_center2_my$w2[,i,j],main=i)
  }
}

for(i in 1:44) ts.plot(output_center2_my$beta1[,i],main=i)
for(i in 1:21) ts.plot(output_center2_my$beta2[,i],main=i)

for(i in 1:100) ts.plot(output_center2_my$theta1[,i],main=i)
for(i in 1:100) ts.plot(output_center2_my$theta2[,i],main=i)


# save.image('lsrm_joint_trial_my_center2.RData')

plot(output_center2_my$w1_estimate[,1], output_center2_my$w1_estimate[,2], pch="", xlim = c(-3.5,3.5), ylim=c(-3.5,3.5))
points(output_center2_my$z_estimate, pch=".")
text(output_center2_my$w1_estimate[,1], output_center2_my$w1_estimate[,2], labels=1:44, cex=0.75, font=2, col='blue')
text(output_center2_my$w2_estimate[,1], output_center2_my$w2_estimate[,2], labels=1:21, cex=0.75, font=2, col='red')
