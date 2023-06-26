// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::export]]
Rcpp::List lsrm_joint_cpp(arma::mat data1, arma::mat data2, const int nsample, const int nitem1, const int nitem2, 
                          const int ndim, const int niter, const int nburn, const int nthin, const int nprint,
                          const double jump_beta1, const double jump_beta2, const double jump_theta, const double jump_z, 
                          const double jump_w1, const double jump_w2, const double pr_mean_beta, const double pr_sd_beta,
                          const double prior_a1, const double prior_a2, const double prior_b1, const double prior_b2,
                          const double pr_mean_theta){
  
  // Rprintf("Error00\n");
  
  int i, j, k, count, accept;
  double num, den, old_like_beta1, old_like_beta2, new_like_beta1, new_like_beta2; 
  double old_like_theta1, old_like_theta2, new_like_theta1, new_like_theta2, pr_sd_theta1 = 1.0, pr_sd_theta2 = 1.0;
  double old_like_z, new_like_z, old_like_w1, old_like_w2, new_like_w1, new_like_w2;
  double ratio, un, post_a1, post_a2, post_b1, post_b2, dist_temp1, dist_temp2, dist_old_temp1, dist_old_temp2, dist_new_temp1, dist_new_temp2;
  double pr_mean_z = 0.0, pr_sd_z = 1.0, pr_mean_w = 0.0, pr_sd_w = 1.0, mle;
  
  arma::dvec oldbeta1(nitem1, fill::randu);
  oldbeta1 = oldbeta1 * 4.0 - 2.0;
  arma::dvec newbeta1 = oldbeta1;
  
  arma::dvec oldbeta2(nitem2, fill::randu);
  oldbeta2 = oldbeta2 * 4.0 - 2.0;
  arma::dvec newbeta2 = oldbeta2;
  
  arma::dvec oldtheta1(nsample, fill::randu);
  oldtheta1 = oldtheta1 * 4.0 - 2.0;
  arma::dvec newtheta1 = oldtheta1;
  
  arma::dvec oldtheta2(nsample, fill::randu);
  oldtheta2 = oldtheta2 * 4.0 - 2.0;
  arma::dvec newtheta2 = oldtheta2;
  
  arma::dmat oldz(nsample,ndim,fill::randu);
  oldz = oldz * 2.0 - 1.0;
  arma::dmat newz = oldz;
  
  arma::dmat oldw1(nitem1,ndim,fill::randu);
  oldw1 = oldw1 * 2.0 - 1.0;
  arma::dmat neww1 = oldw1;
  
  arma::dmat oldw2(nitem2,ndim,fill::randu);
  oldw2 = oldw2 * 2.0 - 1.0;
  arma::dmat neww2 = oldw2;
  
  arma::dmat samp_beta1((niter-nburn)/nthin, nitem1, fill::zeros);
  arma::dmat samp_beta2((niter-nburn)/nthin, nitem2, fill::zeros);
  arma::dmat samp_theta1((niter-nburn)/nthin, nsample, fill::zeros);
  arma::dmat samp_theta2((niter-nburn)/nthin, nsample, fill::zeros);
  arma::dcube samp_z(((niter-nburn)/nthin), nsample, ndim, fill::zeros);
  arma::dcube samp_w1(((niter-nburn)/nthin), nitem1, ndim, fill::zeros);
  arma::dcube samp_w2(((niter-nburn)/nthin), nitem2, ndim, fill::zeros);
  arma::dvec samp_sd_theta1((niter-nburn)/nthin, fill::zeros);
  arma::dvec samp_sd_theta2((niter-nburn)/nthin, fill::zeros);
  arma::dvec sample_mle((niter-nburn)/nthin, fill::zeros);
  
  arma::dvec accept_beta1(nitem1, fill::zeros);
  arma::dvec accept_beta2(nitem2, fill::zeros);
  arma::dvec accept_theta1(nsample, fill::zeros);
  arma::dvec accept_theta2(nsample, fill::zeros);
  arma::dvec accept_z(nsample, fill::zeros);
  arma::dvec accept_w1(nitem1, fill::zeros);
  arma::dvec accept_w2(nitem2, fill::zeros);
  
  accept = count = 0;
  
  arma::dmat dist1(nsample,nitem1,fill::zeros);
  arma::dmat dist2(nsample,nitem2,fill::zeros);
  arma::dvec old_dist_k1(nitem1,fill::zeros);
  arma::dvec new_dist_k1(nitem1,fill::zeros);
  arma::dvec old_dist_k2(nitem2,fill::zeros);
  arma::dvec new_dist_k2(nitem2,fill::zeros);
  arma::dvec old_dist_i1(nsample,fill::zeros);
  arma::dvec new_dist_i1(nsample,fill::zeros);
  arma::dvec old_dist_i2(nsample,fill::zeros);
  arma::dvec new_dist_i2(nsample,fill::zeros);
  
  // Rprintf("Error01\n");
  
  for(int iter = 0; iter < niter; iter++){
    
    //dist1(j,i) is distance of z_j and w1_i
    dist1.fill(0.0);
    for(i = 0; i < nitem1; i++){
      for(k = 0; k < nsample; k++){
        dist_temp1 = 0.0;
        for(j = 0; j < ndim; j++) dist_temp1 += std::pow((oldz(k,j)-oldw1(i,j)), 2.0);
        dist1(k,i) = std::sqrt(dist_temp1);
      }
    }
    
    //dist2(j,i) is distance of z_j and w2_i
    dist2.fill(0.0);
    for(i = 0; i < nitem2; i++){
      for(k = 0; k < nsample; k++){
        dist_temp2 = 0.0;
        for(j = 0; j < ndim; j++) dist_temp2 += std::pow((oldz(k,j)-oldw2(i,j)), 2.0);
        dist2(k,i) = std::sqrt(dist_temp2);
      }
    }
    
    // Rprintf("Error02\n");
    
    // beta1 update
    for(i = 0; i < nitem1; i++){
      newbeta1(i) = R::rnorm(oldbeta1(i), jump_beta1);
      old_like_beta1 = new_like_beta1 = 0.0;
      for(k = 0; k < nsample; k++){
        if(data1(k,i) == 1.0) new_like_beta1 += -std::log(1.0 + std::exp(-(newbeta1(i) + oldtheta1(k) - dist1(k,i))));
        else new_like_beta1 += -std::log(1.0 + std::exp(newbeta1(i) + oldtheta1(k) - dist1(k,i)));
        if(data1(k,i) == 1.0) old_like_beta1 += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - dist1(k,i))));
        else old_like_beta1 += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - dist1(k,i)));
      }
      
      num = new_like_beta1 + R::dnorm4(newbeta1(i), pr_mean_beta, pr_sd_beta, 1);
      den = old_like_beta1 + R::dnorm4(oldbeta1(i), pr_mean_beta, pr_sd_beta, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldbeta1(i) = newbeta1(i);
        accept_beta1(i) += 1.0 / (niter * 1.0);
      }
      else newbeta1(i) = oldbeta1(i);
    }
    
    // Rprintf("Error03\n");
    
    // beta2 update
    for(i = 0; i < nitem2; i++){
      newbeta2(i) = R::rnorm(oldbeta2(i), jump_beta2);
      old_like_beta2 = new_like_beta2 = 0.0;
      for(k = 0; k < nsample; k++){
        if(data2(k,i) == 1.0) new_like_beta2 += -std::log(1.0 + std::exp(-(newbeta2(i) + oldtheta2(k) - dist2(k,i))));
        else new_like_beta2 += -std::log(1.0 + std::exp(newbeta2(i) + oldtheta2(k) - dist2(k,i)));
        if(data2(k,i) == 1.0) old_like_beta2 += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - dist2(k,i))));
        else old_like_beta2 += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - dist2(k,i)));
      }
      
      num = new_like_beta2 + R::dnorm4(newbeta2(i), pr_mean_beta, pr_sd_beta, 1);
      den = old_like_beta2 + R::dnorm4(oldbeta2(i), pr_mean_beta, pr_sd_beta, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldbeta2(i) = newbeta2(i);
        accept_beta2(i) += 1.0 / (niter * 1.0);
      }
      else newbeta2(i) = oldbeta2(i);
    }
    
    // Rprintf("Error04\n");
    
    // theta1 update
    for(k = 0; k < nsample; k++){
      newtheta1(k) = R::rnorm(oldtheta1(k), jump_theta);
      old_like_theta1 = new_like_theta1 = 0.0;
      
      for(i = 0; i < nitem1; i++){
        if(data1(k,i) == 1.0) new_like_theta1 += -std::log(1.0 + std::exp(-(oldbeta1(i) + newtheta1(k) - dist1(k,i))));
        else new_like_theta1 += -std::log(1.0 + std::exp(oldbeta1(i) + newtheta1(k) - dist1(k,i)));
        if(data1(k,i) == 1.0) old_like_theta1 += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - dist1(k,i))));
        else old_like_theta1 += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - dist1(k,i)));
      }
      num = new_like_theta1 + R::dnorm4(newtheta1(k), pr_mean_theta, pr_sd_theta1, 1);
      den = old_like_theta1 + R::dnorm4(oldtheta1(k), pr_mean_theta, pr_sd_theta1, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldtheta1(k) = newtheta1(k);
        accept_theta1(k) += 1.0 / (niter * 1.0);
      }
      else newtheta1(k) = oldtheta1(k);
    }
    
    // Rprintf("Error05\n");
    
    // theta2 update
    for(k = 0; k < nsample; k++){
      newtheta2(k) = R::rnorm(oldtheta2(k), jump_theta);
      old_like_theta2 = new_like_theta2 = 0.0;
      
      for(i = 0; i < nitem2; i++){
        if(data2(k,i) == 1.0) new_like_theta2 += -std::log(1.0 + std::exp(-(oldbeta2(i) + newtheta2(k) - dist2(k,i))));
        else new_like_theta2 += -std::log(1.0 + std::exp(oldbeta2(i) + newtheta2(k) - dist2(k,i)));
        if(data2(k,i) == 1.0) old_like_theta2 += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - dist2(k,i))));
        else old_like_theta2 += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - dist2(k,i)));
      }
      num = new_like_theta2 + R::dnorm4(newtheta2(k), pr_mean_theta, pr_sd_theta2, 1);
      den = old_like_theta2 + R::dnorm4(oldtheta2(k), pr_mean_theta, pr_sd_theta2, 1);
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        oldtheta2(k) = newtheta2(k);
        accept_theta2(k) += 1.0 / (niter * 1.0);
      }
      else newtheta2(k) = oldtheta2(k);
    }
    
    // Rprintf("Error06\n");
    
    
    // zj update
    for(k = 0; k < nsample; k++){
      for(j = 0; j < ndim; j++) newz(k,j) = R::rnorm(oldz(k,j), jump_z);
      old_like_z = new_like_z = 0.0;
      
      //calculate distance of oldw1 and newz
      for(i = 0; i < nitem1; i++){
        dist_old_temp1 = dist_new_temp1 = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp1 += std::pow((newz(k,j)-oldw1(i,j)), 2.0);
          dist_old_temp1 += std::pow((oldz(k,j)-oldw1(i,j)), 2.0);
        }
        new_dist_k1(i) = sqrt(dist_new_temp1);
        old_dist_k1(i) = sqrt(dist_old_temp1);
      }
      
      //calculate distance of oldw2 and newz
      for(i = 0; i < nitem2; i++){
        dist_old_temp2 = dist_new_temp2 = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp2 += std::pow((newz(k,j)-oldw2(i,j)), 2.0);
          dist_old_temp2 += std::pow((oldz(k,j)-oldw2(i,j)), 2.0);
        }
        new_dist_k2(i) = sqrt(dist_new_temp2);
        old_dist_k2(i) = sqrt(dist_old_temp2);
      }
      
      //calculate likelihood
      for(i = 0; i < nitem1; i++){
        if(data1(k,i) == 1.0) new_like_z += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - new_dist_k1(i))));
        else new_like_z += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - new_dist_k1(i)));
        if(data1(k,i) == 1.0) old_like_z += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - old_dist_k1(i))));
        else old_like_z += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - old_dist_k1(i)));
      }
      
      for(i = 0; i < nitem2; i++){
        if(data2(k,i) == 1.0) new_like_z += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - new_dist_k2(i))));
        else new_like_z += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - new_dist_k2(i)));
        if(data2(k,i) == 1.0) old_like_z += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - old_dist_k2(i))));
        else old_like_z += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - old_dist_k2(i)));
      }
      
      num = den = 0.0;
      for(j = 0; j < ndim; j++){
        num += R::dnorm4(newz(k,j),pr_mean_z,pr_sd_z,1);
        den += R::dnorm4(oldz(k,j),pr_mean_z,pr_sd_z,1);
      }
      //Rprintf("%.3f %.3f %.3f %.3f\n", num, den, new_like_z, old_like_z);
      //arma::dvec newzz = dmvnorm(newz.cols(2*j,2*j+1),pr_mean_z,pr_cov_z,TRUE);
      //arma::dvec oldzz = dmvnorm(oldz.cols(2*j,2*j+1),pr_mean_z,pr_cov_z,TRUE);
      
      num += new_like_z;
      den += old_like_z;
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        for(j = 0; j < ndim; j++) oldz(k,j) = newz(k,j);
        accept_z(k) += 1.0 / (niter * 1.0);
      }
      else{
        for(j = 0; j < ndim; j++) newz(k,j) = oldz(k,j);
      }
    }
    
    // w1i update
    for(i = 0; i < nitem1; i++){
      for(j = 0; j < ndim; j++) neww1(i,j) = R::rnorm(oldw1(i,j), jump_w1);
      old_like_w1 = new_like_w1 = 0.0;
      
      //calculate distance of neww1 and oldz
      for(k = 0; k < nsample; k++){
        dist_old_temp1 = dist_new_temp1 = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp1 += std::pow((oldz(k,j)-neww1(i,j)), 2.0);
          dist_old_temp1 += std::pow((oldz(k,j)-oldw1(i,j)), 2.0);
        }
        new_dist_i1(k) = sqrt(dist_new_temp1);
        old_dist_i1(k) = sqrt(dist_old_temp1);
      }
      
      //calculate likelihood
      for(k = 0; k < nsample; k++){
        if(data1(k,i) == 1.0) new_like_w1 += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - new_dist_i1(k))));
        else new_like_w1 += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - new_dist_i1(k)));
        if(data1(k,i) == 1.0) old_like_w1 += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - old_dist_i1(k))));
        else old_like_w1 += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k) - old_dist_i1(k)));
      }
      
      num = den = 0.0;
      for(j = 0; j < ndim; j++){
        num += R::dnorm4(neww1(i,j),pr_mean_w,pr_sd_w,1);
        den += R::dnorm4(oldw1(i,j),pr_mean_w,pr_sd_w,1);
      }
      
      num += new_like_w1;
      den += old_like_w1;
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        for(j = 0; j < ndim; j++) oldw1(i,j) = neww1(i,j);
        accept_w1(i) += 1.0 / (niter * 1.0);
      }
      else{
        for(j = 0; j < ndim; j++) neww1(i,j) = oldw1(i,j);
      } 
    }
    
    // w2i update
    for(i = 0; i < nitem2; i++){
      for(j = 0; j < ndim; j++) neww2(i,j) = R::rnorm(oldw2(i,j), jump_w2);
      old_like_w2 = new_like_w2 = 0.0;
      
      //calculate distance of neww2 and oldz
      for(k = 0; k < nsample; k++){
        dist_old_temp2 = dist_new_temp2 = 0.0;
        for(j = 0; j < ndim; j++){
          dist_new_temp2 += std::pow((oldz(k,j)-neww2(i,j)), 2.0);
          dist_old_temp2 += std::pow((oldz(k,j)-oldw2(i,j)), 2.0);
        }
        new_dist_i2(k) = sqrt(dist_new_temp2);
        old_dist_i2(k) = sqrt(dist_old_temp2);
      }
      
      //calculate likelihood
      for(k = 0; k < nsample; k++){
        if(data2(k,i) == 1.0) new_like_w2 += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - new_dist_i2(k))));
        else new_like_w2 += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - new_dist_i2(k)));
        if(data2(k,i) == 1.0) old_like_w2 += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - old_dist_i2(k))));
        else old_like_w2 += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k) - old_dist_i2(k)));
      }
      
      num = den = 0.0;
      for(j = 0; j < ndim; j++){
        num += R::dnorm4(neww2(i,j),pr_mean_w,pr_sd_w,1);
        den += R::dnorm4(oldw2(i,j),pr_mean_w,pr_sd_w,1);
      }
      
      num += new_like_w2;
      den += old_like_w2;
      ratio = num - den;
      
      if(ratio > 0.0) accept = 1;
      else{
        un = R::runif(0,1);
        if(std::log(un) < ratio) accept = 1;
        else accept = 0;
      }
      
      if(accept == 1){
        for(j = 0; j < ndim; j++) oldw2(i,j) = neww2(i,j);
        accept_w2(i) += 1.0 / (niter * 1.0);
      }
      else{
        for(j = 0; j < ndim; j++) neww2(i,j) = oldw2(i,j);
      } 
    }
    
    //sigma_theta1 update with gibbs
    post_a1 = 2 * prior_a1  + nsample;
    post_b1 = prior_b1;
    for(j = 0; j < nsample; j++) post_b1 += std::pow((oldtheta1(j) - pr_mean_theta), 2.0);
    pr_sd_theta1 = std::sqrt(2 * post_b1 *(1.0 /  R::rchisq(post_a1)));
    
    //sigma_theta2 update with gibbs
    post_a2 = 2 * prior_a2 + nsample;
    post_b2 = prior_b2;
    for(j = 0; j < nsample; j++) post_b2 += std::pow((oldtheta2(j) - pr_mean_theta), 2.0);
    pr_sd_theta2 = std::sqrt(2 * post_b2 *(1.0 / R::rchisq(post_a2)));
    
    if(iter >= nburn && iter % nthin == 0){
      for(i = 0; i < nitem1; i++) samp_beta1(count,i) = oldbeta1(i);
      for(k = 0; k < nsample; k++) samp_theta1(count,k) = oldtheta1(k);
      for(i = 0; i < nitem1; i++){
        for(j = 0; j < ndim; j++){
          samp_w1(count,i,j) = oldw1(i,j);
        }
      }
      for(k = 0; k < nsample; k++){
        for(j = 0; j < ndim; j++){
          samp_z(count,k,j) = oldz(k,j);
        }
      }
      
      samp_sd_theta1(count) = pr_sd_theta1;
      
      for(i = 0; i < nitem2; i++) samp_beta2(count,i) = oldbeta2(i);
      for(k = 0; k < nsample; k++) samp_theta2(count,k) = oldtheta2(k);
      for(i = 0; i < nitem2; i++){
        for(j = 0; j < ndim; j++){
          samp_w2(count,i,j) = oldw2(i,j);
        }
      }
      
      samp_sd_theta2(count) = pr_sd_theta2;
      
      //dist1(j,i) is distance of z_j and w1_i
      dist1.fill(0.0);
      for(i = 0; i < nitem1; i++){
        for(k = 0; k < nsample; k++){
          dist_temp1 = 0.0;
          for(j = 0; j < ndim; j++) dist_temp1 += std::pow((oldz(k,j)-oldw1(i,j)), 2.0);
          dist1(k,i) = std::sqrt(dist_temp1);
        }
      }
      
      //dist2(j,i) is distance of z_j and w2_i
      dist2.fill(0.0);
      for(i = 0; i < nitem2; i++){
        for(k = 0; k < nsample; k++){
          dist_temp2 = 0.0;
          for(j = 0; j < ndim; j++) dist_temp2 += std::pow((oldz(k,j)-oldw2(i,j)), 2.0);
          dist2(k,i) = std::sqrt(dist_temp2);
        }
      }
      
      mle = 0.0;
      for(i = 0; i < nitem1; i++) mle += R::dnorm4(oldbeta1(i), pr_mean_beta, pr_sd_beta, 1);
      for(i = 0; i < nitem2; i++) mle += R::dnorm4(oldbeta2(i), pr_mean_beta, pr_sd_beta, 1);
      for(k = 0; k < nsample; k++) mle += R::dnorm4(oldtheta1(k), pr_mean_theta, pr_sd_theta1, 1);
      for(k = 0; k < nsample; k++) mle += R::dnorm4(oldtheta2(k), pr_mean_theta, pr_sd_theta2, 1);
      for(i = 0; i < nitem1; i++)
        for(j = 0; j < ndim; j++) mle += R::dnorm4(oldw1(i,j),pr_mean_w,pr_sd_w,1);
      for(i = 0; i < nitem2; i++)
        for(j = 0; j < ndim; j++) mle += R::dnorm4(oldw2(i,j),pr_mean_w,pr_sd_w,1);
      for(k = 0; k < nsample; k++)
        for(j = 0; j < ndim; j++) mle += R::dnorm4(oldz(k,j),pr_mean_z,pr_sd_z,1);
      for(k = 0; k < nsample; k++){
        for(i = 0; i < nitem1; i++){
          if(data1(k,i) == 1.0) mle += -std::log(1.0 + std::exp(-(oldbeta1(i) + oldtheta1(k) - dist1(k,i))));
          else mle += -std::log(1.0 + std::exp(oldbeta1(i) + oldtheta1(k)- dist1(k,i)));
        }
      }
      for(k = 0; k < nsample; k++){
        for(i = 0; i < nitem2; i++){
          if(data2(k,i) == 1.0) mle += -std::log(1.0 + std::exp(-(oldbeta2(i) + oldtheta2(k) - dist2(k,i))));
          else mle += -std::log(1.0 + std::exp(oldbeta2(i) + oldtheta2(k)- dist2(k,i)));
        }
      }
      sample_mle(count) = mle;
      count++;
    } // burn, thin
    
    if(iter % nprint == 0){
      Rprintf("Iteration: %.5u ", iter); 
      for(i = 0 ; i < nitem1 ; i++ ) {
        Rprintf("% .3f ", oldbeta1(i));
      }
      for(i = 0 ; i < nitem2 ; i++ ) {
        Rprintf("% .3f ", oldbeta2(i));
      }
      Rprintf(" %.3f ", pr_sd_theta1);
      Rprintf(" %.3f\n", pr_sd_theta2);
    }
    
  } //for end
  
  Rcpp::List output;
  output["beta1"] = samp_beta1;
  output["beta2"] = samp_beta2;
  output["theta1"] = samp_theta1;
  output["theta2"] = samp_theta2;
  output["z"] = samp_z;
  output["w1"] = samp_w1;
  output["w2"] = samp_w2;
  output["sigma_theta1"] = samp_sd_theta1;
  output["sigma_theta2"] = samp_sd_theta2;
  output["map"] = sample_mle;
  output["accept_beta1"] = accept_beta1;
  output["accept_beta2"] = accept_beta2;
  output["accept_theta1"] = accept_theta1;
  output["accept_theta2"] = accept_theta2;
  output["accept_z"] = accept_z;
  output["accept_w1"] = accept_w1;
  output["accept_w2"] = accept_w2;
  
  return(output);
  
} // function end


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

-