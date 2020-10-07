#include <RcppArmadillo.h>

//' Confidence sets for Causal Discovery
//' 
//' Test whether the goodness of fit for a linear model
//'
//' @param dat an n by p matrix of covariates
//' @param Y a vector of length n containin the dependent variable
//' @param K a vector containing integer values greater than 2 which will be used for the test statistic
//' @param aggType 0 is L-inf, 1 is L1, 2 is L1, 3 is all
//' @param bs the number of bootstrap resamples for the null distribution
//' @param intercept 
//' @return Product of v1 and v2

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List singleTestcppMulti(const arma::mat & dat, const arma::colvec & Y, arma::ivec K, int aggType, int bs, int intercept) {
  
  
  int n = dat.n_rows;
  int p = dat.n_cols;
  arma::mat X = dat;
  int k_len = K.n_elem;
  
  
  arma::uvec randInd;
  
  
  arma::vec nullDistInf = arma::zeros(bs);
  arma::vec nullDistOne = arma::zeros(bs);
  arma::vec nullDistTwo = arma::zeros(bs);
  
  double testStatInf = 0;
  double testStatOne = 0;
  double testStatTwo = 0;
  
  
  // If including an intercept
  if(intercept){
    X = arma::join_horiz(arma::ones(n), dat);
    p = p + 1;
  }
  

  
  
  // Form hat matrix and compute residuals
  arma::mat hatMat = (arma::eye(n, n) - X * solve(X.t() * X, X.t()));
  arma::colvec res = hatMat * Y;
  
  
  
  arma::cube Q(p, n, k_len);
  for(int k_ind = 0; k_ind < k_len; k_ind++){
    
    Q.slice(k_ind) = pow(X, K(k_ind) ).t() *  hatMat ;
    
  }
  
  

  arma::colvec statCounter(p); 

  statCounter.zeros();  
  for(int k_ind = 0; k_ind < k_len; k_ind++){

    statCounter = abs(pow(X, K(k_ind)).t() * res) / sqrt(n);

  
    if(aggType == 0){
      
      // if aggType is max
      testStatInf += max(statCounter) ;
      
    } else if (aggType == 1){
      
      testStatOne += accu(statCounter);
      
    } else  if (aggType == 2) {
      
      testStatTwo += norm(statCounter, "fro") ;
      
    } else {
  
      testStatInf += max(statCounter) ;
      testStatOne += accu(statCounter) ;
      testStatTwo += norm(statCounter , "fro") ;
      
    }
  }
  
  
  for(int i = 0; i < bs; i++ ){
    
    randInd =  arma::conv_to<arma::uvec>::from(randi( n, arma::distr_param(0,n-1) ));
    
    statCounter.zeros();
    
    
    for(int k_ind = 0; k_ind < k_len; k_ind++){
      statCounter = abs(Q.slice(k_ind) * res.elem(randInd)) / sqrt(n - p);
      
      if(aggType == 0){
        
        // if aggType is max
        nullDistInf(i) += max(statCounter);
        
      } else if (aggType == 1){
        
        // if aggType is sum
        nullDistOne(i) += accu(statCounter);
        
      } else if(aggType == 2){
        
        // if aggType is sum
        nullDistTwo(i) += norm(statCounter, "fro");
        
      } else {
        
        // if aggType is max
        nullDistInf(i) += max(statCounter);
        nullDistOne(i) += accu(statCounter);
        nullDistTwo(i) += norm(statCounter, "fro");
        
      }
    }
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("nullDistInf")=nullDistInf,
                            Rcpp::Named("nullDistOne")=nullDistOne,
                            Rcpp::Named("nullDistTwo")=nullDistTwo,

                            Rcpp::Named("testStatInf") = testStatInf,
                            Rcpp::Named("testStatOne") = testStatOne,
                            Rcpp::Named("testStatTwo") = testStatTwo);
}





// aggType 0 is L-Inf
// aggType 1 is  L-1
// aggType 2 is just L-2
// aggType 3 is all

/*
 * Given the indices `ordered' tests whether the last element could be a child of the other indices
 * 
 */

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List exhaustiveHelperMulti(const arma::uvec ordered, const arma::mat & dat,
                                 arma::ivec K, int aggType, int bs, int intercept) {
  
  int n = dat.n_rows;
  int p = ordered.n_elem - 1;
  arma::mat X = dat.cols(ordered.subvec(0, p -1));
  arma::colvec Y = dat.col(ordered(p));
  int k_len = K.n_elem;
  arma::colvec statCounter(p); 
  arma::uvec randInd;
  
  arma::vec nullDistInf = arma::zeros(bs);
  arma::vec nullDistOne = arma::zeros(bs);
  arma::vec nullDistTwo = arma::zeros(bs);
  
  double testStatInf = 0;
  double testStatOne = 0;
  double testStatTwo = 0;
  
  
  // add intercept to X if needed
  if(intercept){
    X = arma::join_horiz(arma::ones(n), X);
    p = p + 1; 
  }
  
  
  // Hat matrix and residuals
  


  arma::mat hatMat = (arma::eye(n, n) - X * solve(X.t() * X, X.t()));
  arma::colvec res = hatMat * Y;
  
  
  // Q mat used later to to compute bootstrap stats 
  arma::cube Q(p, n, k_len);
  for(int k_ind = 0; k_ind < k_len; k_ind++){
    
    Q.slice(k_ind) = pow(X, K(k_ind) ).t() *  hatMat;
    
  }
  
  
  
  
  // Calculate observed testStat
  statCounter.zeros();  
  for(int k_ind = 0; k_ind < k_len; k_ind++){
    
    statCounter = abs(pow(X, K(k_ind)).t() * res) / sqrt(n);
    
    
    if(aggType == 0){
      
      // if aggType is max
      testStatInf += max(statCounter);
      
    } else if (aggType == 1){
      
      testStatOne += accu(statCounter);
      
    } else  if (aggType == 2) {
      
      testStatTwo += norm(statCounter, "fro") ;
      
    } else {
      
      testStatInf += max(statCounter);
      testStatOne += accu(statCounter);
      testStatTwo += norm(statCounter, "fro") ;
      
    }
  }
  
  
  // Each bootstrap sample
  for(int i = 0; i < bs; i++ ){
    
    // Sample residuals
    randInd =  arma::conv_to<arma::uvec>::from(randi( n, arma::distr_param(0,n-1) ));
    statCounter.zeros();
    
    
    for(int k_ind = 0; k_ind < k_len; k_ind++){
      statCounter = abs(Q.slice(k_ind) * res.elem(randInd)) / sqrt(n - p);
      
      if(aggType == 0){
        
        // if aggType is max
        nullDistInf(i) += max(statCounter);
        
      } else if (aggType == 1){
        
        // if aggType is sum
        nullDistOne(i) += accu(statCounter);
        
      } else if(aggType == 2){
        
        // if aggType is sum
        nullDistTwo(i) += norm(statCounter, "fro");
        
      } else {
        
        // if aggType is max
        nullDistInf(i) += max(statCounter);
        nullDistOne(i) += accu(statCounter);
        nullDistTwo(i) += norm(statCounter, "fro");
        
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("nullDistInf")=nullDistInf,
                            Rcpp::Named("nullDistOne")=nullDistOne,
                            Rcpp::Named("nullDistTwo")=nullDistTwo,
                            
                            Rcpp::Named("testStatInf") = testStatInf,
                            Rcpp::Named("testStatOne") = testStatOne,
                            Rcpp::Named("testStatTwo") = testStatTwo);
}


