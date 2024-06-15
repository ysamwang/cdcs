#include <RcppArmadillo.h>

/*
 * Helper for R function 'branchAndBound'
 * Takes in a set of ancestors and set of children to test
 * 
 * Inputs:
 * ancest: indices of ancestors
 * children: indices of possible children (of the ancestors)
 * dat: n x p matrix of data where v column is Y_v
 * K: a vector of degrees for test statistic
 * aggType: indicates which norm to take for test statistic
 * bs: number of bootstrap draws
 * intercept: should an additional intercept be included?
 * 
 * Returns:
 * pval: p-value for each child in children
 */

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List bnbHelperMulti(const arma::uvec ancest, const arma::uvec children, const arma::mat & dat,
                          arma::vec K, int aggType,  int bs, int intercept) {
  
  
  int n = dat.n_rows;
  int p = ancest.n_elem;
  
  // pulls out the relevant rows for X
  arma::mat X = dat.cols(ancest);  
  
  // adds an intercept and updates p if necessary
  if(intercept){
    X = arma::join_horiz(arma::ones(n), X);
    p = p + 1;
  }
  
  
  
  arma::colvec res;
  arma::vec nullDist = arma::zeros(bs);  
  arma::vec pval = arma::zeros(children.n_elem);
  arma::uvec randInd;
  arma::mat Z;
  double testStat;
  
  
  // computes the hat matrix
  arma::mat hatMat = (arma::eye(n, n) - X * solve(X.t() * X, X.t()));
  
  
  // If there are multiple K (the degree of the moment to check)
  int k_len = K.n_elem;

  // Q is the part that doesn't need to be recomputed each time
  // each slice corresponds to a specific value of K 
  arma::cube Q(p, n, k_len); 
  for(int k_ind = 0; k_ind < k_len; k_ind++){
    
    if(std::floor(K(k_ind)) == K(k_ind)){
      Q.slice(k_ind) = pow(X, K(k_ind) ).t() *  hatMat ;
      
    } else {
      Q.slice(k_ind) = pow(abs(X), K(k_ind) ).t() *  hatMat ;
      
    }
    
  
    
  }
  
  
  // for loop to check each child  
  for(int j = 0; j < children.n_elem; j++){
    
    
    // Reset statistic and nullDist
    testStat = 0;
    nullDist.zeros();
    
    // computes residuals
    res = hatMat * dat.col( children(j) );
    
    // computes test statistic
    for(int k_ind = 0; k_ind < k_len; k_ind++){
      
      if(std::floor(K(k_ind)) == K(k_ind)){
         Z = pow(X, K(k_ind)).t();
        
      } else {
        Z = pow(abs(X), K(k_ind)).t();
        
      }
      
      if ( aggType == 0){
        
        testStat += max( Z * res) / sqrt(n);
        
      } else if (aggType == 1){
        
        testStat += accu( abs(Z * res) ) / sqrt(n) ;
        
      } else {
        
        testStat += norm( Z * res / sqrt(n), "fro") ;
        
      }
      
    }
    
    
    for(int i = 0; i < bs; i++ ){
      
      // Draw indices for bootstrap sample 
      randInd =  arma::conv_to<arma::uvec>::from(randi( n, arma::distr_param(0,n-1) ));
      
      
      // compute the statistic for the bootstrap draw
      for(int k_ind = 0; k_ind < k_len; k_ind++){
        
        if (aggType == 0){
          
          nullDist(i) += max( abs(Q.slice(k_ind) * res.elem(randInd))) / sqrt((n-p))  ;
          
        } else if (aggType == 1){
          
          nullDist(i) += accu( abs(Q.slice(k_ind) * res.elem(randInd))) / sqrt((n-p))  ;
          
        } else {
          
          nullDist(i) += norm( Q.slice(k_ind) * res.elem(randInd)  / sqrt(n-p)  , "fro") ;
          
        }
      }
    }
    
    // compute p-value for jth child
    pval(j) = (double) arma::sum(nullDist >= testStat) / bs;
  }
  
  // return all p-values
  return Rcpp::List::create(Rcpp::Named("pVals")=pval);
}



