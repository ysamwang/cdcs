#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

//' bnbHelperGof
//' 
//' This function is a helper for the branchAndBound function
//'
//' @param ancest indices of ancestor in the dat matrix
//' @param children indices of children in the dat matrix
//' @param dat n x p matrix with observations
//' @param G an n x k x |ancest| array of where each variable (corresponding to slice) has k test functions 
//' @param  withinAgg indicates which norm to take when combining test statistics within
//'    a variable but across test functions
//' @param aggType indicates which norm to take when combining test statistics across variables
//' @param bs number of bootstrap draws
//' @param intercept should an additional intercept be included?
//' @return 
//' pval: p-value for each child in children
//' @export
//[[Rcpp::export]]
Rcpp::List bnbHelperGof(const arma::uvec ancest, const arma::uvec children, const arma::mat & dat,
                          const arma::cube & G, int withinAgg, int aggType,  int bs, int intercept) {
  
  
  int n = dat.n_rows;
  // Number of test functions for each variable
  int k = G.n_cols;
  
  // pulls out the relevant rows for X
  arma::mat X = dat.cols(ancest);  
  // adds an intercept and updates p if necessary
  if(intercept){
    X = arma::join_horiz(arma::ones(n), X);
  }
  int p = X.n_cols;  
  
  
  arma::colvec res(n);
  arma::vec nullDist(bs);  
  arma::vec pval = arma::zeros(children.n_elem);
  arma::uvec randInd(n);
  
  double testStat;
  
  
  // computes the hat matrix
  arma::mat hatMat = (arma::eye(n, n) - X * solve(X.t() * X, X.t()));
  
  
  // Q mat used later to to compute bootstrap stats 
  arma::cube Q(k, n, ancest.n_elem);
  
  for(int j = 0; j < ancest.n_elem; j++){
    Q.slice(j) = G.slice( ancest(j) ).t() * hatMat ;
  }
  
  arma::mat statMat(k, ancest.n_elem);
  arma::vec testStatisticVec(ancest.n_elem);
  arma::vec nullStatVec(ancest.n_elem);
  
  
  // for loop to check each child  
  for(int j = 0; j < children.n_elem; j++){
    
    // Reset statistic and nullDist
    testStat = 0;
    nullDist.zeros();
    
    // computes residuals
    res = hatMat * dat.col( children(j) );
    // K x p matrix: each column corresponds 
    
    statMat = abs(Q.each_slice() * res / sqrt(n)); 

    
    if(withinAgg == 1){
      
      for(int z = 0; z < ancest.n_elem; z++){
        testStatisticVec(z) = norm(statMat.col(z),1);
      }
      
      
    } else if (withinAgg == 2){
      
      for(int z = 0; z < ancest.n_elem; z++){
        testStatisticVec(z) = norm(statMat.col(z), "fro");
      }
      
    } else if (withinAgg == 3){
      
      for(int z = 0; z < ancest.n_elem; z++){
        testStatisticVec(z) = norm(statMat.col(z),"inf");
      }
      
    }
    
    if(aggType == 1){
      testStat = norm(testStatisticVec, 1);
    } else if (aggType == 2){
      testStat = norm(testStatisticVec, "fro");
    } else if (aggType == 3){
      testStat = norm(testStatisticVec, "inf");
    }
    
    
    
    // Compute null distribution for 
    for(int i = 0; i < bs; i++ ){
      
      // Draw indices for bootstrap sample 
      randInd =  arma::conv_to<arma::uvec>::from(randi( n, arma::distr_param(0,n-1) ));
      
      
      // compute the statistic for the bootstrap draw
      nullStatVec.zeros();
      
      statMat = abs(Q.each_slice() * res.elem(randInd) / sqrt(n - p)); 
      
      
      if(withinAgg == 1){
        
        for(int z = 0; z < ancest.n_elem; z++){
          nullStatVec(z) = norm(statMat.col(z),1);
        }
        
        
      } else if (withinAgg == 2){
        
        for(int z = 0; z < ancest.n_elem; z++){
          nullStatVec(z) = norm(statMat.col(z),2);
        }
        
      } else if (withinAgg == 3){
        
        for(int z = 0; z < ancest.n_elem; z++){
          nullStatVec(z) = norm(statMat.col(z),"inf");
        }
      }
      
      if(aggType == 1){
        nullDist(i) = norm(nullStatVec, 1);
      } else if (aggType == 2){
        nullDist(i) = norm(nullStatVec, "fro");
      } else {
        nullDist(i) = norm(nullStatVec, "inf");
      }
    }

    
    // compute p-value for jth child
    pval(j) = (double) (arma::sum(nullDist >= testStat) + 1) / (bs + 1);
  }
  
  // return all p-values
  return Rcpp::List::create(Rcpp::Named("pVals")=pval);
}



