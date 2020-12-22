// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y Nx1 Numeric vector.
//' @param X NxP Numeric matrix.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
// [[Rcpp::export]]
SEXP fitOLS(const arma::colvec y, const arma::mat X){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = X.n_cols;
  // Information
  const arma::mat A = X.t()*X;
  // Estimate beta
  const arma::vec b = arma::solve(A,X.t()*y,arma::solve_opts::likely_sympd);
  // Calculate residuals
  const arma::vec eps = (y-X*b);
  // Scale
  const double v = arma::as_scalar(eps.t()*eps/(n-p));
  // Information
  const arma::mat Ibb = A/v;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("V")=v,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("Resid")=eps);
}