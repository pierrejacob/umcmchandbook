#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Fast crossproduct of single matrix
// Calculates t(X)*X using RcppEigen
// [[Rcpp::export]]
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X){
  const int n = X.cols();
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}

// Fast crossproduct of two matrices
// [[Rcpp::export]]
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
  return Eigen::MatrixXd(X*Y);
}