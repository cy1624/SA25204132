// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat flexible_reader(Rcpp::RObject data) {
  arma::mat result;
  
  if (data.isS4()) {
    Rcpp::warning("S4 object detected, attempting conversion");
    Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(data);
    result = Rcpp::as<arma::mat>(mat);
  }
  else if (Rf_isMatrix(data)) {
    Rcpp::NumericMatrix mat = Rcpp::as<Rcpp::NumericMatrix>(data);
    result = Rcpp::as<arma::mat>(mat);
  }
  else if (Rf_isNewList(data) || TYPEOF(data) == VECSXP) {
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(data);
    
    int n_cols = df.size();
    int n_rows = df.nrow();
    result.set_size(n_rows, n_cols);
    
    for (int i = 0; i < n_cols; i++) {
      Rcpp::NumericVector col = df[i];
      result.col(i) = Rcpp::as<arma::vec>(col);
    }
  }
  else {
    Rcpp::stop("Input must be a matrix or dataframe");
  }
  return result;
}

// [[Rcpp::export]]
arma::mat Sigma_k_cpp(const arma::mat& data, int k) {
  int n = data.n_rows;
  int p = data.n_cols;
  
  if (k >= n) {
    Rcpp::stop("k must be less than n");
  }
  
  arma::mat data1 = data.submat(0, 0, n-k-1, p-1);
  arma::mat data2 = data.submat(k, 0, n-1, p-1);
  
  return arma::cov(data1, data2);
}

//' @title Calculating the numbers of factors with Rcpp
//' @description Rcpp version of \code{FactorNum}. See more in \code{FactorNum}.
//' @param data a data frame or matrix 
//' @param k0 a pre-specified integer used in calculating covariance matrix. The exact value won't affect the result much and it could be proper to choose it as 5(the default value) in practice.
//' @param j0 a pre-specified integer representing the maximum of the number of factors. It must be less than the number of variables. The default value is one forth of the number of variables(rounded down to the nearest integer).
//' @return a list consists of \code{r0}, \code{r}, \code{R} and \code{tao}. \code{r0} is the number of common factors. \code{r} is the number of specific factors. \code{R} is all the ratios of eigenvalue sums. \code{tao} is the indices of local maximums of R, which is sorted into descending order by the corresponding value of R.
//' @examples
//' \dontrun{
//' data(CSI300_2)
//' nums <- FactorNum_cpp(CSI300_2)   
//' cat(c(nums$r0, nums$r))
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List FactorNum_cpp(Rcpp::RObject data, int k0 = 5, int j0 = -1) {
  arma::mat data0 = flexible_reader(data);
  int p = data0.n_cols;
  
  if (j0 == -1) {
    j0 = std::ceil(p / 4.0);
  }
  
  arma::mat EigenValues(j0 + 1, k0 + 2, arma::fill::zeros);
  
  for (int k = 0; k <= k0 + 1; k++) {
    arma::mat Sigma = Sigma_k_cpp(data0, k);
    arma::mat M = Sigma * Sigma.t();
    
    arma::vec eigval;
    arma::mat eigvec;
    
    if (!arma::eig_sym(eigval, eigvec, M)) {
      Rcpp::stop("Eigenvalue decomposition failed");
    }
    
    eigval = arma::reverse(eigval);
    
    int num_vals = std::min(j0 + 1, (int)eigval.n_elem);
    EigenValues.submat(0, k, num_vals - 1, k) = eigval.subvec(0, num_vals - 1);
  }
  
  arma::vec EigenSums(j0 + 1, arma::fill::zeros);
  for (int j = 0; j <= j0; j++) {
    double sum = 0.0;
    for (int k = 0; k <= k0 + 1; k++) {
      sum += EigenValues(j, k);
    }
    EigenSums(j) = sum;
  }
  
  arma::vec R(j0 + 1);
  R(0) = 1.0;
  for (int i = 1; i <= j0; i++) {
    if (EigenSums(i) > 0) {
      R(i) = EigenSums(i - 1) / EigenSums(i);
    } else {
      R(i) = 0.0;
    }
  }
  
  std::vector<int> LocalMax_index;
  for (int i = 1; i < j0; i++) {
    if (R(i) > R(i - 1) && R(i) > R(i + 1)) {
      LocalMax_index.push_back(i);
    }
  }
  
  std::sort(LocalMax_index.begin(), LocalMax_index.end(),
            [&R](int a, int b) { return R(a) > R(b); });
  
  std::vector<int> tao;
  for (int idx : LocalMax_index) {
    tao.push_back(idx);
  }
  
  int tao1 = (tao.size() > 0) ? tao[0] : 0;
  int tao2 = (tao.size() > 1) ? tao[1] : 0;
  
  // 计算r0和r
  int r0 = std::min(tao1, tao2);
  int r = std::max(tao1, tao2) - std::min(tao1, tao2);
  arma::vec R_out = R.subvec(1, j0);  
  
  return Rcpp::List::create(
    Rcpp::Named("r0") = r0,
    Rcpp::Named("r") = r,
    Rcpp::Named("R") = R_out,
    Rcpp::Named("tao") = tao
  );
}

//' @title Calculating the loading matrices of factors with Rcpp
//' @description Rcpp version of \code{Loadings_Factor}. See more in \code{Loadings_Factor}.
//' @param data a data frame or matrix 
//' @param k0 a pre-specified integer used in calculating the numbers of factors. It works only when r is \code{NULL}.
//' @param j0 a pre-specified integer used in calculating the numbers of factors. It works only when r is \code{NULL}.
//' @param r the numbers of factors. The first value is the numbers of common factors and the second value is the numbers of specific factors. When r is \code{NULL}, the function \code{FactorNum} will be used to calculate the numbers of the factors.
//' @return a list consists of \code{Loadings_CommonFactor} and \code{Loadings_SpecificFactor}, which represent the loadings of common factors and specific factors, respectively.
//' @examples
//' \dontrun{
//' data(CSI300_2)
//' Loadings_Factor_cpp(CSI300_2)   
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List Loadings_Factor_cpp(Rcpp::RObject data, int k0 = 5, int j0 = -1, 
                               Rcpp::Nullable<Rcpp::IntegerVector> r = R_NilValue) {
  arma::mat data0 = flexible_reader(data);
  int p = data0.n_cols;
  
  if (j0 == -1) {
    j0 = std::ceil(p / 4.0);
  }
  
  arma::mat M1(p, p, arma::fill::zeros);
  for (int i = 0; i <= k0; i++) {
    arma::mat Sigma = Sigma_k_cpp(data0, i);
    M1 += Sigma * Sigma.t();
  }
  
  int r0, r1;
  if (r.isNull()) {
    Rcpp::List factorNum = FactorNum_cpp(data, k0, j0);
    r0 = factorNum["r0"];
    r1 = factorNum["r"];
  } else {
    Rcpp::IntegerVector r_vec(r);
    r0 = r_vec[0];
    r1 = r_vec[1];
  }
  
  r0 = std::min(r0, p);
  r1 = std::min(r1, p);
  
  arma::vec eigval1;
  arma::mat eigvec1;
  
  if (!arma::eig_sym(eigval1, eigvec1, M1)) {
    Rcpp::stop("Eigenvalue decomposition failed for M1");
  }
  
  arma::uvec indices = arma::sort_index(eigval1, "descend");
  eigval1 = eigval1.elem(indices);
  eigvec1 = eigvec1.cols(indices);
  
  arma::mat A = eigvec1.cols(0, r0 - 1);
  
  arma::mat I = arma::eye(p, p);
  arma::mat data1 = data0 * (I - A * A.t());
  
  arma::mat M2(p, p, arma::fill::zeros);
  for (int i = 0; i <= k0; i++) {
    arma::mat Sigma = Sigma_k_cpp(data1, i);
    M2 += Sigma * Sigma.t();
  }
  
  arma::vec eigval2;
  arma::mat eigvec2;
  
  if (!arma::eig_sym(eigval2, eigvec2, M2)) {
    Rcpp::stop("Eigenvalue decomposition failed for M2");
  }
  
  indices = arma::sort_index(eigval2, "descend");
  eigval2 = eigval2.elem(indices);
  eigvec2 = eigvec2.cols(indices);
  
  arma::mat B = eigvec2.cols(0, r1 - 1);
  
  return Rcpp::List::create(
    Rcpp::Named("Loadings_CommonFactor") = A,
    Rcpp::Named("Loadings_SpecificFactor") = B
  );
}
