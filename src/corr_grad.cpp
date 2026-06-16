// Aggregated correlation AND its first/second derivatives w.r.t. the scale phi,
// at a single phi (for continuous-phi "direct" MCML, task #6). Exponential kernel
// (Matern kappa = 1/2) only, where the phi-derivatives have clean closed forms:
//
//   rho(d)       = exp(-d/phi)
//   d rho/dphi   = exp(-d/phi) * (d / phi^2)
//   d2 rho/dphi2 = exp(-d/phi) * (d^2/phi^4 - 2 d/phi^3)
//
// and the aggregated quantities are the w_ik w_jl - weighted double sums of these
// over candidate points (the discretised double integral and its derivatives).

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

//' Aggregated correlation and its phi-derivatives at one phi (C++, exponential)
//'
//' @param coords list of N candidate-point matrices (n_i x 2).
//' @param weights list of N weight vectors (each summing to 1), or empty for the
//'   unweighted (mean) case.
//' @param phi spatial scale (> 0).
//' @param weighted logical; population-weighted aggregation.
//' @param nthreads OpenMP threads (<= 0 = default).
//' @return list with N x N matrices \code{R}, \code{dR} (dR/dphi) and \code{d2R}
//'   (d2R/dphi2).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List corr_and_grad_cpp(const Rcpp::List& coords,
                             const Rcpp::List& weights,
                             double phi,
                             bool weighted,
                             int nthreads = 0) {
  const int N = coords.size();
  std::vector<arma::mat> X(N);
  std::vector<arma::vec> W(N);
  for (int i = 0; i < N; ++i) {
    X[i] = as<arma::mat>(coords[i]);
    if (weighted) W[i] = as<arma::vec>(weights[i]);
  }

  arma::mat R(N, N, arma::fill::zeros);
  arma::mat dR(N, N, arma::fill::zeros);
  arma::mat d2R(N, N, arma::fill::zeros);

  std::vector<std::pair<int,int> > pairs;
  pairs.reserve((size_t)N * (N + 1) / 2);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      pairs.push_back(std::make_pair(i, j));
  const long long npairs = (long long)pairs.size();

  const double p2 = phi * phi, p3 = p2 * phi, p4 = p2 * p2;

#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
  #pragma omp parallel for schedule(dynamic)
#endif
  for (long long idx = 0; idx < npairs; ++idx) {
    const int i = pairs[idx].first, j = pairs[idx].second;
    const arma::mat& Xi = X[i]; const arma::mat& Xj = X[j];
    const int ni = Xi.n_rows, nj = Xj.n_rows;
    const double inv_npair = 1.0 / ((double)ni * (double)nj);

    double sR = 0.0, sD = 0.0, sD2 = 0.0;
    for (int k = 0; k < ni; ++k) {
      const double xk = Xi(k, 0), yk = Xi(k, 1);
      const double wk = weighted ? W[i](k) : 0.0;
      for (int l = 0; l < nj; ++l) {
        const double dx = xk - Xj(l, 0), dy = yk - Xj(l, 1);
        const double d  = std::sqrt(dx * dx + dy * dy);
        const double e  = std::exp(-d / phi);
        const double w  = weighted ? wk * W[j](l) : inv_npair;
        sR  += w * e;
        sD  += w * e * (d / p2);
        sD2 += w * e * (d * d / p4 - 2.0 * d / p3);
      }
    }
    R(i, j) = R(j, i) = sR;
    dR(i, j) = dR(j, i) = sD;
    d2R(i, j) = d2R(j, i) = sD2;
  }

  return Rcpp::List::create(Rcpp::Named("R") = R,
                            Rcpp::Named("dR") = dR,
                            Rcpp::Named("d2R") = d2R);
}
