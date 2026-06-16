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

// Matern correlation and its first/second phi-derivatives at distance d, with
// a = sqrt(2 kappa) d / phi (a = d/phi for kappa = 1/2). Supports kappa in
// {1/2, 3/2, 5/2}; all are smooth at d = 0 (a = 0 -> rho = 1, derivatives 0).
static inline void matern_dphi(double d, double phi, double kappa,
                               double& rho, double& drho, double& d2rho) {
  const double phi2 = phi * phi;
  if (kappa == 0.5) {
    double a = d / phi, e = std::exp(-a);
    rho = e; drho = a * e / phi; d2rho = a * (a - 2.0) * e / phi2;
  } else if (kappa == 1.5) {
    double a = std::sqrt(3.0) * d / phi, e = std::exp(-a);
    rho = (1.0 + a) * e;
    drho = a * a * e / phi;
    d2rho = a * a * (a - 3.0) * e / phi2;
  } else {  // kappa == 2.5
    double a = std::sqrt(5.0) * d / phi, e = std::exp(-a);
    rho = (1.0 + a + a * a / 3.0) * e;
    drho = a * a * (1.0 + a) * e / (3.0 * phi);
    d2rho = a * a * (a * a - 3.0 * a - 3.0) * e / (3.0 * phi2);
  }
}

//' Aggregated correlation and its phi-derivatives at one phi (C++, Matern)
//'
//' @param coords list of N candidate-point matrices (n_i x 2).
//' @param weights list of N weight vectors (each summing to 1), or empty for the
//'   unweighted (mean) case.
//' @param phi spatial scale (> 0).
//' @param kappa Matern smoothness; one of 0.5, 1.5, 2.5.
//' @param weighted logical; population-weighted aggregation.
//' @param nthreads OpenMP threads (<= 0 = default).
//' @return list with N x N matrices \code{R}, \code{dR} (dR/dphi) and \code{d2R}
//'   (d2R/dphi2).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List corr_and_grad_cpp(const Rcpp::List& coords,
                             const Rcpp::List& weights,
                             double phi,
                             double kappa,
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
        const double w  = weighted ? wk * W[j](l) : inv_npair;
        double rho, drho, d2rho;
        matern_dphi(d, phi, kappa, rho, drho, d2rho);
        sR  += w * rho;
        sD  += w * drho;
        sD2 += w * d2rho;
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
