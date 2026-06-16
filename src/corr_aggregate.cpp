// Aggregated correlation assembly for SDALGCP2 (bottleneck B1).
//
// Builds the N x N x n_phi array of region-level correlations
//
//   Sigma(phi)_{ij} = sum_{k,l} w_ik w_jl * matern(||x_ik - x_jl||, phi, kappa)
//
// (weighted), or the unweighted mean over point pairs. Each pairwise distance is
// computed once; all phi values are reduced in a tight inner loop with no
// intermediate 3-D array. Parallelised over region pairs with OpenMP.
//
// This mirrors precomputeCorrMatrix() in SDALGCP but avoids R-level outer()/array
// allocation and the pdist dependency.

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Matern correlation for a single (non-negative) distance d with unit variance.
// kappa = 0.5 -> exponential exp(-d/phi); 1.5 and 2.5 have closed forms; other
// kappa fall back to the Bessel-K form via R's bessel function is avoided here,
// so callers should restrict to {0.5, 1.5, 2.5} for the C++ fast path.
static inline double matern_corr(double d, double phi, double kappa) {
  if (d <= 0.0) return 1.0;
  double u = d / phi;
  if (kappa == 0.5) {
    return std::exp(-u);
  } else if (kappa == 1.5) {
    double a = std::sqrt(3.0) * u;
    return (1.0 + a) * std::exp(-a);
  } else if (kappa == 2.5) {
    double a = std::sqrt(5.0) * u;
    return (1.0 + a + a * a / 3.0) * std::exp(-a);
  }
  // Fallback: exponential (callers needing general kappa use the R path).
  return std::exp(-u);
}

//' Aggregated correlation array (C++)
//'
//' @param coords list of length N; element i is an n_i x 2 numeric matrix of
//'   candidate-point coordinates inside region i.
//' @param weights list of length N of weight vectors (each summing to 1), or an
//'   empty list for the unweighted (mean) case.
//' @param phi numeric vector of spatial scale parameters.
//' @param kappa Matern smoothness (0.5, 1.5 or 2.5 use closed forms).
//' @param weighted logical; TRUE for population-weighted aggregation.
//' @param nthreads number of OpenMP threads (<=0 uses the OpenMP default).
//' @return a numeric array of dimension N x N x length(phi).
//' @keywords internal
// [[Rcpp::export]]
arma::cube corr_aggregate_cpp(const Rcpp::List& coords,
                              const Rcpp::List& weights,
                              const arma::vec& phi,
                              double kappa,
                              bool weighted,
                              int nthreads = 0) {
  const int N = coords.size();
  const int P = phi.n_elem;

  // Pull coordinates (and weights) into std::vectors of arma objects once.
  std::vector<arma::mat> X(N);
  std::vector<arma::vec> W(N);
  for (int i = 0; i < N; ++i) {
    X[i] = as<arma::mat>(coords[i]);
    if (weighted) W[i] = as<arma::vec>(weights[i]);
  }

  arma::cube R(N, N, P, arma::fill::zeros);

  // Build the list of upper-triangle region pairs (i <= j) to balance threads.
  std::vector<std::pair<int,int> > pairs;
  pairs.reserve((size_t)N * (N + 1) / 2);
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j)
      pairs.push_back(std::make_pair(i, j));
  const long long npairs = (long long)pairs.size();

#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
#endif

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (long long idx = 0; idx < npairs; ++idx) {
    const int i = pairs[idx].first;
    const int j = pairs[idx].second;

    const arma::mat& Xi = X[i];
    const arma::mat& Xj = X[j];
    const int ni = Xi.n_rows;
    const int nj = Xj.n_rows;

    // Accumulate the reduction for every phi.
    arma::vec acc(P, arma::fill::zeros);
    const double inv_npair = 1.0 / ((double)ni * (double)nj);

    for (int k = 0; k < ni; ++k) {
      const double xk = Xi(k, 0), yk = Xi(k, 1);
      const double wk = weighted ? W[i](k) : 0.0;
      for (int l = 0; l < nj; ++l) {
        const double dx = xk - Xj(l, 0);
        const double dy = yk - Xj(l, 1);
        const double d  = std::sqrt(dx * dx + dy * dy);
        const double wkl = weighted ? wk * W[j](l) : inv_npair;
        for (int p = 0; p < P; ++p) {
          acc(p) += wkl * matern_corr(d, phi(p), kappa);
        }
      }
    }

    for (int p = 0; p < P; ++p) {
      R(i, j, p) = acc(p);
      R(j, i, p) = acc(p);
    }
  }

  return R;
}

//' Aggregated cross-covariance between prediction points and regions (C++)
//'
//' Mirrors compute_cross_cov() for continuous prediction (bottleneck B5):
//' returns an n_pred x N matrix with entries
//'   sum_l w_jl * matern(||x_pred - x_jl||, phi, kappa)  (weighted) or the mean.
//' @param pred n_pred x 2 matrix of prediction coordinates.
//' @param coords list of N region point matrices.
//' @param weights list of N weight vectors, or empty for unweighted.
//' @param phi single spatial scale parameter.
//' @param kappa Matern smoothness.
//' @param weighted logical.
//' @param nthreads OpenMP threads (<=0 default).
//' @keywords internal
// [[Rcpp::export]]
arma::mat cross_cov_cpp(const arma::mat& pred,
                        const Rcpp::List& coords,
                        const Rcpp::List& weights,
                        double phi,
                        double kappa,
                        bool weighted,
                        int nthreads = 0) {
  const int N = coords.size();
  const int npred = pred.n_rows;

  std::vector<arma::mat> X(N);
  std::vector<arma::vec> W(N);
  for (int j = 0; j < N; ++j) {
    X[j] = as<arma::mat>(coords[j]);
    if (weighted) W[j] = as<arma::vec>(weights[j]);
  }

  arma::mat Rc(npred, N, arma::fill::zeros);

#ifdef _OPENMP
  if (nthreads > 0) omp_set_num_threads(nthreads);
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < npred; ++i) {
    const double xp = pred(i, 0), yp = pred(i, 1);
    for (int j = 0; j < N; ++j) {
      const arma::mat& Xj = X[j];
      const int nj = Xj.n_rows;
      double s = 0.0;
      for (int l = 0; l < nj; ++l) {
        const double dx = xp - Xj(l, 0);
        const double dy = yp - Xj(l, 1);
        const double d  = std::sqrt(dx * dx + dy * dy);
        const double w  = weighted ? W[j](l) : (1.0 / (double)nj);
        s += w * matern_corr(d, phi, kappa);
      }
      Rc(i, j) = s;
    }
  }
  return Rc;
}
