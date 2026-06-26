// Laplace mode + Langevin-Hastings (MALA) sampler for the SDA-LGCP latent field.
// Poisson, non-nested case — the path the spatial fit uses.
//
// Two pieces:
//   laplace_mode_poisson_cpp(): conditional mode of [S | Y] by damped Newton, plus
//     the Laplace covariance Sigma.tilde = (Sigma^{-1} + diag(m e^S))^{-1}.
//   mala_poisson_cpp(): the adaptive MALA loop. Given the same mode/Sigma.tilde
//     and RNG seed the draws are reproducible (draw order: d normals then one
//     uniform per iteration).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Conditional mode and Laplace covariance for [S | Y], Poisson, non-nested (C++)
//'
//' @param y count vector.
//' @param m offset vector (e.g. expected counts / population).
//' @param mu prior mean vector of the latent field.
//' @param Sigma prior covariance matrix.
//' @param tol convergence tolerance on the gradient infinity-norm.
//' @param maxit maximum Newton iterations.
//' @return list with \code{mode} and \code{Sigma_tilde}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List laplace_mode_poisson_cpp(const arma::vec& y,
                                    const arma::vec& m,
                                    const arma::vec& mu,
                                    const arma::mat& Sigma,
                                    double tol = 1e-8,
                                    int maxit = 100) {
  arma::mat Sigma_inv = arma::inv_sympd(Sigma);
  arma::vec S = mu;

  double f_old = -arma::datum::inf;
  for (int it = 0; it < maxit; ++it) {
    arma::vec eS = arma::exp(S);
    arma::vec h = m % eS;                      // m * exp(S)
    arma::vec diff = S - mu;
    arma::vec grad = (y - h) - Sigma_inv * diff;

    // Negative Hessian of the (concave) objective: Sigma_inv + diag(h), SPD.
    arma::mat negH = Sigma_inv;
    negH.diag() += h;

    arma::vec step = arma::solve(negH, grad, arma::solve_opts::likely_sympd);

    // Damped Newton with simple backtracking to guarantee ascent.
    double t = 1.0;
    double f_cur = arma::accu(y % S - h) - 0.5 * arma::dot(diff, Sigma_inv * diff);
    arma::vec S_new = S;
    for (int ls = 0; ls < 30; ++ls) {
      S_new = S + t * step;
      arma::vec eSn = arma::exp(S_new);
      arma::vec dn = S_new - mu;
      double f_new = arma::accu(y % S_new - m % eSn) - 0.5 * arma::dot(dn, Sigma_inv * dn);
      if (f_new >= f_cur) { f_cur = f_new; break; }
      t *= 0.5;
    }
    S = S_new;

    if (arma::norm(grad, "inf") < tol || std::abs(f_cur - f_old) < tol * (std::abs(f_cur) + tol)) {
      f_old = f_cur;
      break;
    }
    f_old = f_cur;
  }

  arma::vec h = m % arma::exp(S);
  arma::mat negH = Sigma_inv;
  negH.diag() += h;
  arma::mat Sigma_tilde = arma::inv_sympd(negH);

  return Rcpp::List::create(Rcpp::Named("mode") = S,
                            Rcpp::Named("Sigma_tilde") = Sigma_tilde);
}

//' Adaptive MALA sampler for [S | Y], Poisson, non-nested (C++)
//'
//' Draw order per iteration is \code{d} normals then one uniform, giving
//' reproducible results under a common seed and the same mode/Sigma.tilde.
//'
//' @param y,m,mu,Sigma data and prior as in \code{laplace_mode_poisson_cpp}.
//' @param mode,Sigma_tilde Laplace mode and covariance (preconditioner).
//' @param n_sim,burnin,thin MCMC length controls.
//' @param h_init initial step size; if not finite, \code{1.65 / d^(1/6)} is used.
//' @param c1,c2 step-size adaptation constants.
//' @return list with \code{samples} (kept x d matrix of S draws) and \code{h}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mala_poisson_cpp(const arma::vec& y,
                            const arma::vec& m,
                            const arma::vec& mu,
                            const arma::mat& Sigma,
                            const arma::vec& mode,
                            const arma::mat& Sigma_tilde,
                            int n_sim, int burnin, int thin,
                            double h_init, double c1, double c2) {
  const arma::uword d = y.n_elem;

  // Preconditioning matrices (mirror the R code exactly).
  arma::mat L = arma::chol(Sigma_tilde, "lower");      // L L' = Sigma_tilde
  arma::mat Sigma_inv = arma::inv_sympd(Sigma);
  arma::mat Sigma_W_inv = L.t() * Sigma_inv * L;       // = solve(A Sigma A'), A = L^{-1}
  arma::vec mu_W = arma::solve(arma::trimatl(L), mu - mode);  // A (mu - mode)

  double h = std::isfinite(h_init) ? h_init : 1.65 / std::pow((double)d, 1.0 / 6.0);

  // Log conditional density in W-space and its Langevin gradient.
  auto cond_dens = [&](const arma::vec& W, const arma::vec& S) {
    arma::vec dW = W - mu_W;
    double llik = arma::accu(y % S - m % arma::exp(S));
    return -0.5 * arma::dot(dW, Sigma_W_inv * dW) + llik;
  };
  auto lang_grad = [&](const arma::vec& W, const arma::vec& S) {
    arma::vec dW = W - mu_W;
    arma::vec gS = y - m % arma::exp(S);
    return arma::vec(-Sigma_W_inv * dW + L.t() * gS);
  };

  arma::vec W_curr(d, arma::fill::zeros);
  arma::vec S_curr = L * W_curr + mode;
  arma::vec mean_curr = W_curr + (h * h / 2.0) * lang_grad(W_curr, S_curr);
  double lp_curr = cond_dens(W_curr, S_curr);

  const int n_keep = (n_sim - burnin) / thin;
  arma::mat sim(n_keep, d, arma::fill::zeros);
  arma::vec h_vec(n_sim, arma::fill::zeros);
  long acc = 0;

  Rcpp::RNGScope scope;  // ties into R's RNG stream

  for (int i = 1; i <= n_sim; ++i) {
    arma::vec z(d);
    for (arma::uword k = 0; k < d; ++k) z(k) = R::norm_rand();  // d normals, in order
    arma::vec W_prop = mean_curr + h * z;
    arma::vec S_prop = L * W_prop + mode;
    arma::vec mean_prop = W_prop + (h * h / 2.0) * lang_grad(W_prop, S_prop);
    double lp_prop = cond_dens(W_prop, S_prop);

    double dprop_curr = -arma::accu(arma::square(W_prop - mean_curr)) / (2.0 * h * h);
    double dprop_prop = -arma::accu(arma::square(W_curr - mean_prop)) / (2.0 * h * h);
    double log_prob = lp_prop + dprop_prop - lp_curr - dprop_curr;

    if (std::log(R::unif_rand()) < log_prob) {     // one uniform
      acc += 1;
      W_curr = W_prop; S_curr = S_prop; lp_curr = lp_prop; mean_curr = mean_prop;
    }

    if (i > burnin && ((i - burnin) % thin) == 0) {
      sim.row((i - burnin) / thin - 1) = S_curr.t();
    }
    h = std::max(0.0, h + c1 * std::pow((double)i, -c2) * ((double)acc / i - 0.57));
    h_vec(i - 1) = h;
  }

  return Rcpp::List::create(Rcpp::Named("samples") = sim,
                            Rcpp::Named("h") = h_vec,
                            Rcpp::Named("acc_rate") = (double)acc / n_sim);
}
