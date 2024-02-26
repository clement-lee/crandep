// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
using namespace arma;





// 00) Prelim
const double lnan(const double l) {
  return (l != l) ? -INFINITY : l;
}

const double lr1() {
  return log(runif(1)[0]);
}

template <class T>
const NumericVector tv(T x) {
  return NumericVector::create(x);
}

const IntegerVector ti(const int x) {
  return IntegerVector::create(x);
}

const LogicalVector tl(const bool x) {
  return LogicalVector::create(x);
}

const double ldnorm(const double x, const double mean, const double sd) {
  return dnorm(tv(x), mean, sd, true)[0];
}

const double ldgamma(const double x, const double shape, const double rate) {
  return dgamma(tv(x), shape, 1.0 / rate, true)[0];
}

const NumericVector lpgamma(const NumericVector x, const double shape, const double rate) {
  return pgamma(x, shape, 1.0 / rate, true, true);
}

const double ldbeta(const double x, const double a, const double b) {
  return dbeta(tv(x), a, b, true)[0];
}

const double ldunif(const double x, const double a, const double b) {
  return dunif(tv(x), a, b, true)[0];
}

template <class T>
void update(T & par_curr, T par_prop, double & lpost_curr, const double lpost_prop, double & s, const int i, const int burnin, bool & accept_reject, const double invt = 1.0, const double factor = 10.0) {
  const double lalpha = invt * (lpost_prop - lpost_curr);
  accept_reject = lalpha > lr1();
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  if (i < burnin) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
}

const double sd_curr(const vec x, const int i) {
  return sqrt(as_scalar(cov(x.head(i))));
}

const double cor_curr(const vec x, const vec y, const int i) {
  return as_scalar(cor(x.head(i), y.head(i)));
}

const int sample_w(const IntegerVector seq, const NumericVector weights) {
  // sample 1 value from seq with weights
  return Rcpp::RcppArmadillo::sample(seq, 1, true, weights)[0];
}

const double lse(const NumericVector x) {
  return log(sum(exp(x)));
}

DataFrame df_scalars(const int iter,
                     const int thin,
                     const int burn,
                     const int freq,
                     const bool mc3) {
  return
    DataFrame::create(Named("iter") = ti(iter),
                      Named("thin") = ti(thin),
                      Named("burn") = ti(burn),
                      Named("freq") = ti(freq),
                      Named("mc3") = tl(mc3));
}

const IntegerVector for_Smix(const int u, const int xmin) {
  const uvec x0 = regspace<uvec>(1, 9);
  uvec x(x0);
  x = join_cols(x, x0 * pow(10, 1));
  x = join_cols(x, x0 * pow(10, 2));
  x = join_cols(x, x0 * pow(10, 3));
  x = join_cols(x, x0 * pow(10, 4));
  x = join_cols(x, x0 * pow(10, 5));
  IntegerVector x1 = wrap(x);
  x1 = x1[(x1 <= u) & (x1 >= xmin)];
  x1.insert(x1.end(), u);
  return x1;
}

const bool ispm1(const double x, const double precision = 1.0e-10) {
  // TRUE if x=+/-1.0 up to floating-point prec
  return abs(abs(x) - 1.0) < precision;
}

const double sqrt1mx2(const double x) {
  // sqrt(1 - x^2)
  return ispm1(x) ? 0.0 : sqrt(1.0 - x * x);
}

const double intdiv(const int a, const int b) {
  return (double) a / (double) b;
}

const double odds(const double p) {
  return p / (1.0 - p);
}




// 01) polylogarithm only
const double lnc_pol(const double alpha,
                     const double gamma,
                     const int xmin,
                     const int xmax) {
  // log of normalising constant for polylog
  // gamma & xmin have to be +ve
  const IntegerVector y0 = tail(seq_len(xmax), xmax - xmin + 1); // xmin:xmax
  const NumericVector y(y0), ly = log(y);
  const double
    xs(xmin), // x_star: (approx) x w/ max value
    lxs = log(xs),
    lnc = lse(- alpha * (ly - lxs) - gamma * (y - xs))
    - alpha * lxs - gamma * xs;
  return lnc;
}

//' Probability mass function (PMF) of Zipf-polylog distribution
//'
//' \code{dpol} returns the PMF at x for the Zipf-polylog distribution with parameters (alpha, theta). The distribution is reduced to the discrete power law when theta = 1.
//'
//' The PMF is proportional to x^(-alpha) * theta^x. It is normalised in order to be a proper PMF.
//' @param x Vector of positive integers
//' @param alpha Real number greater than 1
//' @param theta Real number in (0, 1]
//' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
//' @return A numeric vector of the same length as x
//' @examples
//' dpol(c(1,2,3,4,5), 1.2, 0.5)
//' @seealso \code{\link{Spol}} for the corresponding survival function, \code{\link{dmix2}} and \code{\link{dmix3}} for the PMFs of the 2-component and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector dpol(const IntegerVector x,
                         const double alpha,
                         const double theta,
                         const int xmax = 100000) {
  // density for polylogarithm
  // x are desired values, not data
  if (is_true(any(x <= 0))) {
    stop("dpol: all of x has to be +ve integers.");
  }
  if (theta <= 0.0 || theta > 1.0) {
    stop("dpol: theta has to be in (0.0, 1.0].");
  }
  else if (theta == 1.0 && alpha <= 1.0) {
    stop("dpol: alpha has to be greater than 1.0 when theta is 1.0.");
  }
  const double gamma = -log(theta);
  const NumericVector x0(x),
    l = - alpha * log(x0) - x0 * gamma
    - lnc_pol(alpha, gamma, min(x), xmax);
  return exp(l);
}

//' Survival function of Zipf-polylog distribution
//'
//' \code{Spol} returns the survival function at x for the Zipf-polylog distribution with parameters (alpha, theta). The distribution is reduced to the discrete power law when theta = 1.
//'
//' @param x Vector of positive integers
//' @param alpha Real number greater than 1
//' @param theta Real number in (0, 1]
//' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
//' @return A numeric vector of the same length as x
//' @examples
//' Spol(c(1,2,3,4,5), 1.2, 0.5)
//' @seealso \code{\link{dpol}} for the corresponding probability mass function, \code{\link{Smix2}} and \code{\link{Smix3}} for the survival functions of the 2-component and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector Spol(const IntegerVector x,
                         const double alpha,
                         const double theta,
                         const int xmax = 100000) {
  // survival for polylogarithm
  // x are the desired values, not data
  if (is_true(any(x <= 0))) {
    stop("Spol: all of x has to be +ve integers.");
  }
  if (theta <= 0.0 || theta > 1.0) {
    stop("Spol: theta has to be in (0.0, 1.0].");
  }
  else if (theta == 1.0 && alpha <= 1.0) {
    stop("Spol: alpha has to be greater than 1.0 when theta is 1.0.");
  }
  const int xmin = min(x);
  const double gamma = -log(theta);
  const NumericVector x0(x);
  IntegerVector seq1 = seq_len(xmax); // 1 to xmax
  seq1 = tail(seq1, xmax - xmin + 1); // xmin to xmax
  const NumericVector
    x1(seq1), // xmin to xmax
    l1 = - alpha * log(x1) - x1 * gamma
    - lnc_pol(alpha, gamma, xmin, xmax),
    cf1 = cumsum(exp(l1));
  NumericVector S(x.size());
  for (int i = 0; i < x.size(); i++) {
    S[i] = 1.0 - cf1[x[i] - xmin] / cf1[xmax - xmin];
  }
  return S;
}

// [[Rcpp::export]]
const double llik_pol(const NumericVector par,
                      const IntegerVector x,
                      const IntegerVector count,
                      const bool powerlaw,
                      const int xmax) {
  // polylogarithm for the whole data
  // x are the unique values (> 1) w/ freq count
  if (x.size() != count.size()) {
    stop("llik_pol: lengths of x & count have to be equal.");
  }
  if (is_true(any(x <= 0))) {
    stop("llik_pol: all of x has to be +ve integers.");
  }
  if (par.size() != 2) {
    stop("llik_pol: length of par has to be 2.");
  }
  const double
    alpha = par[0], theta = powerlaw ? 1.0 : par[1],
    gamma = -log(theta), n(sum(count));
  const int xmin = min(x);
  const NumericVector x0(x), c0(count);
  double l;
  if (theta <= 0.0 || theta > 1.0 ||
      (powerlaw && (alpha <= 1.0))) {
    l = -INFINITY;
  }
  else {
    l = - alpha * sum(c0 * log(x0))
    - gamma * sum(c0 * x0)
    - n * lnc_pol(alpha, gamma, xmin, xmax);
  }
  return lnan(l);
}

// [[Rcpp::export]]
const double lpost_pol(const IntegerVector x,
                       const IntegerVector count,
                       const double alpha,
                       const double theta,
                       const double a_alpha,
                       const double b_alpha,
                       const double a_theta,
                       const double b_theta,
                       const double powerlaw,
                       const int xmax,
                       double & llik,
                       const double invt = 1.0) {
  // checks in llik_pol
  double l;
  if (theta <= 0.0 || theta > 1.0 ||
      (powerlaw && (alpha <= 1.0))) {
    l = -INFINITY;
  }
  else {
    const NumericVector par = NumericVector::create(alpha, theta);
    llik = llik_pol(par, x, count, powerlaw, xmax);
    // invt for marginal likelihood, not parallel tempering
    l = llik * invt +
      ldnorm(alpha, a_alpha, b_alpha) +
      (powerlaw ? 0.0 : ldbeta(theta, a_theta, b_theta));
  }
  return lnan(l);
}

//' Markov chain Monte Carlo for Zipf-polylog distribution
//'
//' \code{mcmc_pol} returns the samples from the posterior of alpha and theta, for fitting the Zipf-polylog distribution to the data x. The samples are obtained using Markov chain Monte Carlo (MCMC). In the MCMC, a Metropolis-Hastings algorithm is used.
//' @param x Vector of the unique values (positive integers) of the data
//' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
//' @param alpha Real number greater than 1, initial value of the parameter
//' @param theta Real number in (0, 1], initial value of the parameter
//' @param a_alpha Real number, mean of the prior normal distribution for alpha
//' @param b_alpha Positive real number, standard deviation of the prior normal distribution for alpha
//' @param a_theta Positive real number, first parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
//' @param b_theta Positive real number, second parameter of the prior beta distribution for theta; ignored if pr_power = 1.0
//' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
//' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
//' @param pr_power Real number in [0, 1], prior probability of the discrete power law
//' @param iter Positive integer representing the length of the MCMC output
//' @param thin Positive integer representing the thinning in the MCMC
//' @param burn Non-negative integer representing the burn-in of the MCMC
//' @param freq Positive integer representing the frequency of the sampled values being printed
//' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC; default c(1.0) i.e. no Metropolis-coupling
//' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
//' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
//' @seealso \code{\link{mcmc_mix2}} and \code{\link{mcmc_mix3}} for MCMC for the 2-component and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
List mcmc_pol(const IntegerVector x,
              const IntegerVector count,
              double alpha,
              double theta,
              const double a_alpha,
              const double b_alpha,
              const double a_theta,
              const double b_theta,
              const double a_pseudo,
              const double b_pseudo,
              const double pr_power, // prior prob of power law
              const int iter,
              const int thin,
              const int burn,
              const int freq,
              const NumericVector invt,
              const int xmax = 100000) { // should be enough
  // 01) save
  DataFrame
    data =
    DataFrame::create(Named("x") = x,
                      Named("count") = count),
    init =
    DataFrame::create(Named("alpha") = tv(alpha),
                      Named("theta") = tv(theta)),
    hyperpars =
    DataFrame::create(Named("a_alpha") = tv(a_alpha),
                      Named("b_alpha") = tv(b_alpha),
                      Named("a_theta") = tv(a_theta),
                      Named("b_theta") = tv(b_theta)),
    gvs_quants =
    DataFrame::create(Named("a_pseudo") = tv(a_pseudo),
                      Named("b_pseudo") = tv(b_pseudo),
                      Named("pr_power") = tv(pr_power)),
    scalars = df_scalars(iter, thin, burn, freq, invt.size() != 1);
  // 02) checks
  // x are the unique values (> 1) w/ freq count
  if (is_true(any(x <= 0))) {
    stop("mcmc_pol: all of x has to be +ve integers.");
  }
  if (invt.at(0) != 1.0) {
    stop("mcmc_pol: 1st element of inverse temperatures (invt) has to be 1.0.");
  }
  if (is_true(any(invt > 1.0)) || is_true(any(invt <= 0.0))) {
    stop("mcmc_pol: all elements of invt must be in (0.0, 1.0].");
  }
  NumericVector invt0 = clone(invt);
  std::reverse(invt0.begin(), invt0.end());
  if (!std::is_sorted(invt0.begin(), invt0.end())) {
    stop("mcmc_pol: invt has to be in decreasing order.");
  }
  // 03) dimensions const to invt.size()
  int k;
  const int K = invt.size();
  //  const double sd0 = 0.001 / sqrt(2.0), sd1 = 2.38 / sqrt(2.0); // see Roberts & Rosenthal (2008)
  double ldiff;
  vec swap_accept(K-1, fill::zeros), swap_count(K-1, fill::zeros);
  bool powl, accept_reject;
  if (pr_power == 1.0) { // power law
    powl = true; // & stay true
    theta = 1.0; // & stay 1.0
  }
  else if (pr_power == 0.0) { // polylog
    powl = false; // & stay false
  }
  else { // model selection
    powl = false;
  }
  // for model selection
  double theta_pseudo, ak, bk, logA_powl, logA_poly;
  running_stat<double> powl_stat;
  // 04) dimensions increase w/ K
  NumericVector
    alpha_curr(K, alpha), alpha_prop(K),
    theta_curr(K, theta), theta_prop(K),
    llik_curr(K), llik_prop(K),
    lpost_curr(K), lpost_prop(K);
  LogicalVector powl_curr(K, powl);
  auto lpost =
    [x, count, a_alpha, b_alpha, a_theta, b_theta, xmax]
    (const double alpha, const double theta, const bool powerlaw, double & llik) {
      return
        lpost_pol(x, count, alpha, theta,
                  a_alpha, b_alpha,
                  a_theta, b_theta,
                  powerlaw, xmax,
                  llik, 1.0);
    };
  for (int k = 0; k < K; k++) {
    lpost_curr.at(k) = lpost(alpha_curr.at(k), theta_curr.at(k), powl_curr.at(k), llik_curr.at(k));
  }
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  mat alpha_burn(burn, K), theta_burn(burn, K);
  NumericVector sd_alpha(K, 0.1), sd_theta(K, 0.1), cor1(K, 0.1);
  const IntegerVector seqKm1 = seq_len(K - 1) - 1;
  // 05) for saving, dimensions const to K
  IntegerVector powl_vec(iter);
  NumericVector alpha_vec(iter), theta_vec(iter), lpost_vec(iter);
  const int xmin = min(x);
  IntegerVector x0 = for_Smix(max(x) - 1, xmin);
  const int n0 = x0.size();
  NumericVector f0(n0), S0(n0); // fitted
  mat f0_mat(iter, n0), S0_mat(iter, n0);
  // 06) run
  int s, t, l;
  for (t = 0; t < iter * thin + burn; t++) {
    for (k = 0; k < K; k++) {
      // alpha & theta
      alpha_prop.at(k) = rnorm(1, alpha_curr.at(k), sd_alpha.at(k))[0];
      lpost_prop.at(k) = lpost(alpha_prop.at(k), theta_curr.at(k), powl_curr.at(k), llik_prop.at(k));
      update(alpha_curr.at(k), alpha_prop.at(k), lpost_curr.at(k), lpost_prop.at(k), sd_alpha.at(k), t, burn, accept_reject, invt.at(k));
      if (accept_reject) {
        llik_curr.at(k) = llik_prop.at(k);
      }
      if (!powl_curr.at(k)) {
        theta_prop.at(k) = rnorm(1, theta_curr.at(k), sd_theta.at(k))[0];
        lpost_prop.at(k) = lpost(alpha_curr.at(k), theta_prop.at(k), powl_curr.at(k), llik_prop.at(k));
        update(theta_curr.at(k), theta_prop.at(k), lpost_curr.at(k), lpost_prop.at(k), sd_theta.at(k), t, burn, accept_reject, invt.at(k));
        if (accept_reject) {
          llik_curr.at(k) = llik_prop.at(k);
        }
      }
      // model selection
      if (pr_power != 0.0 && pr_power != 1.0 && t >= burn) {
        // only model select after burn-in
        ak = invt.at(k) * (a_pseudo - 1.0) + 1.0;
        bk = invt.at(k) * (b_pseudo - 1.0) + 1.0;
        if (powl_curr.at(k)) { // currently power law
          theta_pseudo = rbeta(1, ak, bk)[0]; // sim from pseudoprior
          logA_powl =
            lpost_curr.at(k) +
            ldbeta(theta_pseudo, ak, bk) +
            log(pr_power);
          logA_poly =
            lpost(alpha_curr.at(k), theta_pseudo, false, llik_prop.at(k)) +
            log(1.0 - pr_power);
          if (lr1() > -log(1.0 + exp(logA_poly - logA_powl))) { // switch to polylog
            powl_curr.at(k) = false;
            theta_curr.at(k) = theta_pseudo;
            lpost_curr.at(k) = logA_poly - log(1.0 - pr_power);
            llik_curr.at(k) = llik_prop.at(k);
          }
        }
        else { // currently polylog
          logA_powl =
            lpost(alpha_curr.at(k), theta_curr.at(k), true, llik_prop.at(k)) +
            ldbeta(theta_curr.at(k), ak, bk) +
            log(pr_power);
          logA_poly =
            lpost_curr.at(k) + log(1.0 - pr_power);
          if (lr1() < -log(1.0 + exp(logA_poly - logA_powl))) { // switch to power law
            powl_curr.at(k) = true;
            theta_curr.at(k) = 1.0;
            lpost_curr.at(k) = logA_powl - ldbeta(theta_curr.at(k), ak, bk) - log(pr_power);
            llik_curr.at(k) = llik_prop.at(k);
          }
        }
        if (k == 0) {
          powl_stat((double) powl_curr.at(k));
        }
      }
    } // loop over k completes
    // 07) Metropolis coupling
    if (K > 1) {
      k = Rcpp::RcppArmadillo::sample(seqKm1, 1, false)[0];
      l = k + 1;
      ldiff = (lpost_curr.at(k) - lpost_curr.at(l)) * (invt.at(l) - invt.at(k));
      if (lr1() < ldiff) {
        swap(alpha_curr.at(k), alpha_curr.at(l));
        swap(theta_curr.at(k), theta_curr.at(l));
        swap(powl_curr.at(k), powl_curr.at(l));
        swap(lpost_curr.at(k), llik_curr.at(l));
        swap(llik_curr.at(k), llik_curr.at(l));
        swap_accept.at(k) += 1.0;
      }
      swap_count.at(k) += 1.0;
    }
    // 08) update cor for printing
    // in case we decide to change to update simultaneously
    for (k = 0; k < K; k++) {
      if (t < burn) {
        if (!powl_curr.at(k)) {
          alpha_burn(t, k) = alpha_curr.at(k);
          theta_burn(t, k) = theta_curr.at(k);
          cor1 = cor_curr(alpha_burn.col(k), theta_burn.col(k), t+1);
        }
      }
    }
    // 09) print
    if ((t + 1) % freq == 0) {
      Rcout << "Iteration " << t + 1;
      Rcout << ": Log-posterior = " << std::setprecision(6) << lpost_curr.at(0) << endl;
      if (t < burn) {
        Rcout << "alpha = " << alpha_curr.at(0) << " (" << sd_alpha.at(0) << ")" << endl;
        Rcout << "theta = " << theta_curr.at(0);
        if (!powl_curr.at(0)) {
          Rcout << " (" << sd_theta.at(0) << ")" << endl;
          Rcout << "cor(alpha, theta) = " << cor1.at(0) << endl;
        }
        else {
          Rcout << endl;
        }
      }
      else {
        Rcout << "alpha = " << alpha_curr.at(0) << endl;
        Rcout << "theta = " << theta_curr.at(0) << endl;
        if (pr_power != 0.0 && pr_power != 1.0) {
          // only print after burn-in & w/ model selection
          Rcout << "power law = " << (powl_curr.at(0) ? "true" : "false") << endl;
          Rcout << "average   = " << powl_stat.mean() << endl;
        }
      }
      if (K > 1) {
        Rcout << "swap rates: " << endl;
        for (k = 0; k < K-1; k++) {
          Rcout << "  b/w inv temp " << std::setprecision(3) << invt.at(k);
          Rcout << " & " << invt.at(k+1) << ": " << swap_accept.at(k) / swap_count.at(k) << endl;
        }
      }
      Rcout << endl;
    }
    // 10) save
    if (t >= burn) {
      s = t - burn + 1;
      if (s % thin == 0) {
        s = s / thin - 1;
        alpha = alpha_curr.at(0);
        theta = theta_curr.at(0);
        alpha_vec[s] = alpha;
        theta_vec[s] = theta;
        powl_vec[s] = (int) powl_curr.at(0);
        lpost_vec[s] = lpost_curr.at(0);
        // prob mass & survival functions
        f0 = dpol(x0, alpha, theta, xmax);
        f0_mat.row(s) = as<rowvec>(f0);
        S0 = Spol(x0, alpha, theta, xmax);
        S0_mat.row(s) = as<rowvec>(S0);
      }
    }
  }
  // 11) output
  const vec p1 = { 0.025, 0.50, 0.975 };
  const mat
    f0_q = quantile(f0_mat, p1, 0),
    S0_q = quantile(S0_mat, p1, 0);
  const double
    po_power = mean((NumericVector) powl_vec),
    bf = odds(po_power) / odds(pr_power);
  gvs_quants["po_power"] = tv(po_power);
  gvs_quants["bf"] = tv(bf);
  DataFrame
    pars =
    DataFrame::create(Named("alpha") = alpha_vec,
                      Named("theta") = theta_vec,
                      Named("powl") = powl_vec,
                      Named("lpost") = lpost_vec),
    fitted =
    DataFrame::create(Named("x") = x0,
                      Named("f_025") = wrap(f0_q.row(0)),
                      Named("f_med") = wrap(f0_q.row(1)),
                      Named("f_975") = wrap(f0_q.row(2)),
                      Named("S_025") = wrap(S0_q.row(0)),
                      Named("S_med") = wrap(S0_q.row(1)),
                      Named("S_975") = wrap(S0_q.row(2)));
  List output =
    List::create(Named("pars") = pars,
                 Named("fitted") = fitted,
                 Named("data") = data,
                 Named("init") = init,
                 Named("hyperpars") = hyperpars,
                 Named("gvs_quants") = gvs_quants,
                 Named("scalars") = scalars,
                 Named("swap_rates") = swap_accept / swap_count,
                 Named("invt") = invt);
  return output;
}





// 02) bulk: polylogarithm
// [[Rcpp::export]]
const double llik_bulk(const NumericVector par,
                       const IntegerVector x,
                       const IntegerVector count,
                       const int v,
                       const int u,
                       const double phil,
                       const bool powerlaw,
                       const bool positive) {
  // polylogarithm b/w v (exclusive) & u (inclusive)
  // x are the unique values (> 1) w/ freq count
  if (x.size() != count.size()) {
    stop("llik_bulk: lengths of x & count have to be equal.");
  }
  if (is_true(any(x <= 0))) {
    stop("llik_bulk: all of x has to be +ve integers.");
  }
  if (par.size() != 2) {
    stop("llik_bulk: length of par has to be 2.");
  }
  const double
    alpha = par[0], theta = powerlaw ? 1.0 : par[1],
    gamma = -log(theta);
  const LogicalVector
    between = (x > v) & (x <= u);
  const NumericVector
    xl(x[between]), cl(count[between]);
  const double nl(sum(cl));
  double l;
  if (v >= u || u >= max(x) ||
      phil <= 0.0 || phil >= 1.0 ||
      //      (powerlaw && (alpha <= 1.0)) ||
      (positive && (alpha <= 0.0)) ||
      theta <= 0.0 || theta > 1.0) {
    // v >= min(x) not included for 3-comp mix
    l = -INFINITY;
  }
  else {
    l =
      nl * log(phil) -
      alpha * sum(cl * log(xl)) - gamma * sum(cl * xl)
      - nl * lnc_pol(alpha, gamma, v + 1, u);
  }
  return lnan(l);
}

// [[Rcpp::export]]
const double lpost_bulk(const NumericVector par,
                        const IntegerVector x,
                        const IntegerVector count,
                        const int v,
                        const int u,
                        const double phil,
                        const double a_alpha,
                        const double b_alpha,
                        const double a_theta,
                        const double b_theta,
                        const bool powerlaw,
                        const bool positive) {
  // checks in llik_bulk
  const double
    alpha = par[0], theta = powerlaw ? 1.0 : par[1];
  double l;
  if (v >= u || u >= max(x) ||
      phil <= 0.0 || phil >= 1.0 ||
      //      (powerlaw && (alpha <= 1.0)) ||
      (positive && (alpha <= 0.0)) ||
      theta <= 0.0 || theta > 1.0) {
    // v >= min(x) not included for 3-comp mix
    // necessary for standalone optim() use
    l = -INFINITY;
  }
  else {
    l =
      llik_bulk(par, x, count, v, u, phil, powerlaw, positive) +
      (powerlaw ? 0.0 : ldbeta(theta, a_theta, b_theta)) +
      ldnorm(alpha, a_alpha, b_alpha);
  }
  return lnan(l);
}





// 03) IGPD only
// [[Rcpp::export]]
const double llik_igpd(const NumericVector par,
                       const IntegerVector x,
                       const IntegerVector count,
                       const int u,
                       const double phiu) {
  // x are the unique values (> 1) w/ freq count
  if (x.size() != count.size()) {
    stop("llik_igp: lengths of x & count have to be equal.");
  }
  if (is_true(any(x <= 0))) {
    stop("llik_igpd: all of x has to be +ve integers.");
  }
  if (par.size() != 2) {
    stop("llik_igpd: length of par has to be 2.");
  }
  const double shape = par[0],
    sigma = par[1],
    sigmau = sigma + shape * u;
  const LogicalVector above = x > u;
  const NumericVector
    xu(x[above]), cu(count[above]);
  const double nu(sum(cu));
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x) ||
      sigma <= 0.0 || sigmau <= 0.0) {
    l = -INFINITY;
  }
  else {
    NumericVector yu, zu;
    if (shape != 0.0) {
      yu = 1.0 + shape / sigmau * (xu - 1.0 - u);
      zu = 1.0 + shape / (sigmau + shape * (xu - 1.0 - u));
      if (is_true(any(yu <= 0.0))) {
        l = -INFINITY;
      }
      else {
        l = sum(cu * log(1.0 - pow(zu, -1.0 / shape))) - 1.0 / shape * sum(cu * log(yu));
      }
    }
    else {
      yu = (xu - 1.0 - u) / sigmau;
      l = nu * log(1.0 - exp(-1.0 / sigmau)) - sum(cu * yu);
    }
    l += nu * log(phiu);
  }
  return lnan(l);
}

// [[Rcpp::export]]
const double lpost_igpd(const NumericVector par,
                        const IntegerVector x,
                        const IntegerVector count,
                        const int u,
                        const double m_shape,
                        const double s_shape,
                        const double a_sigma,
                        const double b_sigma,
                        const double phiu) {
  if (x.size() != count.size()) {
    stop("lpost_igpd: lengths of x & count have to be equal.");
  }
  const double shape = par[0], sigma = par[1];
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x) ||
      sigma <= 0.0) {
    l = -INFINITY;
  }
  else {
    l =
      llik_igpd(par, x, count, u, phiu) +
      ldnorm(shape, m_shape, s_shape) +
      ldgamma(sigma, a_sigma, b_sigma);
  }
  return lnan(l);
}





// 04) 2-component mixture distribution

// [[Rcpp::export]]
const double lpost_mix1(const IntegerVector x,
                        const IntegerVector count,
                        const int u,
                        const double alpha1,
                        const double theta1,
                        const double alpha2,
                        const double a_psiu,
                        const double b_psiu,
                        const double a_alpha1,
                        const double b_alpha1,
                        const double a_theta1,
                        const double b_theta1,
                        const double a_alpha2,
                        const double b_alpha2,
                        const bool positive, // alpha1 bound to be positive?
                        const int xmax = 100000) {
  // T(Z)P below u, power law above u
  if (x.size() != count.size()) {
    stop("lpost_mix1: lengths of x & count have to be equal.");
  }
  const int v = min(x) - 1;
  const LogicalVector above = x > u;
  const IntegerVector xu = x[above], cu = count[above];
  const double nu(sum(cu));
  const NumericVector
    par1 = NumericVector::create(alpha1, theta1),
    par2 = NumericVector::create(alpha2, 1.0);
  const double
    phiu = intdiv(sum(cu), sum(count)),
    phil = 1.0 - phiu,
    psiu = intdiv(cu.size(), count.size());
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x)) {
    l = -INFINITY;
  }
  else {
    l =
      llik_bulk(par1, x, count, v, u, phil, false, positive) +
      llik_pol(par2, xu, cu, true, xmax) +
      nu * log(phiu) + // not in llik_pol()
      ldunif(psiu, a_psiu, b_psiu) +
      ldnorm(alpha1, a_alpha1, b_alpha1) +
      ldbeta(theta1, a_theta1, b_theta1) +
      ldnorm(alpha2, a_alpha2, b_alpha2);
  }
  return lnan(l);
}

//' Markov chain Monte Carlo for TZP-power-law mixture
//'
//' \code{mcmc_mix1} returns the posterior samples of the parameters, for fitting the TZP-power-law mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
//'
//' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The threshold u is treated as a parameter and therefore sampled. The hyperparameters are used in the following priors: u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha1 ~ Normal(mean = a_alpha1, sd = b_alpha1); theta1 ~ Beta(a_theta1, b_theta1); alpha2 ~ Normal(mean = a_alpha2, sd = b_alpha2)
//' @param x Vector of the unique values (positive integers) of the data
//' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
//' @param u_set Positive integer vector of the values u will be sampled from
//' @param u Positive integer, initial value of the threshold
//' @param alpha1 Real number, initial value of the parameter
//' @param theta1 Real number in (0, 1], initial value of the parameter
//' @param alpha2 Real number greater than 1, initial value of the parameter
//' @param a_psiu,b_psiu,a_alpha1,b_alpha1,a_theta1,b_theta1,a_alpha2,b_alpha2 Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
//' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE)?
//' @param iter Positive integer representing the length of the MCMC output
//' @param thin Positive integer representing the thinning in the MCMC
//' @param burn Non-negative integer representing the burn-in of the MCMC
//' @param freq Positive integer representing the frequency of the sampled values being printed
//' @param xmax Scalar (default 100000), positive integer limit for computing the normalising constant
//' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
//' @seealso \code{\link{mcmc_pol}}, \code{\link{mcmc_mix2}} and \code{\link{mcmc_mix3}} for MCMC for the Zipf-polylog, and 2-component and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
List mcmc_mix1(const IntegerVector x,
               const IntegerVector count,
               const IntegerVector u_set,
               int u,
               double alpha1,
               double theta1,
               double alpha2,
               const double a_psiu,
               const double b_psiu,
               const double a_alpha1,
               const double b_alpha1,
               const double a_theta1,
               const double b_theta1,
               const double a_alpha2,
               const double b_alpha2,
               const bool positive,
               const int iter,
               const int thin,
               const int burn,
               const int freq,
               const int xmax = 100000) {
  // 01) save
  DataFrame
    data =
    DataFrame::create(Named("x") = x,
                      Named("count") = count),
    init =
    DataFrame::create(Named("u") = ti(u),
                      Named("alpha1") = tv(alpha1),
                      Named("theta1") = tv(theta1),
                      Named("alpha2") = tv(alpha2)),
    hyperpars =
    DataFrame::create(Named("a_psiu") = tv(a_psiu),
                      Named("b_psiu") = tv(b_psiu),
                      Named("a_alpha1") = tv(a_alpha1),
                      Named("b_alpha1") = tv(b_alpha1),
                      Named("a_theta1") = tv(a_theta1),
                      Named("b_theta1") = tv(b_theta1),
                      Named("a_alpha2") = tv(a_alpha2),
                      Named("b_alpha2") = tv(b_alpha2)),
    scalars = df_scalars(iter, thin, burn, freq, false);
  // 02) checks
  // x are the unique values (> 1) w/ freq count
  if (is_true(any(x <= 0))) {
    stop("mcmc_mix1: all of x has to be +ve integers.");
  }
  if (is_true(all(u_set != u))) {
    stop("mcmc_mix1: u must be in u_set.");
  }
  // 03)
  int k;
  const int K = 1, // for future dev
    m = u_set.size();
  const IntegerVector seqm = seq_len(m) - 1;
  NumericVector w(m), w1(m);
  const double sd0 = 0.001 / sqrt(2.0), sd1 = 2.38 / sqrt(2.0); // see Roberts & Rosenthal (2008)
  bool accept_reject;
  // 04)
  IntegerVector u_curr(K, u);
  NumericVector
    alpha1_curr(K, alpha1), alpha1_prop(K),
    theta1_curr(K, theta1), theta1_prop(K),
    alpha2_curr(K, alpha2), alpha2_prop(K),
    lpost_curr(K), lpost_prop(K);
  auto lpost =
    [x, count,
     a_psiu, b_psiu, a_alpha1, b_alpha1,
     a_theta1, b_theta1, a_alpha2, b_alpha2,
     xmax, positive]
    (const int u, const double alpha1,
     const double theta1, const double alpha2) {
      return
        lpost_mix1(x, count, u, alpha1, theta1, alpha2,
                   a_psiu, b_psiu, a_alpha1, b_alpha1,
                   a_theta1, b_theta1, a_alpha2, b_alpha2,
                   positive, xmax);
    };
  for (k = 0; k < K; k++) {
    lpost_curr.at(k) = lpost(u_curr.at(k), alpha1_curr.at(k), theta1_curr.at(k), alpha2_curr.at(k));
  }
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  mat alpha1_burn(burn, K), theta1_burn(burn, K), alpha2_burn(burn, K);
  NumericVector sd_alpha1(K, 0.1), sd_theta1(K, 0.1), sd_alpha2(K, 0.1), cor1(K, 0.1);
  const IntegerVector seqKm1 = seq_len(K - 1) - 1;
  // 05) for saving
  IntegerVector u_vec(iter);
  NumericVector alpha1_vec(iter), theta1_vec(iter), alpha2_vec(iter), phiu_vec(iter), lpost_vec(iter);
  const int n = sum(count);
  double phiu;
  IntegerVector cu;
  // 06) run
  int s, t;
  for (t = 0; t < iter * thin + burn; t++) {
    for (k = 0; k < K; k++) {
      // u
      std::fill(w.begin(), w.end(), NA_REAL);
      for (s = 0; s < u_set.size(); s++) {
        w[s] = lpost(u_set[s], alpha1_curr.at(k), theta1_curr.at(k), alpha2_curr.at(k));
      }
      w1 = exp(1.0 * (w - max(w)));
      s = sample_w(seqm, w1);
      u_curr.at(k) = u_set[s];
      lpost_curr.at(k) = w[s];
      // alpha1 & theta1
      if (t < burn &&
          (isnan(cor1.at(k)) || ispm1(cor1.at(k)) ||
           isnan(sd_alpha1.at(k)) || isnan(sd_theta1.at(k)) || lr1() < log(0.05))) {
        alpha1_prop.at(k) = rnorm(1, alpha1_curr.at(k), sd0)[0];
        theta1_prop.at(k) = rnorm(1, theta1_curr.at(k), sd0)[0];
      }
      else {
        alpha1_prop.at(k) = rnorm(1, alpha1_curr.at(k), sd1 * sd_alpha1.at(k))[0];
        theta1_prop.at(k) = rnorm(1, theta1_curr.at(k) + sd_theta1.at(k) / sd_alpha1.at(k) * cor1.at(k) * (alpha1_prop.at(k) - alpha1_curr.at(k)), sd1 * sqrt1mx2(cor1.at(k)) * sd_theta1.at(k))[0];
      }
      lpost_prop.at(k) = lpost(u_curr.at(k), alpha1_prop.at(k), theta1_prop.at(k), alpha2_curr.at(k));
      if (1.0 * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
        alpha1_curr.at(k) = alpha1_prop.at(k);
        theta1_curr.at(k) = theta1_prop.at(k);
        lpost_curr.at(k) = lpost_prop.at(k);
      }
      // alpha2
      alpha2_prop.at(k) = rnorm(1, alpha2_curr.at(k), sd_alpha2.at(k))[0];
      lpost_prop.at(k) = lpost(u_curr.at(k), alpha1_curr.at(k), theta1_curr.at(k), alpha2_prop.at(k));
      update(alpha2_curr.at(k), alpha2_prop.at(k), lpost_curr.at(k), lpost_prop.at(k), sd_alpha2.at(k), t, burn, accept_reject, 1.0);
    } // loop over k completes
    // 07)
    // 08) adaptive update sds
    for (k = 0; k < K; k++) {
      if (t < burn) {
        alpha1_burn(t, k) = alpha1_curr.at(k);
        theta1_burn(t, k) = theta1_curr.at(k);
        sd_alpha1.at(k) = sd_curr(alpha1_burn.col(k), t+1);
        sd_theta1.at(k) = sd_curr(theta1_burn.col(k), t+1);
        cor1.at(k) = cor_curr(alpha1_burn.col(k), theta1_burn.col(k), t+1);
      }
    }
    // 09) print
    if ((t + 1) % freq == 0) {
      Rcout << "Iteration " << t + 1;
      Rcout << ": Log-posterior = " << std::setprecision(6) << lpost_curr.at(0) << endl;
      Rcout << "  u   = " << u_curr.at(0) << endl;
      if (t < burn) {
        Rcout << "alpha1 = " << alpha1_curr.at(0) << " (" << sd_alpha1.at(0) << ")" << endl;
        Rcout << "theta1 = " << theta1_curr.at(0);
        Rcout << " (" << sd_theta1.at(0) << ")" << endl;
        Rcout << "cor(alpha1, theta1) = " << cor1.at(0) << endl;
        Rcout << "alpha2 = " << alpha2_curr.at(0) << " (" << sd_alpha2.at(0) << ")" << endl;
      }
      else {
        Rcout << "alpha1 = " << alpha1_curr.at(0) << endl;
        Rcout << "theta1 = " << theta1_curr.at(0) << endl;
        Rcout << "alpha2 = " << alpha2_curr.at(0) << endl;
      }
      Rcout << endl;
    }
    // 10) save
    if (t >= burn) {
      s = t - burn + 1;
      if (s % thin == 0) {
        s = s / thin - 1;
        u = u_curr.at(0);
        alpha1 = alpha1_curr.at(0);
        theta1 = theta1_curr.at(0);
        alpha2 = alpha2_curr.at(0);
        cu = count[x > u];
        phiu = intdiv(sum(cu), n);
        u_vec[s] = u;
        alpha1_vec[s] = alpha1;
        theta1_vec[s] = theta1;
        alpha2_vec[s] = alpha2;
        phiu_vec[s] = phiu;
        lpost_vec[s] = lpost_curr.at(0);
      }
    }
  }
  // 11) output
  DataFrame
    pars =
    DataFrame::create(Named("u") = u_vec,
                      Named("alpha1") = alpha1_vec,
                      Named("theta1") = theta1_vec,
                      Named("alpha2") = alpha2_vec,
                      Named("phiu") = phiu_vec,
                      Named("lpost") = lpost_vec);
  List output =
    List::create(Named("pars") = pars,
                 Named("data") = data,
                 Named("u_set") = u_set,
                 Named("init") = init,
                 Named("hyperpars") = hyperpars,
                 Named("scalars") = scalars);
  return output;
}

//' Probability mass function (PMF) of 2-component discrete extreme value mixture distribution
//'
//' \code{dmix2} returns the PMF at x for the 2-component discrete extreme value mixture distribution. The components below and above the threshold u are the (truncated) Zipf-polylog(alpha,theta) and the generalised Pareto(shape, sigma) distributions, respectively.
//' @param x Vector of positive integers
//' @param u Positive integer representing the threshold
//' @param alpha Real number, first parameter of the Zipf-polylog component
//' @param theta Real number in (0, 1], second parameter of the Zipf-polylog component
//' @param shape Real number, shape parameter of the generalised Pareto component
//' @param sigma Real number, scale parameter of the generalised Pareto component
//' @param phiu Real number in (0, 1), exceedance rate of the threshold u
//' @return A numeric vector of the same length as x
//' @seealso \code{\link{Smix2}} for the corresponding survival function, \code{\link{dpol}} and \code{\link{dmix3}} for the PMFs of the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector dmix2(const IntegerVector x,
                          const int u,
                          const double alpha,
                          const double theta,
                          const double shape,
                          const double sigma,
                          const double phiu) {
  // density for 2-component mixture
  // x are desired values, not data
  // v is required for the dist.
  if (is_true(any(x <= 0))) {
    stop("dmix2: all of x has to be +ve integers.");
  }
  if (theta <= 0.0 || theta > 1.0) {
    stop("dmix2: theta has to be in (0.0, 1.0].");
  }
  if (sigma <= 0.0) {
    stop("dmix2: sigma has to be positive.");
  }
  if (phiu <= 0.0 || phiu >= 1.0) {
    stop("dmix2: phiu has to be in (0.0, 1.0).");
  }
  const int v = min(x) - 1;
  const double sigu = sigma + shape * u;
  const NumericVector
    x0(x),
    y1 = 1.0 + shape / sigu * (x0 - 1.0 - u),
    z1 = 1.0 + shape / (sigu + shape * (x0 - 1.0 - u)),
    fu = exp(log(1.0 - pow(z1, -1.0 / shape)) - 1.0 / shape * log(y1) + log(phiu));
  IntegerVector seq2 = seq_len(u); // 1 to u
  seq2 = tail(seq2, u - v); // (v+1) to u
  const NumericVector x2(seq2); // (v+1) to u
  NumericVector l2 = - alpha * log(x2) + x2 * log(theta);
  l2 = l2 - max(l2);
  const NumericVector cf2 = cumsum(exp(l2));
  NumericVector fl(x.size(), NA_REAL);
  for (int i = 0; i < x.size(); i++) {
    if (x[i] >= (v+1) && x[i] <= u) {
      fl[i] = exp(log(1.0 - phiu) + l2[x[i]-(v+1)] - log(cf2[u-(v+1)]));
    }
  }
  const NumericVector f = ifelse(x <= u, fl, fu);
  return f;
}

//' Survival function of 2-component discrete extreme value mixture distribution
//'
//' \code{Smix2} returns the survival function at x for the 2-component discrete extreme value mixture distribution. The components below and above the threshold u are the (truncated) Zipf-polylog(alpha,theta) and the generalised Pareto(shape, sigma) distributions, respectively.
//' @param x Vector of positive integers
//' @param u Positive integer representing the threshold
//' @param alpha Real number, first parameter of the Zipf-polylog component
//' @param theta Real number in (0, 1], second parameter of the Zipf-polylog component
//' @param shape Real number, shape parameter of the generalised Pareto component
//' @param sigma Real number, scale parameter of the generalised Pareto component
//' @param phiu Real number in (0, 1), exceedance rate of the threshold u
//' @return A numeric vector of the same length as x
//' @seealso \code{\link{dmix2}} for the corresponding probability mass function, \code{\link{Spol}} and \code{\link{Smix3}} for the survival functions of the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector Smix2(const IntegerVector x,
                          const int u,
                          const double alpha,
                          const double theta,
                          const double shape,
                          const double sigma,
                          const double phiu) {
  // survival for 2-component mixture
  // x are the desired values, not data
  // v is required for the dist.
  if (is_true(any(x <= 0))) {
    stop("Smix2: all of x has to be +ve integers.");
  }
  if (theta <= 0.0 || theta > 1.0) {
    stop("Smix2: theta has to be in (0.0, 1.0].");
  }
  if (sigma <= 0.0) {
    stop("Smix2: sigma has to be positive.");
  }
  if (phiu <= 0.0 || phiu >= 1.0) {
    stop("Smix2: phiu has to be in (0.0, 1.0).");
  }
  const int v = min(x) - 1;
  const double sigu = sigma + shape * u;
  const NumericVector x0(x);
  NumericVector y = 1.0 + shape / sigu * (x0 - u);
  y = ifelse(y < 0.0, 0.0, y);
  const NumericVector Su = pow(y, -1.0 / shape) * phiu;
  IntegerVector seq2 = seq_len(u); // 1 to u
  seq2 = tail(seq2, u - v); // (v+1) to u
  const NumericVector x2(seq2); // (v+1) to u
  NumericVector l2 = - alpha * log(x2) + x2 * log(theta);
  l2 = l2 - max(l2);
  const NumericVector cf2 = cumsum(exp(l2));
  NumericVector Sl(x.size());
  for (int i = 0; i < x.size(); i++) {
    if (x[i] >= (v+1) && x[i] < u) {
      Sl[i] = phiu + (1.0 - phiu) * (1.0 - cf2[x[i]-(v+1)] / cf2[u-(v+1)]);
    }
    else if (x[i] == u) {
      Sl[i] = phiu;
    }
    else {
      Sl[i] = NA_REAL;
    }
  }
  const NumericVector S = ifelse(x <= u, Sl, Su);
  return S;
}

// [[Rcpp::export]]
const double lpost_mix2(const IntegerVector x,
                        const IntegerVector count,
                        const int u,
                        const double alpha,
                        const double theta,
                        const double shape,
                        const double sigma,
                        const double a_psiu,
                        const double b_psiu,
                        const double a_alpha,
                        const double b_alpha,
                        const double a_theta,
                        const double b_theta,
                        const double m_shape,
                        const double s_shape,
                        const double a_sigma,
                        const double b_sigma,
                        const bool powerlaw,
                        const bool positive) { // alpha bound to be positive?
  if (x.size() != count.size()) {
    stop("lpost_mix2: lengths of x & count have to be equal.");
  }
  const int v = min(x) - 1;
  const LogicalVector above = x > u;
  const NumericVector
    xu(x[above]), cu(count[above]),
    par1 = NumericVector::create(alpha, theta),
    par2 = NumericVector::create(shape, sigma);
  const double
    phiu = intdiv(sum(cu), sum(count)),
    phil = 1.0 - phiu,
    psiu = intdiv(cu.size(), count.size());
  double l;
  if (u <= 1 || u <= min(x) || u >= max(x)) {
    l = -INFINITY;
  }
  else {
    l =
      llik_bulk(par1, x, count, v, u, phil, powerlaw, positive) +
      llik_igpd(par2, x, count, u, phiu) +
      ldunif(psiu, a_psiu, b_psiu) +
      (powerlaw ? 0.0 : ldbeta(theta, a_theta, b_theta)) +
      ldnorm(shape, m_shape, s_shape) +
      ldgamma(sigma, a_sigma, b_sigma) +
      ldnorm(alpha, a_alpha, b_alpha);
  }
  return lnan(l);
}

//' Markov chain Monte Carlo for 2-component discrete extreme value mixture distribution
//'
//' \code{mcmc_mix2} returns the posterior samples of the parameters, for fitting the 2-component discrete extreme value mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
//'
//' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The threshold u is treated as a parameter and therefore sampled. The hyperparameters are used in the following priors: u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha ~ Normal(mean = a_alpha, sd = b_alpha); theta ~ Beta(a_theta, b_theta); shape ~ Normal(mean = m_shape, sd = s_shape); sigma ~ Gamma(a_sigma, scale = b_sigma). If pr_power = 1.0, the discrete power law (below u) is assumed, and the samples of theta will be all 1.0. If pr_power is in (0.0, 1.0), model selection between the polylog distribution and the discrete power law will be performed within the MCMC.
//' @param x Vector of the unique values (positive integers) of the data
//' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
//' @param u_set Positive integer vector of the values u will be sampled from
//' @param u Positive integer, initial value of the threshold
//' @param alpha Real number greater than 1, initial value of the parameter
//' @param theta Real number in (0, 1], initial value of the parameter
//' @param shape Real number, initial value of the parameter
//' @param sigma Positive real number, initial value of the parameter
//' @param a_psiu,b_psiu,a_alpha,b_alpha,a_theta,b_theta,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
//' @param positive Boolean, is alpha positive (TRUE) or unbounded (FALSE)?
//' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
//' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta in model selection; ignored if pr_power = 1.0
//' @param pr_power Real number in [0, 1], prior probability of the discrete power law (below u)
//' @param iter Positive integer representing the length of the MCMC output
//' @param thin Positive integer representing the thinning in the MCMC
//' @param burn Non-negative integer representing the burn-in of the MCMC
//' @param freq Positive integer representing the frequency of the sampled values being printed
//' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC; default c(1.0) i.e. no Metropolis-coupling
//' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
//' @seealso \code{\link{mcmc_pol}} and \code{\link{mcmc_mix3}} for MCMC for the Zipf-polylog and 3-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
List mcmc_mix2(const IntegerVector x,
               const IntegerVector count,
               const IntegerVector u_set,
               int u,
               double alpha,
               double theta,
               double shape,
               double sigma,
               const double a_psiu,
               const double b_psiu,
               const double a_alpha,
               const double b_alpha,
               const double a_theta,
               const double b_theta,
               const double m_shape,
               const double s_shape,
               const double a_sigma,
               const double b_sigma,
               const bool positive, // alpha +ve / unbound?
               const double a_pseudo,
               const double b_pseudo,
               const double pr_power, // prior prob of power law
               const int iter,
               const int thin,
               const int burn,
               const int freq,
               const NumericVector invt) {
  // 01) save
  DataFrame
    data =
    DataFrame::create(Named("x") = x,
                      Named("count") = count),
    init =
    DataFrame::create(Named("u") = ti(u),
                      Named("alpha") = tv(alpha),
                      Named("theta") = tv(theta),
                      Named("shape") = tv(shape),
                      Named("sigma") = tv(sigma)),
    hyperpars =
    DataFrame::create(Named("a_psiu") = tv(a_psiu),
                      Named("b_psiu") = tv(b_psiu),
                      Named("a_alpha") = tv(a_alpha),
                      Named("b_alpha") = tv(b_alpha),
                      Named("a_theta") = tv(a_theta),
                      Named("b_theta") = tv(b_theta),
                      Named("m_shape") = tv(m_shape),
                      Named("s_shape") = tv(s_shape),
                      Named("a_sigma") = tv(a_sigma),
                      Named("b_sigma") = tv(b_sigma),
                      Named("positive") = tl(positive)),
    gvs_quants =
    DataFrame::create(Named("a_pseudo") = tv(a_pseudo),
                      Named("b_pseudo") = tv(b_pseudo),
                      Named("pr_power") = tv(pr_power)),
    scalars = df_scalars(iter, thin, burn, freq, invt.size() != 1);
  // 02) checks
  // x are the unique values (> 1) w/ freq count
  if (is_true(any(x <= 0))) {
    stop("mcmc_mix2: all of x has to be +ve integers.");
  }
  if (is_true(all(u_set != u))) {
    stop("mcmc_mix2: u must be in u_set.");
  }
  if (invt.at(0) != 1.0) {
    stop("mcmc_mix2: 1st element of inverse temperatures (invt) has to be 1.0.");
  }
  if (is_true(any(invt > 1.0)) || is_true(any(invt <= 0.0))) {
    stop("mcmc_mix2: all elements of invt must be in (0.0, 1.0].");
  }
  NumericVector invt0 = clone(invt);
  std::reverse(invt0.begin(), invt0.end());
  if (!std::is_sorted(invt0.begin(), invt0.end())) {
    stop("mcmc_mix2: invt has to be in decreasing order.");
  }
  // 03) dimensions const to invt.size()
  int k;
  const int K = invt.size(), m = u_set.size();
  const IntegerVector seqm = seq_len(m) - 1;
  NumericVector w(m), w1(m);
  const double sd0 = 0.001 / sqrt(2.0), sd1 = 2.38 / sqrt(2.0); // see Roberts & Rosenthal (2008)
  double ldiff;
  vec swap_accept(K-1, fill::zeros), swap_count(K-1, fill::zeros);
  bool powl, accept_reject;
  if (pr_power == 1.0) { // power law
    powl = true; // & stay true
    theta = 1.0; // & stay 1.0
  }
  else if (pr_power == 0.0) { // polylog
    powl = false; // & stay false
  }
  else { // model selection
    powl = false;
  }
  // for model selection
  double theta_pseudo, ak, bk, logA_powl, logA_poly;
  running_stat<double> powl_stat;
  // 04) dimensions increase w/ K
  IntegerVector u_curr(K, u);
  NumericVector
    alpha_curr(K, alpha), alpha_prop(K),
    theta_curr(K, theta), theta_prop(K),
    shape_curr(K, shape), shape_prop(K),
    sigma_curr(K, sigma), sigma_prop(K),
    lpost_curr(K), lpost_prop(K);
  LogicalVector powl_curr(K, powl);
  auto lpost =
    [x, count, a_psiu, b_psiu,
     a_alpha, b_alpha, a_theta, b_theta,
     m_shape, s_shape, a_sigma, b_sigma,
     positive]
    (const int u,
     const double alpha, const double theta,
     const double shape, const double sigma,
     const bool powerlaw) {
      return lpost_mix2(x, count, u,
                        alpha, theta, shape, sigma,
                        a_psiu, b_psiu,
                        a_alpha, b_alpha, a_theta, b_theta,
                        m_shape, s_shape, a_sigma, b_sigma,
                        powerlaw, positive);
    };
  for (int k = 0; k < K; k++) {
    lpost_curr.at(k) = lpost(u_curr.at(k), alpha_curr.at(k), theta_curr.at(k), shape_curr.at(k), sigma_curr.at(k), powl_curr.at(k));
  }
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  mat alpha_burn(burn, K), theta_burn(burn, K), shape_burn(burn, K), sigma_burn(burn, K);
  NumericVector sd_alpha(K, 0.1), sd_theta(K, 0.1), sd_shape(K, 0.1), sd_sigma(K, 0.1), cor1(K, 0.1), cor2(K, 0.1);
  // these sds & cors not as influential as sd0
  const IntegerVector seqKm1 = seq_len(K - 1) - 1;
  // 05) for saving, dimensions const to K
  IntegerVector u_vec(iter), powl_vec(iter);
  NumericVector alpha_vec(iter), theta_vec(iter), shape_vec(iter), sigma_vec(iter), phiu_vec(iter), lpost_vec(iter);
  const int n = sum(count);
  double phiu;
  const int xmin = min(x);
  IntegerVector x0 = for_Smix(max(x) - 1, xmin), cu;
  const int n0 = x0.size();
  NumericVector f0(n0), S0(n0); // fitted
  mat f0_mat(iter, n0), S0_mat(iter, n0);
  // 06) run
  int s, t, l;
  for (t = 0; t < iter * thin + burn; t++) {
    for (k = 0; k < K; k++) {
      // u
      std::fill(w.begin(), w.end(), NA_REAL);
      for (s = 0; s < u_set.size(); s++) {
        w[s] = lpost(u_set[s], alpha_curr.at(k), theta_curr.at(k), shape_curr.at(k), sigma_curr.at(k), powl_curr.at(k));
      }
      w1 = exp(invt.at(k) * (w - max(w)));
      s = sample_w(seqm, w1); // reused as loop done
      u_curr.at(k) = u_set[s];
      lpost_curr.at(k) = w[s];
      // alpha & theta
      if (powl_curr.at(k)) {
        alpha_prop.at(k) = rnorm(1, alpha_curr.at(k), sd_alpha.at(k))[0];
        lpost_prop.at(k) = lpost(u_curr.at(k), alpha_prop.at(k), theta_curr.at(k), shape_curr.at(k), sigma_curr.at(k), powl_curr.at(k));
        update(alpha_curr.at(k), alpha_prop.at(k), lpost_curr.at(k), lpost_prop.at(k), sd_alpha.at(k), t, burn, accept_reject, invt.at(k));
      }
      else {
        if (t < burn &&
            (isnan(cor1.at(k)) || ispm1(cor1.at(k)) ||
             isnan(sd_alpha.at(k)) || isnan(sd_theta.at(k)) || lr1() < log(0.05))) {
          alpha_prop.at(k) = rnorm(1, alpha_curr.at(k), sd0)[0];
          theta_prop.at(k) = rnorm(1, theta_curr.at(k), sd0)[0];
        }
        else {
          alpha_prop.at(k) = rnorm(1, alpha_curr.at(k), sd1 * sd_alpha.at(k))[0];
          theta_prop.at(k) = rnorm(1, theta_curr.at(k) + sd_theta.at(k) / sd_alpha.at(k) * cor1.at(k) * (alpha_prop.at(k) - alpha_curr.at(k)), sd1 * sqrt1mx2(cor1.at(k)) * sd_theta.at(k))[0];
        }
        lpost_prop.at(k) = lpost(u_curr.at(k), alpha_prop.at(k), theta_prop.at(k), shape_curr.at(k), sigma_curr.at(k), powl_curr.at(k));
        if (invt.at(k) * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
          alpha_curr.at(k) = alpha_prop.at(k);
          theta_curr.at(k) = theta_prop.at(k);
          lpost_curr.at(k) = lpost_prop.at(k);
        }
      }
      // shape & sigma
      if (t < burn &&
          (isnan(cor2.at(k)) || ispm1(cor2.at(k)) ||
           isnan(sd_shape.at(k)) || isnan(sd_sigma.at(k)) || lr1() < log(0.05))) {
        shape_prop.at(k) = rnorm(1, shape_curr.at(k), sd0)[0];
        sigma_prop.at(k) = rnorm(1, sigma_curr.at(k), sd0)[0];
      }
      else {
        shape_prop.at(k) = rnorm(1, shape_curr.at(k), sd1 * sd_shape.at(k))[0];
        sigma_prop.at(k) = rnorm(1, sigma_curr.at(k) + sd_sigma.at(k) / sd_shape.at(k) * cor2.at(k) * (shape_prop.at(k) - shape_curr.at(k)), sd1 * sqrt1mx2(cor2.at(k)) * sd_sigma.at(k))[0];
      }
      lpost_prop.at(k) = lpost(u_curr.at(k), alpha_curr.at(k), theta_curr.at(k), shape_prop.at(k), sigma_prop.at(k), powl_curr.at(k));
      if (invt.at(k) * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
        shape_curr.at(k) = shape_prop.at(k);
        sigma_curr.at(k) = sigma_prop.at(k);
        lpost_curr.at(k) = lpost_prop.at(k);
      }
      // model selection
      if (pr_power != 0.0 && pr_power != 1.0 && t >= burn) {
        // only model select after burn-in
        ak = invt.at(k) * (a_pseudo - 1.0) + 1.0;
        bk = invt.at(k) * (b_pseudo - 1.0) + 1.0;
        if (powl_curr.at(k)) { // currently power law
          theta_pseudo = rbeta(1, ak, bk)[0]; // sim from pseudoprior
          logA_powl =
            lpost_curr.at(k) +
            ldbeta(theta_pseudo, ak, bk) +
            log(pr_power);
          logA_poly =
            lpost(u_curr.at(k), alpha_curr.at(k), theta_pseudo, shape_curr.at(k), sigma_curr.at(k), false) +
            log(1.0 - pr_power);
          if (lr1() > -log(1.0 + exp(logA_poly - logA_powl))) { // switch to polylog
            powl_curr.at(k) = false;
            theta_curr.at(k) = theta_pseudo;
            lpost_curr.at(k) = logA_poly - log(1.0 - pr_power);
          }
        }
        else { // currently polylog
          logA_powl =
            lpost(u_curr.at(k), alpha_curr.at(k), theta_curr.at(k), shape_curr.at(k), sigma_curr.at(k), true) +
            ldbeta(theta_curr.at(k), ak, bk) +
            log(pr_power);
          logA_poly =
            lpost_curr.at(k) +
            log(1.0 - pr_power);
          if (lr1() < -log(1.0 + exp(logA_poly - logA_powl))) { // switch to power law
            powl_curr.at(k) = true;
            theta_curr.at(k) = 1.0;
            lpost_curr.at(k) = logA_powl - ldbeta(theta_curr.at(k), ak, bk) - log(pr_power);
          }
        }
        if (k == 0) {
          powl_stat((double) powl_curr.at(k)); // only track cold chain
        }
      }
    } // loop over k completes
    // 07) Metropolis coupling
    if (K > 1) {
      k = Rcpp::RcppArmadillo::sample(seqKm1, 1, false)[0];
      l = k + 1;
      ldiff = (lpost_curr.at(k) - lpost_curr.at(l)) * (invt.at(l) - invt.at(k));
      if (lr1() < ldiff) {
        swap(u_curr.at(k), u_curr.at(l));
        swap(alpha_curr.at(k), alpha_curr.at(l));
        swap(theta_curr.at(k), theta_curr.at(l));
        swap(shape_curr.at(k), shape_curr.at(l));
        swap(sigma_curr.at(k), sigma_curr.at(l));
        swap(lpost_curr.at(k), lpost_curr.at(l));
        swap(powl_curr.at(k), powl_curr.at(l));
        swap_accept.at(k) += 1.0;
      }
      swap_count.at(k) += 1.0;
    }
    // 08) adaptive update sds after chain swap
    for (k = 0; k < K; k++) {
      if (t < burn) {
        if (!powl_curr.at(k)) {
          alpha_burn(t, k) = alpha_curr.at(k);
          theta_burn(t, k) = theta_curr.at(k);
          sd_alpha.at(k) = sd_curr(alpha_burn.col(k), t+1);
          sd_theta.at(k) = sd_curr(theta_burn.col(k), t+1);
          cor1.at(k) = cor_curr(alpha_burn.col(k), theta_burn.col(k), t+1);
        }
        shape_burn(t, k) = shape_curr.at(k);
        sigma_burn(t, k) = sigma_curr.at(k);
        sd_shape.at(k) = sd_curr(shape_burn.col(k), t+1);
        sd_sigma.at(k) = sd_curr(sigma_burn.col(k), t+1);
        cor2.at(k) = cor_curr(shape_burn.col(k), sigma_burn.col(k), t+1);
      }
    }
    // 09) print
    if ((t + 1) % freq == 0) {
      Rcout << "Iteration " << t + 1;
      Rcout << ": Log-posterior = " << std::setprecision(6) << lpost_curr.at(0) << endl;
      Rcout << "  u   = " << u_curr.at(0) << endl;
      if (t < burn) {
        Rcout << "alpha = " << alpha_curr.at(0) << " (" << sd_alpha.at(0) << ")" << endl;
        Rcout << "theta = " << theta_curr.at(0);
        if (!powl_curr.at(0)) {
          Rcout << " (" << sd_theta.at(0) << ")" << endl;
          Rcout << "cor(alpha, theta) = " << cor1.at(0) << endl;
        }
        else {
          Rcout << endl;
        }
        Rcout << "shape = " << shape_curr.at(0) << " (" << sd_shape.at(0) << ")" << endl;
        Rcout << "sigma = " << sigma_curr.at(0) << " (" << sd_sigma.at(0) << ")" << endl;
        Rcout << "cor(shape, sigma) = " << cor2.at(0) << endl;
      }
      else {
        Rcout << "alpha = " << alpha_curr.at(0) << endl;
        Rcout << "theta = " << theta_curr.at(0) << endl;
        Rcout << "shape = " << shape_curr.at(0) << endl;
        Rcout << "sigma = " << sigma_curr.at(0) << endl;
        if (pr_power != 0.0 && pr_power != 1.0) {
          // only print after burn-in & w/ model selection
          Rcout << "power law = " << (powl_curr.at(0) ? "true" : "false") << endl;
          Rcout << "average   = " << powl_stat.mean() << endl;
        }
      }
      if (K > 1) {
        Rcout << "swap rates: " << endl;
        for (k = 0; k < K-1; k++) {
          Rcout << "  b/w inv temp " << std::setprecision(3) << invt.at(k);
          Rcout << " & " << invt.at(k+1) << ": " << swap_accept.at(k) / swap_count.at(k) << endl;
        }
      }
      Rcout << endl;
    }
    // 10) save
    if (t >= burn) {
      s = t - burn + 1;
      if (s % thin == 0) {
        s = s / thin - 1;
        u = u_curr.at(0);
        alpha = alpha_curr.at(0);
        theta = theta_curr.at(0);
        shape = shape_curr.at(0);
        sigma = sigma_curr.at(0);
        cu = count[x > u];
        phiu = intdiv(sum(cu), n);
        u_vec[s] = u;
        alpha_vec[s] = alpha;
        theta_vec[s] = theta;
        shape_vec[s] = shape;
        sigma_vec[s] = sigma;
        powl_vec[s] = (int) powl_curr.at(0);
        phiu_vec[s] = phiu;
        lpost_vec[s] = lpost_curr.at(0);
        // prob mass & survival functions
        f0 = dmix2(x0, u, alpha, theta, shape, sigma, phiu);
        f0_mat.row(s) = as<rowvec>(f0);
        S0 = Smix2(x0, u, alpha, theta, shape, sigma, phiu);
        S0_mat.row(s) = as<rowvec>(S0);
      }
    }
  }
  // 11) output
  const vec p1 = { 0.025, 0.50, 0.975 };
  const mat
    f0_q = quantile(f0_mat, p1, 0), // 3 x n0
    S0_q = quantile(S0_mat, p1, 0); // 3 x n0
  const double
    po_power = mean((NumericVector) powl_vec),
    bf = odds(po_power) / odds(pr_power);
  gvs_quants["po_power"] = tv(po_power);
  gvs_quants["bf"] = tv(bf);
  DataFrame
    pars =
    DataFrame::create(Named("u") = u_vec,
                      Named("alpha") = alpha_vec,
                      Named("theta") = theta_vec,
                      Named("shape") = shape_vec,
                      Named("sigma") = sigma_vec,
                      Named("powl") = powl_vec,
                      Named("phiu") = phiu_vec,
                      Named("lpost") = lpost_vec),
    fitted =
    DataFrame::create(Named("x") = x0,
                      Named("f_025") = wrap(f0_q.row(0)),
                      Named("f_med") = wrap(f0_q.row(1)),
                      Named("f_975") = wrap(f0_q.row(2)),
                      Named("S_025") = wrap(S0_q.row(0)),
                      Named("S_med") = wrap(S0_q.row(1)),
                      Named("S_975") = wrap(S0_q.row(2)));
  List output =
    List::create(Named("pars") = pars,
                 Named("fitted") = fitted,
                 Named("data") = data,
                 Named("u_set") = u_set,
                 Named("init") = init,
                 Named("hyperpars") = hyperpars,
                 Named("gvs_quants") = gvs_quants,
                 Named("scalars") = scalars,
                 Named("swap_rates") = swap_accept / swap_count,
                 Named("invt") = invt);
  return output;
}





// 05) 3-component mixture distribution
//' Probability mass function (PMF) of 3-component discrete extreme value mixture distribution
//'
//' \code{dmix3} returns the PMF at x for the 3-component discrete extreme value mixture distribution. The component below v is the (truncated) Zipf-polylog(alpha1,theta1) distribution, between v & u the (truncated) Zipf-polylog(alpha2,theta2) distribution, and above u the generalised Pareto(shape, sigma) distribution.
//' @param x Vector of positive integers
//' @param v Positive integer representing the lower threshold
//' @param u Positive integer representing the upper threshold
//' @param alpha1 Real number, first parameter of the Zipf-polylog component below v
//' @param theta1 Real number in (0, 1], second parameter of the Zipf-polylog component below v
//' @param alpha2 Real number, first parameter of the Zipf-polylog component between v & u
//' @param theta2 Real number in (0, 1], second parameter of the Zipf-polylog component between v & u
//' @param shape Real number, shape parameter of the generalised Pareto component
//' @param sigma Real number, scale parameter of the generalised Pareto component
//' @param phi1 Real number in (0, 1), proportion of values below v
//' @param phi2 Real number in (0, 1), proportion of values between v & u
//' @param phiu Real number in (0, 1), exceedance rate of the threshold u
//' @return A numeric vector of the same length as x
//' @seealso \code{\link{Smix3}} for the corresponding survival function, \code{\link{dpol}} and \code{\link{dmix2}} for the PMFs of the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector dmix3(const IntegerVector x,
                          const int v,
                          const int u,
                          const double alpha1,
                          const double theta1,
                          const double alpha2,
                          const double theta2,
                          const double shape,
                          const double sigma,
                          const double phi1,
                          const double phi2,
                          const double phiu) {
  // density for 3-component mixture
  const int w = min(x) - 1;
  const double sigu = sigma + shape * u;
  const NumericVector
    x0(x),
    y1 = 1.0 + shape / sigu * (x0 - 1.0 - u),
    z1 = 1.0 + shape / (sigu + shape * (x0 - 1.0 - u)),
    fu = exp(log(1.0 - pow(z1, -1.0 / shape)) - 1.0 / shape * log(y1) + log(phiu));
  const IntegerVector
    seqv = seq_len(v), sequ = seq_len(u), // 1 to v, 1 to u
    seq1 = tail(seqv, v - w), // (w+1) to v
    seq2 = tail(sequ, u - v); // (v+1) to u
  const NumericVector
    x1(seq1), x2(seq2); // (w+1) to v, (v+1) to u
  NumericVector
    l1 = - alpha1 * log(x1) + x1 * log(theta1),
    l2 = - alpha2 * log(x2) + x2 * log(theta2);
  l1 = l1 - max(l1);
  l2 = l2 - max(l2);
  const NumericVector cf1 = cumsum(exp(l1)), cf2 = cumsum(exp(l2));
  NumericVector fl(x.size(), NA_REAL);
  for (int i = 0; i < x.size(); i++) {
    if (x[i] <= v) {
      fl[i] = exp(log(phi1) + l1[x[i]-(w+1)] - log(cf1[v-(w+1)]));
    }
    else if (x[i] <= u) {
      fl[i] = exp(log(phi2) + l2[x[i]-(v+1)] - log(cf2[u-(v+1)]));
    }
  }
  const NumericVector f = ifelse(x <= u, fl, fu);
  return f;
}

//' Survival function of 3-component discrete extreme value mixture distribution
//'
//' \code{Smix3} returns the survival function at x for the 3-component discrete extreme value mixture distribution. The component below v is the (truncated) Zipf-polylog(alpha1,theta1) distribution, between v & u the (truncated) Zipf-polylog(alpha2,theta2) distribution, and above u the generalised Pareto(shape, sigma) distribution.
//' @param x Vector of positive integers
//' @param v Positive integer representing the lower threshold
//' @param u Positive integer representing the upper threshold
//' @param alpha1 Real number, first parameter of the Zipf-polylog component below v
//' @param theta1 Real number in (0, 1], second parameter of the Zipf-polylog component below v
//' @param alpha2 Real number, first parameter of the Zipf-polylog component between v & u
//' @param theta2 Real number in (0, 1], second parameter of the Zipf-polylog component between v & u
//' @param shape Real number, shape parameter of the generalised Pareto component
//' @param sigma Real number, scale parameter of the generalised Pareto component
//' @param phi1 Real number in (0, 1), proportion of values below v
//' @param phi2 Real number in (0, 1), proportion of values between v & u
//' @param phiu Real number in (0, 1), exceedance rate of the threshold u
//' @return A numeric vector of the same length as x
//' @seealso \code{\link{dmix3}} for the corresponding probability mass function, \code{\link{Spol}} and \code{\link{Smix2}} for the survival functions of the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
const NumericVector Smix3(const IntegerVector x,
                          const int v,
                          const int u,
                          const double alpha1,
                          const double theta1,
                          const double alpha2,
                          const double theta2,
                          const double shape,
                          const double sigma,
                          const double phi1,
                          const double phi2,
                          const double phiu) {
  // survival for 3-component mixture
  const int w = min(x) - 1;
  const double sigu = sigma + shape * u;
  const NumericVector x0(x);
  NumericVector y = 1.0 + shape / sigu * (x0 - u);
  y = ifelse(y < 0.0, 0.0, y);
  const NumericVector Su = pow(y, -1.0 / shape) * phiu;
  const IntegerVector
    seqv = seq_len(v), sequ = seq_len(u), // 1 to v, 1 to u
    seq1 = tail(seqv, v - w), // (w+1) to v
    seq2 = tail(sequ, u - v); // (v+1) to u
  const NumericVector
    x1(seq1), x2(seq2); // (w+1) to v, (v+1) to u
  NumericVector
    l1 = - alpha1 * log(x1) + x1 * log(theta1),
    l2 = - alpha2 * log(x2) + x2 * log(theta2);
  l1 = l1 - max(l1);
  l2 = l2 - max(l2);
  const NumericVector cf1 = cumsum(exp(l1)), cf2 = cumsum(exp(l2));
  NumericVector Sl(x.size());
  for (int i = 0; i < x.size(); i++) {
    if (x[i] >= (w+1) && x[i] < v) {
      Sl[i] = phiu + phi2 + phi1 * (1.0 - cf1[x[i]-(w+1)] / cf1[v-(w+1)]);
    }
    else if (x[i] >= v && x[i] < u) {
      Sl[i] = phiu + phi2 * (1.0 - cf2[x[i]-(v+1)] / cf2[u-(v+1)]);
    }
    else if (x[i] == u) {
      Sl[i] = phiu;
    }
    else {
      Sl[i] = NA_REAL;
    }
  }
  const NumericVector S = ifelse(x <= u, Sl, Su);
  return S;
}

// [[Rcpp::export]]
const double lpost_mix3(const IntegerVector x,
                        const IntegerVector count,
                        const int v,
                        const int u,
                        const double alpha1,
                        const double theta1,
                        const double alpha2,
                        const double theta2,
                        const double shape,
                        const double sigma,
                        const double a_psi1,
                        const double a_psi2,
                        const double a_psiu,
                        const double b_psiu,
                        const double a_alpha1,
                        const double b_alpha1,
                        const double a_theta1,
                        const double b_theta1,
                        const double a_alpha2,
                        const double b_alpha2,
                        const double a_theta2,
                        const double b_theta2,
                        const double m_shape,
                        const double s_shape,
                        const double a_sigma,
                        const double b_sigma,
                        const bool powerlaw1,
                        const bool powerlaw2,
                        const bool positive1,
                        const bool positive2) {
  if (x.size() != count.size()) {
    stop("lpost_mix3: lengths of x & count have to be equal.");
  }
  const int w = min(x) - 1;
  const LogicalVector
    below = x <= v, between = (x > v) & (x <= u), above = x > u;
  const NumericVector
    x1(x[below]), c1(count[below]),
    x2(x[between]), c2(count[between]),
    xu(x[above]), cu(count[above]),
    par1 = NumericVector::create(alpha1, theta1),
    par2 = NumericVector::create(alpha2, theta2),
    paru = NumericVector::create(shape, sigma);
  const int n = sum(count), m = x.size(); // m = count.size()
  const double
    phi1 = intdiv(sum(c1), n),
    phi2 = intdiv(sum(c2), n),
    phiu = intdiv(sum(cu), n),
    psi1 = intdiv(x1.size(), m),
    //    psi2 = intdiv(x2.size(), m),
    psiu = intdiv(xu.size(), m);
  double l;
  if (v <= min(x) || v >= u || u >= max(x)) {
    l = -INFINITY;
  }
  else {
    l =
      llik_bulk(par1, x, count, w, v, phi1, powerlaw1, positive1) +
      llik_bulk(par2, x, count, v, u, phi2, powerlaw2, positive2) +
      llik_igpd(paru, x, count, u, phiu) +
      ldbeta(psi1 / (1.0 - psiu), a_psi1, a_psi2) +
      ldunif(psiu, a_psiu, b_psiu) +
      (powerlaw1 ? 0.0 : ldbeta(theta1, a_theta1, b_theta1)) +
      (powerlaw2 ? 0.0 : ldbeta(theta2, a_theta2, b_theta2)) +
      ldnorm(shape, m_shape, s_shape) +
      ldgamma(sigma, a_sigma, b_sigma) +
      ldnorm(alpha1, a_alpha1, b_alpha1) +
      ldnorm(alpha2, a_alpha2, b_alpha2);
  }
  return lnan(l);
}

//' Markov chain Monte Carlo for 3-component discrete extreme value mixture distribution
//'
//' \code{mcmc_mix3} returns the posterior samples of the parameters, for fitting the 3-component discrete extreme value mixture distribution. The samples are obtained using Markov chain Monte Carlo (MCMC).
//'
//' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. The thresholds v and u are treated as parameters and therefore sampled. The hyperparameters are used in the following priors: psi1 / (1.0 - psiu) ~ Beta(a_psi1, a_psi2); u is such that the implied unique exceedance probability psiu ~ Uniform(a_psi, b_psi); alpha1 ~ Normal(mean = a_alpha1, sd = b_alpha1); theta1 ~ Beta(a_theta1, b_theta1); alpha2 ~ Normal(mean = a_alpha2, sd = b_alpha2); theta2 ~ Beta(a_theta2, b_theta2); shape ~ Normal(mean = m_shape, sd = s_shape); sigma ~ Gamma(a_sigma, scale = b_sigma). If pr_power2 = 1.0, the discrete power law (between v and u) is assumed, and the samples of theta2 will be all 1.0. If pr_power2 is in (0.0, 1.0), model selection between the polylog distribution and the discrete power law will be performed within the MCMC.
//' @param x Vector of the unique values (positive integers) of the data
//' @param count Vector of the same length as x that contains the counts of each unique value in the full data, which is essentially rep(x, count)
//' @param v_set Positive integer vector of the values v will be sampled from
//' @param u_set Positive integer vector of the values u will be sampled from
//' @param v Positive integer, initial value of the lower threshold
//' @param u Positive integer, initial value of the upper threshold
//' @param alpha1 Real number greater than 1, initial value of the parameter
//' @param theta1 Real number in (0, 1], initial value of the parameter
//' @param alpha2 Real number greater than 1, initial value of the parameter
//' @param theta2 Real number in (0, 1], initial value of the parameter
//' @param shape Real number, initial value of the parameter
//' @param sigma Positive real number, initial value of the parameter
//' @param a_psi1,a_psi2,a_psiu,b_psiu,a_alpha1,b_alpha1,a_theta1,b_theta1,a_alpha2,b_alpha2,a_theta2,b_theta2,m_shape,s_shape,a_sigma,b_sigma Scalars, real numbers representing the hyperparameters of the prior distributions for the respective parameters. See details for the specification of the priors.
//' @param powerlaw1 Boolean, is the discrete power law assumed for below v?
//' @param positive1 Boolean, is alpha1 positive (TRUE) or unbounded (FALSE)?
//' @param positive2 Boolean, is alpha2 positive (TRUE) or unbounded (FALSE)?
//' @param a_pseudo Positive real number, first parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
//' @param b_pseudo Positive real number, second parameter of the pseudoprior beta distribution for theta2 in model selection; ignored if pr_power2 = 1.0
//' @param pr_power2 Real number in [0, 1], prior probability of the discrete power law (between v and u)
//' @param iter Positive integer representing the length of the MCMC output
//' @param thin Positive integer representing the thinning in the MCMC
//' @param burn Non-negative integer representing the burn-in of the MCMC
//' @param freq Positive integer representing the frequency of the sampled values being printed
//' @param invt Vector of the inverse temperatures for Metropolis-coupled MCMC; default c(1.0) i.e. no Metropolis-coupling
//' @return A list: $pars is a data frame of iter rows of the MCMC samples, $fitted is a data frame of length(x) rows with the fitted values, amongst other quantities related to the MCMC
//' @seealso \code{\link{mcmc_pol}} and \code{\link{mcmc_mix2}} for MCMC for the Zipf-polylog and 2-component discrete extreme value mixture distributions, respectively.
//' @export
// [[Rcpp::export]]
List mcmc_mix3(const IntegerVector x,
               const IntegerVector count,
               const IntegerVector v_set,
               const IntegerVector u_set,
               int v,
               int u,
               double alpha1,
               double theta1,
               double alpha2,
               double theta2,
               double shape,
               double sigma,
               const double a_psi1,
               const double a_psi2,
               const double a_psiu,
               const double b_psiu,
               const double a_alpha1,
               const double b_alpha1,
               const double a_theta1,
               const double b_theta1,
               const double a_alpha2,
               const double b_alpha2,
               const double a_theta2,
               const double b_theta2,
               const double m_shape,
               const double s_shape,
               const double a_sigma,
               const double b_sigma,
               const bool powerlaw1, // power law for left component
               const bool positive1, // alpha1 +ve / unbounded
               const bool positive2, // alpha2 +ve / unbounded
               const double a_pseudo,
               const double b_pseudo,
               const double pr_power2, // prior prob of power law for middle component
               const int iter,
               const int thin,
               const int burn,
               const int freq,
               const NumericVector invt) {
  // 01) save
  DataFrame
    data =
    DataFrame::create(Named("x") = x,
                      Named("count") = count),
    init =
    DataFrame::create(Named("v") = ti(v),
                      Named("u") = ti(u),
                      Named("alpha1") = tv(alpha1),
                      Named("theta1") = tv(theta1),
                      Named("alpha2") = tv(alpha2),
                      Named("theta2") = tv(theta2),
                      Named("shape") = tv(shape),
                      Named("sigma") = tv(sigma)),
    hyperpars =
    DataFrame::create(Named("a_psi1") = tv(a_psi1),
                      Named("a_psi2") = tv(a_psi2),
                      Named("a_psiu") = tv(a_psiu),
                      Named("b_psiu") = tv(b_psiu),
                      Named("a_alpha1") = tv(a_alpha1),
                      Named("b_alpha1") = tv(b_alpha1),
                      Named("a_theta1") = tv(a_theta1),
                      Named("b_theta1") = tv(b_theta1),
                      Named("a_alpha2") = tv(a_alpha2),
                      Named("b_alpha2") = tv(b_alpha2),
                      Named("a_theta2") = tv(a_theta2),
                      Named("b_theta2") = tv(b_theta2),
                      Named("m_shape") = tv(m_shape),
                      Named("s_shape") = tv(s_shape),
                      Named("a_sigma") = tv(a_sigma),
                      Named("b_sigma") = tv(b_sigma),
                      Named("powerlaw1") = tl(powerlaw1),
                      Named("positive1") = tl(positive1),
                      Named("positive2") = tl(positive2)),
    gvs_quants =
    DataFrame::create(Named("a_pseudo") = tv(a_pseudo),
                      Named("b_pseudo") = tv(b_pseudo),
                      Named("pr_power2") = tv(pr_power2)),
    scalars = df_scalars(iter, thin, burn, freq, invt.size() != 1);
  // 02) checks
  // x are the unique values (> 1) w/ freq count
  if (is_true(any(x <= 0))) {
    stop("mcmc_mix3: all of x has to be +ve integers.");
  }
  if (is_true(all(v_set != v))) {
    stop("mcmc_mix3: v must be in v_set.");
  }
  if (is_true(all(u_set != u))) {
    stop("mcmc_mix3: u must be in u_set.");
  }
  if (invt.at(0) != 1.0) {
    stop("mcmc_mix3: 1st element of inverse temperatures (invt) has to be 1.0.");
  }
  if (is_true(any(invt > 1.0)) || is_true(any(invt <= 0.0))) {
    stop("mcmc_mix3: all elements of invt must be in (0.0, 1.0].");
  }
  NumericVector invt0 = clone(invt);
  std::reverse(invt0.begin(), invt0.end());
  if (!std::is_sorted(invt0.begin(), invt0.end())) {
    stop("mcmc_mix3: invt has to be in decreasing order.");
  }
  // 03) dimensions const to invt.size()
  int k;
  const int K = invt.size(), m = u_set.size();
  const IntegerVector seqm = seq_len(m) - 1;
  NumericVector w(m), w1(m);
  const double sd0 = 0.001 / sqrt(2.0), sd1 = 2.38 / sqrt(2.0); // see Roberts & Rosenthal (2008)
  double ldiff;
  vec swap_accept(K-1, fill::zeros), swap_count(K-1, fill::zeros);
  bool powl, accept_reject;
  if (pr_power2 == 1.0) { // power law
    powl = true; // & stay true
    theta2 = 1.0; // & stay 1.0
  }
  else if (pr_power2 == 0.0) { // polylog
    powl = false; // & stay false
  }
  else { // model selection
    powl = false;
  }
  // for model selection
  double theta2_pseudo, ak, bk, logA_powl, logA_poly;
  running_stat<double> powl_stat;
  // 04) dimensions increase w/ K
  IntegerVector v_curr(K, v), u_curr(K, u);
  NumericVector
    alpha1_curr(K, alpha1), alpha1_prop(K),
    theta1_curr(K, theta1), theta1_prop(K),
    alpha2_curr(K, alpha2), alpha2_prop(K),
    theta2_curr(K, theta2), theta2_prop(K),
    shape_curr(K, shape), shape_prop(K),
    sigma_curr(K, sigma), sigma_prop(K),
    lpost_curr(K), lpost_prop(K);
  LogicalVector powl_curr(K, powl);
  auto lpost =
    [x, count,
     a_psi1, a_psi2, a_psiu, b_psiu,
     a_alpha1, b_alpha1, a_theta1, b_theta1,
     a_alpha2, b_alpha2, a_theta2, b_theta2,
     m_shape, s_shape, a_sigma, b_sigma,
     powerlaw1, positive1, positive2]
    (const int v, const int u,
     const double alpha1, const double theta1,
     const double alpha2, const double theta2,
     const double shape, const double sigma,
     const bool powerlaw2) {
      return lpost_mix3(x, count, v, u,
                        alpha1, theta1, alpha2, theta2,
                        shape, sigma,
                        a_psi1, a_psi2, a_psiu, b_psiu,
                        a_alpha1, b_alpha1, a_theta1, b_theta1,
                        a_alpha2, b_alpha2, a_theta2, b_theta2,
                        m_shape, s_shape, a_sigma, b_sigma,
                        powerlaw1, powerlaw2,
                        positive1, positive2);
    };
  for (int k = 0; k < K; k++) {
    lpost_curr.at(k) = lpost(v_curr.at(k), u_curr.at(k),
                             alpha1_curr.at(k), theta1_curr.at(k),
                             alpha2_curr.at(k), theta2_curr.at(k),
                             shape_curr.at(k), sigma_curr.at(k),
                             powl_curr.at(k));
  }
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  mat
    alpha1_burn(burn, K), theta1_burn(burn, K),
    alpha2_burn(burn, K), theta2_burn(burn, K),
    shape_burn(burn, K), sigma_burn(burn, K);
  NumericVector
    sd_alpha1(K, 0.1), sd_theta1(K, 0.1),
    sd_alpha2(K, 0.1), sd_theta2(K, 0.1),
    sd_shape(K, 0.1), sd_sigma(K, 0.1),
    cor1(K, 0.1), cor2(K, 0.1), cor3(K, 0.1);
  // these sds & cors not as influential as sd0
  const IntegerVector seqKm1 = seq_len(K - 1) - 1;
  // 05) for saving, dimensions const to K
  IntegerVector
    v_vec(iter), u_vec(iter), powl_vec(iter);
  NumericVector
    alpha1_vec(iter), theta1_vec(iter),
    alpha2_vec(iter), theta2_vec(iter),
    shape_vec(iter), sigma_vec(iter),
    phi1_vec(iter), phi2_vec(iter),
    phiu_vec(iter), lpost_vec(iter);
  const int n = sum(count);
  double phi1, phi2, phiu;
  const int xmin = min(x);
  IntegerVector x0 = for_Smix(max(x) - 1, xmin), c1, c2, cu;
  const int n0 = x0.size();
  NumericVector f0(n0), S0(n0); // fitted
  mat f0_mat(iter, n0), S0_mat(iter, n0);
  // 06) run
  int s, t, l;
  for (t = 0; t < iter * thin + burn; t++) {
    for (k = 0; k < K; k++) {
      // v, u
      std::fill(w.begin(), w.end(), NA_REAL);
      for (s = 0; s < u_set.size(); s++) {
        w[s] = lpost(v_set[s], u_set[s],
                     alpha1_curr.at(k), theta1_curr.at(k),
                     alpha2_curr.at(k), theta2_curr.at(k),
                     shape_curr.at(k), sigma_curr.at(k),
                     powl_curr.at(k));
      }
      w1 = exp(invt.at(k) * (w - max(w)));
      s = sample_w(seqm, w1); // reused as loop done
      v_curr.at(k) = v_set[s];
      u_curr.at(k) = u_set[s];
      lpost_curr.at(k) = w[s];
      // alpha1 & theta1
      if (t < burn &&
          (isnan(cor1.at(k)) || ispm1(cor1.at(k)) ||
           isnan(sd_alpha1.at(k)) || isnan(sd_theta1.at(k)) || lr1() < log(0.05))) {
        alpha1_prop.at(k) = rnorm(1, alpha1_curr.at(k), sd0)[0];
        theta1_prop.at(k) = rnorm(1, theta1_curr.at(k), sd0)[0];
      }
      else {
        alpha1_prop.at(k) = rnorm(1, alpha1_curr.at(k), sd1 * sd_alpha1.at(k))[0];
        theta1_prop.at(k) = rnorm(1, theta1_curr.at(k) + sd_theta1.at(k) / sd_alpha1.at(k) * cor1.at(k) * (alpha1_prop.at(k) - alpha1_curr.at(k)), sd1 * sqrt1mx2(cor1.at(k)) * sd_theta1.at(k))[0];
      }
      lpost_prop.at(k) = lpost(v_curr.at(k), u_curr.at(k),
                               alpha1_prop.at(k), theta1_prop.at(k),
                               alpha2_curr.at(k), theta2_curr.at(k),
                               shape_curr.at(k), sigma_curr.at(k),
                               powl_curr.at(k));
      if (invt.at(k) * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
        alpha1_curr.at(k) = alpha1_prop.at(k);
        theta1_curr.at(k) = theta1_prop.at(k);
        lpost_curr.at(k) = lpost_prop.at(k);
      }
      // alpha2 & theta2
      if (powl_curr.at(k)) {
        alpha2_prop.at(k) = rnorm(1, alpha2_curr.at(k), sd_alpha2.at(k))[0];
        lpost_prop.at(k) = lpost(v_curr.at(k), u_curr.at(k),
                                 alpha1_curr.at(k), theta1_curr.at(k),
                                 alpha2_prop.at(k), theta2_curr.at(k),
                                 shape_curr.at(k), sigma_curr.at(k),
                                 powl_curr.at(k));
        update(alpha2_curr.at(k), alpha2_prop.at(k), lpost_curr.at(k), lpost_prop.at(k), sd_alpha2.at(k), t, burn, accept_reject, invt.at(k));
      }
      else {
        if (t < burn &&
            (isnan(cor2.at(k)) || ispm1(cor2.at(k)) ||
             isnan(sd_alpha2.at(k)) || isnan(sd_theta2.at(k)) || lr1() < log(0.05))) {
          alpha2_prop.at(k) = rnorm(1, alpha2_curr.at(k), sd0)[0];
          theta2_prop.at(k) = rnorm(1, theta2_curr.at(k), sd0)[0];
        }
        else {
          alpha2_prop.at(k) = rnorm(1, alpha2_curr.at(k), sd1 * sd_alpha2.at(k))[0];
          theta2_prop.at(k) = rnorm(1, theta2_curr.at(k) + sd_theta2.at(k) / sd_alpha2.at(k) * cor2.at(k) * (alpha2_prop.at(k) - alpha2_curr.at(k)), sd1 * sqrt1mx2(cor2.at(k)) * sd_theta2.at(k))[0];
        }
        lpost_prop.at(k) = lpost(v_curr.at(k), u_curr.at(k),
                                 alpha1_curr.at(k), theta1_curr.at(k),
                                 alpha2_prop.at(k), theta2_prop.at(k),
                                 shape_curr.at(k), sigma_curr.at(k),
                                 powl_curr.at(k));
        if (invt.at(k) * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
          alpha2_curr.at(k) = alpha2_prop.at(k);
          theta2_curr.at(k) = theta2_prop.at(k);
          lpost_curr.at(k) = lpost_prop.at(k);
        }
      }
      // shape & sigma
      if (t < burn &&
          (isnan(cor3.at(k)) || ispm1(cor3.at(k)) ||
           isnan(sd_shape.at(k)) || isnan(sd_sigma.at(k)) || lr1() < log(0.05))) {
        shape_prop.at(k) = rnorm(1, shape_curr.at(k), sd0)[0];
        sigma_prop.at(k) = rnorm(1, sigma_curr.at(k), sd0)[0];
      }
      else {
        shape_prop.at(k) = rnorm(1, shape_curr.at(k), sd1 * sd_shape.at(k))[0];
        sigma_prop.at(k) = rnorm(1, sigma_curr.at(k) + sd_sigma.at(k) / sd_shape.at(k) * cor3.at(k) * (shape_prop.at(k) - shape_curr.at(k)), sd1 * sqrt1mx2(cor3.at(k)) * sd_sigma.at(k))[0];
      }
      lpost_prop.at(k) = lpost(v_curr.at(k), u_curr.at(k),
                               alpha1_curr.at(k), theta1_curr.at(k),
                               alpha2_curr.at(k), theta2_curr.at(k),
                               shape_prop.at(k), sigma_prop.at(k),
                               powl_curr.at(k));
      if (invt.at(k) * (lpost_prop.at(k) - lpost_curr.at(k)) > lr1()) {
        shape_curr.at(k) = shape_prop.at(k);
        sigma_curr.at(k) = sigma_prop.at(k);
        lpost_curr.at(k) = lpost_prop.at(k);
      }
      // model selection
      if (pr_power2 != 0.0 && pr_power2 != 1.0 && t >= burn) {
        // only model select after burn-in
        ak = invt.at(k) * (a_pseudo - 1.0) + 1.0;
        bk = invt.at(k) * (b_pseudo - 1.0) + 1.0;
        if (powl_curr.at(k)) { // currently power law
          theta2_pseudo = rbeta(1, ak, bk)[0]; // sim from pseudoprior
          logA_powl =
            lpost_curr.at(k) +
            ldbeta(theta2_pseudo, ak, bk) +
            log(pr_power2);
          logA_poly =
            lpost(v_curr.at(k), u_curr.at(k),
                  alpha1_curr.at(k), theta1_curr.at(k),
                  alpha2_curr.at(k), theta2_pseudo,
                  shape_curr.at(k), sigma_curr.at(k),
                  false) +
            log(1.0 - pr_power2);
          if (lr1() > -log(1.0 + exp(logA_poly - logA_powl))) { // switch to polylog
            powl_curr.at(k) = false;
            theta2_curr.at(k) = theta2_pseudo;
            lpost_curr.at(k) = logA_poly - log(1.0 - pr_power2);
          }
        }
        else { // currently polylog
          logA_powl =
            lpost(v_curr.at(k), u_curr.at(k),
                  alpha1_curr.at(k), theta1_curr.at(k),
                  alpha2_curr.at(k), theta2_curr.at(k),
                  shape_curr.at(k), sigma_curr.at(k),
                  true) +
            ldbeta(theta2_curr.at(k), ak, bk) +
            log(pr_power2);
          logA_poly =
            lpost_curr.at(k) +
            log(1.0 - pr_power2);
          if (lr1() < -log(1.0 + exp(logA_poly - logA_powl))) { // switch to power law
            powl_curr.at(k) = true;
            theta2_curr.at(k) = 1.0;
            lpost_curr.at(k) = logA_powl - ldbeta(theta2_curr.at(k), ak, bk) - log(pr_power2);
          }
        }
        if (k == 0) {
          powl_stat((double) powl_curr.at(k)); // only track cold chain
        }
      }
    } // loop over k completes
    // 07) Metropolis coupling
    if (K > 1) {
      k = Rcpp::RcppArmadillo::sample(seqKm1, 1, false)[0];
      l = k + 1;
      ldiff = (lpost_curr.at(k) - lpost_curr.at(l)) * (invt.at(l) - invt.at(k));
      if (lr1() < ldiff) {
        swap(v_curr.at(k), v_curr.at(l));
        swap(u_curr.at(k), u_curr.at(l));
        swap(alpha1_curr.at(k), alpha1_curr.at(l));
        swap(theta1_curr.at(k), theta1_curr.at(l));
        swap(alpha2_curr.at(k), alpha2_curr.at(l));
        swap(theta2_curr.at(k), theta2_curr.at(l));
        swap(shape_curr.at(k), shape_curr.at(l));
        swap(sigma_curr.at(k), sigma_curr.at(l));
        swap(lpost_curr.at(k), lpost_curr.at(l));
        swap(powl_curr.at(k), powl_curr.at(l));
        swap_accept.at(k) += 1.0;
      }
      swap_count.at(k) += 1.0;
    }
    // 08) adaptive update sds after chain swap
    for (k = 0; k < K; k++) {
      if (t < burn) {
        alpha1_burn(t, k) = alpha1_curr.at(k);
        theta1_burn(t, k) = theta1_curr.at(k);
        sd_alpha1.at(k) = sd_curr(alpha1_burn.col(k), t+1);
        sd_theta1.at(k) = sd_curr(theta1_burn.col(k), t+1);
        cor1.at(k) = cor_curr(alpha1_burn.col(k), theta1_burn.col(k), t+1);
        if (!powl_curr.at(k)) {
          alpha2_burn(t, k) = alpha2_curr.at(k);
          theta2_burn(t, k) = theta2_curr.at(k);
          sd_alpha2.at(k) = sd_curr(alpha2_burn.col(k), t+1);
          sd_theta2.at(k) = sd_curr(theta2_burn.col(k), t+1);
          cor2.at(k) = cor_curr(alpha2_burn.col(k), theta2_burn.col(k), t+1);
        }
        shape_burn(t, k) = shape_curr.at(k);
        sigma_burn(t, k) = sigma_curr.at(k);
        sd_shape.at(k) = sd_curr(shape_burn.col(k), t+1);
        sd_sigma.at(k) = sd_curr(sigma_burn.col(k), t+1);
        cor3.at(k) = cor_curr(shape_burn.col(k), sigma_burn.col(k), t+1);
      }
    }
    // 09) print
    if ((t + 1) % freq == 0) {
      Rcout << "Iteration " << t + 1;
      Rcout << ": Log-posterior = " << std::setprecision(6) << lpost_curr.at(0) << endl;
      Rcout << "   v   = " << v_curr.at(0) << endl;
      Rcout << "   u   = " << u_curr.at(0) << endl;
      if (t < burn) {
        Rcout << "alpha1 = " << alpha1_curr.at(0) << " (" << sd_alpha1.at(0) << ")" << endl;
        Rcout << "theta1 = " << theta1_curr.at(0);
        Rcout << " (" << sd_theta1.at(0) << ")" << endl;
        Rcout << "cor(alpha1, theta1) = " << cor1.at(0) << endl;
        Rcout << "alpha2 = " << alpha2_curr.at(0) << " (" << sd_alpha2.at(0) << ")" << endl;
        Rcout << "theta2 = " << theta2_curr.at(0);
        if (!powl_curr.at(0)) {
          Rcout << " (" << sd_theta2.at(0) << ")" << endl;
          Rcout << "cor(alpha2, theta2) = " << cor2.at(0) << endl;
        }
        else {
          Rcout << endl;
        }
        Rcout << "shape  = " << shape_curr.at(0) << " (" << sd_shape.at(0) << ")" << endl;
        Rcout << "sigma  = " << sigma_curr.at(0) << " (" << sd_sigma.at(0) << ")" << endl;
        Rcout << "cor(shape, sigma) = " << cor3.at(0) << endl;
      }
      else {
        Rcout << "alpha1 = " << alpha1_curr.at(0) << endl;
        Rcout << "theta1 = " << theta1_curr.at(0) << endl;
        Rcout << "alpha2 = " << alpha2_curr.at(0) << endl;
        Rcout << "theta2 = " << theta2_curr.at(0) << endl;
        Rcout << "shape  = " << shape_curr.at(0) << endl;
        Rcout << "sigma  = " << sigma_curr.at(0) << endl;
        if (pr_power2 != 0.0 && pr_power2 != 1.0) {
          // only print after burn-in & w / model selection
          Rcout << "power law = " << (powl_curr.at(0) ? "true" : "false") << endl;
          Rcout << "average   = " << powl_stat.mean() << endl;
        }
      }
      if (K > 1) {
        Rcout << "swap rates: " << endl;
        for (k = 0; k < K-1; k++) {
          Rcout << "  b/w inv temp " << std::setprecision(3) << invt.at(k);
          Rcout << " & " << invt.at(k+1) << ": " << swap_accept.at(k) / swap_count.at(k) << endl;
        }
      }
      Rcout << endl;
    }
    // 10) save
    if (t >= burn) {
      s = t - burn + 1;
      if (s % thin == 0) {
        s = s / thin - 1;
        // cold chain
        v = v_curr.at(0);
        u = u_curr.at(0);
        alpha1 = alpha1_curr.at(0);
        theta1 = theta1_curr.at(0);
        alpha2 = alpha2_curr.at(0);
        theta2 = theta2_curr.at(0);
        shape = shape_curr.at(0);
        sigma = sigma_curr.at(0);
        // compute phis
        c1 = count[x <= v];
        c2 = count[(x > v) & (x <= u)];
        cu = count[x > u];
        phi1 = intdiv(sum(c1), n);
        phi2 = intdiv(sum(c2), n);
        phiu = intdiv(sum(cu), n);
        // save
        v_vec[s] = v;
        u_vec[s] = u;
        alpha1_vec[s] = alpha1;
        theta1_vec[s] = theta1;
        alpha2_vec[s] = alpha2;
        theta2_vec[s] = theta2;
        shape_vec[s] = shape;
        sigma_vec[s] = sigma;
        powl_vec[s] = (int) powl_curr.at(0);
        phi1_vec[s] = phi1;
        phi2_vec[s] = phi2;
        phiu_vec[s] = phiu;
        lpost_vec[s] = lpost_curr.at(0);
        // prob mass & survival functions
        f0 = dmix3(x0, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu);
        f0_mat.row(s) = as<rowvec>(f0);
        S0 = Smix3(x0, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu);
        S0_mat.row(s) = as<rowvec>(S0);
      }
    }
  }
  // 11) output
  const vec p1 = { 0.025, 0.50, 0.975 };
  const mat
    f0_q = quantile(f0_mat, p1, 0), // 3 x n0
    S0_q = quantile(S0_mat, p1, 0); // 3 x n0
  const double
    po_power2 = mean((NumericVector) powl_vec),
    bf = odds(po_power2) / odds(pr_power2);
  gvs_quants["po_power2"] = tv(po_power2);
  gvs_quants["bf"] = tv(bf);
  DataFrame
    pars =
    DataFrame::create(Named("v") = v_vec,
                      Named("u") = u_vec,
                      Named("alpha1") = alpha1_vec,
                      Named("theta1") = theta1_vec,
                      Named("alpha2") = alpha2_vec,
                      Named("theta2") = theta2_vec,
                      Named("shape") = shape_vec,
                      Named("sigma") = sigma_vec,
                      Named("powl") = powl_vec,
                      Named("phi1") = phi1_vec,
                      Named("phi2") = phi2_vec,
                      Named("phiu") = phiu_vec,
                      Named("lpost") = lpost_vec),
    fitted =
    DataFrame::create(Named("x") = x0,
                      Named("f_025") = wrap(f0_q.row(0)),
                      Named("f_med") = wrap(f0_q.row(1)),
                      Named("f_975") = wrap(f0_q.row(2)),
                      Named("S_025") = wrap(S0_q.row(0)),
                      Named("S_med") = wrap(S0_q.row(1)),
                      Named("S_975") = wrap(S0_q.row(2)));
  List output =
    List::create(Named("pars") = pars,
                 Named("fitted") = fitted,
                 Named("data") = data,
                 Named("v_set") = v_set,
                 Named("u_set") = u_set,
                 Named("init") = init,
                 Named("hyperpars") = hyperpars,
                 Named("gvs_quants") = gvs_quants,
                 Named("scalars") = scalars,
                 Named("swap_rates") = swap_accept / swap_count,
                 Named("invt") = invt);
  return output;
}
