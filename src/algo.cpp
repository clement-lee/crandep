// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_zeta.h>
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

const NumericVector tv(const double x) {
  return NumericVector::create(x);
}

const IntegerVector ti(const int x) {
  return IntegerVector::create(x);
}

template <class T>
void update(T & par_curr, T par_prop, double & lpost_curr, const double lpost_prop, double & s, const int i, const int burnin, const double factor = 10.0) {
  const double lalpha = lpost_prop - lpost_curr;
  const bool accept_reject = lalpha > lr1();
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

const double lq(const int from, const int to, const double s) {
  const double
    lower = log((to + 0.0) / (from + 0.0)),
    upper = log((to + 1.0) / (from + 0.0));
  return log(pnorm(tv(upper), 0.0, s)[0] - pnorm(tv(lower), 0.0, s)[0]);
}





// 01) Component dists: discrete power laws & GPD
void hzeta_check(const double alpha, const int u) {
  if (alpha <= 1.0) {
    stop("hzeta_check: 1st argument of gsl_sf_hzeta() (exponent of power-law) has to be strictly greater than 1.0.");
  }
  if (u <= 0) {
    stop("hzeta_check: 2nd argument of gsl_sf_hzeta() has to be a positive integer.");
  }
}

//' Probability mass function (PMF) of discrete power law
//'
//' \code{dupp} returns the PMF at x for the discrete power law with exponent alpha, for values greater than or equal to u. 
//'
//' The PMF is proportional to x^(-alpha). To be a proper PMF, it is normalised by 1/hzeta(alpha, u), where hzeta is the Hurwitz zeta function i.e. hzeta(y, z) = z^(-y) + (z+1)^(-y) + (z+2)^(-y) + ... Any values below u will have PMF equal to 0.0.
//' @param x Vector of positive integers
//' @param u Scalar, non-negative integer threshold
//' @param alpha Scalar, a numeric value greater than 1.0
//' @param give_log Boolean, whether the PMF should be returned on the log scale. If 'FALSE', the PMF is returned on the original scale.
//' @return A numeric vector of the same length as x
//' @examples
//' dupp(c(10,20,30,40,50), 12, 2.0, FALSE)
//' dupp(c(10,20,30,40,50), 12, 2.0, TRUE)
//' @seealso \code{\link{Supp}} for the corresponding survival function, \code{\link{dmix}} for the PMF of the discrete extreme value mixture distribution.
//' @export
// [[Rcpp::export]]
const NumericVector dupp(const NumericVector x, const int u, const double alpha, const bool give_log = false) {
  // density f() of discrete power law >= u
  const int n = x.size();
  NumericVector v(n);
  hzeta_check(alpha, u);
  const double log_denom = log(gsl_sf_hzeta(alpha, u));
  for (int i = 0; i < n; i++) {
    v[i] = -alpha * log(x[i]) - log_denom;
  }
  v = ifelse(x < u, -INFINITY, v);
  NumericVector output;
  if (give_log) {
	output = v;
  }
  else {
	output = exp(v);
  }
  return output;
}

//' Survival function of discrete power law
//'
//' \code{Supp} returns the survival function at x for the discrete power law with exponent alpha, for values greater than or equal to u.
//'
//' The survival function used is S(x) = Pr(X >= x), where X is a random variable following the discrete power law. The inclusion of x in the sum means S(x) may not necessarily equal to Pr(X > x) as the distribution is discrete. In the case of discrete power law, it can be shown that S(x) = hzeta(alpha, x)/hzeta(alpha, u), where hzeta is the Hurwitz zeta function i.e. hzeta(y, z) = z^(-y) + (z+1)^(-y) + (z+2)^(-y) + ...
//' @param x Vector of positive integers
//' @param u Scalar, non-negative integer threshold
//' @param alpha Scalar, a numeric value greater than 1.0
//' @return A numeric vector of the same length as x
//' @examples
//' Supp(c(10,20,30,40,50), 12, 2.0)
//' @seealso \code{\link{dupp}} for the corresponding probability mass function, \code{\link{Smix}} for the survival function of the discrete extreme value mixture distribution.
//' @export
// [[Rcpp::export]]
const NumericVector Supp(const NumericVector x, const int u, const double alpha) {
  // survival f() of discrete power law >= u
  const int n = x.size();
  NumericVector v(n);
  hzeta_check(alpha, u);
  const double log_denom = log(gsl_sf_hzeta(alpha, u));
  for (int i = 0; i < n; i++) {
    v[i] = exp(log(gsl_sf_hzeta(alpha, x[i])) - log_denom);
  }
  v = ifelse(x < u, 1.0, v);
  return v;
}

const double llik_upp(const NumericVector par, const NumericVector x, const int u) {
  // discrete power law >= u
  if (par.size() != 1) {
    stop("llik_upp: length of par has to be 1.");
  }
  const double xi1 = par[0], alpha = 1.0 / xi1 + 1.0;
  double l;
  const NumericVector xu = x[x >= u];
  const double nu(xu.size());
  if (xi1 <= 0.0 || u <= 0) {
    l = -INFINITY;
  }
  else {
    gsl_set_error_handler_off();
    gsl_sf_result result;
    int status = gsl_sf_hzeta_e(alpha, u, &result);
    if (status != GSL_SUCCESS) {
      l = -INFINITY;
    }
    else {
      l = -alpha * sum(log(xu)) - nu * log(result.val);
    }
  }
  return lnan(l);
}

const double lpost_upp(const NumericVector x, const int u, const double xi1, const double a_xi1, const double b_xi1) {
  const NumericVector par = tv(xi1);
  const double l = llik_upp(par, x, u) + dunif(par, a_xi1, b_xi1, true)[0];
  return lnan(l);
}

//' Markov chain Monte Carlo for discrete power law
//'
//' \code{mcmc_upp} returns the samples from the posterior of xi1, for fitting the discrete power law to the data x. The samples are obtained using Markov chain Monte Carlo (MCMC).
//'
//' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. Unlike \code{mcmc_mix}, the threshold u is treated as fixed in \code{mcmc_upp}.
//' @param x Vector of positive integers, representing the data
//' @param u Scalar, non-negative integer threshold
//' @param xi1 Scalar, initial value of the shape parameter
//' @param a_xi1 Scalar, lower bound of the uniform distribution as the prior of xi1
//' @param b_xi1 Scalar, upper bound of the uniform distribution as the prior of xi1
//' @param N Scalar, positive integer representing the length of the output chain i.e. the number of rows in the returned data frame
//' @param thin Scalar, positive integer representing the thinning in the MCMC
//' @param burnin Scalar, non-negative integer representing the burn-in of the MCMC
//' @param print_freq Scalar, positive integer representing the frequency of printing the sampled values
//' @return A data frame containing N rows and 2 columns which represent xi1 and the log-posterior density (lpost)
//' @seealso \code{\link{mcmc_mix}} for MCMC for the discrete extreme value mixture distribution.
//' @export
// [[Rcpp::export]]
DataFrame mcmc_upp(const NumericVector x, const int u, const double xi1, const double a_xi1, const double b_xi1, const int N = 20000, const int thin = 10, const int burnin = 20000, const int print_freq = 10000) {
  double xi1_curr = xi1, xi1_prop, lpost_curr = lpost_upp(x, u, xi1_curr, a_xi1, b_xi1), lpost_prop, sd_xi1 = 0.1;
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  NumericVector xi1_vec(N), lpost_vec(N);
  int s, t;
  for (t = 0; t < N * thin + burnin; t++){
	xi1_prop = rnorm(1, xi1_curr, sd_xi1)[0];
	lpost_prop = lpost_upp(x, u, xi1_prop, a_xi1, b_xi1);
	update(xi1_curr, xi1_prop, lpost_curr, lpost_prop, sd_xi1, t, burnin);
	if ((t + 1) % print_freq == 0) {
	  Rcout << "Iteration " << t + 1;
	  Rcout << ": Log-posterior = " << lpost_curr;
	  Rcout << endl;
	  if (t < burnin) {
		Rcout << "xi1 = " << xi1_curr << " (" << sd_xi1 << ")" << endl;
	  }
	}
	if (t >= burnin) {
	  s = t - burnin + 1;
	  if (s % thin == 0) {
		s = s / thin - 1;
		xi1_vec[s] = xi1_curr;
		lpost_vec[s] = lpost_curr;
	  }
	}
  }
  Rcout << "Final check: ldiff = " << lpost_upp(x, u, xi1_curr, a_xi1, b_xi1) - lpost_curr << endl << endl;
  DataFrame output = DataFrame::create(Named("xi1") = xi1_vec,
									   Named("lpost") = lpost_vec);
  return output;
}

const double diff_hzeta(const double alpha, const int u) {
  double d;
  gsl_set_error_handler_off();
  gsl_sf_result result1, result2;
  int status1 = gsl_sf_hzeta_e(alpha, 1, &result1),
    status2 = gsl_sf_hzeta_e(alpha, u+1, &result2);
  if (status1 != GSL_SUCCESS || status2 != GSL_SUCCESS) {
    d = nan("");
  }
  else {
    d = result1.val - result2.val;
  }
  return d;
}

const double llik_low(const NumericVector par, const NumericVector x, const int u, const double phi) {
  // discrete power law <= u
  if (par.size() != 1) {
    stop("llik_low: length of par has to be 1.");
  }
  const double xi1 = par[0], alpha = 1.0 / xi1 + 1.0;
  double l;
  const NumericVector xl = x[x <= u];
  const double nl(xl.size());
  if (xi1 <= 0.0 || u <= 0) {
    l = -INFINITY;
  }
  else {
    l = -alpha * sum(log(xl)) - nl * log(diff_hzeta(alpha, u));
  }
  l += nl * log(1.0 - phi);
  return lnan(l);
}

const double llik_geo(const NumericVector par, const NumericVector x, const int u, const double phi) {
  // geometric <= u
  if (par.size() != 1) {
    stop("llik_geo: length of par has to be 1.");
  }
  const double xi1 = par[0],
    p = 1.0 - exp(-1.0 / xi1), q = 1.0 - p;
  double l;
  const NumericVector xl = x[x <= u];
  const double nl(xl.size());
  if (xi1 <= 0.0 || u <= 0) {
    l = -INFINITY;
  }
  else {
    l = sum(xl) * log(q) + nl * (log(p) - log(q) - log(1.0 - pow(q, u))) ;
  }
  l += nl * log(1.0 - phi);
  return lnan(l);
}

const double llik_gpd(const NumericVector par, const NumericVector x, const int u, const double phi) {
  if (par.size() != 2) {
    stop("llik_gpd: length of par has to be 2.");
  }
  const double xi2 = par[0],
    sig = par[1],
    sigu = sig + xi2 * u;
  double l;
  if (sig <= 0.0 || sigu <= 0.0) {
    l = -INFINITY;
  }
  else {
    const NumericVector xu = x[x > u];
    const double nu(xu.size());
    NumericVector yu, zu;
    if (xi2 != 0.0) {
      yu = 1.0 + xi2 / sigu * (xu - 1.0 - u);
      zu = 1.0 + xi2 / (sigu + xi2 * (xu - 1.0 - u));
      if (is_true(any(yu <= 0.0))) {
        l = -INFINITY;
      }
      else {
        l = sum(log(1.0 - pow(zu, -1.0 / xi2))) - 1.0 / xi2 * sum(log(yu));
      }
    }
    else {
      yu = (xu - 1.0 - u) / sigu;
      l = nu * log(1.0 - exp(-1.0 / sigu)) - sum(yu);
    }
    l += nu * log(phi);
  }
  return lnan(l);
}





// 02) Mixture
const double par2phi(const double xi1, const double xi2, const double sig, const int u, const bool geo) {
  const double
    sigu = sig + xi2 * u,
    denu1 = (xi2 == 0.0) ? (1.0 - exp(-1.0 / sigu)) : (1.0 - pow(1.0 + xi2 / sigu, -1.0 / xi2)),
    alpha = 1.0 / xi1 + 1.0,
    p = 1.0 - exp(-1.0 / xi1), q = 1.0 - p;
  double denom, phi;
  if (geo) {
    denom = 1.0 - pow(q, u);
    phi = 1.0 / (1.0 + xi1 * denom / sigu / pow(q, u));
  }
  else {
    denom = diff_hzeta(alpha, u);
    phi = 1.0 / (1.0 + denu1 * denom * pow(u + 1.0, alpha));
  }
  return phi;
}

//' Probability mass function (PMF) of discrete extreme value mixture distribution
//'
//' \code{dmix} returns the PMF at x for the discrete extreme value mixture distribution.
//' @param x Vector of positive integers
//' @param xi1 Scalar, shape parameter for values below or equal to u
//' @param xi2 Scalar, shape parameter of integer generalised Pareto distribution (IGPD), for values above u
//' @param sig Scalar, scale parameter of IGPD, for values above u
//' @param u Scalar, positive integer threshold
//' @param phi Scalar, exceedance probability of u, between 0.0 and 1.0 exclusive
//' @param geo Boolean. If 'TRUE', the geometric distribution is used for the values below u. If 'FALSE', the discrete power law is used.
//' @param give_log Boolean, whether the PMF should be returned on the log scale. If 'FALSE', the PMF is returned on the original scale.
//' @return A numeric vector of the same length as x
//' @examples
//' dmix(10:15, 2.0, 0.5, 1.0, 12, 0.2, TRUE)
//' dmix(10:15, 2.0, 0.5, 1.0, 12, 0.2, FALSE)
//' dmix(10:15, 2.0, 0.5, 1.0, 12, 0.2, FALSE, TRUE)
//' @seealso \code{\link{Smix}} for the corresponding survival function, \code{\link{dupp}} for the probability mass function of the discrete power law.
//' @export
// [[Rcpp::export]]
const NumericVector dmix(const NumericVector x, const double xi1, const double xi2, const double sig, const int u, const double phi, const bool geo, const bool give_log = false) {
  const double
    sigu = sig + xi2 * u,
    alpha = 1.0 / xi1 + 1.0,
    p = 1.0 - exp(-1.0 / xi1), q = 1.0 - p;
  NumericVector fl; double denom;
  if (geo) {
    denom = 1.0 - pow(q, u);
    fl = log(p) + (x - 1.0) * log(q) + log(1.0 - phi);
  }
  else {
    hzeta_check(alpha, u);
    denom = diff_hzeta(alpha, u);
    fl = - alpha * log(x) - log(denom) + log(1.0 - phi);
  }
  const NumericVector
    y = 1.0 + xi2 / sigu * (x - 1.0 - u),
    z = 1.0 + xi2 / (sigu + xi2 * (x - 1.0 - u)),
    fu = log(1.0 - pow(z, -1.0 / xi2)) - 1.0 / xi2 * log(y) + log(phi),
    ld = ifelse(x <= u, fl, fu);
  NumericVector output;
  if (give_log) {
    output = ld;
  }
  else {
    output = exp(ld);
  }
  return output;
}

//' Survival function of discrete extreme value mixture distribution
//'
//' \code{Smix} returns the survival function at x for the discrete extreme value mixture distribution.
//' @param x Vector of positive integers
//' @param xi1 Scalar, shape parameter for values below or equal to u
//' @param xi2 Scalar, shape parameter of integer generalised Pareto distribution (IGPD), for values above u
//' @param sig Scalar, scale parameter of IGPD, for values above u
//' @param u Scalar, positive integer threshold
//' @param phi Scalar, exceedance probability of u, between 0.0 and 1.0 exclusive
//' @param geo Boolean. If 'TRUE', the geometric distribution is used for the values below u. If 'FALSE', the discrete power law is used.
//' @return A numeric vector of the same length as x
//' @examples
//' Smix(10:15, 2.0, 0.5, 1.0, 12, 0.2, TRUE)
//' Smix(10:15, 2.0, 0.5, 1.0, 12, 0.2, FALSE)
//' @seealso \code{\link{dmix}} for the corresponding probability mass function, \code{\link{Supp}} for the survival function of the discrete power law.
//' @export
// [[Rcpp::export]]
const NumericVector Smix(const NumericVector x, const double xi1, const double xi2, const double sig, const int u, const double phi, const bool geo) {
  const double
    sigu = sig + xi2 * u,
    alpha = 1.0 / xi1 + 1.0,
    p = 1.0 - exp(-1.0 / xi1), q = 1.0 - p;
  NumericVector cP(x.size()); double denom;
  if (geo) {
    denom = 1.0 - pow(q, u);
    cP = 1.0 - exp((x - 1.0) * log(q));
  }
  else {
    hzeta_check(alpha, u);
    denom = diff_hzeta(alpha, u);
    for (int i = 0; i < x.size(); i++) {
      cP[i] = diff_hzeta(alpha, x[i] - 1);
    }
  }
  const NumericVector
    y = 1.0 + xi2 / sigu * (x - 1.0 - u),
    Sl = phi + (1.0 - phi) * (1.0 - cP / denom),
    Su = pow(y, -1.0 / xi2) * phi,
    S = ifelse(x <= u, Sl, Su);
  return S;
}

const double lpost_mix(const NumericVector x,
                       const int u,
                       const double xi1,
                       const double xi2,
                       const double sig,
                       const bool cont,
                       const bool geo,
                       const double a_phi,
                       const double b_phi,
                       const double a_xi1,
                       const double b_xi1,
                       const double m_xi2,
                       const double s_xi2,
                       const double a_sig,
                       const double b_sig,
                       const double pcont) {
  //// check x are in fact +ve integers?
  const NumericVector xu = x[x > u],
    par1 = NumericVector::create(xi1),
    par2 = NumericVector::create(xi2, sig);
  const double nu(xu.size()), n(x.size());
  double l;
  if (u <= 1 || u <= min(x) || u > max(x)) {
    l = -INFINITY;
  }
  else {
    const double phi = cont ? par2phi(xi1, xi2, sig, u, geo) : (nu / n);
    l =
      (geo ?
       llik_geo(par1, x, u, phi) :
       llik_low(par1, x, u, phi)) +
      llik_gpd(par2, x, u, phi) +
      dunif(tv(phi), a_phi, b_phi, true)[0] +
      dunif(tv(xi1), a_xi1, b_xi1, true)[0] +
      dnorm(tv(xi2), m_xi2, s_xi2, true)[0] +
      dgamma(tv(sig), a_sig, 1.0 / b_sig, true)[0] +
      dbinom(ti((int) cont), 1, pcont, true)[0];
  }
  return lnan(l);
}

//' Markov chain Monte Carlo for discrete extreme value mixture distribution
//'
//' \code{mcmc_mix} returns the samples from the joint posterior of the parameters (u, xi1, xi2, sig), for fitting the discrete extreme value mixture distribution (DEVMD) to the data x. The samples are obtained using Markov chain Monte Carlo (MCMC).
//'
//' In the MCMC, a componentwise Metropolis-Hastings algorithm is used. Unlike \code{mcmc_upp}, the threshold u is treated as a parameter in \code{mcmc_mix} and therefore inferred. The 8 hyperparameters are used in the following priors: u is such that the implied exceedance probability phi ~ Uniform(a_phi, b_phi); xi1 ~ Uniform(a_xi1, b_xi1); xi2 ~ Normal(mean = m_xi2, sd = s_xi2); sig ~ Gamma(shape = a_sig, rate = b_sig). If pcont = 0.0, only the unconstrained version of the DEVMD is fitted; if pcont = 1.0, only the continuity constrained version is fitted. Setting pcont between 0.0 and 1.0 allows both versions to be fitted, if model selection between the two is of interest.
//' @param x Vector of positive integers, representing the data
//' @param u Scalar, initial value of the positive integer threshold
//' @param xi1 Scalar, initial value of the parameter for values below or equal to u
//' @param xi2 Scalar, initial value of the shape parameter of the integer generalised Pareto distribution (IGPD), for values above u
//' @param sig Scalar, initial value of the scale parameter of IGPD, for values above u
//' @param cont Boolean, whether the continuity constraint is imposed at u
//' @param geo Boolean. If 'TRUE', the geometric distribution is used for the values below u. If 'FALSE', the discrete power law is used.
//' @param a_phi,b_phi,a_xi1,b_xi1,m_xi2,s_xi2,a_sig,b_sig Scalars, representing the hyperparameters of the prior distributions of the respective parameters. See details for the specification of the priors.
//' @param pcont Scalar, between 0.0 and 1.0, representing the prior probability of the continuity constrained version, for model selection.
//' @param N Scalar, positive integer representing the length of the output chain i.e. the number of rows in the returned data frame
//' @param thin Scalar, positive integer representing the thinning in the MCMC
//' @param burnin Scalar, non-negative integer representing the burn-in of the MCMC
//' @param print_freq Scalar, positive integer representing the frequency of printing the sampled values
//' @return A data frame containing N rows and 7 columns which represent (in this order) the 4 parameters (u, xi1, xi2, sig), the implied exceedance probability (phi), the log-posterior density (lpost), and whether the continuity constraint is imposed (cont).
//' @seealso \code{\link{mcmc_upp}} for MCMC for the discrete power law.
//' @export
// [[Rcpp::export]]
DataFrame mcmc_mix(const NumericVector x,
                   const int u,
                   const double xi1,
                   const double xi2,
                   const double sig,
                   const bool cont,
                   const bool geo,
                   const double a_phi,
                   const double b_phi,
                   const double a_xi1,
                   const double b_xi1,
                   const double m_xi2,
                   const double s_xi2,
                   const double a_sig,
                   const double b_sig,
                   const double pcont,
                   const int N = 20000,
                   const int thin = 100,
                   const int burnin = 20000,
                   const int print_freq = 10000) {
  //// check x are in fact +ve integers?
  int u_curr = u, u_prop;
  double
    xi1_curr = xi1, xi1_prop,
    xi2_curr = xi2, xi2_prop,
    sig_curr = sig, sig_prop,
    lpost_curr, lpost_prop, lprob;
  bool cont_curr = (pcont == 0.0) ? false : (pcont == 1.0) ? true : cont;
  auto lpost = [x, geo, a_phi, b_phi, a_xi1, b_xi1, m_xi2, s_xi2, a_sig, b_sig, pcont](const double u, const double xi1, const double xi2, const double sig, const bool cont) {
    return lpost_mix(x, u, xi1, xi2, sig, cont, geo, a_phi, b_phi, a_xi1, b_xi1, m_xi2, s_xi2, a_sig, b_sig, pcont);
  };
  lpost_curr = lpost(u_curr, xi1_curr, xi2_curr, sig_curr, cont_curr);
  Rcout << "Iteration 0: Log-posterior = " << lpost_curr << endl;
  const double n(x.size());
  IntegerVector u_vec(N);
  NumericVector xu, xi1_vec(N), xi2_vec(N), sig_vec(N), phi_vec(N), lpost_vec(N);
  LogicalVector cont_vec(N);
  const double
    sd0 = 0.1 / sqrt(2.0),
    sd1 = 2.34 / sqrt(2.0);
  vec u_burn(burnin), xi1_burn(burnin), xi2_burn(burnin), sig_burn(burnin);
  double sd_u = 1.0, sd_xi1 = 0.1, sd_xi2, sd_sig, cor2, factor = 10.0, lalpha;
  running_stat<double> cont_stat;
  int s, t;
  for (t = 0; t < N * thin + burnin; t++) {
    // u
    u_prop = (int) (u_curr * exp(rnorm(1, 0.0, sd_u)[0])); // (int) rounds down
    lpost_prop = lpost(u_prop, xi1_curr, xi2_curr, sig_curr, cont_curr);
    lalpha = lpost_prop - lpost_curr +
      lq(u_prop, u_curr, sd_u) - lq(u_curr, u_prop, sd_u);
    if (lalpha > lr1()) {
      u_curr = u_prop;
      lpost_curr = lpost_prop;
      if (t < burnin) {
        sd_u = sqrt(sd_u * sd_u + 3.0 * sd_u * sd_u / factor / sqrt(t + 1.0));
      }
    }
    else {
      if (t < burnin) {
        sd_u = sqrt(sd_u * sd_u - 1.0 * sd_u * sd_u / factor / sqrt(t + 1.0));
      }
    }
    // xi1
    xi1_prop = rnorm(1, xi1_curr, sd_xi1)[0];
    lpost_prop = lpost(u_curr, xi1_prop, xi2_curr, sig_curr, cont_curr);
    update(xi1_curr, xi1_prop, lpost_curr, lpost_prop, sd_xi1, t, burnin, factor);
    // xi2 & sigma
    if (t < burnin / 2 || lr1() < log(0.05)) {
      xi2_prop = rnorm(1, xi2_curr, sd0)[0];
      sig_prop = rnorm(1, sig_curr, sd0)[0];
    }
    else {
      xi2_prop = rnorm(1, xi2_curr, sd1 * sd_xi2)[0];
      sig_prop = rnorm(1, sig_curr + sd_sig / sd_xi2 * cor2 * (xi2_prop - xi2_curr), sd1 * sqrt(1.0 - cor2 * cor2) * sd_sig)[0];
    }
    lpost_prop = lpost(u_curr, xi1_curr, xi2_prop, sig_prop, cont_curr);
    if (lpost_prop - lpost_curr > lr1()) {
      xi2_curr = xi2_prop;
      sig_curr = sig_prop;
      lpost_curr = lpost_prop;
    }
    // cont
    lprob = -log(1.0 + exp(lpost(u_curr, xi1_curr, xi2_curr, sig_curr, false) - lpost(u_curr, xi1_curr, xi2_curr, sig_curr, true)));
    cont_curr = (lr1() < lprob);
    cont_stat((double) cont_curr);
    lpost_curr = lpost(u_curr, xi1_curr, xi2_curr, sig_curr, cont_curr);
    // update
    if (t < burnin) {
      xi2_burn[t] = xi2_curr;
      sig_burn[t] = sig_curr;
      sd_xi2 = sd_curr(xi2_burn, t+1);
      sd_sig = sd_curr(sig_burn, t+1);
      cor2 = cor_curr(xi2_burn, sig_burn, t+1);
    }
    if ((t + 1) % print_freq == 0) {
      Rcout << "Iteration " << t + 1;
      Rcout << ": Log-posterior = " << lpost_curr;
      Rcout << endl;
      if (t < burnin) {
        Rcout << "u = " << u_curr << " (" << sd_u << ")" << endl;
        Rcout << "xi1 = " << xi1_curr << " (" << sd_xi1 << ")" << endl;
        Rcout << "xi2 = " << xi2_curr << " (" << sd_xi2 << ")" << endl;
        Rcout << "sig = " << sig_curr << " (" << sd_sig << ")" << endl;
        Rcout << "cont = " << cont_stat.mean() << endl;
      }
    }
    if (t >= burnin) {
      s = t - burnin + 1;
      if (s % thin == 0) {
        s = s / thin - 1;
        u_vec[s] = u_curr;
        xi1_vec[s] = xi1_curr;
        xi2_vec[s] = xi2_curr;
        sig_vec[s] = sig_curr;
        xu = x[x > u_curr];
        phi_vec[s] = cont_curr ? par2phi(xi1_curr, xi2_curr, sig_curr, u_curr, geo) : (xu.size() / n);
        lpost_vec[s] = lpost_curr;
        cont_vec[s] = cont_curr;
      }
    }
  }
  Rcout << "Final check: ldiff = " << lpost(u_curr, xi1_curr, xi2_curr, sig_curr, cont_curr) - lpost_curr << endl << endl;
  DataFrame output = DataFrame::create(Named("u") = u_vec,
                                       Named("xi1") = xi1_vec,
                                       Named("xi2") = xi2_vec,
                                       Named("sig") = sig_vec,
                                       Named("phi") = phi_vec,
                                       Named("lpost") = lpost_vec,
                                       Named("cont") = cont_vec);
  return output;
}
