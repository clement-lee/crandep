// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dpol
const NumericVector dpol(const IntegerVector x, const double alpha, const double theta, const int x_max);
RcppExport SEXP _crandep_dpol(SEXP xSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(dpol(x, alpha, theta, x_max));
    return rcpp_result_gen;
END_RCPP
}
// Spol
const NumericVector Spol(const IntegerVector x, const double alpha, const double theta, const int x_max);
RcppExport SEXP _crandep_Spol(SEXP xSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(Spol(x, alpha, theta, x_max));
    return rcpp_result_gen;
END_RCPP
}
// llik_pol
const double llik_pol(const NumericVector par, const IntegerVector x, const IntegerVector count, const bool powerlaw, const int x_max);
RcppExport SEXP _crandep_llik_pol(SEXP parSEXP, SEXP xSEXP, SEXP countSEXP, SEXP powerlawSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw(powerlawSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_pol(par, x, count, powerlaw, x_max));
    return rcpp_result_gen;
END_RCPP
}
// lpost_pol
const double lpost_pol(const IntegerVector x, const IntegerVector count, const double alpha, const double theta, const double a_alpha, const double b_alpha, const double a_theta, const double b_theta, const double powerlaw, const int x_max, double& llik, const double invt);
RcppExport SEXP _crandep_lpost_pol(SEXP xSEXP, SEXP countSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP, SEXP powerlawSEXP, SEXP x_maxSEXP, SEXP llikSEXP, SEXP invtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type powerlaw(powerlawSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< double& >::type llik(llikSEXP);
    Rcpp::traits::input_parameter< const double >::type invt(invtSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_pol(x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, powerlaw, x_max, llik, invt));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_pol
List mcmc_pol(const IntegerVector x, const IntegerVector count, double alpha, double theta, const double a_alpha, const double b_alpha, const double a_theta, const double b_theta, const double a_pseudo, const double b_pseudo, const double pr_power, const int iter, const int thin, const int burn, const int freq, const NumericVector invt, const bool mc3_or_marg, const int x_max);
RcppExport SEXP _crandep_mcmc_pol(SEXP xSEXP, SEXP countSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP, SEXP a_pseudoSEXP, SEXP b_pseudoSEXP, SEXP pr_powerSEXP, SEXP iterSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP freqSEXP, SEXP invtSEXP, SEXP mc3_or_margSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_pseudo(a_pseudoSEXP);
    Rcpp::traits::input_parameter< const double >::type b_pseudo(b_pseudoSEXP);
    Rcpp::traits::input_parameter< const double >::type pr_power(pr_powerSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const int >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type invt(invtSEXP);
    Rcpp::traits::input_parameter< const bool >::type mc3_or_marg(mc3_or_margSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_pol(x, count, alpha, theta, a_alpha, b_alpha, a_theta, b_theta, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg, x_max));
    return rcpp_result_gen;
END_RCPP
}
// llik_bulk
const double llik_bulk(const NumericVector par, const IntegerVector x, const IntegerVector count, const int v, const int u, const double phil, const bool powerlaw, const bool positive);
RcppExport SEXP _crandep_llik_bulk(SEXP parSEXP, SEXP xSEXP, SEXP countSEXP, SEXP vSEXP, SEXP uSEXP, SEXP philSEXP, SEXP powerlawSEXP, SEXP positiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type phil(philSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw(powerlawSEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_bulk(par, x, count, v, u, phil, powerlaw, positive));
    return rcpp_result_gen;
END_RCPP
}
// lpost_bulk
const double lpost_bulk(const NumericVector par, const IntegerVector x, const IntegerVector count, const int v, const int u, const double phil, const double a_alpha, const double b_alpha, const double a_theta, const double b_theta, const bool powerlaw, const bool positive);
RcppExport SEXP _crandep_lpost_bulk(SEXP parSEXP, SEXP xSEXP, SEXP countSEXP, SEXP vSEXP, SEXP uSEXP, SEXP philSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP, SEXP powerlawSEXP, SEXP positiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type phil(philSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw(powerlawSEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_bulk(par, x, count, v, u, phil, a_alpha, b_alpha, a_theta, b_theta, powerlaw, positive));
    return rcpp_result_gen;
END_RCPP
}
// llik_igpd
const double llik_igpd(const NumericVector par, const IntegerVector x, const IntegerVector count, const int u, const double phiu);
RcppExport SEXP _crandep_llik_igpd(SEXP parSEXP, SEXP xSEXP, SEXP countSEXP, SEXP uSEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_igpd(par, x, count, u, phiu));
    return rcpp_result_gen;
END_RCPP
}
// lpost_igpd
const double lpost_igpd(const NumericVector par, const IntegerVector x, const IntegerVector count, const int u, const double m_shape, const double s_shape, const double a_sigma, const double b_sigma, const double phiu);
RcppExport SEXP _crandep_lpost_igpd(SEXP parSEXP, SEXP xSEXP, SEXP countSEXP, SEXP uSEXP, SEXP m_shapeSEXP, SEXP s_shapeSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type m_shape(m_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type s_shape(s_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigma(b_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_igpd(par, x, count, u, m_shape, s_shape, a_sigma, b_sigma, phiu));
    return rcpp_result_gen;
END_RCPP
}
// lpost_mix1
const double lpost_mix1(const IntegerVector x, const IntegerVector count, const int u, const double alpha1, const double theta1, const double alpha2, const double a_psiu, const double b_psiu, const double a_alpha1, const double b_alpha1, const double a_theta1, const double b_theta1, const double a_alpha2, const double b_alpha2, const bool positive, const int x_max, double& llik, const double invt);
RcppExport SEXP _crandep_lpost_mix1(SEXP xSEXP, SEXP countSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alpha1SEXP, SEXP b_alpha1SEXP, SEXP a_theta1SEXP, SEXP b_theta1SEXP, SEXP a_alpha2SEXP, SEXP b_alpha2SEXP, SEXP positiveSEXP, SEXP x_maxSEXP, SEXP llikSEXP, SEXP invtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha1(a_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha1(b_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta1(a_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta1(b_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha2(a_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha2(b_alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< double& >::type llik(llikSEXP);
    Rcpp::traits::input_parameter< const double >::type invt(invtSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_mix1(x, count, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, x_max, llik, invt));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_mix1
List mcmc_mix1(const IntegerVector x, const IntegerVector count, const IntegerVector u_set, int u, double alpha1, double theta1, double alpha2, const double a_psiu, const double b_psiu, const double a_alpha1, const double b_alpha1, const double a_theta1, const double b_theta1, const double a_alpha2, const double b_alpha2, const bool positive, const int iter, const int thin, const int burn, const int freq, const NumericVector invt, const bool mc3_or_marg, const int x_max);
RcppExport SEXP _crandep_mcmc_mix1(SEXP xSEXP, SEXP countSEXP, SEXP u_setSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alpha1SEXP, SEXP b_alpha1SEXP, SEXP a_theta1SEXP, SEXP b_theta1SEXP, SEXP a_alpha2SEXP, SEXP b_alpha2SEXP, SEXP positiveSEXP, SEXP iterSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP freqSEXP, SEXP invtSEXP, SEXP mc3_or_margSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type u_set(u_setSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha1(a_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha1(b_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta1(a_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta1(b_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha2(a_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha2(b_alpha2SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const int >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type invt(invtSEXP);
    Rcpp::traits::input_parameter< const bool >::type mc3_or_marg(mc3_or_margSEXP);
    Rcpp::traits::input_parameter< const int >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_mix1(x, count, u_set, u, alpha1, theta1, alpha2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, positive, iter, thin, burn, freq, invt, mc3_or_marg, x_max));
    return rcpp_result_gen;
END_RCPP
}
// dmix2
const NumericVector dmix2(const IntegerVector x, const int u, const double alpha, const double theta, const double shape, const double sigma, const double phiu);
RcppExport SEXP _crandep_dmix2(SEXP xSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmix2(x, u, alpha, theta, shape, sigma, phiu));
    return rcpp_result_gen;
END_RCPP
}
// Smix2
const NumericVector Smix2(const IntegerVector x, const int u, const double alpha, const double theta, const double shape, const double sigma, const double phiu);
RcppExport SEXP _crandep_Smix2(SEXP xSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(Smix2(x, u, alpha, theta, shape, sigma, phiu));
    return rcpp_result_gen;
END_RCPP
}
// lpost_mix2
const double lpost_mix2(const IntegerVector x, const IntegerVector count, const int u, const double alpha, const double theta, const double shape, const double sigma, const double a_psiu, const double b_psiu, const double a_alpha, const double b_alpha, const double a_theta, const double b_theta, const double m_shape, const double s_shape, const double a_sigma, const double b_sigma, const bool powerlaw, const bool positive, double& llik, const double invt, const bool constrained);
RcppExport SEXP _crandep_lpost_mix2(SEXP xSEXP, SEXP countSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP, SEXP m_shapeSEXP, SEXP s_shapeSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP, SEXP powerlawSEXP, SEXP positiveSEXP, SEXP llikSEXP, SEXP invtSEXP, SEXP constrainedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type m_shape(m_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type s_shape(s_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigma(b_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw(powerlawSEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< double& >::type llik(llikSEXP);
    Rcpp::traits::input_parameter< const double >::type invt(invtSEXP);
    Rcpp::traits::input_parameter< const bool >::type constrained(constrainedSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_mix2(x, count, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, powerlaw, positive, llik, invt, constrained));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_mix2
List mcmc_mix2(const IntegerVector x, const IntegerVector count, const IntegerVector u_set, int u, double alpha, double theta, double shape, double sigma, const double a_psiu, const double b_psiu, const double a_alpha, const double b_alpha, const double a_theta, const double b_theta, const double m_shape, const double s_shape, const double a_sigma, const double b_sigma, const bool positive, const double a_pseudo, const double b_pseudo, double pr_power, const int iter, const int thin, const int burn, const int freq, const NumericVector invt, const bool mc3_or_marg, const bool constrained);
RcppExport SEXP _crandep_mcmc_mix2(SEXP xSEXP, SEXP countSEXP, SEXP u_setSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alphaSEXP, SEXP b_alphaSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP, SEXP m_shapeSEXP, SEXP s_shapeSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP, SEXP positiveSEXP, SEXP a_pseudoSEXP, SEXP b_pseudoSEXP, SEXP pr_powerSEXP, SEXP iterSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP freqSEXP, SEXP invtSEXP, SEXP mc3_or_margSEXP, SEXP constrainedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type u_set(u_setSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha(a_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha(b_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type m_shape(m_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type s_shape(s_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigma(b_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const double >::type a_pseudo(a_pseudoSEXP);
    Rcpp::traits::input_parameter< const double >::type b_pseudo(b_pseudoSEXP);
    Rcpp::traits::input_parameter< double >::type pr_power(pr_powerSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const int >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type invt(invtSEXP);
    Rcpp::traits::input_parameter< const bool >::type mc3_or_marg(mc3_or_margSEXP);
    Rcpp::traits::input_parameter< const bool >::type constrained(constrainedSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_mix2(x, count, u_set, u, alpha, theta, shape, sigma, a_psiu, b_psiu, a_alpha, b_alpha, a_theta, b_theta, m_shape, s_shape, a_sigma, b_sigma, positive, a_pseudo, b_pseudo, pr_power, iter, thin, burn, freq, invt, mc3_or_marg, constrained));
    return rcpp_result_gen;
END_RCPP
}
// dmix3
const NumericVector dmix3(const IntegerVector x, const int v, const int u, const double alpha1, const double theta1, const double alpha2, const double theta2, const double shape, const double sigma, const double phi1, const double phi2, const double phiu);
RcppExport SEXP _crandep_dmix3(SEXP xSEXP, SEXP vSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP theta2SEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP phi1SEXP, SEXP phi2SEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< const double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(dmix3(x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu));
    return rcpp_result_gen;
END_RCPP
}
// Smix3
const NumericVector Smix3(const IntegerVector x, const int v, const int u, const double alpha1, const double theta1, const double alpha2, const double theta2, const double shape, const double sigma, const double phi1, const double phi2, const double phiu);
RcppExport SEXP _crandep_Smix3(SEXP xSEXP, SEXP vSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP theta2SEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP phi1SEXP, SEXP phi2SEXP, SEXP phiuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< const double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< const double >::type phiu(phiuSEXP);
    rcpp_result_gen = Rcpp::wrap(Smix3(x, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, phi1, phi2, phiu));
    return rcpp_result_gen;
END_RCPP
}
// lpost_mix3
const double lpost_mix3(const IntegerVector x, const IntegerVector count, const int v, const int u, const double alpha1, const double theta1, const double alpha2, const double theta2, const double shape, const double sigma, const double a_psi1, const double a_psi2, const double a_psiu, const double b_psiu, const double a_alpha1, const double b_alpha1, const double a_theta1, const double b_theta1, const double a_alpha2, const double b_alpha2, const double a_theta2, const double b_theta2, const double m_shape, const double s_shape, const double a_sigma, const double b_sigma, const bool powerlaw1, const bool powerlaw2, const bool positive1, const bool positive2, double& llik, const double invt);
RcppExport SEXP _crandep_lpost_mix3(SEXP xSEXP, SEXP countSEXP, SEXP vSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP theta2SEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP a_psi1SEXP, SEXP a_psi2SEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alpha1SEXP, SEXP b_alpha1SEXP, SEXP a_theta1SEXP, SEXP b_theta1SEXP, SEXP a_alpha2SEXP, SEXP b_alpha2SEXP, SEXP a_theta2SEXP, SEXP b_theta2SEXP, SEXP m_shapeSEXP, SEXP s_shapeSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP, SEXP powerlaw1SEXP, SEXP powerlaw2SEXP, SEXP positive1SEXP, SEXP positive2SEXP, SEXP llikSEXP, SEXP invtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi1(a_psi1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi2(a_psi2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha1(a_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha1(b_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta1(a_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta1(b_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha2(a_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha2(b_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta2(a_theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta2(b_theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type m_shape(m_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type s_shape(s_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigma(b_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw1(powerlaw1SEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw2(powerlaw2SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive1(positive1SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive2(positive2SEXP);
    Rcpp::traits::input_parameter< double& >::type llik(llikSEXP);
    Rcpp::traits::input_parameter< const double >::type invt(invtSEXP);
    rcpp_result_gen = Rcpp::wrap(lpost_mix3(x, count, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, powerlaw2, positive1, positive2, llik, invt));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_mix3
List mcmc_mix3(const IntegerVector x, const IntegerVector count, const IntegerVector v_set, const IntegerVector u_set, int v, int u, double alpha1, double theta1, double alpha2, double theta2, double shape, double sigma, const double a_psi1, const double a_psi2, const double a_psiu, const double b_psiu, const double a_alpha1, const double b_alpha1, const double a_theta1, const double b_theta1, const double a_alpha2, const double b_alpha2, const double a_theta2, const double b_theta2, const double m_shape, const double s_shape, const double a_sigma, const double b_sigma, const bool powerlaw1, const bool positive1, const bool positive2, const double a_pseudo, const double b_pseudo, const double pr_power2, const int iter, const int thin, const int burn, const int freq, const NumericVector invt, const bool mc3_or_marg);
RcppExport SEXP _crandep_mcmc_mix3(SEXP xSEXP, SEXP countSEXP, SEXP v_setSEXP, SEXP u_setSEXP, SEXP vSEXP, SEXP uSEXP, SEXP alpha1SEXP, SEXP theta1SEXP, SEXP alpha2SEXP, SEXP theta2SEXP, SEXP shapeSEXP, SEXP sigmaSEXP, SEXP a_psi1SEXP, SEXP a_psi2SEXP, SEXP a_psiuSEXP, SEXP b_psiuSEXP, SEXP a_alpha1SEXP, SEXP b_alpha1SEXP, SEXP a_theta1SEXP, SEXP b_theta1SEXP, SEXP a_alpha2SEXP, SEXP b_alpha2SEXP, SEXP a_theta2SEXP, SEXP b_theta2SEXP, SEXP m_shapeSEXP, SEXP s_shapeSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP, SEXP powerlaw1SEXP, SEXP positive1SEXP, SEXP positive2SEXP, SEXP a_pseudoSEXP, SEXP b_pseudoSEXP, SEXP pr_power2SEXP, SEXP iterSEXP, SEXP thinSEXP, SEXP burnSEXP, SEXP freqSEXP, SEXP invtSEXP, SEXP mc3_or_margSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type v_set(v_setSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type u_set(u_setSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< double >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi1(a_psi1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi2(a_psi2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_psiu(a_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psiu(b_psiuSEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha1(a_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha1(b_alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta1(a_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta1(b_theta1SEXP);
    Rcpp::traits::input_parameter< const double >::type a_alpha2(a_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_alpha2(b_alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta2(a_theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta2(b_theta2SEXP);
    Rcpp::traits::input_parameter< const double >::type m_shape(m_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type s_shape(s_shapeSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigma(b_sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type powerlaw1(powerlaw1SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive1(positive1SEXP);
    Rcpp::traits::input_parameter< const bool >::type positive2(positive2SEXP);
    Rcpp::traits::input_parameter< const double >::type a_pseudo(a_pseudoSEXP);
    Rcpp::traits::input_parameter< const double >::type b_pseudo(b_pseudoSEXP);
    Rcpp::traits::input_parameter< const double >::type pr_power2(pr_power2SEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const int >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type invt(invtSEXP);
    Rcpp::traits::input_parameter< const bool >::type mc3_or_marg(mc3_or_margSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_mix3(x, count, v_set, u_set, v, u, alpha1, theta1, alpha2, theta2, shape, sigma, a_psi1, a_psi2, a_psiu, b_psiu, a_alpha1, b_alpha1, a_theta1, b_theta1, a_alpha2, b_alpha2, a_theta2, b_theta2, m_shape, s_shape, a_sigma, b_sigma, powerlaw1, positive1, positive2, a_pseudo, b_pseudo, pr_power2, iter, thin, burn, freq, invt, mc3_or_marg));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_crandep_dpol", (DL_FUNC) &_crandep_dpol, 4},
    {"_crandep_Spol", (DL_FUNC) &_crandep_Spol, 4},
    {"_crandep_llik_pol", (DL_FUNC) &_crandep_llik_pol, 5},
    {"_crandep_lpost_pol", (DL_FUNC) &_crandep_lpost_pol, 12},
    {"_crandep_mcmc_pol", (DL_FUNC) &_crandep_mcmc_pol, 18},
    {"_crandep_llik_bulk", (DL_FUNC) &_crandep_llik_bulk, 8},
    {"_crandep_lpost_bulk", (DL_FUNC) &_crandep_lpost_bulk, 12},
    {"_crandep_llik_igpd", (DL_FUNC) &_crandep_llik_igpd, 5},
    {"_crandep_lpost_igpd", (DL_FUNC) &_crandep_lpost_igpd, 9},
    {"_crandep_lpost_mix1", (DL_FUNC) &_crandep_lpost_mix1, 18},
    {"_crandep_mcmc_mix1", (DL_FUNC) &_crandep_mcmc_mix1, 23},
    {"_crandep_dmix2", (DL_FUNC) &_crandep_dmix2, 7},
    {"_crandep_Smix2", (DL_FUNC) &_crandep_Smix2, 7},
    {"_crandep_lpost_mix2", (DL_FUNC) &_crandep_lpost_mix2, 22},
    {"_crandep_mcmc_mix2", (DL_FUNC) &_crandep_mcmc_mix2, 29},
    {"_crandep_dmix3", (DL_FUNC) &_crandep_dmix3, 12},
    {"_crandep_Smix3", (DL_FUNC) &_crandep_Smix3, 12},
    {"_crandep_lpost_mix3", (DL_FUNC) &_crandep_lpost_mix3, 32},
    {"_crandep_mcmc_mix3", (DL_FUNC) &_crandep_mcmc_mix3, 40},
    {NULL, NULL, 0}
};

RcppExport void R_init_crandep(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
