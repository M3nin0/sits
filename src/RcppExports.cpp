// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bayes_smoother
arma::mat bayes_smoother(const arma::mat& m, const arma::uword m_nrow, const arma::uword m_ncol, const arma::mat& w, const arma::mat& sigma, bool covar_sigma0);
RcppExport SEXP _sits_bayes_smoother(SEXP mSEXP, SEXP m_nrowSEXP, SEXP m_ncolSEXP, SEXP wSEXP, SEXP sigmaSEXP, SEXP covar_sigma0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_nrow(m_nrowSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_ncol(m_ncolSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type covar_sigma0(covar_sigma0SEXP);
    rcpp_result_gen = Rcpp::wrap(bayes_smoother(m, m_nrow, m_ncol, w, sigma, covar_sigma0));
    return rcpp_result_gen;
END_RCPP
}
// kernel_smoother
arma::mat kernel_smoother(const arma::mat& m, const arma::uword m_nrow, const arma::uword m_ncol, const arma::mat& w, const bool normalised);
RcppExport SEXP _sits_kernel_smoother(SEXP mSEXP, SEXP m_nrowSEXP, SEXP m_ncolSEXP, SEXP wSEXP, SEXP normalisedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_nrow(m_nrowSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_ncol(m_ncolSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const bool >::type normalised(normalisedSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_smoother(m, m_nrow, m_ncol, w, normalised));
    return rcpp_result_gen;
END_RCPP
}
// bilateral_smoother
arma::mat bilateral_smoother(const arma::mat& m, const arma::uword m_nrow, const arma::uword m_ncol, const arma::mat& w, double tau);
RcppExport SEXP _sits_bilateral_smoother(SEXP mSEXP, SEXP m_nrowSEXP, SEXP m_ncolSEXP, SEXP wSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_nrow(m_nrowSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_ncol(m_ncolSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(bilateral_smoother(m, m_nrow, m_ncol, w, tau));
    return rcpp_result_gen;
END_RCPP
}
// average_probs
NumericMatrix average_probs(const List& data_lst);
RcppExport SEXP _sits_average_probs(SEXP data_lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    rcpp_result_gen = Rcpp::wrap(average_probs(data_lst));
    return rcpp_result_gen;
END_RCPP
}
// weighted_probs
NumericMatrix weighted_probs(const List& data_lst, const NumericVector& weights);
RcppExport SEXP _sits_weighted_probs(SEXP data_lstSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_probs(data_lst, weights));
    return rcpp_result_gen;
END_RCPP
}
// weighted_uncert_probs
NumericMatrix weighted_uncert_probs(const List& data_lst, const List& unc_lst);
RcppExport SEXP _sits_weighted_uncert_probs(SEXP data_lstSEXP, SEXP unc_lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    Rcpp::traits::input_parameter< const List& >::type unc_lst(unc_lstSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_uncert_probs(data_lst, unc_lst));
    return rcpp_result_gen;
END_RCPP
}
// entropy_probs
IntegerVector entropy_probs(const IntegerMatrix& mtx, const int& n);
RcppExport SEXP _sits_entropy_probs(SEXP mtxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(entropy_probs(mtx, n));
    return rcpp_result_gen;
END_RCPP
}
// kernel_fun
NumericVector kernel_fun(const NumericMatrix& data, const int band, const int img_nrow, const int img_ncol, const int window_size, const int fun);
RcppExport SEXP _sits_kernel_fun(SEXP dataSEXP, SEXP bandSEXP, SEXP img_nrowSEXP, SEXP img_ncolSEXP, SEXP window_sizeSEXP, SEXP funSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< const int >::type img_nrow(img_nrowSEXP);
    Rcpp::traits::input_parameter< const int >::type img_ncol(img_ncolSEXP);
    Rcpp::traits::input_parameter< const int >::type window_size(window_sizeSEXP);
    Rcpp::traits::input_parameter< const int >::type fun(funSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_fun(data, band, img_nrow, img_ncol, window_size, fun));
    return rcpp_result_gen;
END_RCPP
}
// label_max_prob
arma::ucolvec label_max_prob(const arma::mat& x);
RcppExport SEXP _sits_label_max_prob(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(label_max_prob(x));
    return rcpp_result_gen;
END_RCPP
}
// least_probs
IntegerVector least_probs(const IntegerMatrix& mtx, const int& n);
RcppExport SEXP _sits_least_probs(SEXP mtxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(least_probs(mtx, n));
    return rcpp_result_gen;
END_RCPP
}
// linear_interp
IntegerMatrix linear_interp(IntegerMatrix& mtx);
RcppExport SEXP _sits_linear_interp(SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix& >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_interp(mtx));
    return rcpp_result_gen;
END_RCPP
}
// linear_interp_vec
IntegerVector linear_interp_vec(IntegerVector& vec);
RcppExport SEXP _sits_linear_interp_vec(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_interp_vec(vec));
    return rcpp_result_gen;
END_RCPP
}
// margin_probs
IntegerVector margin_probs(const IntegerMatrix& mtx, const int& n);
RcppExport SEXP _sits_margin_probs(SEXP mtxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(margin_probs(mtx, n));
    return rcpp_result_gen;
END_RCPP
}
// nnls_solver
arma::mat nnls_solver(const arma::mat x, const arma::mat A, const bool rmse, const int iterate, const float tolerance);
RcppExport SEXP _sits_nnls_solver(SEXP xSEXP, SEXP ASEXP, SEXP rmseSEXP, SEXP iterateSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const bool >::type rmse(rmseSEXP);
    Rcpp::traits::input_parameter< const int >::type iterate(iterateSEXP);
    Rcpp::traits::input_parameter< const float >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(nnls_solver(x, A, rmse, iterate, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// normalize_data
NumericMatrix normalize_data(const NumericMatrix& data, const double& min, const double& max);
RcppExport SEXP _sits_normalize_data(SEXP dataSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_data(data, min, max));
    return rcpp_result_gen;
END_RCPP
}
// ratio_probs
IntegerVector ratio_probs(const IntegerMatrix& mtx, const int& n);
RcppExport SEXP _sits_ratio_probs(SEXP mtxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type mtx(mtxSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(ratio_probs(mtx, n));
    return rcpp_result_gen;
END_RCPP
}
// max_sampling
DataFrame max_sampling(const IntegerMatrix& data, const int band, const int img_nrow, const int img_ncol, const int window_size);
RcppExport SEXP _sits_max_sampling(SEXP dataSEXP, SEXP bandSEXP, SEXP img_nrowSEXP, SEXP img_ncolSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< const int >::type img_nrow(img_nrowSEXP);
    Rcpp::traits::input_parameter< const int >::type img_ncol(img_ncolSEXP);
    Rcpp::traits::input_parameter< const int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(max_sampling(data, band, img_nrow, img_ncol, window_size));
    return rcpp_result_gen;
END_RCPP
}
// smooth_sg
arma::vec smooth_sg(const arma::vec& data, const arma::mat& f_res, const int& p, const int& n);
RcppExport SEXP _sits_smooth_sg(SEXP dataSEXP, SEXP f_resSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_res(f_resSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_sg(data, f_res, p, n));
    return rcpp_result_gen;
END_RCPP
}
// smooth_sg_mtx
arma::mat smooth_sg_mtx(const arma::mat& data, const arma::mat& f_res, const int& p, const int& n);
RcppExport SEXP _sits_smooth_sg_mtx(SEXP dataSEXP, SEXP f_resSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_res(f_resSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_sg_mtx(data, f_res, p, n));
    return rcpp_result_gen;
END_RCPP
}
// smooth_whit
NumericVector smooth_whit(const NumericVector& data, const double& lambda, const int& length);
RcppExport SEXP _sits_smooth_whit(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_whit(data, lambda, length));
    return rcpp_result_gen;
END_RCPP
}
// smooth_whit_mtx
NumericMatrix smooth_whit_mtx(const NumericMatrix& data, const double& lambda, const int& length);
RcppExport SEXP _sits_smooth_whit_mtx(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_whit_mtx(data, lambda, length));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sits_bayes_smoother", (DL_FUNC) &_sits_bayes_smoother, 6},
    {"_sits_kernel_smoother", (DL_FUNC) &_sits_kernel_smoother, 5},
    {"_sits_bilateral_smoother", (DL_FUNC) &_sits_bilateral_smoother, 5},
    {"_sits_average_probs", (DL_FUNC) &_sits_average_probs, 1},
    {"_sits_weighted_probs", (DL_FUNC) &_sits_weighted_probs, 2},
    {"_sits_weighted_uncert_probs", (DL_FUNC) &_sits_weighted_uncert_probs, 2},
    {"_sits_entropy_probs", (DL_FUNC) &_sits_entropy_probs, 2},
    {"_sits_kernel_fun", (DL_FUNC) &_sits_kernel_fun, 6},
    {"_sits_label_max_prob", (DL_FUNC) &_sits_label_max_prob, 1},
    {"_sits_least_probs", (DL_FUNC) &_sits_least_probs, 2},
    {"_sits_linear_interp", (DL_FUNC) &_sits_linear_interp, 1},
    {"_sits_linear_interp_vec", (DL_FUNC) &_sits_linear_interp_vec, 1},
    {"_sits_margin_probs", (DL_FUNC) &_sits_margin_probs, 2},
    {"_sits_nnls_solver", (DL_FUNC) &_sits_nnls_solver, 5},
    {"_sits_normalize_data", (DL_FUNC) &_sits_normalize_data, 3},
    {"_sits_ratio_probs", (DL_FUNC) &_sits_ratio_probs, 2},
    {"_sits_max_sampling", (DL_FUNC) &_sits_max_sampling, 5},
    {"_sits_smooth_sg", (DL_FUNC) &_sits_smooth_sg, 4},
    {"_sits_smooth_sg_mtx", (DL_FUNC) &_sits_smooth_sg_mtx, 4},
    {"_sits_smooth_whit", (DL_FUNC) &_sits_smooth_whit, 3},
    {"_sits_smooth_whit_mtx", (DL_FUNC) &_sits_smooth_whit_mtx, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_sits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
