// Rcpp interface and native registration. Equivalent to compileAttributes().
#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
Rcpp::List egcar_native_solve(Rcpp::List context, Rcpp::List state,
                             Rcpp::List controls, bool group, bool verbose);
Rcpp::NumericVector egcar_native_block_product(Rcpp::List C,
    Rcpp::IntegerVector edge_k, Rcpp::IntegerVector edge_l,
    Rcpp::List local, Rcpp::List position, Rcpp::NumericVector v);
RcppExport SEXP _egcar_egcar_native_block_product(SEXP C, SEXP k, SEXP l,
                                                SEXP local, SEXP position, SEXP v) {
  BEGIN_RCPP
  return Rcpp::wrap(egcar_native_block_product(Rcpp::as<Rcpp::List>(C),
    Rcpp::as<Rcpp::IntegerVector>(k), Rcpp::as<Rcpp::IntegerVector>(l),
    Rcpp::as<Rcpp::List>(local), Rcpp::as<Rcpp::List>(position),
    Rcpp::as<Rcpp::NumericVector>(v)));
  END_RCPP
}
RcppExport SEXP _egcar_egcar_native_solve(SEXP contextSEXP, SEXP stateSEXP,
                                       SEXP controlsSEXP, SEXP groupSEXP,
                                       SEXP verboseSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter<Rcpp::List>::type context(contextSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type state(stateSEXP);
  Rcpp::traits::input_parameter<Rcpp::List>::type controls(controlsSEXP);
  Rcpp::traits::input_parameter<bool>::type group(groupSEXP);
  Rcpp::traits::input_parameter<bool>::type verbose(verboseSEXP);
  rcpp_result_gen = Rcpp::wrap(egcar_native_solve(context, state, controls, group, verbose));
  return rcpp_result_gen;
  END_RCPP
}
static const R_CallMethodDef CallEntries[] = {
  {"_egcar_egcar_native_block_product", (DL_FUNC) &_egcar_egcar_native_block_product, 6},
  {"_egcar_egcar_native_solve", (DL_FUNC) &_egcar_egcar_native_solve, 5},
  {NULL, NULL, 0}
};
RcppExport void R_init_egcar(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
