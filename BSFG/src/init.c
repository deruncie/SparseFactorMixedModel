#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BSFG_find_candidate_states(SEXP, SEXP, SEXP);
extern SEXP BSFG_get_fitted_set_c(SEXP, SEXP, SEXP);
extern SEXP BSFG_LDLt_notSparse(SEXP);
extern SEXP BSFG_LDLt_sparse(SEXP);
extern SEXP BSFG_log_p_h2s(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_rgig_multiple(SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_rstdnorm_mat(SEXP, SEXP);
extern SEXP BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_set_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_delta_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_factors_scores_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s(SEXP, SEXP);
extern SEXP BSFG_sample_h2s_discrete_MH_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_fixedEffects_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagK(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagR(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_ZKZts_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_SxD(SEXP, SEXP);
extern SEXP BSFG_SxS(SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rcpp_module_boot_Pre_calculations(void);

static const R_CallMethodDef CallEntries[] = {
    {"BSFG_find_candidate_states",                             (DL_FUNC) &BSFG_find_candidate_states,                             3},
    {"BSFG_get_fitted_set_c",                                  (DL_FUNC) &BSFG_get_fitted_set_c,                                  3},
    {"BSFG_LDLt_notSparse",                                    (DL_FUNC) &BSFG_LDLt_notSparse,                                    1},
    {"BSFG_LDLt_sparse",                                       (DL_FUNC) &BSFG_LDLt_sparse,                                       1},
    {"BSFG_log_p_h2s",                                         (DL_FUNC) &BSFG_log_p_h2s,                                         5},
    {"BSFG_rgig_multiple",                                     (DL_FUNC) &BSFG_rgig_multiple,                                     4},
    {"BSFG_rstdnorm_mat",                                      (DL_FUNC) &BSFG_rstdnorm_mat,                                      2},
    {"BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen", (DL_FUNC) &BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen, 8},
    {"BSFG_sample_coefs_set_c",                                (DL_FUNC) &BSFG_sample_coefs_set_c,                                5},
    {"BSFG_sample_delta_c_Eigen",                              (DL_FUNC) &BSFG_sample_delta_c_Eigen,                              6},
    {"BSFG_sample_factors_scores_c",                           (DL_FUNC) &BSFG_sample_factors_scores_c,                           5},
    {"BSFG_sample_h2s",                                        (DL_FUNC) &BSFG_sample_h2s,                                        2},
    {"BSFG_sample_h2s_discrete_MH_c",                          (DL_FUNC) &BSFG_sample_h2s_discrete_MH_c,                          8},
    {"BSFG_sample_MME_fixedEffects_c",                         (DL_FUNC) &BSFG_sample_MME_fixedEffects_c,                         8},
    {"BSFG_sample_MME_single_diagK",                           (DL_FUNC) &BSFG_sample_MME_single_diagK,                           7},
    {"BSFG_sample_MME_single_diagR",                           (DL_FUNC) &BSFG_sample_MME_single_diagR,                           5},
    {"BSFG_sample_MME_ZKZts_c",                                (DL_FUNC) &BSFG_sample_MME_ZKZts_c,                                7},
    {"BSFG_SxD",                                               (DL_FUNC) &BSFG_SxD,                                               2},
    {"BSFG_SxS",                                               (DL_FUNC) &BSFG_SxS,                                               2},
    {"BSFG_tot_prec_scores",                                   (DL_FUNC) &BSFG_tot_prec_scores,                                   4},
    {"_rcpp_module_boot_Pre_calculations",                     (DL_FUNC) &_rcpp_module_boot_Pre_calculations,                     0},
    {NULL, NULL, 0}
};

void R_init_BSFG(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
