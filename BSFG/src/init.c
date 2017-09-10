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
extern SEXP BSFG_log_p_h2s_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_log_p_h2s_fast_missing(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_rgig_multiple(SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_cis_coefs_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_parallel_sparse_missing_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_set_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_delta_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_factors_scores_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_factors_scores_sparse_missing_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s(SEXP, SEXP);
extern SEXP BSFG_sample_h2s_discrete_MH_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s_discrete_MH_fast_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s_discrete_MH_fast_missing_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_fixedEffects_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_fixedEffects_cis_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagK_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_ZKZts_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_randomEffects_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_randomEffects_parallel_sparse_missing_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_SxD(SEXP, SEXP);
extern SEXP BSFG_SxS(SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores(SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores_c(SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores_missing_c(SEXP, SEXP, SEXP);
extern SEXP BSFG_uncorrelated_prec_mat(SEXP, SEXP, SEXP);
extern SEXP _rcpp_module_boot_Pre_calculations(void);

static const R_CallMethodDef CallEntries[] = {
    {"BSFG_find_candidate_states",                                (DL_FUNC) &BSFG_find_candidate_states,                                 3},
    {"BSFG_get_fitted_set_c",                                     (DL_FUNC) &BSFG_get_fitted_set_c,                                      3},
    {"BSFG_LDLt_notSparse",                                       (DL_FUNC) &BSFG_LDLt_notSparse,                                        1},
    {"BSFG_LDLt_sparse",                                          (DL_FUNC) &BSFG_LDLt_sparse,                                           1},
    {"BSFG_log_p_h2s",                                            (DL_FUNC) &BSFG_log_p_h2s,                                             5},
    {"BSFG_log_p_h2s_fast",                                       (DL_FUNC) &BSFG_log_p_h2s_fast,                                        5},
    {"BSFG_log_p_h2s_fast_missing",                               (DL_FUNC) &BSFG_log_p_h2s_fast_missing,                                5},
    {"BSFG_rgig_multiple",                                        (DL_FUNC) &BSFG_rgig_multiple,                                         4},
    {"BSFG_sample_cis_coefs_parallel_sparse_c_Eigen",             (DL_FUNC) &BSFG_sample_cis_coefs_parallel_sparse_c_Eigen,              8},
    {"BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen",    (DL_FUNC) &BSFG_sample_coefs_hierarchical_parallel_sparse_c_Eigen,    12},
    {"BSFG_sample_coefs_parallel_sparse_c_Eigen",                 (DL_FUNC) &BSFG_sample_coefs_parallel_sparse_c_Eigen,                  6},
    {"BSFG_sample_coefs_parallel_sparse_missing_c_Eigen",         (DL_FUNC) &BSFG_sample_coefs_parallel_sparse_missing_c_Eigen,          8},
    {"BSFG_sample_coefs_set_c",                                   (DL_FUNC) &BSFG_sample_coefs_set_c,                                    7},
    {"BSFG_sample_delta_c_Eigen",                                 (DL_FUNC) &BSFG_sample_delta_c_Eigen,                                  6},
    {"BSFG_sample_factors_scores_sparse_c_Eigen",                 (DL_FUNC) &BSFG_sample_factors_scores_sparse_c_Eigen,                  5},
    {"BSFG_sample_factors_scores_sparse_missing_c_Eigen",         (DL_FUNC) &BSFG_sample_factors_scores_sparse_missing_c_Eigen,          7},
    {"BSFG_sample_h2s",                                           (DL_FUNC) &BSFG_sample_h2s,                                            2},
    {"BSFG_sample_h2s_discrete_MH_c",                             (DL_FUNC) &BSFG_sample_h2s_discrete_MH_c,                             10},
    {"BSFG_sample_h2s_discrete_MH_fast_c",                        (DL_FUNC) &BSFG_sample_h2s_discrete_MH_fast_c,                         8},
    {"BSFG_sample_h2s_discrete_MH_fast_missing_c",                (DL_FUNC) &BSFG_sample_h2s_discrete_MH_fast_missing_c,                 8},
    {"BSFG_sample_MME_fixedEffects_c",                            (DL_FUNC) &BSFG_sample_MME_fixedEffects_c,                            10},
    {"BSFG_sample_MME_fixedEffects_cis_c",                        (DL_FUNC) &BSFG_sample_MME_fixedEffects_cis_c,                        13},
    {"BSFG_sample_MME_single_diagK_c",                            (DL_FUNC) &BSFG_sample_MME_single_diagK_c,                             8},
    {"BSFG_sample_MME_single_diagR",                              (DL_FUNC) &BSFG_sample_MME_single_diagR,                               8},
    {"BSFG_sample_MME_ZKZts_c",                                   (DL_FUNC) &BSFG_sample_MME_ZKZts_c,                                    9},
    {"BSFG_sample_randomEffects_parallel_sparse_c_Eigen",         (DL_FUNC) &BSFG_sample_randomEffects_parallel_sparse_c_Eigen,          6},
    {"BSFG_sample_randomEffects_parallel_sparse_missing_c_Eigen", (DL_FUNC) &BSFG_sample_randomEffects_parallel_sparse_missing_c_Eigen,  6},
    {"BSFG_SxD",                                                  (DL_FUNC) &BSFG_SxD,                                                   2},
    {"BSFG_SxS",                                                  (DL_FUNC) &BSFG_SxS,                                                   2},
    {"BSFG_tot_prec_scores",                                      (DL_FUNC) &BSFG_tot_prec_scores,                                       4},
    {"BSFG_tot_prec_scores_c",                                    (DL_FUNC) &BSFG_tot_prec_scores_c,                                     2},
    {"BSFG_tot_prec_scores_missing_c",                            (DL_FUNC) &BSFG_tot_prec_scores_missing_c,                             3},
    {"BSFG_uncorrelated_prec_mat",                                (DL_FUNC) &BSFG_uncorrelated_prec_mat,                                 3},
    {"_rcpp_module_boot_Pre_calculations",                        (DL_FUNC) &_rcpp_module_boot_Pre_calculations,                         0},
    {NULL, NULL, 0}
};

void R_init_BSFG(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
