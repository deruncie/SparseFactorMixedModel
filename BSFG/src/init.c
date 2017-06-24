#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BSFG_LDLt(SEXP);
extern SEXP BSFG_find_candidate_states(SEXP, SEXP, SEXP);
extern SEXP BSFG_log_p_h2s(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_log_p_h2s_fast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_cis_coefs_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_parallel_sparse_c_Eigen_group(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_coefs_set_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_delta_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_factors_scores_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s(SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_h2s_discrete_MH_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_fixedEffects_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagK_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_single_diagR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_MME_ZKZts_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_sample_randomEffects_parallel_sparse_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores(SEXP, SEXP, SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores_c(SEXP, SEXP, SEXP);
extern SEXP BSFG_tot_prec_scores_withX_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rcpp_module_boot_Pre_calculations(void);

static const R_CallMethodDef CallEntries[] = {
    {"BSFG_LDLt",                                         (DL_FUNC) &BSFG_LDLt,                                          1},
    {"BSFG_find_candidate_states",                        (DL_FUNC) &BSFG_find_candidate_states,                         3},
    {"BSFG_log_p_h2s",                                    (DL_FUNC) &BSFG_log_p_h2s,                                     5},
    {"BSFG_log_p_h2s_fast",                               (DL_FUNC) &BSFG_log_p_h2s_fast,                                5},
    {"BSFG_sample_coefs_parallel_sparse_c_Eigen",         (DL_FUNC) &BSFG_sample_coefs_parallel_sparse_c_Eigen,         11},
    {"BSFG_sample_cis_coefs_parallel_sparse_c_Eigen",     (DL_FUNC) &BSFG_sample_cis_coefs_parallel_sparse_c_Eigen,     14},
    {"BSFG_sample_coefs_parallel_sparse_c_Eigen_group",   (DL_FUNC) &BSFG_sample_coefs_parallel_sparse_c_Eigen_group,   13},
    {"BSFG_sample_coefs_set_c",                           (DL_FUNC) &BSFG_sample_coefs_set_c,                            9},
    {"BSFG_sample_delta_c_Eigen",                         (DL_FUNC) &BSFG_sample_delta_c_Eigen,                          9},
    {"BSFG_sample_factors_scores_sparse_c_Eigen",         (DL_FUNC) &BSFG_sample_factors_scores_sparse_c_Eigen,          6},
    {"BSFG_sample_h2s",                                   (DL_FUNC) &BSFG_sample_h2s,                                    3},
    {"BSFG_sample_h2s_discrete_MH_c",                     (DL_FUNC) &BSFG_sample_h2s_discrete_MH_c,                     10},
    {"BSFG_sample_MME_fixedEffects_c",                    (DL_FUNC) &BSFG_sample_MME_fixedEffects_c,                    10},
    {"BSFG_sample_MME_single_diagK_c",                    (DL_FUNC) &BSFG_sample_MME_single_diagK_c,                     8},
    {"BSFG_sample_MME_single_diagR",                      (DL_FUNC) &BSFG_sample_MME_single_diagR,                       8},
    {"BSFG_sample_MME_ZKZts_c",                           (DL_FUNC) &BSFG_sample_MME_ZKZts_c,                            9},
    {"BSFG_sample_randomEffects_parallel_sparse_c_Eigen", (DL_FUNC) &BSFG_sample_randomEffects_parallel_sparse_c_Eigen,  7},
    {"BSFG_tot_prec_scores",                              (DL_FUNC) &BSFG_tot_prec_scores,                               4},
    {"BSFG_tot_prec_scores_c",                            (DL_FUNC) &BSFG_tot_prec_scores_c,                             3},
    {"BSFG_tot_prec_scores_withX_c",                      (DL_FUNC) &BSFG_tot_prec_scores_withX_c,                       5},
    {"_rcpp_module_boot_Pre_calculations",                (DL_FUNC) &_rcpp_module_boot_Pre_calculations,                 0},
    {NULL, NULL, 0}
};

void R_init_BSFG(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
