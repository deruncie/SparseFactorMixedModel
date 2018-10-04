#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BSFG_find_candidate_states(SEXP, SEXP, SEXP);
extern SEXP _BSFG_LDLt(SEXP);
extern SEXP _BSFG_log_p_h2s(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_make_chol_V_list(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_make_chol_ZtZ_Kinv_list(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_matrix_multiply_toDense(SEXP, SEXP);
extern SEXP _BSFG_record_sample_Posterior_array(SEXP, SEXP, SEXP);
extern SEXP _BSFG_regression_sampler_parallel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_rstdnorm_mat(SEXP, SEXP);
extern SEXP _BSFG_sample_factors_scores_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_sample_h2s(SEXP, SEXP);
extern SEXP _BSFG_sample_h2s_discrete_MH_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_sample_MME_single_diagR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_sample_MME_ZKZts_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BSFG_sample_tau2_delta_c_Eigen_v2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BSFG_find_candidate_states",         (DL_FUNC) &_BSFG_find_candidate_states,          3},
    {"_BSFG_LDLt",                          (DL_FUNC) &_BSFG_LDLt,                           1},
    {"_BSFG_log_p_h2s",                     (DL_FUNC) &_BSFG_log_p_h2s,                      5},
    {"_BSFG_make_chol_V_list",              (DL_FUNC) &_BSFG_make_chol_V_list,               7},
    {"_BSFG_make_chol_ZtZ_Kinv_list",       (DL_FUNC) &_BSFG_make_chol_ZtZ_Kinv_list,        8},
    {"_BSFG_matrix_multiply_toDense",       (DL_FUNC) &_BSFG_matrix_multiply_toDense,        2},
    {"_BSFG_record_sample_Posterior_array", (DL_FUNC) &_BSFG_record_sample_Posterior_array,  3},
    {"_BSFG_regression_sampler_parallel",   (DL_FUNC) &_BSFG_regression_sampler_parallel,   15},
    {"_BSFG_rstdnorm_mat",                  (DL_FUNC) &_BSFG_rstdnorm_mat,                   2},
    {"_BSFG_sample_factors_scores_c",       (DL_FUNC) &_BSFG_sample_factors_scores_c,        5},
    {"_BSFG_sample_h2s",                    (DL_FUNC) &_BSFG_sample_h2s,                     2},
    {"_BSFG_sample_h2s_discrete_MH_c",      (DL_FUNC) &_BSFG_sample_h2s_discrete_MH_c,       8},
    {"_BSFG_sample_MME_single_diagR",       (DL_FUNC) &_BSFG_sample_MME_single_diagR,        6},
    {"_BSFG_sample_MME_ZKZts_c",            (DL_FUNC) &_BSFG_sample_MME_ZKZts_c,             7},
    {"_BSFG_sample_tau2_delta_c_Eigen_v2",  (DL_FUNC) &_BSFG_sample_tau2_delta_c_Eigen_v2,   9},
    {NULL, NULL, 0}
};

void R_init_BSFG(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
