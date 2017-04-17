#include "BSFG_types.h"

using namespace Rcpp;

struct randomEffect_C_Cholesky{
  SpMat chol_Ci,chol_K_inv;
  MatrixXd P;
};

class randomEffect_C_Cholesky_database {
public:
  randomEffect_C_Cholesky_database(List chol_Ki_mats_list, MatrixXd h2s_matrix, MSpMat ZtZ,int grainSize);
  SpMat get_chol_Ci(int i) {return randomEffect_C_Cholesky_list[i-1].chol_Ci;}
  SpMat get_chol_K_inv_i(int i) {return randomEffect_C_Cholesky_list[i-1].chol_K_inv;}
private:
  VectorXd h2s;
  std::vector<randomEffect_C_Cholesky> randomEffect_C_Cholesky_list;
};

struct Sigma_Cholesky{
  SpMat chol_Sigma;
  double log_det;
};


class Sigma_Cholesky_database {
public:
  Sigma_Cholesky_database(List chol_Ki_mats_list, MatrixXd h2s_matrix, int grainSize);
  SpMat get_chol_Sigma(int i) {return Sigma_Cholesky_list[i-1].chol_Sigma;}
  double get_log_det(int i) {return Sigma_Cholesky_list[i-1].log_det;}
private:
  VectorXd h2s;
  std::vector<Sigma_Cholesky> Sigma_Cholesky_list;
};
