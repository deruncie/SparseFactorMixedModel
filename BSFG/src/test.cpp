// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
//
// struct R_matrix2 {
//   Map<MatrixXd> dense;
//   MSpMat sparse;
//   bool isDense;
//   R_matrix2(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// };
//
//
// void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix2>& X_vector){
//   // null_matrices
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   Map<MatrixXd> M_null_d(null_d.data(),0,0);
//   SpMat null_s = null_d.sparseView();
//   MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     if(Rf_isMatrix(Xi_)){
//       Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
//       R_matrix2 Xim(Xi,M_null_s,true);
//       X_vector.push_back(Xim);
//     } else{
//       MSpMat Xi = as<MSpMat>(Xi_);
//       R_matrix2 Xim(M_null_d,Xi,false);
//       X_vector.push_back(Xim);
//     }
//   }
// }
//
// SpMat make_chol_R2(const std::vector<R_matrix2>& ZKZts, const VectorXd h2s, const double tol){  //std::vector<Map<MatrixXd> > ZKZts
//   // Map<MatrixXd> ZKZts_0 = as<Map<MatrixXd> >(ZKZts[0]);
//   int n;
//   bool dense = false;
//   if(ZKZts[0].isDense) {
//     n = ZKZts[0].dense.rows();
//     dense = true;
//   } else{
//     n = ZKZts[0].sparse.rows();
//   }
//   int h = h2s.size();
//   MatrixXd Rd(n,n);
//   SpMat Rs(n,n);
//   Rs.setZero();
//   if(dense) {
//     Rd.setZero();
//   }
//   for(int i = 0; i < h; i++){
//     if(ZKZts[i].isDense) {
//       if(!dense) {
//         Rd = Rs.toDense();
//         dense = true;
//       }
//       Rd += h2s[i] * ZKZts[i].dense;
//     } else{
//       if(dense) {
//         Rd += h2s[i] * ZKZts[i].sparse;
//       } else{
//         Rs += h2s[i] * ZKZts[i].sparse;
//       }
//     }
//   }
//   if(dense) {
//     Rd.diagonal().array() += (1.0-h2s.sum());
//     // return Rd.sparseView(0,tol);
//     Eigen::LLT<MatrixXd> chol_R(Rd);
//     MatrixXd chol_R_U = chol_R.matrixU();
//     return chol_R_U.sparseView(0,tol);
//   } else{
//     for(int i = 0; i < n; i++){
//       Rs.coeffRef(i,i) += (1.0-h2s.sum());
//     }
//     // return Rs;
//     // diagonal().array() += (1.0-h2s.sum());
//     Eigen::SimplicialLLT<SpMat,Lower,NaturalOrdering<int> > chol_R(Rs);
//     MatrixXd chol_R_U = chol_R.matrixL();
//     Rcout << chol_R_U.diagonal() << std::endl;
//     return chol_R_U.sparseView(0,tol);
//   }
// }
//
//
// struct make_chol_R_worker2 : public RcppParallel::Worker {
//   const std::vector<R_matrix2>& ZKZts;
//   const Map<MatrixXd> h2s_matrix;
//   const double tol;
//   std::vector<SpMat>& chol_R_list;
//
//   make_chol_R_worker2(
//     const std::vector<R_matrix2>& ZKZts, //const std::vector<Map<MatrixXd> > ZKZts,
//     const Map<MatrixXd> h2s_matrix,
//     const double tol,
//     std::vector<SpMat>& chol_R_list
//   ):
//     ZKZts(ZKZts), h2s_matrix(h2s_matrix), tol(tol),
//     chol_R_list(chol_R_list) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       chol_R_list[j] = make_chol_R2(ZKZts, h2s_matrix.col(j), tol);
//     }
//   }
// };
//
// // [[Rcpp::export()]]
// Rcpp::List make_chol_V_list2(Rcpp::List ZKZts_,
//                             Map<MatrixXd> h2s_matrix,
//                             double drop0_tol,
//                             SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
//                             int ncores) {
//   int s = h2s_matrix.cols();
//
//   std::vector<R_matrix2> ZKZts;
//   load_R_matrices_list(ZKZts_,ZKZts);
//
//   std::vector<SpMat> chol_R_list;
//   chol_R_list.reserve(s);
//   for(int i = 0; i < s; i++){
//     chol_R_list.push_back(SpMat(0,0));
//   }
//
//   // int n_groups = s/ncores;
//   // for(int i = 0; i <= n_groups; i++){
//   //   make_chol_R_worker2 make_chols(ZKZts,h2s_matrix,drop0_tol,chol_R_list);
//   //   int start = ncores*i;
//   //   int end = std::min(static_cast<double>(s),ncores*(i+1.0));
//   //   RcppParallel::parallelFor(start,end,make_chols,1);
//   //   int pb_state = as<int>(getTxtProgressBar(pb));
//   //   setTxtProgressBar(pb,pb_state+(end-start));
//   // }
//   for(int j = 0; j < s; j++) {
//     chol_R_list[j] = make_chol_R2(ZKZts, h2s_matrix.col(j), drop0_tol);
//   }
//
//   Rcpp::List chol_R_list_out(s);
//   for(int i = 0; i < s; i++){
//     chol_R_list_out[i] = chol_R_list[i];
//   }
//
//   return(chol_R_list_out);
// }
