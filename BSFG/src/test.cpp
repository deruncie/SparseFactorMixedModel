// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
// // // // // // // // // #include<Eigen/SparseCholesky>
//
// // [[Rcpp::export()]]
// VectorXi sample_h2s2(
//     Map<ArrayXXd> log_ps,
//     int grainSize
// )
// {
//
//   int p = log_ps.cols();
//   VectorXd rs = as<VectorXd>(runif(p));
//
//   VectorXi h2s_index(p);
//
//
//   struct sampleColumn : public RcppParallel::Worker {
//     ArrayXXd log_ps;
//     VectorXd rs;
//     VectorXi &h2s_index;
//
//     sampleColumn(ArrayXXd log_ps,
//                  VectorXd rs,
//                  VectorXi &h2s_index):
//       log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       int b = log_ps.rows();
//       for(std::size_t j = begin; j < end; j++){
//         // for(int j = 0; j < p; j++){
//         double max_col = log_ps.col(j).maxCoeff();
//         double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
//         VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
//         h2s_index[j] = 1;
//         double cumsum = 0;
//         for(int i = 0; i < b; i++){
//           cumsum += ps_j[i];
//           if(rs[j] > cumsum) {
//             h2s_index[j] ++;
//           }
//         }
//       }
//     }
//   };
//   sampleColumn sampler(log_ps,rs,h2s_index);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return h2s_index;
// }
