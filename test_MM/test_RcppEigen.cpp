// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

// using namespace std;
// using namespace Rcpp;
// using namespace RcppParallel;

// typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> SparseCholesky;
// typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> SparseCholesky;
// typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
// typedef CPivQR::PermutationType Permutation;

// // [[Rcpp::export]]
// Eigen::VectorXd my_C_solve(SparseCholesky C, Eigen::VectorXd b){
// 	Eigen::VectorXd x = C.solve(b);
// 	return(x);
// }


typedef Eigen::MappedSparseMatrix< double > mappedSparseMatrix ;
typedef Eigen::Map< Eigen::VectorXd > mappedVector ;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd cgSparse(
    const mappedSparseMatrix A,
    const mappedVector b
) {
    Eigen::ConjugateGradient< mappedSparseMatrix, Eigen::Lower > cg( A ) ;
    return cg.solve( b ) ;
}

// NumericVector sample_MME_single_diagA_c(Eigen::VectorXd y, Eigen::MatrixXd W,  Eigen::SparseMatrix C,
// 				Eigen::SparseMatrix RinvSqW, Eigen::VectorXd prior_mean, Eigen::VectorXd prior_prec,
// 				Cholesky Cholesky_R, Eigen::SparseMatrix chol_R, Eigen::SparseMatrix R_Perm, double tot_Y_prec, 
// 				Eigen::VectorXd randn_theta, Eigen::VectorXd randn_e)
// 				) {
// 	int n = y.size();
// 	int n_theta = prior_prec.size();

// 	Eigen::VectorXd theta_star = prior_mean + randn_theta ./ sqrt(prior_prec);
// 	Eigen::VectorXd e_star = chol_R * (randn_e / sqrt(tot_Y_prec));
// 	Eigen::MatrixXd W_theta_star = W * theta_star;
// 	Eigen::VectorXd y_resid = R_perm * (y - W_theta_star - e_star);

// 	y_resid = Cholesky_R.solve(y_resid);
// 	Eigen::MatrixXd WtRinvy = crossprod(RinvSqW, y_resid) * tot_Y_prec



// 	Cholesky Cholesky_C, Eigen::VectorXd pe, Eigen::SparseMatrix chol_A_inv,
// 				double tot_Y_prec, Eigen::VectorXd randn_theta, Eigen::VectorXd randn_e)
// 				) {
// 	int n = Y.rows();
// 	int p = Y.cols();
// 	n_theta = chol_A_inv.rows();

// 	ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
	
// 	Eigen::VectorXd theta_star = solver.compute(chol_A_inv).solve(randn_theta);
// 	Eigen::VectorXd e_star = randn_e ./ sqrt(pe);
// 	Eigen::SparseMatrix W_theta_star = W * theta_star;
// 	Eigen


// 	x = solver.compute(A).solve(b);
// 	typedef Eigen::SimplicialLDLT<SpMat> Cholesky;
// }
// sample_MME_single_diagA = function(y, W, C, RinvSqW, prior_mean,prior_prec,Cholesky_R,chol_R,R_Perm,tot_Y_prec,randn_theta = NULL, randn_e = NULL) {
// 	# R is aZAZ + bI
// 	# 	 then form chol_R
// 	# 	 checkh whether solve(t(chol_Rinv),theta) or chol_R %*% theta is better.
// 	# G is diagonal (fixed effect) - prior prec
// 	# prior_prec must be > 0
// 	n = length(y)
// 	n_theta = length(prior_prec)
// 	if(is.null(randn_theta)) {
// 	  randn_theta = rnorm(n_theta)
// 	}
// 	if(is.null(randn_e)) {
// 	  randn_e = rnorm(n)
// 	}
// 	theta_star = prior_mean + randn_theta/sqrt(prior_prec)
// 	e_star = (chol_R) %*% randn_e / sqrt(tot_Y_prec)
// 	W_theta_star = W %*% theta_star
// 	y_resid = y - W_theta_star - e_star@x
// 	if(!is.null(R_Perm)) {
// 		y_resid_p = R_Perm %*% y_resid
// 	} else{
// 		y_resid_p = y_resid
// 	}

// 	WtRinvy = crossprod(RinvSqW, solve(Cholesky_R,y_resid_p,'L')) * tot_Y_prec

// 	theta_tilda = solve(C,WtRinvy)

// 	theta = theta_tilda@x + theta_star
// 	theta
// }