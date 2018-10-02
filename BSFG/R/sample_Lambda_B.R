# sample_Lambda_B = function(BSFG_state,grainSize = 1,...) {
#   data_matrices  = BSFG_state$data_matrices
#   priors         = BSFG_state$priors
#   run_parameters = BSFG_state$run_parameters
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state_names = names(current_state)
#   current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
#     k = ncol(Lambda)
#
#     # -----Sample Lambda and B ------------------ #
#     #conditioning on W, F, marginalizing over random effects (conditional on resid_h2)
#     n_coefs = b + k
#     prior_mean = matrix(0,n_coefs,p)
#     if(b > 0) {
#       prior_prec = rbind(B_prec,t(Plam))
#     } else{ # b == 0
#       prior_prec = t(Plam)
#     }
#     Eta_tilde = Eta
#     if(b_QTL > 0){
#       Eta_tilde = Eta - QTL_resid_Z %**% (QTL_resid_X %*% B_QTL)
#     }
#     for(set in seq_along(Missing_data_map)){
#       cols = Missing_data_map[[set]]$Y_cols
#       rows = Missing_data_map[[set]]$Y_obs
#       if(length(cols) == 0 || length(rows) == 0) next
#
#       if(is.null(cis_genotypes)){
#         library(microbenchmark)
#         QtE = Qt_list[[set]] %**% Eta_tilde[rows,cols,drop=FALSE]
#         QtX = cbind(QtX_list[[set]], Qt_list[[set]] %**% F[rows,,drop=FALSE])
#         Sigma_Choleskys_list_dense = list()
#         Sigma_Choleskys_list_dense[[1]] = lapply(Sigma_Choleskys_list[[1]],function(x) within(x,{chol_Sigma = as.matrix(chol_Sigma)}))
#         x = chol(tcrossprod(rstdnorm_mat(n,n)) + diag(1,n))
#         Sigma_Choleskys_list_dense[[1]][[1]]$chol_Sigma = x
#         Sigma_Choleskys_list[[1]][[1]]$chol_Sigma = as(x,'CsparseMatrix')
#         resid_h2_index[] = 1
#         microbenchmark(
#         # coefs = sample_MME_fixedEffects_c(QtE,
#         #                                   QtX,
#         #                                   Sigma_Choleskys_list[[set]],
#         #                                   resid_h2_index[cols],
#         #                                   tot_Eta_prec[cols],
#         #                                   prior_mean[,cols,drop=FALSE],
#         #                                   prior_prec[,cols,drop=FALSE],
#         #                                   grainSize
#         # ) ,
#         # coefs2 = sample_MME_fixedEffects_c2(QtE,
#         #                                     QtX,
#         #                                     Sigma_Choleskys_list_dense[[set]],
#         #                                     resid_h2_index[cols],
#         #                                     tot_Eta_prec[cols],
#         #                                     prior_mean[,cols,drop=FALSE],
#         #                                     prior_prec[,cols,drop=FALSE],
#         #                                     grainSize
#         # ),
#         coefs3 = sample_MME_fixedEffects_c3(QtE,
#                                             QtX,
#                                             Sigma_Choleskys_list_dense[[set]],
#                                             resid_h2_index[cols],
#                                             tot_Eta_prec[cols],
#                                             prior_mean[,cols,drop=FALSE],
#                                             prior_prec[,cols,drop=FALSE],
#                                             grainSize
#         ),
#         coefs4 = sample_MME_fixedEffects_c4(QtE,
#                                             QtX,
#                                             Sigma_Choleskys_list_dense[[set]],
#                                             resid_h2_index[cols],
#                                             tot_Eta_prec[cols],
#                                             prior_mean[,cols,drop=FALSE],
#                                             prior_prec[,cols,drop=FALSE],
#                                             grainSize
#         ),
#         coefs5 = sample_MME_fixedEffects_c5(QtE,
#                                             QtX,
#                                             Sigma_Choleskys_list_dense[[set]],
#                                             resid_h2_index[cols],
#                                             tot_Eta_prec[cols],
#                                             prior_mean[,cols,drop=FALSE],
#                                             prior_prec[,cols,drop=FALSE],
#                                             grainSize
#         )
#         ,times=4)
#         #
#         # microbenchmark({
#         # set.seed(1)
#         # x2=rstdnorm_mat2(n,p)
#         # },{
#         #   set.seed(1)
#         #   x3=rstdnorm_mat3(n,p)
#         # },{
#         #   set.seed(1)
#         #   x4=rstdnorm_mat4(n,p)
#         # },times=100)
#         # microbenchmark(
#         #   test_rstdnorm2(n,p),
#         #   test_rstdnorm3(n,p)
#         # )
#         if(b > 0) {
#           B[,cols] = coefs[1:b,,drop=FALSE]
#         }
#         Lambda[cols,] = t(coefs[b + 1:k,,drop=FALSE])
#       } else{
#         result = sample_MME_fixedEffects_cis_c(Qt_list[[set]] %**% Eta_tilde[rows,cols,drop=FALSE],
#                                                cbind(QtX_list[[set]], Qt_list[[set]] %**% F[rows,,drop=FALSE]),
#                                                Qt_cis_genotypes_list[[set]],
#                                                Sigma_Choleskys_list[[set]],
#                                                resid_h2_index[cols],
#                                                tot_Eta_prec[cols],
#                                                prior_mean[,cols,drop=FALSE],
#                                                prior_prec[,cols,drop=FALSE],
#                                                cis_effects_index[cols],
#                                                sum(n_cis_effects),
#                                                grainSize
#         )
#         coefs = result[[1]]
#         if(b > 0) {
#           B[,cols] = coefs[1:b,,drop=FALSE]
#         }
#         Lambda[cols,] = t(coefs[b + 1:k,,drop=FALSE])
#         cis_index = do.call(c,lapply(cols,function(x) if(n_cis_effects[x] > 0) cis_effects_index[x] + 1:n_cis_effects[x]-1))
#         cis_effects[cis_index] = result[[2]][cis_index]
#       }
#     }
#
#     XB = X %**% B
#     if(!is.null(cis_genotypes)){
#       for(j in 1:p){
#         if(n_cis_effects[j] > 0){
#           cis_X_j = cis_genotypes[[j]]
#           XB[,j] = XB[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]:(cis_effects_index[j+1]-1)]
#         }
#       }
#     }
#
#     # sample B_QTL
#     if(b_QTL > 0){
#       Eta_tilde = Eta - XB - F %**% t(Lambda)
#       prior_mean = matrix(0,b_QTL,p)
#       prior_prec = B_QTL_prec * tot_Eta_prec[rep(1,b_QTL),,drop=FALSE]  # prior for B_QTL includes tot_Eta_prec
#
#       for(set in seq_along(Missing_data_map)){
#         cols = Missing_data_map[[set]]$Y_cols
#         rows = Missing_data_map[[set]]$Y_obs
#         if(length(cols) == 0 || length(rows) == 0) next
#
#         B_QTL[,cols] = sample_MME_fixedEffects_hierarchical_c(Qt_list[[set]] %**% Eta_tilde[rows,cols,drop=FALSE],
#                                                              Qt_QTL_resid_Z_list[[set]],
#                                                              QTL_resid_X,
#                                                              Sigma_Choleskys_list[[set]],
#                                                              resid_h2_index[cols],
#                                                              tot_Eta_prec[cols],
#                                                              prior_mean[,cols,drop=FALSE],
#                                                              prior_prec[,cols,drop=FALSE],
#                                                              grainSize)
#       }
#       XB = XB + QTL_resid_Z %**% (QTL_resid_X %*% B_QTL)
#     }
#
#   }))
#   current_state = current_state[current_state_names]
#
#   return(current_state)
# }
#
