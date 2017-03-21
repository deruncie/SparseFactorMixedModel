initialize_BSFG.fast_BSFG = function(BSFG_state, K_mats = NULL, chol_Ki_mats = NULL,verbose=T,...){

    Y          = BSFG_state$data_matrices$Y
    X_F        = BSFG_state$data_matrices$X_F
    Z_matrices = BSFG_state$data_matrices$Z_matrices
    Z          = BSFG_state$data_matrices$Z
    h2s_matrix = BSFG_state$data_matrices$h2s_matrix

    RE_names   = rownames(h2s_matrix)
    n_RE       = length(RE_names)
    stopifnot(n_RE == 1)  # must be only 1 random effect

    traitnames     = BSFG_state$traitnames
    priors         = BSFG_state$priors
    run_parameters = BSFG_state$run_parameters
    run_variables  = BSFG_state$run_variables

    p    = BSFG_state$run_variables$p
    n    = BSFG_state$run_variables$n
    r    = BSFG_state$run_variables$r_RE
    b    = BSFG_state$run_variables$b
    b_F  = BSFG_state$run_variables$b_F

# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- #

   # --- transcript-level model
    # p-vector of gene precisions after taking removing effect of factors.
	#  Prior: Gamma distribution for each element
    #       shape = tot_Eta_prec_shape
    #       rate = tot_Eta_prec_rate
    tot_Eta_prec = with(priors,matrix(rgamma(p,shape = tot_Eta_prec_shape,rate = tot_Eta_prec_rate),nrow = 1))
    colnames(tot_Eta_prec) = traitnames

    # Factors:
    #  initial number of factors
    k = run_parameters$k_init

    # Factor loading precisions (except column penalty tauh).
	 #  Prior: Gamma distribution for each element.
     #       shape = Lambda_df/2
     #       rate = Lambda_df/2
     #    Marginilizes to t-distribution with Lambda_df degrees of freedom
     #    on each factor loading, conditional on tauh
    Lambda_prec = with(priors,matrix(rgamma(p*k,shape = Lambda_df/2,rate = Lambda_df/2),nr = p,nc = k))

    # Factor penalty. tauh(h) = \prod_{i=1}^h \delta_i
	 #  Prior: Gamma distribution for each element of delta
     #     delta_1:
     #       shape = delta_1_shape
     #       rate = delta_1_rate
     #     delta_2 ... delta_m:
     #       shape = delta_2_shape
     #       rate = delta_2_rate
    delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
    tauh  = matrix(cumprod(delta),nrow=1)

    # Total Factor loading precisions Lambda_prec * tauh
    Plam = sweep(Lambda_prec,2,tauh,'*')

    # Lambda - factor loadings
     #   Prior: Normal distribution for each element.
     #       mu = 0
     #       sd = sqrt(1/Plam)
    Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)
    rownames(Lambda) = traitnames

    # p-vector of heritabilities for residuals on factors.
	#  Prior: discrete prior on interval [0,1)
    # h2s = with(run_parameters,seq(0,1,length = h2_divisions+1)[1:h2_divisions])
    # resid_h2 = with(priors,sample(h2s,p,prob = h2_priors_resids, replace = T))
    resid_h2_index = sample(1:ncol(h2s_matrix),p,replace=T)
    resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]

    # Genetic effects not accounted for by factors.
    #   Prior: Normal distribution on each element.
    #       mean = 0
    #       sd = 1./sqrt(U_R_prec)' on each row
    U_R = matrix(rnorm(p*r,0,sqrt(resid_h2 * tot_Eta_prec)),nr = r,nc = p, byrow = T)
    colnames(U_R) = traitnames
    rownames(U_R) = colnames(Z)

    # Latent factor variances
    tot_F_prec = with(priors,matrix(rgamma(k,shape = tot_F_prec_shape,rate = tot_F_prec_rate),nrow=1))

    F_h2_index = sample(1:ncol(h2s_matrix),k,replace=T)
    F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]
    # F_h2 = with(priors,sample(h2s,k,prob = h2_priors_factors, replace = T))

    # Genetic effects on the factor scores.
	#  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(F_h2') for each row.
    U_F = matrix(rnorm(k*r,0,sqrt(F_h2)),nr = r,nc = k, byrow = T)
    rownames(U_F) = colnames(Z)

    # Factor fixed effects
    B_F = matrix(rnorm(b_F * k),b_F,k)

    # Full Factor scores. Combination of genetic and residual variation on
    # each factor.
	#  Prior: Normal distribution for each element
    #       mean = Z   * U_F
    #       sd = sqrt(1-F_h2') for each row.
    F = as.matrix(X_F %*% B_F + Z %*% U_F + matrix(rnorm(k*n,0,sqrt(1-F_h2)),nr = n,nc = k, byrow = T))

    # Fixed effect coefficients.
	#  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(1/fixed_effects_prec)
    B = matrix(rnorm(b*p),nr = b, nc = p)
    colnames(B) = traitnames

    if(b > 0) {
      tau_B = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_prec_shape, rate = priors$fixed_prec_rate)),nrow=1)
      tau_B_F = tau_B[1,-1,drop=FALSE]
    } else{
      tau_B = matrix(0,ncol=0,nrow=1)
      tau_B_F = tau_B
    }
    prec_B = matrix(tau_B,nrow = b, ncol = p)
    prec_B_F = matrix(tau_B_F,nrow = b_F, ncol = k)


# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
            k              = k,
            Lambda_prec    = Lambda_prec,
            delta          = delta,
            tauh           = tauh,
            Plam           = Plam,
            Lambda         = Lambda,
            tot_Eta_prec   = tot_Eta_prec,
            resid_h2_index = resid_h2_index,
            resid_h2       = resid_h2,
            tot_F_prec     = tot_F_prec,
            F_h2_index     = F_h2_index,
            F_h2           = F_h2,
            U_F            = U_F,
            F              = F,
            U_R            = U_R,
            B              = B,
            B_F            = B_F,
            tau_B          = tau_B,
            tau_B_F        = tau_B_F,
            prec_B         = prec_B,
            prec_B_F       = prec_B_F,
            traitnames     = traitnames,
            nrun           = 0,
            total_time     = 0
    	)
    BSFG_state$current_state = current_state

# ------------------------------------ #
# ----Precalculate some matrices------ #
# ------------------------------------ #

    # recover()
    #invert the random effect covariance matrices
    K = K_mats[[1]]
    chol_Kinv = chol_Ki_mats[[1]]

    #pre-calculate transformation parameters to diagonalize aI + bZKZ for fast
    #inversion: inv(aI + bZKZ) = 1/b*U*diag(1./(s+a/b))*U'
    #uses singular value decomposition of ZKZ for stability when ZKZ is low
    #rank
#     XZ = [X_f Z  ]
#     [U,S,~]          = svd(XZ*blkdiag(1e6*eye(b_f),K)*XZ')

    result = svd(Z %*% K %*% t(Z))
    invert_aI_bZKZ = list(
        U = Matrix(result$u,sparse=T),
        s = result$d
    )

    #genetic effect variances of factor traits
    # diagonalizing a*Z  '*Z   + b*Kinv for fast inversion
    #diagonalize mixed model equations for fast inversion:
    # inv(a*Z  '*Z   + b*Kinv) = U*diag(1./(a.*s1+b.*s2))*U'
    #similar to fixed effects + random effects 1 above, but no fixed effects.

    ZZt = crossprod(Z)
    svd_ZZt = svd(ZZt)
    ZZt_sqrt = t(sweep(svd_ZZt$u,2,sqrt(svd_ZZt$d),'*'))
    result = GSVD_2_c(ZZt_sqrt,as.matrix(chol_Kinv))

    invert_aZZt_Kinv = list(
        U = drop0(Matrix(t(solve(result$X)),sparse=T),tol = 1e-14),
        # U = t(solve(result$X)),
			s1 = diag(result$C)^2,
			s2 = diag(result$S)^2
		)


# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

    run_variables = c(run_variables,list(
      chol_Kinv        = chol_Kinv,
			invert_aI_bZKZ   = invert_aI_bZKZ,
			invert_aZZt_Kinv = invert_aZZt_Kinv
    ))

    RNG = list(
    	Random.seed = .Random.seed,
    	RNGkind = RNGkind()
    )

    BSFG_state$run_variables = run_variables
    BSFG_state$RNG = RNG

    return(BSFG_state)
}
