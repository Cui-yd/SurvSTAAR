exactScore = function(objNull = objNull, G_mat = G_mat, use_SPA = FALSE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  t_P_mat = objNull$t_P_mat
  W_mat = objNull$W_mat
  X_mat = objNull$X_mat
  PX_mat = objNull$PX_mat
  XVX_inv = objNull$XVX_inv
  res_null = objNull$res_null
  t_X_mat = objNull$t_X_mat
  X_XXinv = objNull$X_XXinv


  PG_mat = t_P_mat %*% G_mat
  GWG_mat = t(W_mat %*% G_mat) %*% G_mat
  GWX_mat = t(W_mat %*% G_mat) %*% X_mat
  GVG_mat = GWG_mat - t(PG_mat) %*% PG_mat
  GVX_mat = GWX_mat - t(PG_mat) %*% PX_mat
  Var_mat = GVG_mat - GVX_mat %*% XVX_inv %*% t(GVX_mat)     ## GQG


  Score = t(res_null) %*% G_mat
  Score = as.vector(Score)
  Variance = diag(Var_mat)
  Stest = as.vector(Score^2 / Variance)
  p_value = as.vector(pchisq(Stest, df = 1, lower.tail = F))

  result = data.frame("Score" = Score, "Variance" = Variance, "Stest" = Stest, "Pvalue" = p_value)


  if (SPA_filter) {
    SPA_index = which(p_value < SPA_filter_cutoff)   ## SPA adjust for those p-value smaller than SPA_filter_cutoff
  } else {
    SPA_index = 1:length(p_value)
  }


  if(use_SPA) {

    result$Pvalue_SPA = result$Pvalue


    if (length(SPA_index) > 0) {

      if (verbose) print(paste0("Single variants analysis: apply SPA adjustment to ", length(SPA_index), " markers"))

      G_SPA = G_mat[, SPA_index]
      G_SPA_Score = Score[SPA_index]
      G_SPA_Var = Variance[SPA_index]
      G_tilde = G_tilde_forSPA(G_SPA, objNull = objNull)


      for(i in 1:length(SPA_index)){

        G_i = G_tilde[, i]
        G_i_Score = G_SPA_Score[i]
        G_i_Var = G_SPA_Var[i]

        pval_i = single_SPA(G_i, G_i_Score, G_i_Var, objNull)

        result$Pvalue_SPA[SPA_index[i]] = pval_i
      }

    } else {

      if (verbose) print(paste0("Under SPA filter, no markers' p-value is smaller than ", SPA_filter_cutoff))

    }

  }

  Uscore_se = sqrt(result$Variance)
  Est = result$Score / result$Variance
  Est_se = 1 / Uscore_se

  if (use_SPA) {

    pvalue_log10 = -log10(result$Pvalue_SPA)
    result$pvalue_log10 = pvalue_log10
    result = result[, c(4:6, 1:3)]
    result$Est = Est
    result$Est_se = Est_se

  } else {

    pvalue_log10 = -log10(result$Pvalue)
    result$pvalue_log10 = pvalue_log10
    result = result[, c(3:5, 1:3)]
    result$Est = Est
    result$Est_se = Est_se

  }

  result = list("result" = result, "Score" = Score, "Covariance" = Var_mat)

  return(result)
}


G_tilde_forSPA = function(G, objNull) {

  t_X_mat = objNull$t_X_mat
  X_XXinv = objNull$X_XXinv

  if (inherits(G, "sparseMatrix")) {G = suppressWarnings(as(G, "matrix"))}
  G_tilde = G - X_XXinv %*% (t(t_X_mat) %*% G)

  return(G_tilde = G_tilde)
}



genoMatrixPlink = function(Geno, bim_data = bim_data, markerIndex,
                           geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                           min_mac_cutoff = NULL, min_maf_cutoff = 0.01) {

  Geno[Geno == -9 | Geno == 9] = NA
  colnames(Geno) = bim_data$V2[markerIndex]
  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }

  G_summary = data.table::rbindlist(G_summary)


  ## variants filter
  if (is.null(min_mac_cutoff)) {
    include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
  } else {
    maf_mac_cutoff = Geno_row * min_maf_cutoff * 2
    if (maf_mac_cutoff < min_mac_cutoff) {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
    } else {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAC >= min_mac_cutoff)
    }
  }

  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  G_summary = cbind(bim_data[markerIndex[include_index], c(1,4:6,2)], G_summary)
  colnames(G_summary)[1:5] = c("CHR", "POS", "REF1", "REF2", "rsID")

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary)

  return(result)
}


genoMatrixGDS = function(Geno, G_info, geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                         min_mac_cutoff = NULL, min_maf_cutoff = 0.01) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    AF = mean(Geno_p, na.rm = T)/2
    ALT_AF = 1-AF
    MAF = AF
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }

  G_summary = data.table::rbindlist(G_summary)
  G_summary = cbind(G_info, G_summary)

  ## variants filter
  if (is.null(min_mac_cutoff)) {
    include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
  } else {
    maf_mac_cutoff = Geno_row * min_maf_cutoff * 2
    if (maf_mac_cutoff < min_mac_cutoff) {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff)
    } else {
      include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAC >= min_mac_cutoff)
    }
  }

  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary, include_index = include_index)

  return(result)
}



genoFlipRV = function(Geno, geno_missing_imputation = c("mean", "minor"),
                      geno_missing_cutoff = 0.1, min_maf_cutoff = 1e-4, rare_maf_cutoff = 0.01, rare_num_cutoff = 2) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF, MAC, missing, and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC, missing_rate)
    G_na[[p]] = missing_index
  }
  G_summary = data.table::rbindlist(G_summary)


  ## variants filter
  include_index = which(G_summary$missing_rate <= geno_missing_cutoff & G_summary$MAF >= min_maf_cutoff & G_summary$MAF <= rare_maf_cutoff)
  Geno = Geno[, include_index, drop = FALSE]
  G_summary = G_summary[include_index, ]
  G_na = G_na[include_index]


  ## check for rare variants number
  if(length(include_index) < rare_num_cutoff) {
    stop('Number of variants in this gene is less than ', rare_num_cutoff, ", will skip this category...", call. = FALSE)
  }


  ## imputation
  switch (geno_missing_imputation,
          "mean" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              imput_p = G_summary$MAF[p]*2
              Geno[na_index, p] = imput_p
            }
          },
          "minor" = {
            for (p in 1:length(include_index)) {
              na_index = G_na[[p]]
              Geno[na_index, p] = 0
            }
          }
  )

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary, include_index = include_index)

  return(result)

}



genoFlip = function(Geno) {

  Geno_col = ncol(Geno)
  Geno_row = nrow(Geno)

  G_summary = list()
  G_na = list()

  ## calculate MAF MAC missing and flip
  for (p in 1:Geno_col) {
    Geno_p = Geno[, p]

    MAF = mean(Geno_p, na.rm = T)/2
    if (MAF > 0.5) {
      MAF = 1-MAF
      Geno_p = 2-Geno_p
      Geno[, p] = Geno_p
    }
    MAC = sum(Geno_p, na.rm = T)
    missing_index = which(is.na(Geno_p))
    missing_rate = length(missing_index)/Geno_row
    G_summary[[p]] = data.frame(MAF, MAC)
    G_na[[p]] = missing_index
  }
  G_summary = data.table::rbindlist(G_summary)

  ## imputation  minor
  for (p in 1:ncol(Geno)) {
    na_index = G_na[[p]]
    Geno[na_index, p] = 0
  }

  ## update imputed MAC
  G_summary$MAC = colSums(Geno)

  result = list(Geno = Geno, G_summary = G_summary)

  return(result)

}


#' @export
rap_load_as = function(file, name) {
  # file <- basename(file)
  data_list <- load(file)
  assign(name, get(data_list[1]), envir = .GlobalEnv)
}


#' @export
argsReshape = function(default_args, args, num_args, log_args) {

  args_list = default_args

  for (arg in args) {

    key_value <- strsplit(arg, "=")[[1]]

    if (length(key_value) == 2) {
      arg_name = sub("--", "", key_value[1])

      if (arg_name %in% names(default_args)) {

        if (arg_name %in% num_args) {

          value <- suppressWarnings(as.numeric(key_value[2]))
          if (is.numeric(value)) {
            args_list[[arg_name]] <- value
          } else {
            warning(paste("Invalid value: ", arg_name, ". A numeric value is required."))
          }

        } else if (arg_name %in% log_args) {

          value <- suppressWarnings(as.logical(key_value[2]))
          if (is.logical(value)) {
            args_list[[arg_name]] <- value
          } else {
            warning(paste("Invalid value: ", arg_name, ". A logical value is required."))
          }

        } else {
          value <- key_value[2]
          args_list[[arg_name]] <- value
        }

      } else {
        warning(paste("Unknown argument: ", arg_name))
      }
    }
  }
  return(args_list)
}


single_SPA = function(G, Score, Variance, objNull) {

  # X_XXinv = objNull$X_XXinv
  # t_X_mat = objNull$t_X_mat

  # G_tilde = G - X_XXinv %*% (t(t_X_mat) %*% G)

  Stest = Score / sqrt(Variance)
  # G_norm = G_tilde/sqrt(Variance)
  G_norm = G/sqrt(Variance)

  N = length(G)
  N1set = 1:N
  N0 = 0
  G2N1 = G_norm
  G2N0 = 0

  pval1 = GetProb_SPA(objNull = objNull, G2N1 = G2N1, G2N0 = G2N0, N1set = N1set, N0 = N0, ztest = abs(Stest), lower.tail = FALSE)
  pval2 = GetProb_SPA(objNull = objNull, G2N1 = G2N1, G2N0 = G2N0, N1set = N1set, N0 = N0, ztest = -abs(Stest), lower.tail = TRUE)
  pval = pval1 + pval2

  return(pval)
}



### Following functions are sourced from the GitHub repository SPACox (https://github.com/WenjianBI/SPACox)

CGF4MartingaleRes = function(fit_null = fit_null, range = c(-100,100), length.out = 1e4, verbose = TRUE) {

  fit_null$use_SPA = TRUE

  ### calculate empirical CGF for martingale residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  res_null = fit_null$res_null

  cumul = NULL
  if(verbose) print("Start calculating empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    exp_res = exp(res_null*t)
    M0 = mean(exp_res)
    M1 = mean(res_null*exp_res)
    M2 = mean(res_null^2*exp_res)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(verbose){if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))}
  }

  K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)


  fit_null = c(fit_null, K_org_emp = K_org_emp, K_1_emp = K_1_emp, K_2_emp = K_2_emp)

  return(fit_null)

}

GetProb_SPA = function(objNull = objNull, G2N1 = G2N1, G2N0 = G2N0, N1set = N1set, N0 = N0, ztest = ztest, lower.tail) {

  out = uniroot(K1_adj, c(-20,20), extendInt = "upX",
                G2N1=G2N1, G2N0=G2N0, N1set=N1set, N0=N0, ztest=ztest, objNull=objNull)
  zeta = out$root

  k1 = K_org(zeta, G2N1=G2N1, G2N0=G2N0, N1set=N1set, N0=N0, objNull=objNull)
  k2 = K2(zeta,  G2N1=G2N1, G2N0=G2N0, N1set=N1set, N0=N0, objNull=objNull)

  temp1 = zeta * ztest - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)

  return(pval)
}

K_org = function(t, G2N1, G2N0, N1set, N0, objNull) {

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2N0 = t1*G2N0
    t2N1 = t1*G2N1
    out[i] = N0*objNull$K_org_emp(t2N0) + sum(objNull$K_org_emp(t2N1))
  }
  return(out)
}

K1_adj = function(t, G2N1, G2N0, N1set, N0, ztest, objNull) {
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2N0 = t1*G2N0
    t2N1 = t1*G2N1
    out[i] = N0*G2N0*objNull$K_1_emp(t2N0) + sum(G2N1*objNull$K_1_emp(t2N1)) - ztest
  }
  return(out)
}

K2 = function(t, G2N1, G2N0, N1set, N0, objNull) {
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2N0 = t1*G2N0
    t2N1 = t1*G2N1
    out[i] = N0*G2N0^2*objNull$K_2_emp(t2N0) + sum(G2N1^2*objNull$K_2_emp(t2N1))
  }
  return(out)
}




### Following function is sourced from the GitHub repository STAAR (https://github.com/xihaoli/STAAR/blob/master/R/CCT.R)

CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}




Burden = function(Geno, Score, Covariance, weight, objNull, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("genotype is not a matrix!")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (length(Score) != ncol(Geno) | ncol(Covariance) != ncol(Geno)) stop("must inpute corresponding score and covariance")

  if (nrow(weight) != ncol(Geno)) stop("dimentions don't match for genotype and weight")

  pval_B = NULL
  G_tilde = G_tilde_forSPA(Geno, objNull)
  G_tilde_w = G_tilde %*% weight

  for (k in 1:ncol(weight)) {

    weight_k = weight[, k]

    Score_k = sum(Score * weight_k)
    Variance_k = as.vector(weight_k %*% Covariance %*% weight_k)

    Stest_k = Score_k^2 / Variance_k
    pval_k = pchisq(Stest_k, df = 1, lower.tail = FALSE)

    if (use_SPA && (!SPA_filter || (SPA_filter && pval_k < SPA_filter_cutoff))) {

      G_w_k = G_tilde_w[,k,drop = F]
      pval_k = single_SPA(G_w_k, Score_k, Variance_k, objNull)

    }

    pval_B = c(pval_B, pval_k)
  }

  pval_B = matrix(pval_B, nrow = 2, byrow = T)

  return(pval_B)

}



ACAT = function(Geno, Score, Covariance, Pvalue, MAC = NULL, weight_A, weight_B,
                objNull, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 60, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("genotype is not a matrix!")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (length(Score) != ncol(Geno) | ncol(Covariance) != ncol(Geno)) stop("must inpute corresponding score and covariance")

  if (nrow(weight_A) != ncol(Geno)) stop("dimentions don't match for genotype and weight")

  if (ncol(weight_A) != ncol(weight_B)) stop("dimentions don't match for two weight matrix")

  if (is.null(MAC)) MAC = colSums(Geno)



  ultra_rare_index = which(MAC <= ultra_rare_mac_cutoff)

  pval_A = NULL

  if(!combine_ultra_rare | length(ultra_rare_index) <= 1) {

    ### Cauchy ###
    if (verbose) print("             Apply Cauchy combination in ACAT test")
    weight = weight_A

    for (k in 1:ncol(weight)) {

      weight_k = weight[, k]
      pval_k = CCT(Pvalue, weight_k)
      pval_A = c(pval_A, pval_k)

    }

  } else if (combine_ultra_rare & length(ultra_rare_index) == ncol(Geno)) {

    ### Burden ###
    if (verbose) print("             All the markers in this gene are ultra rare, apply burden test to replace ACAT test")

    weight = weight_B
    pval_A = Burden(Geno, Score, Covariance, weight, objNull, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = FALSE)

  } else {

    ### Cauchy + Burden ###
    if (verbose) print(paste0("             Combine ", length(ultra_rare_index), " ultra rare variants in ACAT test"))

    ## Burden test
    weight_B_rare = matrix(weight_B[ultra_rare_index, ], nrow = length(ultra_rare_index))
    weight = rbind(colMeans(weight_A[ultra_rare_index, ]), weight_A[-ultra_rare_index, ])

    G_ultra_rare = as.matrix(Geno[, ultra_rare_index])
    Pvalue_sub = Pvalue[-ultra_rare_index]
    Score_rare = Score[ultra_rare_index]
    Covariance_rare = as.matrix(Covariance[ultra_rare_index, ultra_rare_index])

    Pval_rare_B = Burden(G_ultra_rare, Score = Score_rare, Covariance = Covariance_rare, weight = weight_B_rare,
                         objNull, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = FALSE)


    for (k in 1:ncol(weight)) {

      ## ACAT test
      Pvalue_k = c(Pval_rare_B[k], Pvalue_sub)
      weight_k = weight[, k]

      pval_k = CCT(Pvalue_k, weight_k)

      pval_A = c(pval_A, pval_k)

    }

  }

  pval_A = matrix(pval_A, nrow = 2, byrow = T)

  return(pval_A)

}



SKAT = function(Geno, Score, Covariance, Pvalue, MAC = NULL, weight_S, weight_B,
                objNull, use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 60, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("genotype is not a matrix!")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (length(Score) != ncol(Geno) | ncol(Covariance) != ncol(Geno)) stop("must inpute corresponding score and covariance")

  if (nrow(weight_S) != ncol(Geno)) stop("dimentions don't match for genotype and weight")

  if (ncol(weight_S) != ncol(weight_B)) stop("dimentions don't match for two weight matrix")

  if (is.null(MAC)) MAC = colSums(Geno)

  pval_S = NULL

  ultra_rare_index = which(MAC <= ultra_rare_mac_cutoff)


  if (use_SPA) {

    ### use SPA ###

    if (combine_ultra_rare & length(ultra_rare_index) > 1) {

      ## combine ultra-rare variants before SKAT test
      if (verbose) print(paste0("             Combine ", length(ultra_rare_index), " ultra rare variants in SKAT test"))


      if (length(ultra_rare_index) == ncol(Geno)) {

        ## if all the variants are ultra-rare, use Burden test ##

        pval_S = Burden(Geno = Geno, Score = Score, Covariance = Covariance, weight = weight_B, objNull = objNull,
                        use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = FALSE)

      } else {

        ## combine the ultra rare variants subset first ##

        weight_B_rare = weight_B[ultra_rare_index, ]
        for(i in 1:ncol(weight_B_rare)) {
          weight_B_rare[,i] = weight_B_rare[,i] / sqrt(sum(weight_B_rare[,i]^2))
        }

        Geno_rare = Geno[, ultra_rare_index]
        Score_rare = Score[ultra_rare_index]
        Covariance_rare = Covariance[ultra_rare_index, ultra_rare_index]

        pval_B_rare = Burden(Geno = Geno_rare, Score = Score_rare, Covariance = Covariance_rare, weight = weight_B_rare, objNull = objNull,
                             use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = FALSE)

        Pvalue_sub = Pvalue[-ultra_rare_index]
        weight_S_new = rbind(colMeans(weight_S[ultra_rare_index, ]), weight_S[-ultra_rare_index, ])

        Covariance_sub = Covariance[-ultra_rare_index, -ultra_rare_index]
        Covariance_rare_cov = Covariance[-ultra_rare_index, ultra_rare_index]

        for (k in 1:ncol(weight_S)) {

          weight_B_rare_k = weight_B_rare[, k]
          Geno_rare_wk = Geno_rare %*% weight_B_rare_k
          Geno_new_wk = cbind(Geno_rare_wk, Geno[, -ultra_rare_index])

          Covariance_rare_Bk = as.vector(weight_B_rare_k %*% Covariance_rare %*% weight_B_rare_k)
          Covariance_rare_cov_Bk = as.vector(Covariance_rare_cov %*% weight_B_rare_k)

          Covariance_wk = Matrix(NA, nrow = ncol(Geno_new_wk), ncol = ncol(Geno_new_wk))
          Covariance_wk[1,1] = Covariance_rare_Bk
          Covariance_wk[2:ncol(Geno_new_wk), 1] = Covariance_rare_cov_Bk
          Covariance_wk[1, 2:ncol(Geno_new_wk)] = Covariance_rare_cov_Bk
          Covariance_wk[2:ncol(Geno_new_wk), 2:ncol(Geno_new_wk)] = Covariance_sub


          Pvalue_wk = c(pval_B_rare[k], Pvalue_sub)
          z_tilde = qnorm(Pvalue_wk/2, mean = 0, sd = 1, lower.tail = F, log.p = F)
          weight_k = weight_S_new[, k]
          Qtest = sum(z_tilde^2 * weight_k^2 * diag(Covariance_wk))

          w_k_mat = matrix(rep(weight_k, ncol(Geno_new_wk)), nrow = ncol(Geno_new_wk), byrow = T)
          Cov_weight = t(w_k_mat) * Covariance_wk * w_k_mat

          if (sum(Cov_weight) == 0) {

            pval_k = NA

          } else {

            lambda = eigen(Cov_weight, only.values = T, symmetric = T)$values
            pval_k = suppressWarnings(CompQuadForm::davies(Qtest, lambda)$Qq)   # p[Q>q]

            if (length(pval_k) == 0 | pval_k > 1 | pval_k <= 0) pval_k = CompQuadForm::liu(Qtest, lambda)

          }

          pval_S = c(pval_S, pval_k)

        }
      }


    } else {

      ## all the variants are not ultra-rare ##

      z_tilde = qnorm(Pvalue/2, mean = 0, sd = 1, lower.tail = F, log.p = F)

      for (k in 1:ncol(weight_S)) {

        weight_k = weight_S[, k]
        Qtest = sum(z_tilde^2 * weight_k^2 * diag(Covariance))   ## equal to sum(z_tilde^2 %*% diag(Cov_weight))

        w_k_mat = matrix(rep(weight_k, ncol(Geno)), nrow = ncol(Geno), byrow = T)
        Cov_weight = t(w_k_mat) * Covariance * w_k_mat

        if (sum(Cov_weight) == 0) {

          pval_k = NA

        } else {

          lambda = eigen(Cov_weight, only.values = T, symmetric = T)$values
          pval_k = suppressWarnings(CompQuadForm::davies(Qtest, lambda)$Qq)   # p[Q>q]

          if (length(pval_k) == 0 | pval_k > 1 | pval_k <= 0) pval_k = CompQuadForm::liu(Qtest, lambda)

        }

        pval_S = c(pval_S, pval_k)

      }

    }

  } else {

    ## not use SPA ##

    for (k in 1:ncol(weight_S)) {

      weight_k = weight_S[, k]
      Qtest = sum(Score^2 * weight_k^2)   ## equal to sum(Stest %*% diag(Cov_weight))

      w_k_mat = matrix(rep(weight_k, ncol(Geno)), nrow = ncol(Geno), byrow = T)
      Cov_weight = t(w_k_mat) * Covariance * w_k_mat

      if (sum(Cov_weight) == 0) {
        pval_k = NA
      } else {

        lambda = eigen(Cov_weight, only.values = T, symmetric = T)$values
        pval_k = suppressWarnings(CompQuadForm::davies(Qtest, lambda)$Qq)   # p[Q>q]

        if (length(pval_k) == 0 | pval_k > 1 | pval_k <= 0) pval_k = CompQuadForm::liu(Qtest, lambda)

      }

      pval_S = c(pval_S, pval_k)

    }

  }

  pval_S = matrix(pval_S, nrow = 2, byrow = T)

  return(pval_S)

}



SurvSTAAR_O = function(Geno, objNull, annotation_rank = NULL, MAC = NULL,
                       use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                       weight_A, weight_B, weight_S,
                       combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("Genotype is not a matrix")

  if (ncol(Geno) != nrow(weight_A) | ncol(Geno) != nrow(weight_B) | ncol(Geno) != nrow(weight_S)) stop("Dimensions don't match for genotype and weights")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (is.null(MAC)) MAC = colSums(Geno)


  ### individual test
  single_test = exactScore(objNull = objNull, G_mat = Geno, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = verbose)
  Score = single_test$Score
  Covariance = single_test$Covariance
  if(use_SPA) {
    Pvalue = single_test$result$Pvalue_SPA
  } else {
    Pvalue = single_test$result$Pvalue
  }


  ### set based test
  if (verbose) print(paste0("SKAT test:   begin at ", Sys.time()))
  Pvalue_S = SKAT(Geno, Score, Covariance, Pvalue, MAC, weight_S, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("ACAT test:   begin at ", Sys.time()))
  Pvalue_A = ACAT(Geno, Score, Covariance, Pvalue, MAC, weight_A, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("Burden test: begin at ", Sys.time()))
  Pvalue_B = Burden(Geno, Score, Covariance, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, verbose)


  ### combine results
  results_pvalue = rbind(Pvalue_S, Pvalue_A, Pvalue_B)
  rownames(results_pvalue) = c("SKAT-(1,25)", "SKAT-(1,1)", "ACAT-(1,25)", "ACAT-(1,1)", "Burden-(1,25)", "Burden-(1,1)")

  if (!is.null(annotation_rank)) {
    colnames(results_pvalue) = c("Beta", colnames(annotation_rank))
  } else {
    colnames(results_pvalue) = "results"
  }


  results_STAAR_O = CCT(c(results_pvalue))
  results_ACAT_O = CCT(c(results_pvalue[, 1]))
  results_STAAR_S = CCT(results_pvalue[1:2, ])
  results_STAAR_A = CCT(results_pvalue[3:4, ])
  results_STAAR_B = CCT(results_pvalue[5:6, ])

  STAAR_S_1_25 = CCT(Pvalue_S[1, ])
  STAAR_S_1_1  = CCT(Pvalue_S[2, ])
  STAAR_A_1_25 = CCT(Pvalue_A[1, ])
  STAAR_A_1_1  = CCT(Pvalue_A[2, ])
  STAAR_B_1_25 = CCT(Pvalue_B[1, ])
  STAAR_B_1_1  = CCT(Pvalue_B[2, ])

  results_SurvSTAAR = c(STAAR_S_1_25, STAAR_S_1_1, STAAR_A_1_25, STAAR_A_1_1, STAAR_B_1_25, STAAR_B_1_1)
  results_pvalue = cbind(results_pvalue, results_SurvSTAAR)


  results = list("SurvSTAAR_O" = results_STAAR_O,
                 "ACAT_O" = results_ACAT_O,
                 "SurvSTAAR_SKAT" = results_STAAR_S,
                 "SurvSTAAR_ACAT" = results_STAAR_A,
                 "SurvSTAAR_Burden" = results_STAAR_B,
                 "SurvSTAAR_pvalue" = results_pvalue)


  return(results)
}



#' @export
SurvSTAAR_Manhattan = function(result_data, chr = "Chr", pos = "pos", name = "Gene",
                               pval = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous"),
                               pch = c(0, 1, 2, 3, 4), max_y,
                               sugline = 5e-7, genoline = 5e-7, annotateP = NULL,
                               col=c("gray30", "gray80"), legned_title = "Functional Categories") {

  result_data = result_data[order(result_data[, chr], result_data[, pos]), ]

  if (length(pval) == 1) {
    result_data[, pval] = -log10(result_data[, pval])
  }
  if (length(pval) > 1) {
    for (j in pval) {
      result_data[, j] = -log10(result_data[, j])
    }
  }

  result_data$index = rep.int(seq_along(unique(result_data[, chr])), times = tapply(result_data[, name],result_data[, chr],length))

  lastbase = 0
  ticks = NULL
  result_data$pos_new = NULL
  for (i in unique(result_data$index)) {
    if (i==1) {
      result_data[result_data$index == i, "pos_new"] = result_data[result_data$index == i, pos]
    } else {
      lastbase = lastbase + max(result_data[result_data$index == (i-1), pos])
      result_data[result_data$index == i,pos] = result_data[result_data$index == i,pos]-min(result_data[result_data$index==i,pos]) +1
      result_data[result_data$index == i,"pos_new"] = result_data[result_data$index == i,pos] + lastbase
    }
  }
  ticks <-tapply(result_data$pos_new, result_data$index, quantile, probs=0.5)
  xlabel = 'Chromosome'
  labs <- unique(result_data[, chr])

  xmax = ceiling(max(result_data$pos_new) * 1.01)
  xmin = floor(max(result_data$pos_new) * -0.01)
  col = rep_len(col, max(result_data$index))


  ### begin plotting
  if (is.null(max_y)) {
    def_args <- list(xaxt='n', bty='o', xaxs='i', yaxs='i', las=1, pch=20,
                     cex.axis=1.5, cex.lab=1.5, cex=1.5,
                     xlim=c(xmin,xmax), ylim=c(0, ceiling(max(result_data[, pval], na.rm = T))),
                     xlab=xlabel, ylab=expression(-log[10](italic(p))))
  } else {
    def_args <- list(xaxt='n', bty='o', xaxs='i', yaxs='i', las=1, pch=20,
                     cex.axis=1.5, cex.lab=1.5, cex=1.5,
                     xlim=c(xmin,xmax), ylim=c(0,max_y),
                     xlab=xlabel, ylab=expression(-log[10](italic(p))))
  }
  par(mar = c(5, 4, 6, 2) + 1.5)
  do.call("plot", c(NA, def_args))
  axis(1, at=ticks, labels=labs, tck = 0, cex.axis = 1.5)

  if (length(pval) == 1) {

    icol=1
    for (i in unique(result_data$index)) {
      points(result_data[result_data$index==i,"pos_new"], result_data[result_data$index==i,pval], col=col[icol], pch=20,)
      icol=icol+1
    }

  } else if (length(pval) > 1) {

    for (j in 1:length(pval)) {
      pch_j = pch[j]
      icol=1
      for (i in unique(result_data$index)) {
        points(result_data[result_data$index==i,"pos_new"], result_data[result_data$index==i,pval[j]], col=col[icol], pch=pch_j)
        icol=icol+1
      }
    }

  }


  sugline = -log10(sugline)
  genoline = -log10(genoline)
  # Add suggestive and genomewide lines
  usr <- par("usr")
  clip(usr[1], usr[2], usr[3], usr[4])
  if (sugline) abline(h=(sugline), col="blue", lwd = 2)
  if (genoline) abline(h=genoline, col="red", lwd = 2)


  # Highlight top SNPs

  if (!is.null(annotateP)) {
    annotateP = -log10(annotateP)

    if (length(pval) == 1) {

      topHits = result_data[which(result_data[, pval] >= annotateP), ]
      par(xpd = TRUE)
      if (nrow(topHits) != 0) {
        with(topHits, textxy(pos_new, topHits[, pval], offset = 0.625,
                             labs = topHits[, name], cex = 1.2))
      }

    } else if (length(pval) > 1) {

      for (j in pval) {
        topHits = result_data[which(result_data[, j] >= annotateP), ]
        par(xpd = TRUE)
        if (nrow(topHits) != 0) {
          with(topHits, textxy(pos_new, topHits[, j], offset = 0.625,
                               labs = topHits[, name], cex = 1.2))
        }
      }
    }
  }


  ### add legend
  if (length(pval) <= 5) {
    legend("top", legend = pval, pch = pch, title = legned_title, bty = "n",
           text.font = 1, cex = 1.4, ncol = length(pval), xpd = TRUE, inset = -0.2)
  } else if (length(pval) > 5) {
    temp_len = floor(length(pval)/2)
    legend("top", legend = pval[1:temp_len], pch = pch[1:temp_len], title = legned_title, bty = "n",
           text.font = 1, cex = 1.4, ncol = temp_len, xpd = TRUE, inset = -0.2)
    legend("top", legend = pval[(temp_len+1):length(pval)], pch = pch[(temp_len+1):length(pval)], bty = "n",
           text.font = 1, cex = 1.4, ncol = (length(pval)-temp_len), xpd = TRUE, inset = -0.1)
  }
  par(xpd = FALSE)
}





#' @export
SurvSTAAR_QQplot = function(pval_result, main = NULL, pch = NULL, max_x = NULL, max_y = NULL,
                            legend_lable = c("pLoF", "pLoF+D", "Missense", "Disruptive Missense", "Synonymous"),
                            legned_title = "Functional Categories") {

  if (inherits(pval_result, "numeric")) {

    pval_result = sort(pval_result)
    if (is.null(max_y)) {
      max_y = ceiling(-log10(min(pval_result, na.rm = T)))
    }
    if (is.null(main)) {
      if (is.null(max_x)) {
        plot(-log10(ppoints(pval_result)), -log10(pval_result), pch = 20, cex = 0.5, col="grey40",
             xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
             ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3)
      } else {
        plot(-log10(ppoints(pval_result)), -log10(pval_result), pch = 20, cex = 0.5, col="grey40",
             xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
             ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3, xlim = c(0, max_x))
      }
    } else {
      if (is.null(max_x)) {
        plot(-log10(ppoints(pval_result)), -log10(pval_result), pch = 20, cex = 0.5, col="grey40",
             xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
             ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3, main = main)
      } else {
        plot(-log10(ppoints(pval_result)), -log10(pval_result), pch = 20, cex = 0.5, col="grey40",
             xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
             ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3, xlim = c(0, max_x), main = main)
      }
    }
    abline(a = 0, b = 1, lty = "longdash", lwd = 1.5)

  } else if (inherits(pval_result, c("data.frame", "matrix"))) {

    if (is.null(max_y)) {
      max_y = ceiling(-log10(min(pval_result, na.rm = T)))
    }
    par(mar = c(5, 4, 2, 2) + 1)
    if (is.null(pch)) pch = c(0:(ncol(pval_result)-1))

    for (j in 1:ncol(pval_result)) {

      pval_result_temp = pval_result[, j]
      pval_result_temp = sort(pval_result_temp)

      if (j == 1) {
        if (is.null(main)) {
          if (is.null(max_x)) {
            plot(-log10(ppoints(pval_result_temp)), -log10(pval_result_temp), pch = pch[j], cex = 0.5, col="grey40",
                 xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
                 ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3)
          } else {
            plot(-log10(ppoints(pval_result_temp)), -log10(pval_result_temp), pch = pch[j], cex = 0.5, col="grey40",
                 xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
                 ylim = c(0, max_y), xlim = c(0, max_x), cex.lab = 1.3, cex.axis = 1.3)
          }
        } else {
          if (is.null(max_x)) {
            plot(-log10(ppoints(pval_result_temp)), -log10(pval_result_temp), pch = pch[j], cex = 0.5, col="grey40",
                 xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
                 ylim = c(0, max_y), cex.lab = 1.3, cex.axis = 1.3, main = main)
          } else {
            plot(-log10(ppoints(pval_result_temp)), -log10(pval_result_temp), pch = pch[j], cex = 0.5, col="grey40",
                 xlab = expression(Expected~~-log[10](p)), ylab = expression(Observed~~-log[10](p)),
                 ylim = c(0, max_y), xlim = c(0, max_x), cex.lab = 1.3, cex.axis = 1.3, main = main)
          }
        }
      } else {
        points(-log10(ppoints(pval_result_temp)), -log10(pval_result_temp), pch = pch[j], cex = 0.5, col="grey40")
      }

    }
    abline(a = 0, b = 1, lty = "longdash", lwd = 1.5)
    if (is.null(legned_title) | isFALSE(legned_title)){
      legend(x = "topleft", legend = legend_lable, pch = pch, # bty = "n",
             text.font = 2, cex = 1.3, ncol = 1, xpd = TRUE)
    } else {
      legend(x = "topleft", legend = legend_lable, pch = pch, title = legned_title, # bty = "n",
             text.font = 2, cex = 1.3, ncol = 1, xpd = TRUE)
    }

  }
}


