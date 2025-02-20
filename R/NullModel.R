#' Fit Cox proportional hazards model under the null hypothesis for related samples.
#'
#' This function fit Cox proportional hazards null model using the \code{\link{coxph}}
#' function from the \code{\link{survival}} package. Tied event times are handled
#' using Breslow’s approximation. Additionally, the function calculates the empirical
#' cumulative generating function (CGF) of the martingale residuals using the
#' \code{CGF4MartingaleRes} function from the \code{SPACox} package.
#'
#' @param genofile a genome file, which can be NULL. If not NULL, it can either
#' be a character string specifying the file name in PLINK format or a SeqVarGDSClass
#' object for aGDS format data. When this parameter is provided (non-NULL), the
#' function will match sample IDs with genome IDs; if it is NULL, no matching
#' will be performed.
#' @param phenofile a character string specifying the file name of the phenotype
#' data, which includes all variables used in the null model.
#' @param LOCO a logical value indicating whether to perform leave-one-chromosome-out (LOCO)
#' in null model fitting (default = TRUE).
#' @param chr a numeric value representing the chromosome. If `LOCO = FALSE`, this
#' can be NULL. If `LOCO = TRUE`, `chr` must be specified as a valid chromosome number.
#' @param statusCol the column name of survival status.
#' @param timeCol the column name of survival time.
#' @param sampleCol the column name of sample IDs.
#' @param covCol a vector of column names of covariates.
#' @param PRSCol the column name of polygenic risk scores (PRS). If `LOCO = FALSE`,
#' it is recommended to use "chr23".
#' @param use_SPA a logical value indicating whether to use SPA (SaddlePoint
#' Approximation) in the subsequent analysis, it is highly recommended under
#' imbalanced phenotypes (default = TRUE).
#' @param range a two-element numeric vector to specify the domain of the
#' empirical CGF (default = c(-100,100)).
#' @param length.out a positive integer for empirical CGF. Larger length.out
#' corresponds to longer calculation time and more accurate estimated empirical
#' CGF (default = 1e4).
#' @param verbose a logical value indicating whether to provides additional detailed
#' information or messages during the execution of this function (default = FALSE).
#'
#' @returns A list containing the model fit from \code{\link{coxph}}.
#' If \code{use_SPA = TRUE}, the list will also include the empirical CGF.
#' @export

NullModel = function(genofile = NULL, phenofile, LOCO = TRUE, chr = NULL,
                     statusCol, timeCol, sampleCol, covCol = NULL, PRSCol = NULL,
                     use_SPA = FALSE, range=c(-100,100), length.out = 1e4, verbose = FALSE) {

  if (!is.null(phenofile)) {
    if (is.character(phenofile)) {
      if (!file.exists(phenofile)) stop("The phenotype file is not exist!")
      use_data = fread(phenofile, data.table = FALSE)
      use_data = na.omit(use_data)
    } else {
      use_data = phenofile
      use_data = na.omit(use_data)
    }
  } else {
    stop("You should input the phenotype data!")
  }

  ### Check files if the genotype is input
  if (!is.null(genofile)) {
    if (inherits(genofile, "SeqVarGDSClass")) {   ## input is an agds file

      sample.geno <- seqGetData(genofile, "sample.id")
      sample.pheno <- use_data[, sampleCol]

      if (class(sample.geno) != class(sample.pheno)) {
        if (is.integer(sample.geno)) {
          sample.pheno = as.integer(sample.pheno)
        } else if (is.character(sample.geno)) {
          sample.pheno = as.character(sample.pheno)
        }
      }

      use_data = use_data[sample.pheno %in% sample.geno, ]

    } else if (is.character(genofile)) {   ## input is a plink file name

      fam_file = paste0(genofile, ".fam")
      if (!file.exists(fam_file)) stop("The plink fam file is not exist!")
      fam_data = fread(fam_file, data.table = F)
      use_data = use_data[use_data[, sampleCol] %in% fam_data$V2, ]

    } else {   ## input is a genotype data matrix or data.frame
      use_data = use_data[use_data[, sampleCol] %in% genofile$IID, ]
    }
  }


  ### sort survival time
  use_data = use_data[order(use_data[, timeCol], decreasing = FALSE), ]
  fail_time = unique(use_data[, timeCol][which(use_data[, statusCol] == 1)])


  if(verbose) print("Fitting null model...")

  #### fitting null model
  if (is.null(PRSCol)) {
    formula_null = formula(paste0("Surv(", timeCol, ",", statusCol, ") ~ ", paste0(covCol, collapse = "+")))
  } else {
    formula_null = formula(paste0("Surv(", timeCol, ",", statusCol, ") ~ ", paste0(covCol, collapse = "+"), " + offset(", PRSCol, ")"))
  }

  null_model = coxph(formula_null, data = use_data, x = TRUE, ties = "breslow")
  if (verbose) print(summary(null_model)$coef)
  if (verbose) print(null_model$concordance[6])
  # formula_null = null_model$formula
  coefficient_null = summary(null_model)$coef
  concordance_null = null_model$concordance[6:7]

  res_null = null_model$residuals
  res_null_var = var(res_null)
  mu_null = use_data$status - res_null
  eta_null = null_model$linear.predictors
  eta_exp = exp(eta_null)  # hazard for each person


  if(verbose) print("Start contrust P matrix...")
  ## generate transposed P matrix (use t_P_mat can accelerate following calculation)
  t_P_nrow = as.numeric(length(fail_time))
  t_P_ncol = as.numeric(nrow(use_data))
  t_P_length = t_P_nrow * t_P_ncol
  t_P_mat = matrix(0, nrow = length(fail_time), ncol = nrow(use_data))
  if(verbose) print(paste("Dimension of the P matrix:", nrow(use_data), length(fail_time)))
  for (j in 1:length(fail_time)) {
    riskset_j = which(use_data$time >= fail_time[j])
    failset_j = which(use_data$time == fail_time[j] & use_data$status == 1)
    d_j = length(failset_j)

    t_P_mat[j, riskset_j] = sqrt(d_j) * eta_exp[riskset_j]/ sum(eta_exp[riskset_j])
    if(verbose){if(j %% 1000 == 0) print(paste0("Complete ", j, "/", length(fail_time)))}
  }


  W_mat = Diagonal(x = mu_null)
  X_mat = as.matrix(use_data[, covCol])   # covariate matrix X
  PX_mat = t_P_mat %*% X_mat              # P^TX

  XWX_mat = t(W_mat %*% X_mat) %*% X_mat  # X^TWX
  XVX_inv = solve(XWX_mat - (t(PX_mat) %*% PX_mat))      # (XVX)^-1 = (X(W-PP)X)^-1 = (XWX-XPPX)^-1

  t_X_mat = cbind(1, X_mat)
  X_XXinv = t_X_mat %*% solve(crossprod(t_X_mat, t_X_mat))


  if (is.null(LOCO) || isFALSE(LOCO)) {
    fit_null = list(res_null = res_null, res_null_var = res_null_var, use_SPA = use_SPA, LOCO = FALSE,
                    mu_null = mu_null, eta_null = eta_null, eta_exp = eta_exp,
                    formula_null = formula_null, coefficient_null = coefficient_null, concordance_null = concordance_null,
                    W_mat = W_mat, X_mat = X_mat, t_P_mat = t_P_mat,
                    PX_mat = PX_mat, XVX_inv = XVX_inv,
                    t_X_mat = t_X_mat, X_XXinv = X_XXinv, use_data = use_data)
  } else {
    fit_null = list(res_null = res_null, res_null_var = res_null_var, use_SPA = use_SPA, LOCO = LOCO, chr = chr,
                    mu_null = mu_null, eta_null = eta_null, eta_exp = eta_exp,
                    formula_null = formula_null, coefficient_null = coefficient_null, concordance_null = concordance_null,
                    W_mat = W_mat, X_mat = X_mat, t_P_mat = t_P_mat,
                    PX_mat = PX_mat, XVX_inv = XVX_inv,
                    t_X_mat = t_X_mat, X_XXinv = X_XXinv, use_data = use_data)
  }


  if(use_SPA){
    fit_null = CGF4MartingaleRes(fit_null = fit_null, range = range, length.out = length.out, verbose = verbose)
  }


  return(fit_null)
}
