#' STAAR procedure using omnibus test
#'
#' The \code{SurvSTAAR} function takes in genotype, the object from fitting the
#' null model, and functional annotation data to analyze the association between
#' a time-to-event phenotype and a variant-set by using SurvSTAAR procedure.
#'
#'
#' @param Geno an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param MAF minor allele frequency, can be NULL.
#' @param MAC minor allele count, can be NULL.
#' @param objNull an object from null model, which is the output of \code{nullModel}.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' @param rare_maf_cutoff a numeric value indicating the maximum minor allele frequency
#' for defining rare variants (default = 0.01).
#' @param rare_num_cutoff a numeric value indicating the minimum number of variants
#' required for analyzing a given variant set (default = 2).
#' @param combine_ultra_rare a logical variable indicating whether to combine ultra-rare
#' variants in the SKAT and ACAT-V test under imbalanced phenotypes (default = TRUE).
#' @param ultra_rare_mac_cutoff a numeric value indicating the maximum MAC for
#' ultra-rare variants (default = 20).
#' @param use_SPA a logical value indicating whether to use SaddlePoint
#' Approximation (SPA) in the analysis. Highly recommended for imbalanced
#' phenotypes (default = TRUE).
#' @param SPA_filter a logical value indicating whether to use SPA to calibrate
#' P-values for imbalanced traits when the P-value is below a pre-specified
#' threshold (default = TRUE).
#' @param SPA_filter_cutoff a numeric value specifying the threshold for recalculating
#' P-values using the SPA method (default = 0.05).
#' @param verbose a logical value indicating whether to provides additional detailed
#' information or messages during the execution of this function (default = FALSE).
#'
#' @returns A list containing the set-based test results from \code{SurvSTAAR.}
#' \code{num_variant} is the number of rare variants in this analyzed set.
#' \code{MAF} is minor allele frequency.
#' \code{MAC} is minor allele count.
#'
#' @export
#'
#'
SurvSTAAR = function(Geno, MAF = NULL, MAC = NULL, objNull, annotation_phred = NULL,
                     rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                     combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                     use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) {
    stop("Genotype is not a matrix!")
  }

  if(inherits(Geno, "sparseMatrix")){
    Geno = as.matrix(Geno)
  }

  if (ncol(Geno) < rare_num_cutoff) {
    stop(paste0("Number of variants in the set is less than ", rare_num_cutoff, ", will skip this category..."), call. = FALSE)
    return(list("results_STAAR_O" = NA))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if (nrow(annotation_phred) != 0 & ncol(Geno) != nrow(annotation_phred)) {
    stop("Dimensions don't match for genotype and annotation!")
  }

  if (inherits(Geno, "sparseMatrix")) {
    Geno <- as.matrix(Geno)
  }

  if (is.null(MAF) | is.null(MAC)) {

    genotype = genoFlip(Geno = Geno)
    MAF = genotype$G_summary$MAF
    MAC = floor(MAF * nrow(Geno) * 2)

    RV_label = as.vector((MAF < rare_maf_cutoff)&(MAF > 0))
    Geno = genotype$Geno[ ,RV_label]

    rm(genotype)

  } else {

    RV_label = as.vector((MAF < rare_maf_cutoff)&(MAF > 0))
    Geno = Geno[ ,RV_label]

  }

  MAF = MAF[RV_label]
  MAC = MAC[RV_label]
  Geno = as(Geno, "CsparseMatrix")
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]


  if(sum(RV_label) < rare_num_cutoff) {

    stop(paste0("Number of rare variant in the set is less than ", rare_num_cutoff, " will skip this category..."), call. = FALSE)
    return(list("results_STAAR_O" = NA))

  } else {

    if (combine_ultra_rare) {
      ultra_rare_length = length(which(MAC < ultra_rare_mac_cutoff))
      if (verbose) print(paste0("Apply STAAR-O test to ", sum(RV_label), " rare variants, with ", ultra_rare_length, " ultra rare variants"))
    } else {
      if (verbose) print(paste0("Apply STAAR-O test to ", sum(RV_label), " rare variants"))
    }

    annotation_rank <- 1 - 10^(-annotation_phred/10)


    w_1 <- dbeta(MAF, 1, 25)
    w_2 <- dbeta(MAF, 1, 1)

    w_a_1 <- w_1^2/dbeta(MAF,0.5,0.5)^2
    w_a_2 <- w_2^2/dbeta(MAF,0.5,0.5)^2

    if(dim(annotation_phred)[2] == 0){

      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1, w_2))
      w_A <- as.matrix(cbind(w_a_1, w_a_2))

    }else{

      ## Burden
      w_B = as.matrix(cbind(w_1, annotation_rank*w_1, w_2, annotation_rank*w_2))
      ## SKAT
      w_S = as.matrix(cbind(w_1, sqrt(annotation_rank)*w_1, w_2, sqrt(annotation_rank)*w_2))
      ## ACAT-V
      w_A = as.matrix(cbind(w_a_1, annotation_rank*w_a_1, w_a_2, annotation_rank*w_a_2))

    }


    pvalues <- SurvSTAAR_O(Geno = Geno, objNull = objNull, annotation_rank, MAC = MAC,
                           use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                           weight_A = w_A, weight_B = w_B, weight_S = w_S,
                           combine_ultra_rare, ultra_rare_mac_cutoff, verbose = verbose)


    return(c(pvalues,
             list(num_variant = sum(RV_label),
                  MAC = MAC,
                  MAF = MAF)))

  }
}
