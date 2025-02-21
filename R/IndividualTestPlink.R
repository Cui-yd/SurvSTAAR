#' Score test for individual variants in a given plink format file
#'
#' The \code{IndividualTestPlink} performs association analysis between a
#' time-to-event (TTE) phenotype and individual variants in a PLINK file using
#' the exact score test.
#'
#' @param objNull an object from null model, which is the output of \code{nullModel}.
#' @param start_loc starting location (position) of the genetic region for each
#' individual variant to be analyzed using score test (default = NULL).
#' @param end_loc ending location (position) of the genetic region for each
#' individual variant to be analyzed using score test (default = NULL).
#' @param genofile a character string specifying the file name of the genotype file.
#' @param chr a numeric value representing the chromosome.
#' @param sampleCol a character string specifying the column name of sample IDs.
#' @param use_SPA a logical value indicating whether to use SaddlePoint
#' Approximation (SPA) in the analysis. Highly recommended for imbalanced
#' phenotypes (default = TRUE).
#' @param SPA_filter a logical value indicating whether to use SPA to calibrate
#' P-values for imbalanced traits when the P-value is below a pre-specified
#' threshold (default = TRUE).
#' @param SPA_filter_cutoff A numeric value specifying the threshold for recalculating
#' P-values using the SPA method (default = 0.05).
#' @param geno_missing_cutoff a numeric value specifying the cutoff for missing rates.
#' Variants with a missing rate higher than this cutoff are excluded (default = 1,
#' do not remove any variants).
#' @param geno_missing_imputation a character string specifying the method for
#' handling missing genotypes. Options are “mean” or “minor” (default = “mean”).
#' @param min_mac_cutoff a numeric variable representing the minimum minor allele
#' count (MAC) to include in the individual test. Can be NULL.
#' @param min_maf_cutoff a numeric variable representing the minimum minor allele
#' frequency (MAF) to include in the individual test (default = 0.01).
#' @param chunk_size a numeric variable representing the number of variants analyzed
#' in a single chunk (default = 1000).
#' @param verbose a logical value indicating whether to provides additional detailed
#' information or messages during the execution of this function (default = FALSE).
#'
#' @import Matrix
#' @import data.table
#' @import seqminer
#' @import methods
#'
#' @returns A data frame containing the score test results and the estimated effect
#' size of the minor allele for each individual variant in the specified genotype file.
#' The columns include: CHR (chromosome), MarkerID (rsID), POS (position),
#' REF1 (alternative allele), REF2 (reference allele), Pvalue (P-value estimated
#' from Chi-squared distribution),	Pvalue_SPA (P-value approximated using the SPA),
#' pvalue_log10 (-log10(P-value). If \code{use_SPA = TRUE}, this is calculated
#' using Pvalue_SPA; otherwise, it is based on Pvalue.), Score (estimated score
#' statistic), Variance (variance of the score statistic), Stest (Score^2 / Variance),
#' Est (estimated effect), Est_se (standard error of the estimated effect size),
#' MAF (minor allele frequency), MAC (minor allele count), missing_rate (missing rate).
#'
#' @export

IndividualTestPlink = function(objNull, start_loc = NULL, end_loc = NULL, genofile, chr, sampleCol,
                               use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                               geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                               min_mac_cutoff = NULL, min_maf_cutoff = 0.01,
                               chunk_size = 1e3, verbose = FALSE) {

  use_data = objNull$use_data
  phenotype.id = objNull$use_data$IID

  bim_file = paste0(genofile, ".bim")
  fam_file = paste0(genofile, ".fam")
  if (!file.exists(bim_file)) stop("The plink bim file is not exist!")
  if (!file.exists(fam_file)) stop("The plink fam file is not exist!")

  fam_data = fread(fam_file, data.table = F)
  bim_data = fread(bim_file, data.table = F)

  if (sum(phenotype.id %in% fam_data$V2) != length(phenotype.id)) {
    stop(paste0("The subjects used for fitting the null model are not included in the provided genotype data.", "\n",
                "       Please refit the null model using this genotype data!"))
  }

  sampleIndex = match(use_data[, sampleCol], fam_data$V2)

  N = length(sampleIndex)

  if (is.null(start_loc) & is.null(end_loc)) {
    if (length(unique(bim_data[,1])) > 1) {
      chr_index = which(bim_data[,1] == chr)
      M = length(chr_index)
    } else {
      M = nrow(bim_data)
    }
  } else {
    if (length(unique(bim_data[,1])) > 1) {
      chr_index = which(bim_data[,1] == chr)
      bim_data_chr = bim_data[chr_index,]
      chr_index = chr_index[which(bim_data_chr[,4] >= start_loc & bim_data_chr[,4] <= end_loc)]
      M = length(chr_index)
    } else {
      M = length(which(bim_data[,4] >= start_loc & bim_data[,4] <= end_loc))
    }
  }

  N_chunk = ceiling(M / chunk_size)
  if (verbose) print(paste0("There are ", M, " variants, devided into ", N_chunk, " chunks"))


  indivResult = NULL

  begin_loop = Sys.time()

  for (chunk_i in 1:N_chunk) {


    if (length(unique(bim_data[,1])) > 1) {
      if (chunk_i == N_chunk) {
        markerIndex = chr_index[((N_chunk-1) * chunk_size + 1): M]
      } else {
        markerIndex = chr_index[1:chunk_size + (chunk_i-1) * chunk_size]
      }
    } else {
      if (chunk_i == N_chunk) {
        markerIndex = ((N_chunk-1) * chunk_size + 1): M
      } else {
        markerIndex = 1:chunk_size + (chunk_i-1) * chunk_size
      }
    }


    if (verbose) {

      print(paste0('Analyzing chunk[', chunk_i, '/', N_chunk, ']   ', round(chunk_i/N_chunk*100, 2), '%   Time: ', Sys.time()))
      print(paste0("Markers: ", length(markerIndex)))
      print(paste0("Samples: ", length(sampleIndex)))
      message(paste0("-------- Reading plink files [", chunk_i, "/", N_chunk, "] ---------------------------------------------------"))

    }


    begin = Sys.time()

    ## read and impute genotype data
    G_mat = readPlinkToMatrixByIndex(genofile, sampleIndex, markerIndex)
    getGeno = genoMatrixPlink(Geno = G_mat, bim_data = bim_data, markerIndex = markerIndex,
                              geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                              min_mac_cutoff = min_mac_cutoff, min_maf_cutoff = min_maf_cutoff)

    G_mat = getGeno$Geno
    G_summary = getGeno$G_summary
    G_maf = G_summary$MAF

    preprocess_time = as.numeric(difftime(Sys.time(), begin, units = "secs"))

    if (verbose) print(paste0("# markers:", ncol(G_mat)))
    if (verbose) print(paste0("Preprocess time: ", round(preprocess_time, 2), " secs"))


    if (verbose) print("Estimating exact score variance...")
    ## SPA for unbalanced phenotype
    if(is.null(use_SPA)) use_SPA = objNull$use_SPA
    begin = Sys.time()
    indivResult_i = exactScore(objNull = objNull, G_mat = G_mat,
                               use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = verbose)
    indivResult_i = indivResult_i$result
    indivResult_i = cbind(G_summary[, c(1:5)], indivResult_i, G_summary[, c(6:8)])
    exactTestTime = as.numeric(difftime(Sys.time(), begin, units = "secs"))
    if(verbose) print(paste0("Exact score variance, time of this chunk: ", round(exactTestTime, 2), " secs"))

    indivResult = rbind(indivResult, indivResult_i)

  }

  if (verbose) {
    print(paste0("The total elapsed time is ", round(difftime(Sys.time(), begin_loop, units = "hour"), 2), " hours"))
  }

  return(indivResult)

}
