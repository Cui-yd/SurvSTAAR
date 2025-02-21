#' Score test for individual variants in a given region of an annotated GDS (aGDS) format file
#'
#' The \code{IndividualTestGDS} performs association analysis between a
#' time-to-event (TTE) phenotype and individual variants in a aGDS file using
#' the exact score test.
#'
#' @param objNull an object from null model, which is the output of \code{nullModel}.
#' @param start_loc starting location (position) of the genetic region for each
#' individual variant to be analyzed using score test.
#' @param end_loc ending location (position) of the genetic region for each
#' individual variant to be analyzed using score test.
#' @param genofile an SeqVarGDSClass object of opened annotated GDS (aGDS) file.
#' @param chr a numeric value representing the chromosome.
#' @param rs_channel channel name of the rsID in the aGDS file. Can be NULL.
#' @param QC_label channel name of the QC label in the aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include
#' "variant", "SNV", or "Indel" (default = "variant").
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
#' @import SeqArray
#' @importFrom dplyr left_join
#' @importFrom SeqVarTools isSNV
#'
#' @returns A data frame containing the score test results and the estimated effect
#' size of the minor allele for each individual variant in the specified genotype file.
#' The columns include: CHR (chromosome), MarkerID (rsID), POS (position),
#' REF (reference allele), ALT (alternative allele), Pvalue (P-value estimated
#' from Chi-squared distribution),	Pvalue_SPA (P-value approximated using the SPA),
#' pvalue_log10 (-log10(P-value). If \code{use_SPA = TRUE}, this is calculated
#' using Pvalue_SPA; otherwise, it is based on Pvalue.), Score (estimated score
#' statistic), Variance (variance of the score statistic), Stest (Score^2 / Variance),
#' Est (estimated effect), Est_se (standard error of the estimated effect size),
#' MAF (minor allele frequency), MAC (minor allele count), missing_rate (missing rate).
#'
#' @export

IndividualTestGDS = function(objNull, start_loc, end_loc, genofile, chr, rs_channel = NULL,
                             QC_label="annotation/filter", variant_type = c("variant","SNV","Indel"),
                             use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                             geno_missing_cutoff = 1, geno_missing_imputation = c("mean", "minor"),
                             min_mac_cutoff = NULL, min_maf_cutoff = 0.01,
                             chunk_size = 5e3, verbose = FALSE) {

  phenotype.id = objNull$use_data$IID
  use_data = objNull$use_data

  id.genotype <- seqGetData(genofile,"sample.id")
  if (class(id.genotype) != class(phenotype.id)) {
    if (is.integer(id.genotype)) {
      phenotype.id = as.integer(phenotype.id)
    } else if (is.character(id.genotype)) {
      phenotype.id = as.character(phenotype.id)
    }
  }
  if (sum(phenotype.id %in% id.genotype) != length(phenotype.id)) {
    stop(paste0("The subjects used for fitting the null model are not included in the provided genotype data.", "\n",
                "       Please refit the null model using this genotype data!"))
  }


  seqSetFilter(genofile, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile,"sample.id")
  id.genotype.merge <- data.frame(id.genotype, index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  filter <- seqGetData(genofile, QC_label)
  if (variant_type=="variant") {
    SNVlist <- filter == "PASS"
  } else if (variant_type=="SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  } else if (variant_type=="Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  position <- as.numeric(seqGetData(genofile, "position"))

  variant.id <- seqGetData(genofile, "variant.id")
  is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
  SNV.id <- variant.id[is.in]

  N = nrow(use_data)
  M = length(SNV.id)

  N_chunk = ceiling(M / chunk_size)
  if (verbose) print(paste0("There are ", M, " variants, devided into ", N_chunk, " chunks"))


  indivResult = NULL

  begin_loop = Sys.time()

  for (chunk_i in 1:N_chunk) {

    if (chunk_i == N_chunk) {
      is.in = ((N_chunk-1) * chunk_size + 1): M
    } else {
      is.in = 1:chunk_size + (chunk_i-1) * chunk_size
    }
    seqSetFilter(genofile, variant.id = SNV.id[is.in], sample.id = phenotype.id)



    if (verbose) {
      print(paste0('Analyzing chunk[', chunk_i, '/', N_chunk, ']   ', round(chunk_i/N_chunk*100, 2), '%   Time: ', Sys.time()))
      print(paste0("Markers: ", length(is.in)))
      print(paste0("Samples: ", N))
    }

    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    CHR <- as.numeric(seqGetData(genofile, "chromosome"))
    POS <- as.numeric(seqGetData(genofile, "position"))
    REF <- as.character(seqGetData(genofile, "$ref"))
    ALT <- as.character(seqGetData(genofile, "$alt"))
    if (!is.null(rs_channel)) {
      rsID = as.character(seqGetData(genofile, rs_channel))
      G_info = data.frame(CHR, POS, REF, ALT, rsID)
    } else {
      G_info = data.frame(CHR, POS, REF, ALT)
    }

    begin = Sys.time()

    ## read and impute genotype data
    getGeno = genoMatrixGDS(Geno = Geno, G_info = G_info,
                            geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                            min_mac_cutoff = min_mac_cutoff, min_maf_cutoff = min_maf_cutoff)
    G_mat = getGeno$Geno
    G_summary = getGeno$G_summary
    G_maf = G_summary$MAF

    if (!all(G_summary$CHR==chr)) {
      warning("chr does not match the chromosome of genofile (the opened aGDS)!")
    }

    preprocess_time = as.numeric(difftime(Sys.time(), begin, units = "secs"))

    if (verbose) print(paste0("# markers:", ncol(G_mat)))
    if (verbose) print(paste0("Preprocess time: ", round(preprocess_time, 2), " secs"))


    if (verbose) print("Estimating exact score variance...")
    ## SPA for unbalanced phenotype
    if(is.null(use_SPA)) use_SPA = objNull$use_SPA
    begin = Sys.time()
    indivResult_i = exactScore(objNull = objNull, G_mat = G_mat, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = verbose)
    indivResult_i = indivResult_i$result
    indivResult_i = cbind(G_summary[,1:5], indivResult_i, G_summary[,6:8])
    exactTestTime = as.numeric(difftime(Sys.time(), begin, units = "secs"))
    if(verbose) print(paste0("Exact score variance, time of this chunk: ", round(exactTestTime, 2), " secs"))

    indivResult = rbind(indivResult, indivResult_i)

    seqResetFilter(genofile)
    gc()
  }

  if (verbose) {
    print(paste0("The total elapsed time is ", round(difftime(Sys.time(), begin_loop, units = "hour"), 2), " hours"))
  }

  return(indivResult)

}
