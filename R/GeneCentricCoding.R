#' Title Gene-centric analysis of coding functional categories using SurvSTAAR procedure
#'
#' The \code{GeneCentricCoding} function takes in chromosome, gene name, functional
#' category, the object of opened annotated GDS file, and the object from null model
#' to analyze the association between a time-to-event phenotype and coding functional
#' categories of a gene by using SurvSTAAR procedure.
#'
#' For each coding functional category, the SurvSTAAR-O p-value is derived from an
#' omnibus test that combines SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1),
#' ACAT-V(1,25), and ACAT-V(1,1), with p-values weighted by each annotation using
#' the Cauchy combination method.
#'
#' @param gene_name a character string specifying the name of the gene to be analyzed.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param objNull an object from null model, which is the output of \code{nullModel}.
#' @param genes_info a matrix containing information on gene names and chromosomes.
#' @param variant_type the type of variant included in the analysis.
#' Options are \code{“SNV”}, \code{“Indel”}, or \code{“variant”} (default = “SNV”).
#' @param categories the coding functional category to be analyzed using SurvSTAAR procedure.
#' Choices include \code{all}, \code{all_ptv}, \code{plof}, \code{plof_ds}, \code{ptv},
#' \code{ptv_ds}, \code{synonymous}, \code{missense}, \code{dmissense} (default = "all").
#' @param rare_maf_cutoff a numeric value indicating the maximum minor allele frequency
#' for defining rare variants (default = 0.01).
#' @param rare_num_cutoff a numeric value indicating the minimum number of variants
#' required for analyzing a given variant set (default = 2).
#' @param geno_missing_cutoff a numeric value specifying the cutoff for missing rates.
#' Variants with a missing rate higher than this cutoff are excluded (default = 1,
#' do not remove any variants).
#' @param geno_missing_imputation a character string specifying the method for
#' handling missing genotypes. Options are “mean” or “minor” (default = “mean”).
#' @param min_maf_cutoff a numeric variable representing the minimum minor allele
#' frequency (MAF) to include in the individual test (default = 0).
#' @param combine_ultra_rare a logical variable indicating whether to combine ultra-rare
#' variants in the SKAT and ACAT-V test under imbalanced phenotypes (default = TRUE).
#' @param ultra_rare_mac_cutoff a numeric value indicating the maximum MAC for
#' ultra-rare variants (default = 20).
#' @param QC_label the channel name of the QC label in the GDS/aGDS file
#' (default = "annotation/filter").
#' @param Annotation_dir the channel name of the annotations in the aGDS file
#' (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing annotation names and their
#' corresponding channel names in the aGDS file.
#' @param Use_annotation_weights a logical variable indicating whether to use annotations
#' as weights (default = TRUE).
#' @param Annotation_name a vector of annotation names used in SurvSTAAR (default = NULL).
#' @param use_SPA a logical value indicating whether to use SaddlePoint
#' Approximation (SPA) in the analysis. Highly recommended for imbalanced
#' phenotypes (default = TRUE).
#' @param SPA_filter a logical value indicating whether to use SPA to calibrate
#' P-values for imbalanced traits when the P-value is below a pre-specified
#' threshold (default = TRUE).
#' @param SPA_filter_cutoff a numeric value specifying the threshold for recalculating
#' P-values using the SPA method (default = 0.05).
#' @param rm_long a logical value indicating whether to remove the long masks
#' during the analysis (default = TRUE).
#' @param rm_long_cutoff a numeric value specifying the cutoff for the number of
#' variants in a mask to remove (default = 3000).
#' @param verbose a logical value indicating whether to provides additional detailed
#' information or messages during the execution of this function (default = FALSE).
#'
#' @returns a list of results containing the SurvSTAAR results.
#' \code{gene_info} is the information about the analyzed gene name and chromosome.
#' \code{SurvSTAAR_O_all} is the p-values of the functional categories from the SurvSTAAR-O test.
#' Other data frames are the detailed results for each functional category.
#'
#'
GeneCentricCoding = function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                             categories = c("all", "all_ptv", "plof", "plof_ds", "ptv", "ptv_ds", "synonymous", "missense", "dmissense"),
                             rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                             geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                             min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                             QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                             Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                             use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                             rm_long = TRUE, rm_long_cutoff = 3000, verbose = FALSE) {


  if (categories == "all") {

    results = coding(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                     rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                     min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                     QC_label = QC_label, Annotation_dir = Annotation_dir,
                     Annotation_name_catalog = Annotation_name_catalog,
                     Use_annotation_weights = Use_annotation_weights,
                     Annotation_name = Annotation_name,
                     use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                     rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "all_ptv") {

    results = coding_incl_ptv(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                              rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                              min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                              QC_label = QC_label, Annotation_dir = Annotation_dir,
                              Annotation_name_catalog = Annotation_name_catalog,
                              Use_annotation_weights = Use_annotation_weights,
                              Annotation_name = Annotation_name,
                              use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                              rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "plof") {

    results = plof(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                   rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                   min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                   QC_label = QC_label, Annotation_dir = Annotation_dir,
                   Annotation_name_catalog = Annotation_name_catalog,
                   Use_annotation_weights = Use_annotation_weights,
                   Annotation_name = Annotation_name,
                   use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                   rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "plof_ds") {

    results = plof_ds(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                      rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                      min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                      QC_label = QC_label, Annotation_dir = Annotation_dir,
                      Annotation_name_catalog = Annotation_name_catalog,
                      Use_annotation_weights = Use_annotation_weights,
                      Annotation_name = Annotation_name,
                      use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                      rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "ptv") {

    results = ptv(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                  rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                  min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                  QC_label = QC_label, Annotation_dir = Annotation_dir,
                  Annotation_name_catalog = Annotation_name_catalog,
                  Use_annotation_weights = Use_annotation_weights,
                  Annotation_name = Annotation_name,
                  use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                  rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "ptv_ds") {

    results = ptv_ds(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                     rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                     min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                     QC_label = QC_label, Annotation_dir = Annotation_dir,
                     Annotation_name_catalog = Annotation_name_catalog,
                     Use_annotation_weights = Use_annotation_weights,
                     Annotation_name = Annotation_name,
                     use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                     rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "synonymous") {

    results = synonymous(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                         rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                         min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                         QC_label = QC_label, Annotation_dir = Annotation_dir,
                         Annotation_name_catalog = Annotation_name_catalog,
                         Use_annotation_weights = Use_annotation_weights,
                         Annotation_name = Annotation_name,
                         use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                         rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "missense") {

    results = missense(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                       rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                       min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                       QC_label = QC_label, Annotation_dir = Annotation_dir,
                       Annotation_name_catalog = Annotation_name_catalog,
                       Use_annotation_weights = Use_annotation_weights,
                       Annotation_name = Annotation_name,
                       use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                       rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  } else if (categories == "dmissense") {

    results = dmissense(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                        rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                        min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                        QC_label = QC_label, Annotation_dir = Annotation_dir,
                        Annotation_name_catalog = Annotation_name_catalog,
                        Use_annotation_weights = Use_annotation_weights,
                        Annotation_name = Annotation_name,
                        use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                        rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)

  }

  seqResetFilter(genofile)
  gc()

  return(results)
}
