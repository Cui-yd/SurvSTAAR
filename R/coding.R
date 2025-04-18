coding <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                   rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                   geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                   min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                   QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                   Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                   use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                   rm_long = TRUE, rm_long_cutoff = 3000, verbose = FALSE) {

  ## individuals
  phenotype.id = objNull$use_data$IID

  ## SPA status
  if(is.null(use_SPA)) use_SPA = objNull$use_SPA

  ## get SNV id, position, REF, ALT (whole genome)
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")

  rm(filter)
  gc()

  ### Gene
  kk <- which(genes_info[,2]==gene_name)
  gene_info_kk = genes_info[kk, 1:2]

  sub_start_loc <- genes_info[kk,3]
  sub_end_loc <- genes_info[kk,4]

  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]

  rm(position)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))


  ## Coding ######

  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
  variant.id.gene <- variant.id.gene[lof.in.coding]

  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL

  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }

          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }

      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }


  ## plof_ds ######

  if (verbose) print(paste0(">>> Analyzing plof_ds ...           ", Sys.time(), " <<<"))

  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  variant.id.gene.category <- variant.id.gene[lof.in.plof]

  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *plof_ds* is less than 2, will skip this category...")
    result.plof_ds = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *plof_ds* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.plof_ds = list("SurvSTAAR_O" = NA)
  } else {

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]

    seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id = as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id = as.character(phenotype.id)
      }
    }

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute & flip & select rare variants
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *plof_ds* is less than 1, will skip this category...")
      result.plof_ds = list("SurvSTAAR_O" = NA)
    } else {
      result.plof_ds = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.plof_ds, "try-error")) {
        result.plof_ds = list("SurvSTAAR_O" = NA)
      }
    }
  }



  ## plof ######

  if (verbose) print(paste0(">>> Analyzing plof ...           ", Sys.time(), " <<<"))

  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
  variant.id.gene.category <- variant.id.gene[lof.in.plof]

  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *plof* is less than 2, will skip this category...")
    result.plof = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *plof* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.plof = list("SurvSTAAR_O" = NA)
  } else {

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]

    seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id = as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id = as.character(phenotype.id)
      }
    }

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute & flip & select rare variants
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *plof* is less than 1, will skip this category...")
      result.plof = list("SurvSTAAR_O" = NA)
    } else {
      result.plof = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                                  rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                  use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.plof, "try-error")) {
        result.plof = list("SurvSTAAR_O" = NA)
      }
    }
  }



  ## synonymous ######

  if (verbose) print(paste0(">>> Analyzing synonymous ...           ", Sys.time(), " <<<"))

  lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *synonymous* is less than 2, will skip this category...")
    result.synonymous = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *synonymous* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.synonymous = list("SurvSTAAR_O" = NA)
  } else {

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous, ]

    seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id = as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id = as.character(phenotype.id)
      }
    }

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute & flip & select rare variants
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *synonymous* is less than 1, will skip this category...")
      result.synonymous = list("SurvSTAAR_O" = NA)
    } else {
      result.synonymous = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                                        rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                        use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.synonymous, "try-error")) {
        result.synonymous = list("SurvSTAAR_O" = NA)
      }
    }
  }



  ## missense ######

  if (verbose) print(paste0(">>> Analyzing missense ...           ", Sys.time(), " <<<"))

  lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.missense]

  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *missense* is less than 2, will skip this category...")
    result.missense = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *missense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.missense = list("SurvSTAAR_O" = NA)
  } else {

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense, ]

    seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id = as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id = as.character(phenotype.id)
      }
    }

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute & flip & select rare variants
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *missense* is less than 1, will skip this category...")
      result.missense = list("SurvSTAAR_O" = NA)
    } else {
      result.missense = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                                      rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                      use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.missense, "try-error")) {
        result.missense = list("SurvSTAAR_O" = NA)
      }
    }
  }



  ## disruptive missense ######

  if (verbose) print(paste0(">>> Analyzing disruptive missense ...           ", Sys.time(), " <<<"))

  lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
  variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *dmissense* is less than 2, will skip this category...")
    result.dmissense = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *dmissense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.dmissense = list("SurvSTAAR_O" = NA)
  } else {

    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense, ]

    seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id = as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id = as.character(phenotype.id)
      }
    }

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index

    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute & flip & select rare variants
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub.category = Anno.Int.PHRED.sub.category[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *dmissense* is less than 1, will skip this category...")
      result.dmissense = list("SurvSTAAR_O" = NA)
    } else {
      result.dmissense = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                                       rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                       use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.dmissense, "try-error")) {
        result.dmissense = list("SurvSTAAR_O" = NA)
      }
    }
  }



  if (!is.na(result.missense$SurvSTAAR_O) & !is.na(result.dmissense$SurvSTAAR_O)) {

    result.missense$SurvSTAAR_pvalue = result.missense$SurvSTAAR_pvalue[, -ncol(result.missense$SurvSTAAR_pvalue)]
    result.missense$SurvSTAAR_pvalue = cbind(result.missense$SurvSTAAR_pvalue, result.dmissense$SurvSTAAR_pvalue[, 1])
    colnames(result.missense$SurvSTAAR_pvalue)[ncol(result.missense$SurvSTAAR_pvalue)] = "Disruptive"

    result.missense_temp = c(
      CCT(result.missense$SurvSTAAR_pvalue[1, ]),
      CCT(result.missense$SurvSTAAR_pvalue[2, ]),
      CCT(result.missense$SurvSTAAR_pvalue[3, ]),
      CCT(result.missense$SurvSTAAR_pvalue[4, ]),
      CCT(result.missense$SurvSTAAR_pvalue[5, ]),
      CCT(result.missense$SurvSTAAR_pvalue[6, ])
    )

    result.missense$SurvSTAAR_O = CCT(result.missense$SurvSTAAR_pvalue)
    result.missense$ACAT_O = CCT(c(result.missense$SurvSTAAR_pvalue[, 1]))
    result.missense$SurvSTAAR_SKAT = CCT(result.missense$SurvSTAAR_pvalue[1:2, ])
    result.missense$SurvSTAAR_ACAT = CCT(result.missense$SurvSTAAR_pvalue[3:4, ])
    result.missense$SurvSTAAR_Burden = CCT(result.missense$SurvSTAAR_pvalue[5:6, ])

    result.missense$SurvSTAAR_pvalue = cbind(result.missense$SurvSTAAR_pvalue, result.missense_temp)
    colnames(result.missense$SurvSTAAR_pvalue)[ncol(result.missense$SurvSTAAR_pvalue)] = "results_SurvSTAAR"

  }



  ## return results #####

  SurvSTAAR_O_all = c("plof" = result.plof$SurvSTAAR_O, "plof_ds" = result.plof_ds$SurvSTAAR_O,
                      "missense" = result.missense$SurvSTAAR_O, "dmissense" = result.dmissense$SurvSTAAR_O,
                      "synonymous" = result.synonymous$SurvSTAAR_O)

  results_return = list("gene_info" = gene_info_kk, "SurvSTAAR_O_all" = SurvSTAAR_O_all,
                        "plof" = result.plof, "plof_ds" = result.plof_ds,
                        "missense" = result.missense, "dmissense" = result.dmissense,
                        "synonymous" = result.synonymous)

}
