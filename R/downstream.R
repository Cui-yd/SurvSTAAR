downstream <- function(gene_name, genofile, objNull, variant_type=c("SNV","Indel","variant"), genes_info,
                       rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                       geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                       min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                       QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                       Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                       use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                       rm_long = TRUE, rm_long_cutoff = 5000, verbose = FALSE) {

  ## individuals
  phenotype.id = objNull$use_data$IID

  ## SPA status
  if(is.null(use_SPA)) use_SPA = objNull$use_SPA

  ## Gene Info ######
  kk = which(genes_info[,2]==gene_name)
  gene_info_kk = genes_info[kk, 1:2]

  ## get SNV id
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

  variant.id <- seqGetData(genofile, "variant.id")

  rm(filter)
  gc()


  ## Downstream ######

  if (verbose) print(paste0(">>> Analyzing downstream ...        ", Sys.time(), " <<<"))

  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]

  rm(GENCODE.Category)
  gc()

  seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)

  rm(variant.id.downstream)
  gc()

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

  rm(GENCODE.Info)
  gc()

  rm(variant_gene_num)
  gc()

  Gene <- as.character(unlist(GENCODE.Info.split))

  rm(GENCODE.Info.split)
  gc()

  seqResetFilter(genofile)

  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]

  if (length(variant.is.in) < 2) {
    message("Variants number of *Downstream* is less than 2, will skip this category...")
    result.downstream = list("SurvSTAAR_O" = NA)
  } else if (rm_long && length(variant.is.in) > rm_long_cutoff) {
    message(paste0("Variants number of *Downstream* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.downstream = list("SurvSTAAR_O" = NA)
  } else {

    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

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

    ## impute & flip
    getGeno = genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                         min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno = getGeno$Geno
    MAF = getGeno$G_summary$MAF
    MAC = getGeno$G_summary$MAC
    Anno.Int.PHRED.sub = Anno.Int.PHRED.sub[getGeno$include_index, ]

    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub[,2]))) != 0) {
      include_index = which(!is.na(Anno.Int.PHRED.sub[,2]))
      Geno = Geno[, include_index]
      MAF = MAF[include_index]
      MAC = MAC[include_index]
      Anno.Int.PHRED.sub = Anno.Int.PHRED.sub[include_index,]
    }

    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *Downstream* is less than 1, will skip this category...")
      result.downstream = list("SurvSTAAR_O" = NA)
    } else {
      result.downstream = try(SurvSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
                                        rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                                        use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE)

      if (inherits(result.downstream, "try-error")) {
        result.downstream = list("SurvSTAAR_O" = NA)
      }
    }
  }
  seqResetFilter(genofile)

  result = c(list("gene_info" = gene_info_kk, "category" = "downstream"), result.downstream)

  return(result)
}
