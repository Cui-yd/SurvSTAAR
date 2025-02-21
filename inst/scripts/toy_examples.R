library(data.table)
library(survival)
library(seqminer)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(Matrix)
library(stats)
library(CompQuadForm)
library(SurvSTAAR)

### The genotype datas are generated using using the calibration coalescent model
### The covariates are generated using R’s random functions (rbinom, rnorm)
### The phenotype is generated based on covariates and the causal gene ("Gene4"),
###     with 30% of its variants randomly selected as causal
### The functional annotations are generated based on standard normal distribution
### The Gencode are randomly generated using following codes:
# selected_vars = sample(1:5000, 3500, replace = FALSE)
# GENCODE.EXONIC.Category = sample(c("stopgain", "stoploss", "nonsynonymous SNV", "synonymous SNV"), 3500,
#                                  prob = c(0.15,0.15,0.35,0.35), replace = TRUE)
# GENCODE.Category = sample(c("splicing", "exonic;splicing", "ncRNA_splicing",
#                             "downstream", "upstream", "UTR3", "UTR5"), 3500,
#                           prob = c(0.15,0.05,0.1,0.2,0.2,0.15,0.15), replace = TRUE)

# fit null model ----------------------------------------------------------

pheno_dat = fread("inst/extdata/pheno.txt", data.table = FALSE)

## use polygenic effects
## declare: we set LOCO = FALSE here, for the toy data only contains one chromosome
objNull = NullModel(phenofile = pheno_dat, LOCO = FALSE, chr = 1,
                    statusCol = "status", timeCol = "time", sampleCol = "IID",
                    covCol = c("sex", "age"), PRSCol = "PE",
                    use_SPA = TRUE, verbose = TRUE)

## do not use polygenic effects
objNull = NullModel(phenofile = pheno_dat, LOCO = FALSE, chr = 1,
                    statusCol = "status", timeCol = "time", sampleCol = "IID",
                    covCol = c("sex", "age"), PRSCol = NULL,
                    use_SPA = TRUE, verbose = TRUE)


# individual analysis -----------------------------------------------------

### use plink format file ######

# Apply individual analysis to all the single variants in the plink file
# when there are not too many variants
individual_plink = IndividualTestPlink(objNull, genofile = "inst/extdata/toy", chr = 1, sampleCol = "IID",
                                       use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                       geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                       min_mac_cutoff = 20, min_maf_cutoff = 0.01,
                                       chunk_size = 1000, verbose = TRUE)

# Apply individual analysis to the selected single variants in the plink file
individual_plink = IndividualTestPlink(objNull, start_loc = 1, end_loc = 3e5, genofile = "inst/extdata/toy",
                                       chr = 1, sampleCol = "IID",
                                       use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                       geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                       min_mac_cutoff = 20, min_maf_cutoff = 0.01,
                                       chunk_size = 1000, verbose = TRUE)


### use agds format file ######

genofile = seqOpen("inst/extdata/toy.gds")
individual_agds = IndividualTestGDS(objNull, start_loc = 1, end_loc = 1e6, genofile = genofile,
                                    chr = 1, rs_channel = "annotation/id",
                                    QC_label = "annotation/info/QC_label", variant_type = "SNV",
                                    use_SPA = TRUE, SPA_filter = TRUE,SPA_filter_cutoff = 0.05,
                                    geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                    min_mac_cutoff = 20, min_maf_cutoff = 0.01,
                                    chunk_size = 1000, verbose = TRUE)



# gene-centric coding analysis --------------------------------------------

# genofile = seqOpen("inst/extdata/toy.gds")
load("inst/extdata/Annotation_name_catalog_toy.rda")
load("inst/extdata/genes_info_toy.rda")
Annotation_name = c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive",
                    "aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                    "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability",
                    "aPC.TF","aPC.Protein")

coding_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  coding_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                       genes_info = genes_info_toy, variant_type = "SNV",
                                       categories = "all", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                       geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                       min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                       QC_label = "annotation/info/QC_label",
                                       Annotation_dir = "annotation/info/FunctionalAnnotation",
                                       Annotation_name_catalog = Annotation_name_catalog_toy,
                                       Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                       use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                       rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(coding_agds[[j]]$SurvSTAAR_O_all)
}



coding_ptv_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  coding_ptv_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                           genes_info = genes_info_toy, variant_type = "SNV",
                                           categories = "all_ptv", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                           geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                           min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                           QC_label = "annotation/info/QC_label",
                                           Annotation_dir = "annotation/info/FunctionalAnnotation",
                                           Annotation_name_catalog = Annotation_name_catalog_toy,
                                           Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                           use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                           rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(coding_ptv_agds[[j]]$SurvSTAAR_O_all)
}


plof_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  plof_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                     genes_info = genes_info_toy, variant_type = "SNV",
                                     categories = "plof", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                     geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                     min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                     QC_label = "annotation/info/QC_label",
                                     Annotation_dir = "annotation/info/FunctionalAnnotation",
                                     Annotation_name_catalog = Annotation_name_catalog_toy,
                                     Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                     use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                     rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(plof_agds[[j]]$SurvSTAAR_O)
}


plof_ds_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  plof_ds_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                        genes_info = genes_info_toy, variant_type = "SNV",
                                        categories = "plof_ds", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                        geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                        min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                        QC_label = "annotation/info/QC_label",
                                        Annotation_dir = "annotation/info/FunctionalAnnotation",
                                        Annotation_name_catalog = Annotation_name_catalog_toy,
                                        Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                        use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                        rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(plof_ds_agds[[j]]$SurvSTAAR_O)
}


ptv_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  ptv_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                    genes_info = genes_info_toy, variant_type = "SNV",
                                    categories = "ptv", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                    geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                    min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                    QC_label = "annotation/info/QC_label",
                                    Annotation_dir = "annotation/info/FunctionalAnnotation",
                                    Annotation_name_catalog = Annotation_name_catalog_toy,
                                    Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                    use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                    rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(ptv_agds[[j]]$SurvSTAAR_O)
}


ptv_ds_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  ptv_ds_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                       genes_info = genes_info_toy, variant_type = "SNV",
                                       categories = "ptv_ds", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                       geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                       min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                       QC_label = "annotation/info/QC_label",
                                       Annotation_dir = "annotation/info/FunctionalAnnotation",
                                       Annotation_name_catalog = Annotation_name_catalog_toy,
                                       Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                       use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                       rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(ptv_agds[[j]]$SurvSTAAR_O)
}


missense_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  missense_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                         genes_info = genes_info_toy, variant_type = "SNV",
                                         categories = "missense", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                         geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                         min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                         QC_label = "annotation/info/QC_label",
                                         Annotation_dir = "annotation/info/FunctionalAnnotation",
                                         Annotation_name_catalog = Annotation_name_catalog_toy,
                                         Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                         use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                         rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(missense_agds[[j]]$SurvSTAAR_O)
}


dmissense_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  dmissense_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                          genes_info = genes_info_toy, variant_type = "SNV",
                                          categories = "dmissense", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                          geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                          min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                          QC_label = "annotation/info/QC_label",
                                          Annotation_dir = "annotation/info/FunctionalAnnotation",
                                          Annotation_name_catalog = Annotation_name_catalog_toy,
                                          Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                          use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                          rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(dmissense_agds[[j]]$SurvSTAAR_O)
}


synonymous_agds = NULL
for (i in 1:5) {
  gene_name = genes_info_toy$Gene[i]
  synonymous_agds[[i]] = GeneCentricCoding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                                           genes_info = genes_info_toy, variant_type = "SNV",
                                           categories = "synonymous", rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                           geno_missing_cutoff = 1, geno_missing_imputation = "mean",
                                           min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                           QC_label = "annotation/info/QC_label",
                                           Annotation_dir = "annotation/info/FunctionalAnnotation",
                                           Annotation_name_catalog = Annotation_name_catalog_toy,
                                           Use_annotation_weights = TRUE, Annotation_name = Annotation_name,
                                           use_SPA = TRUE, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                           rm_long = FALSE, verbose = TRUE)
}
for (j in 1:5) {
  print(synonymous_agds[[j]]$SurvSTAAR_O)
}


### The simulated fake data currently do not support all gene-centric noncoding
### and ncRNA category analyses; updates are planned soon.
### You can try these analyses by your own sequencing data using the SurvSTAARpipeline scripts,
### which can be found at https://github.com/Cui-yd/SurvSTAARpipeline. Enjoy! :)


