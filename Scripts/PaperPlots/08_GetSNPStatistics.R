#'---
#' title: 1000G variant statistics
#' author: Christian Mertes
#' wb:
#'  input: 
#'   - snpCounts: 'Data/1000G_snp_counts.tsv'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#' 
#' # How to get the data from web
#'
#'
#' 1000G variant number statistics
#' 
#' * Download from: https://public.tableau.com/profile/yali.xue#!/vizhome/1000G_phase3_per_individual_count_all/Sheet1
#' 
#' * convert to UTF8: iconv -f UTF-16LE -t UTF-8 tmp.tsv > Data/1000G_snp_counts.tsv
#' 


#' 
#' # Load data set
#' 
library(data.table)
dt <- fread('./Data/1000G_snp_counts.tsv')

#' 
#' ## Tags for frequency and consequence
#' 
unique(dt$`Freq BIN`)
unique(dt[,`Consequence Type`])


not_used_tags <- c(
    "CLINVAR_clnsig=4", "CLINVAR_clnsig=5", "ctcf_insulator", 
    "FunSeq_score>1.5", "GWAS_3465", "HGMD-DM",
    "Large_Deletion_NON_OVERLAP_GENE", "Large_Deletion_OVERLAP_GENE", 
    "mature_miRNA_variant", "PHOSPHORYLATION_evidence<=5", 
    "PHOSPHORYLATION_evidence>5", "small_indel_all_NON_OVERLAP_GENE", 
    "small_indel_all_OVERLAP_GENE", "small_indel_filtered_NON_OVERLAP_GENE", 
    "small_indel_filtered_OVERLAP_GENE", "synonymous_variant", 
    "total_site", "unannotated_tfbs"
)

not_used_freq_tags <- c(">5%", "0.5%-5%")


#'
#' ### Define tags
#' 
protein_affecting_tags <- c(
    "frameshift_variant", "frameshift_variant_aloft", "inframe", 
    "missense_variant", "splice_acceptor_variant", 
    "splice_acceptor_variant_aloft", "splice_donor_variant",
    "splice_donor_variant_aloft", "stop_gained", "stop_gained_aloft",
    "stop_lost"
)

intergenic_region_tags <- c(
    "intergenic_variant", "distal_enhancer"
)

non_coding_tags <- c(
    "3_prime_UTR_variant", "5_prime_UTR_variant", 
    "intron_variant", "promoter", "proximal_enhancer", 
    "regulatory_region_variant", "TF_binding_site_variant")


rare_freq <- c("singleton-0.5%", "singleton")


#'
#' # Get distributions of variants all medians per group
#' 


#' ## Non coding + intergenic variants
#'
sort(c(non_coding_tags, intergenic_region_tags))
rare_freq
res <- dt[`Freq BIN` %in% rare_freq & 
            `Consequence Type` %in% c(non_coding_tags, intergenic_region_tags), 
    sum(`Median Site Count`),by='Sup Pop']
res
mean(res[,V1])
median(res[,V1])


#' ## Non coding variants
#'
non_coding_tags
rare_freq
res <- dt[`Freq BIN` %in% rare_freq & `Consequence Type` %in% non_coding_tags, 
   sum(`Median Site Count`),by='Sup Pop']
res
mean(res[,V1])
median(res[,V1])


#' 
#' ## Protein affecting variants
#' 
protein_affecting_tags
rare_freq
res <- dt[`Freq BIN` %in% rare_freq & `Consequence Type` %in% protein_affecting_tags, 
   sum(`Median Site Count`),by='Sup Pop']
res
mean(res[,V1])
median(res[,V1])

