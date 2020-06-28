require(vcfR)
require(tibble)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)

### Script to convert a vcf to a bedassle file
# MAKE SURE THE POPULATIONS IN THE POPS FILE ARE IN THE CORRECT ORDER
# The populations variable name needs to be titled "pops" and the individuals need to be titled "individuals"

vcf_to_counts <- function(vcf_path, pops_path) {
  
  vcf <- read.vcfR(vcf_path)
  
  vcf_tidy <- vcfR2tidy(vcf)
  
  pops_df <- read_csv(pops_path)
  
  
  pops_vcf <- dplyr::left_join(pops_df, 
                        vcf_tidy$gt %>% 
                          rename(individuals = Indiv), 
                        by = "individuals") %>% 
    dplyr::left_join(vcf_tidy$fix, by = "ChromKey") %>% 
    dplyr::select(pops, individuals, CHROM, gt_GT) %>% 
    dplyr::mutate(gt_GT = replace_na(gt_GT, "0/0")) %>% #Even though this isn't the true genotype, this coding makes it easier to count the alternative alleles
    dplyr::mutate(a_1 = str_split_fixed(gt_GT, "/", 2)[,1] %>% as.numeric(), 
           a_2 = str_split_fixed(gt_GT, "/", 2)[,2] %>% as.numeric()) %>% 
    dplyr::mutate(alt_count = a_1 + a_2) %>% 
    dplyr::select(-gt_GT, -a_1, -a_2)
  
  allele_counts <- pops_vcf %>% 
    dplyr::group_by(pops, CHROM) %>% 
    dplyr::summarize(alt_count = sum(alt_count)) %>%
    tidyr::pivot_wider(names_from = CHROM, values_from = alt_count) %>% 
    as.matrix()
  
  return(allele_counts[,-1])
}




