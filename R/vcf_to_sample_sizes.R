require(vcfR)
require(tibble)
require(dplyr)
require(tidyr)
require(readr)
require(purrr)

vcf_to_sample_sizes <- function(vcf_path, pops_path){
  d_vcf <- vcfR::read.vcfR(vcf_path)
  
  d_tidy <- vcfR2tidy(d_vcf, format_fields = c("GT"))
  
  
  d_pops_df <- read_csv(pops_path)
  
  ## This made a list, where list items were the individuals contained in a population. Leaving it in in case I ever need it
  #d_pops_list <- d_pops_df %>% 
  #  split(d_pops_df$pops) %>% 
  #  map(function(x) {
  #    s <- dplyr::dplyr::select(x, individuals)
  #    as.vector(s$individuals)
  #  })
  
  pops_vcf <- left_join(d_pops_df, 
                        d_tidy$gt %>% 
                          rename(individuals = Indiv), 
                        by = "individuals") %>% 
    left_join(d_tidy$fix, by = "ChromKey") %>% 
    dplyr::select(pops, individuals, CHROM, gt_GT)
  
  
  locus_counts <- pops_vcf %>%
    mutate(gt_GT = replace_na(gt_GT, -1)) %>% 
    dplyr::group_by(pops, CHROM) %>% 
    dplyr::summarize(n = sum(gt_GT > -1) * 2) %>% # multiplying by two to get the number of chromosomes
    pivot_wider(names_from = CHROM, values_from = n) %>% 
    map_dfc(function(x) replace_na(x, 0)) %>% 
    dplyr::select(-pops) %>% 
    as.matrix()
  
  return(locus_counts) 
}

