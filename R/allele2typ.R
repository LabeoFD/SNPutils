#' Convert Results to TYP Format
#'
#' Converts genotype results from genotype2allele() into TYP format,
#' which provides a structured representation with forward alleles.
#'
#' @param conversion_results Results from genotype2allele()
#'
#' @return A data frame in TYP format with Sample_ID, SNP_Name, Allele1_Forward, and Allele2_Forward columns
#' @export
#'
#' @examples
#' \dontrun{
#' results <- genotype2allele(genotype_data, annotation_data)
#' typ_data <- allele2typ(results)
#' }
allele2typ <- function(conversion_results) {
  pre_typ <- conversion_results$results %>% 
    dplyr::relocate(Sample_ID, probeset_id, allele_genotype) %>% 
    dplyr::rename(SNP_Name = probeset_id) %>% 
    dplyr::mutate(
      Allele1_Forward = stringr::str_extract(allele_genotype, "^[^/]+"),
      Allele2_Forward = stringr::str_extract(allele_genotype, "[^/]+$")
    ) %>% 
    dplyr::select(-allele_genotype)
  
  return(pre_typ)
}
