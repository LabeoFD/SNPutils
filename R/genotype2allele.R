#' Convert Genotype Data to Allele Pairs
#'
#' Converts numerical genotype data (0, 1, 2) to actual allele pairs using annotation data.
#'
#' @param genotype_data A data frame containing genotype data with probeset_id column
#' @param annotation_data A data frame containing annotation data with Probe.Set.ID, Allele.A, and Allele.B columns
#' @param validate_alleles Logical indicating whether to validate allele assignments
#' @param verbose Logical indicating whether to print validation results
#'
#' @return A list containing validation results, converted genotypes
#' @export
#'
#' @examples
#' \dontrun{
#' results <- genotype2allele(genotype_data, annotation_data)
#' }
genotype2allele <- function(genotype_data, annotation_data, validate_alleles = TRUE, verbose = TRUE) {
  # Prepare data tables
  dt_genotype <- data.table::as.data.table(genotype_data)
  dt_annotation <- data.table::as.data.table(annotation_data)
  
  # Set keys for optimized operations
  data.table::setkey(dt_genotype, probeset_id)
  data.table::setkey(dt_annotation, Probe.Set.ID)
  
  # Transform genotype data to long format
  dt_genotype_long <- dt_genotype %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = -probeset_id,
      names_to = "Sample_ID",
      values_to = "genotype"
    ) %>%
    data.table::as.data.table()
  
  # Create allele lookup table
  allele_lookup <- dt_annotation %>%
    tibble::as_tibble() %>%
    dplyr::select(
      probeset_id = Probe.Set.ID,
      allele_a = Allele.A,
      allele_b = Allele.B
    )
  
  # Validation results placeholder
  validation_results <- list(
    allele_validation = NULL,
    rule_violations = NULL
  )
  
  # Perform validation if requested
  if (validate_alleles) {
    validation_results <- validateAlleles(allele_lookup, verbose)
  }
  
  # Convert genotype codes to allele pairs
  results <- dt_genotype_long %>%
    dplyr::left_join(allele_lookup, by = "probeset_id") %>%
    dplyr::mutate(
      allele_genotype = dplyr::case_when(
        
	# First check if alleles exist (handle missing matches)
        is.na(allele_a) | is.na(allele_b) ~ NA_character_,					 
        # Indels
        allele_a == "-" & genotype == 0 ~ "D/D",
        allele_a == "-" & genotype == 1 ~ "D/I",
        allele_a == "-" & genotype == 2 ~ "I/I",
        
        # Standard SNPs
        genotype == 0 ~ paste0(allele_a, "/", allele_a),
        genotype == 1 ~ paste0(allele_a, "/", allele_b),
        genotype == 2 ~ paste0(allele_b, "/", allele_b),
        
        # Catch other cases
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(probeset_id, Sample_ID, allele_genotype)
  
  # Return results including validation
  return(list(
    validation = validation_results,
    results = results
  ))
}

