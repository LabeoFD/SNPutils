#' Validate Allele Assignments
#'
#' Validates allele assignments according to standard rules for genetic data.
#'
#' @param allele_lookup A data frame with probeset_id, allele_a, and allele_b columns
#' @param verbose Logical indicating whether to print validation results
#'
#' @return A list containing validation information and rule violations
#' @keywords internal
validateAlleles <- function(allele_lookup, verbose = TRUE) {
  # Classify variant types and check rules
  allele_validation <- allele_lookup %>%
    dplyr::mutate(
      # Variant classification
      is_indel = (allele_a == "-"),
      is_multibase = (!is_indel & (nchar(allele_a) > 1 | nchar(allele_b) > 1)),
      is_AT_SNP = (!is_indel & !is_multibase & (allele_a %in% c("A", "T") & allele_b %in% c("A", "T"))),
      is_CG_SNP = (!is_indel & !is_multibase & (allele_a %in% c("C", "G") & allele_b %in% c("C", "G"))),
      is_other_SNP = (!is_indel & !is_multibase & !is_AT_SNP & !is_CG_SNP),
      
      # Rule checks
      rule2_check = dplyr::case_when(
        is_AT_SNP | is_CG_SNP ~ (allele_a < allele_b),
        TRUE ~ TRUE
      ),
      rule3_check = dplyr::case_when(
        is_other_SNP ~ (allele_a %in% c("A", "T") & allele_b %in% c("C", "G")),
        TRUE ~ TRUE
      ),
      rule4_check = dplyr::case_when(
        is_indel ~ (allele_b != "-"),
        TRUE ~ TRUE
      ),
      rule5_check = dplyr::case_when(
        is_multibase ~ (allele_a < allele_b),
        TRUE ~ TRUE
      )
    )
  
  # Identify rule violations
  rule_violations <- allele_validation %>%
    dplyr::filter(!rule2_check | !rule3_check | !rule4_check | !rule5_check) %>%
    dplyr::mutate(
      violations = dplyr::case_when(
        !rule2_check ~ "Rule 2 violation: AT/CG SNP not in alphabetical order",
        !rule3_check ~ "Rule 3 violation: Non-AT/CG SNP with incorrect allele assignment",
        !rule4_check ~ "Rule 4 violation: Indel with incorrect allele assignment",
        !rule5_check ~ "Rule 5 violation: Multi-base alleles not in alphabetical order",
        TRUE ~ "Unknown violation"
      )
    )
  
  # Report violations if verbose
  if (verbose && nrow(rule_violations) > 0) {
    message("Found ", nrow(rule_violations), " rule violations in annotation data")
    print(table(rule_violations$violations))
    message("First 10 violations:")
    print(head(rule_violations %>% dplyr::select(probeset_id, allele_a, allele_b, violations), 10))
  } else if (verbose) {
    message("No rule violations found - all rules are respected.")
  }
  
  return(list(
    allele_validation = allele_validation,
    rule_violations = rule_violations
  ))
}
