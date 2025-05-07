# Declare global variables used in package functions to avoid R CMD check notes
utils::globalVariables(c(
  "Allele.A", "Allele.B", "Probe.Set.ID", "Sample_ID", "allele_a", "allele_b",
  "allele_genotype", "is_AT_SNP", "is_CG_SNP", "is_indel", "is_multibase",
  "probeset_id", "rule2_check", "rule3_check", "rule4_check", "rule5_check",
  "violations"
))
