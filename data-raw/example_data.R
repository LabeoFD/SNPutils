# Create example genotype data
example_genotypes <- data.frame(
  probeset_id = paste0("SNP", 1:20),
  Sample1 = sample(c(0, 1, 2), 20, replace = TRUE),
  Sample2 = sample(c(0, 1, 2), 20, replace = TRUE),
  Sample3 = sample(c(0, 1, 2), 20, replace = TRUE),
  stringsAsFactors = FALSE
)

# Create corresponding annotation data for regular SNPs (first 10)
example_annotations <- data.frame(
  Probe.Set.ID = paste0("SNP", 1:20),
  Chromosome = sample(c(1:22, "X", "Y"), 20, replace = TRUE),
  Physical.Position = sample(1:100000000, 20),
  Allele.A = c(sample(c("A", "T"), 10, replace = TRUE), rep("-", 10)),
  Allele.B = c(sample(c("C", "G"), 10, replace = TRUE), character(10)),
  stringsAsFactors = FALSE
)

# Define indel insertions
indel_insertions <- c(
  "A", "G", "T", "C", "AT", "CG", "ATG", "CCAT", "GGATA", "ACGTAC"
)

# Use purrr to add indel insertions to the last 10 SNPs
library(purrr)
example_annotations$Allele.B[11:20] <- indel_insertions

# Save the datasets
usethis::use_data(example_genotypes, overwrite = TRUE)
usethis::use_data(example_annotations, overwrite = TRUE)

# Now create additional data sets derived from these
# (Uncomment these after implementing the functions)

# # Generate conversion results
# example_results <- genotype2allele(example_genotypes, example_annotations)
# usethis::use_data(example_results, overwrite = TRUE)
# 
# # Generate TYP format example
# example_typ_format <- allele2typ(example_results)
# usethis::use_data(example_typ_format, overwrite = TRUE)
