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

# Create header with metadata (as character vector)
header <- c(
  "#%chip_type=Axiom_GW_Hu_SNP",
  "#%lib_set_name=Axiom_GW_Hu_SNP.r2", 
  "#%lib_set_version=r2",
  "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
  "#%create_date=2024-01-15T10:30:45Z",
  "#%guid=1234567890abcdef",
  "#%program_name=apt-probeset-genotype",
  "#%program_version=1.21.0",
  "#%program_company=Affymetrix"
)

# Create sample data with realistic call rates
sample_data <- data.frame(
  cel_files = paste0("Sample_", sprintf("%03d", 1:15), ".CEL"),
  total_call_rate = c(97.5, 94.2, 98.1, 96.8, 93.5, 99.2, 95.7, 92.1, 
                     97.9, 96.3, 94.8, 98.5, 95.2, 93.8, 97.1),
  average_heterozygosity = c(0.230, 0.285, 0.215, 0.250, 0.305, 0.220, 0.265, 0.315,
                            0.240, 0.275, 0.295, 0.200, 0.255, 0.280, 0.235),
  median_intensity = c(2150, 1890, 2340, 2100, 1750, 2500, 2200, 1680,
                      2380, 2050, 1920, 2450, 2180, 1820, 2280),
  stringsAsFactors = FALSE
)

# Combine into example data
example_axiom_report <- list(
  header = header,
  data = sample_data
)

# Save the datasets
usethis::use_data(example_genotypes, overwrite = TRUE)
usethis::use_data(example_annotations, overwrite = TRUE)
usethis::use_data(example_axiom_report, overwrite = TRUE)

# Now create additional data sets derived from these
# (Uncomment these after implementing the functions)

# # Generate conversion results
# example_results <- genotype2allele(example_genotypes, example_annotations)
# usethis::use_data(example_results, overwrite = TRUE)
# 
# # Generate TYP format example
# example_typ_format <- allele2typ(example_results)
# usethis::use_data(example_typ_format, overwrite = TRUE)
