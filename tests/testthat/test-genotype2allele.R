test_that("genotype2allele converts basic genotypes correctly", {
  # Create sample data frames for testing
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2", "SNP3"),
    Sample1 = c(0, 1, 2),
    Sample2 = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2", "SNP3"),
    Allele.A = c("A", "C", "G"),
    Allele.B = c("G", "T", "T"),
    stringsAsFactors = FALSE
  )
  
  # Run the function
  result <- genotype2allele(genotype_data, annotation_data)
  
  # Check structure of result
  expect_type(result, "list")
  expect_named(result, c("validation", "results"))
  
  # Check results dataframe
  expect_s3_class(result$results, "data.frame")
  expect_named(result$results, c("probeset_id", "Sample_ID", "allele_genotype"))
  
  # Number of rows should be number of SNPs * number of samples
  expect_equal(nrow(result$results), nrow(genotype_data) * (ncol(genotype_data) - 1))
  
  # Check specific allele assignments
  sample1_snp1 <- result$results %>% 
    dplyr::filter(probeset_id == "SNP1", Sample_ID == "Sample1")
  expect_equal(sample1_snp1$allele_genotype, "A/A")  # genotype 0 = homozygous for allele A
  
  sample2_snp1 <- result$results %>% 
    dplyr::filter(probeset_id == "SNP1", Sample_ID == "Sample2")
  expect_equal(sample2_snp1$allele_genotype, "A/G")  # genotype 1 = heterozygous
  
  sample1_snp3 <- result$results %>% 
    dplyr::filter(probeset_id == "SNP3", Sample_ID == "Sample1")
  expect_equal(sample1_snp3$allele_genotype, "T/T")  # genotype 2 = homozygous for allele B
})

test_that("genotype2allele handles indels correctly", {
  # Create sample data with indels
  genotype_data <- data.frame(
    probeset_id = c("INS1", "INS2"),
    Sample1 = c(0, 1),
    Sample2 = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("INS1", "INS2"),
    Allele.A = c("-", "-"),  # Deletion marker
    Allele.B = c("AGT", "C"),
    stringsAsFactors = FALSE
  )
  
  # Run the function
  result <- genotype2allele(genotype_data, annotation_data)
  
  # Check indel-specific encodings
  ins1_sample1 <- result$results %>% 
    dplyr::filter(probeset_id == "INS1", Sample_ID == "Sample1")
  expect_equal(ins1_sample1$allele_genotype, "D/D")  # genotype 0 for indel = D/D
  
  ins1_sample2 <- result$results %>% 
    dplyr::filter(probeset_id == "INS1", Sample_ID == "Sample2")
  expect_equal(ins1_sample2$allele_genotype, "D/I")  # genotype 1 for indel = D/I
  
  ins2_sample2 <- result$results %>% 
    dplyr::filter(probeset_id == "INS2", Sample_ID == "Sample2")
  expect_equal(ins2_sample2$allele_genotype, "I/I")  # genotype 2 for indel = I/I
})

test_that("genotype2allele handles missing/NA values correctly", {
  # Create sample data with NA values
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    Sample1 = c(NA, 1),
    Sample2 = c(1, NA),
    stringsAsFactors = FALSE
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2"),
    Allele.A = c("A", "C"),
    Allele.B = c("G", "T"),
    stringsAsFactors = FALSE
  )
  
  # Run the function
  result <- genotype2allele(genotype_data, annotation_data)
  
  # Check NA values are preserved
  na_entries <- result$results %>%
    dplyr::filter(is.na(allele_genotype))
  
  expect_equal(nrow(na_entries), 2)  # Should have 2 NA values
})

test_that("genotype2allele validation parameter works correctly", {
  # Create sample data 
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    Sample1 = c(0, 1),
    stringsAsFactors = FALSE
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2"),
    Allele.A = c("A", "C"),
    Allele.B = c("G", "T"),
    stringsAsFactors = FALSE
  )
  
  # With validation
  result_with <- genotype2allele(genotype_data, annotation_data, validate_alleles = TRUE)
  expect_false(is.null(result_with$validation$allele_validation))
  
  # Without validation
  result_without <- genotype2allele(genotype_data, annotation_data, validate_alleles = FALSE)
  expect_null(result_without$validation$allele_validation)
})

test_that("genotype2allele handles mismatched IDs gracefully", {
  # Create sample data with mismatches
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2", "SNPX"),  # SNPX doesn't exist in annotation
    Sample1 = c(0, 1, 2),
    stringsAsFactors = FALSE
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2", "SNP3"),  # SNP3 doesn't exist in genotype
    Allele.A = c("A", "C", "G"),
    Allele.B = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )
  
  # Run the function
  result <- genotype2allele(genotype_data, annotation_data)
  
  # Check results
  expect_equal(nrow(result$results), 3)  # Still 3 rows (one per genotype)
  
  # SNPX should have NA for allele_genotype
  snpx_row <- result$results %>% dplyr::filter(probeset_id == "SNPX")
  expect_true(all(is.na(snpx_row$allele_genotype)))
})
