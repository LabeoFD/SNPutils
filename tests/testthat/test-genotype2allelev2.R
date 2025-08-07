library(testthat)
library(SNPutils)
library(dplyr)

test_that("genotype2allelev2 works with basic input", {
  # Create test data
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2", "SNP3"),
    Sample1 = c(0, 1, 2),
    Sample2 = c(1, 2, 0)
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2", "SNP3"),
    Allele.A = c("A", "C", "T"),
    Allele.B = c("G", "T", "G")
  )
  
  # Run function
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Test structure
  expect_type(result, "list")
  expect_true(all(c("validation", "results") %in% names(result)))
  
  # Test results structure
  expect_s3_class(result$results, "data.frame")
  expect_true(all(c("probeset_id", "Sample_ID", "allele_genotype") %in% names(result$results)))
  
  # Test specific conversions
  # SNP1, Sample1: genotype 0 -> A/A
  expect_true("A/A" %in% result$results$allele_genotype)
  # SNP2, Sample1: genotype 1 -> C/T
  expect_true("C/T" %in% result$results$allele_genotype)
  # SNP3, Sample1: genotype 2 -> G/G
  expect_true("G/G" %in% result$results$allele_genotype)
})

test_that("genotype2allelev2 handles missing genotypes correctly", {
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    Sample1 = c(-1, 0),
    Sample2 = c(1, -1)
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2"),
    Allele.A = c("A", "C"),
    Allele.B = c("G", "T")
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Test that -1 genotypes become "-"
  missing_genotypes <- result$results %>%
  filter(probeset_id == "SNP1", Sample_ID == "Sample1") %>%
  pull(allele_genotype)
  expect_equal(missing_genotypes, "-")
  
  missing_genotypes2 <- result$results %>%
  filter(probeset_id == "SNP2", Sample_ID == "Sample2") %>%
  pull(allele_genotype)
  expect_equal(missing_genotypes2, "-")
})

test_that("genotype2allelev2 handles indels correctly", {
  # Test with deletion marker
  genotype_data <- data.frame(
    probeset_id = c("INDEL1", "INDEL2"),
    Sample1 = c(0, 1),
    Sample2 = c(1, 2)
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("INDEL1", "INDEL2"),
    Allele.A = c("-", "A"),
    Allele.B = c("ATG", "CGTAA")  # Multi-character allele
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Test indel conversions
  # INDEL1, Sample1: genotype 0 with deletion -> D/D
  indel_result1 <- result$results %>%
  filter(probeset_id == "INDEL1", Sample_ID == "Sample1") %>%
  pull(allele_genotype)
  expect_equal(indel_result1, "D/D")
  
  # INDEL1, Sample2: genotype 1 with deletion -> D/I
  indel_result2 <- result$results %>%
  filter(probeset_id == "INDEL1", Sample_ID == "Sample2") %>%
  pull(allele_genotype)
  expect_equal(indel_result2, "D/I")
  
  # INDEL2, Sample2: genotype 2 with multi-character -> I/I
  indel_result3 <- result$results %>%
  filter(probeset_id == "INDEL2", Sample_ID == "Sample2") %>%
  pull(allele_genotype)
  expect_equal(indel_result3, "I/I")
})

test_that("genotype2allelev2 handles space-to-dash conversion", {
  genotype_data <- data.frame(
    probeset_id = "SNP1",
    Sample1 = 0
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = "SNP1",
    Allele.A = " ",  # Space should be converted to dash
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Check that space was converted to dash in processing
  # Since Allele.A = " " should be converted to "-", and we have genotype 0,
  # this should be treated as an indel (D/D)
  expect_equal(result$results$allele_genotype, "D/D")
})

test_that("genotype2allelev2 handles missing annotation data", {
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    Sample1 = c(0, 1)
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = "SNP1",  # Missing SNP2
    Allele.A = "A",
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # SNP2 should have NA for allele_genotype since no annotation
  snp2_result <- result$results %>%
  filter(probeset_id == "SNP2") %>%
  pull(allele_genotype)
  expect_true(is.na(snp2_result))
})

test_that("genotype2allelev2 validation parameter works", {
  genotype_data <- data.frame(
    probeset_id = "SNP1",
    Sample1 = 0
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = "SNP1",
    Allele.A = "A",
    Allele.B = "G"
  )
  
  # Test with validation = FALSE
  result_no_val <- genotype2allelev2(genotype_data, annotation_data, 
                                     validate_alleles = FALSE, verbose = FALSE)
  expect_null(result_no_val$validation$allele_validation)
  expect_null(result_no_val$validation$rule_violations)
  
  # Test with validation = TRUE (if validateAlleles function exists)
  # This test might need to be skipped if validateAlleles is not implemented
  if (exists("validateAlleles", mode = "function")) {
    result_with_val <- genotype2allelev2(genotype_data, annotation_data, 
                                         validate_alleles = TRUE, verbose = FALSE)
    expect_type(result_with_val$validation, "list")
  }
})

test_that("genotype2allelev2 returns correct data types", {
  genotype_data <- data.frame(
    probeset_id = "SNP1",
    Sample1 = 1
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = "SNP1",
    Allele.A = "A",
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Test column types
  expect_type(result$results$probeset_id, "character")
  expect_type(result$results$Sample_ID, "character")
  expect_type(result$results$allele_genotype, "character")
  
  # Test that result has expected dimensions
  expect_equal(nrow(result$results), 1)  # 1 SNP × 1 sample = 1 row
  expect_equal(ncol(result$results), 3)  # probeset_id, Sample_ID, allele_genotype
})

test_that("genotype2allelev2 handles multiple samples correctly", {
  genotype_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    Sample1 = c(0, 1),
    Sample2 = c(1, 2),
    Sample3 = c(2, 0)
  )
  
  annotation_data <- data.frame(
    Probe.Set.ID = c("SNP1", "SNP2"),
    Allele.A = c("A", "C"),
    Allele.B = c("G", "T")
  )
  
  result <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
  
  # Should have 2 SNPs × 3 samples = 6 rows
  expect_equal(nrow(result$results), 6)
  
  # Check that all samples are represented
  expect_true(all(c("Sample1", "Sample2", "Sample3") %in% result$results$Sample_ID))
  
  # Check that all SNPs are represented
  expect_true(all(c("SNP1", "SNP2") %in% result$results$probeset_id))
})
