test_that("genotype2allelev2 handles missing genotypes (-1)", {
  test_genotypes <- data.frame(
    probeset_id = "MISSING_GENO",
    Sample1 = -1,
    Sample2 = 0,
    Sample3 = -1
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "MISSING_GENO",
    Chromosome = "1",
    Physical.Position = 1000,
    Allele.A = "A",
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Missing genotypes (-1) should become "-"
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "-")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "-")
  
  # Normal genotype should still work
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "A/A")
})

test_that("genotype2allelev2 handles space-to-dash conversion in Allele.A", {
  test_genotypes <- data.frame(
    probeset_id = "SPACE_CONVERSION",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "SPACE_CONVERSION",
    Chromosome = "1",
    Physical.Position = 1000,
    Allele.A = " ",      # Space should be converted to "-"
    Allele.B = "ATG"     # Regular insertion
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels because Allele.A becomes "-"
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 handles both alleles as spaces", {
  test_genotypes <- data.frame(
    probeset_id = "BOTH_SPACES",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "BOTH_SPACES",
    Chromosome = "2",
    Physical.Position = 2000,
    Allele.A = " ",      # Space → "-"
    Allele.B = " "       # Space stays as space (no conversion for Allele.B)
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels because Allele.A becomes "-"
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 distinguishes standard SNPs from multi-character variants", {
  test_genotypes <- data.frame(
    probeset_id = c("STANDARD", "MULTI_CHAR"),
    Sample1 = c(2, 2),
    Sample2 = c(1, 1)
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = c("STANDARD", "MULTI_CHAR"),
    Chromosome = c("1", "2"),
    Physical.Position = c(1000, 2000),
    Allele.A = c("A", "ATG"),        # Single vs multi-character
    Allele.B = c("G", "CCCG")        # Single vs multi-character
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Standard SNP genotype 2 should be B/B (G/G)
  standard_result <- results_df[results_df$probeset_id == "STANDARD" & 
                               results_df$Sample_ID == "Sample1", ]$allele_genotype
  expect_equal(standard_result, "G/G")
  
  # Multi-character genotype 2 should be I/I
  multi_result <- results_df[results_df$probeset_id == "MULTI_CHAR" & 
                            results_df$Sample_ID == "Sample1", ]$allele_genotype
  expect_equal(multi_result, "I/I")
})

test_that("genotype2allelev2 handles empty or unusual input gracefully", {
  # Test with minimal data
  test_genotypes <- data.frame(
    probeset_id = "MINIMAL",
    Sample1 = 0
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "MINIMAL",
    Chromosome = "1",
    Physical.Position = 1000,
    Allele.A = "A",
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  
  # Should handle minimal input without error
  expect_type(result, "list")
  expect_equal(nrow(result$results), 1)
  expect_equal(result$results$allele_genotype, "A/A")
})

test_that("genotype2allelev2 handles genotype codes outside normal range", {
  test_genotypes <- data.frame(
    probeset_id = "UNUSUAL_CODES",
    Sample1 = 3,      # Outside normal 0,1,2 range
    Sample2 = 99,     # Very unusual code
    Sample3 = 0       # Normal code for comparison
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "UNUSUAL_CODES",
    Chromosome = "1",
    Physical.Position = 1000,
    Allele.A = "A",
    Allele.B = "G"
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Unusual codes should result in NA (handled by case_when TRUE condition)
  expect_true(is.na(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype))
  expect_true(is.na(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype))
  
  # Normal code should still work
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "A/A")
})
