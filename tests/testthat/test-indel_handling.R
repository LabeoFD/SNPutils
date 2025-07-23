test_that("genotype2allelev2 handles standard indels with deletion marker", {
  test_genotypes <- data.frame(
    probeset_id = "STANDARD_INDEL",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "STANDARD_INDEL",
    Chromosome = "1",
    Physical.Position = 1000,
    Allele.A = "-",      # Deletion marker
    Allele.B = "ATG"     # Insertion sequence
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 detects indels by multi-character alleles", {
  test_genotypes <- data.frame(
    probeset_id = "MULTI_CHAR_INDEL",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "MULTI_CHAR_INDEL",
    Chromosome = "2",
    Physical.Position = 2000,
    Allele.A = "A",         # Single character
    Allele.B = "CCCG"       # Multi-character = indel
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels because Allele.B has length > 1
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 handles both alleles being multi-character", {
  test_genotypes <- data.frame(
    probeset_id = "BOTH_MULTI",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "BOTH_MULTI",
    Chromosome = "3",
    Physical.Position = 3000,
    Allele.A = "ATG",       # Multi-character
    Allele.B = "CCCG"       # Multi-character
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 handles identical long sequences", {
  test_genotypes <- data.frame(
    probeset_id = "IDENTICAL_LONG",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "IDENTICAL_LONG",
    Chromosome = "4",
    Physical.Position = 4000,
    Allele.A = "ATTCTGACTACCACAACTAAACATCTATGC",
    Allele.B = "ATTCTGACTACCACAACTAAACATCTATGC"
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels because both alleles have length > 1
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})

test_that("genotype2allelev2 handles short identical multi-character alleles", {
  test_genotypes <- data.frame(
    probeset_id = "IDENTICAL_SHORT",
    Sample1 = 0,
    Sample2 = 1,
    Sample3 = 2
  )
  
  test_annotations <- data.frame(
    Probe.Set.ID = "IDENTICAL_SHORT",
    Chromosome = "5",
    Physical.Position = 5000,
    Allele.A = "CA",
    Allele.B = "CA"
  )
  
  result <- genotype2allelev2(test_genotypes, test_annotations, verbose = FALSE)
  results_df <- result$results
  
  # Should be treated as indels because both alleles have length > 1
  expect_equal(results_df[results_df$Sample_ID == "Sample1", ]$allele_genotype, "D/D")
  expect_equal(results_df[results_df$Sample_ID == "Sample2", ]$allele_genotype, "D/I")
  expect_equal(results_df[results_df$Sample_ID == "Sample3", ]$allele_genotype, "I/I")
})
