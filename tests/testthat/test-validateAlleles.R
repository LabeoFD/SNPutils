test_that("validateAlleles identifies AT/CG SNP order violations", {
  # Create test data with AT SNP violation
  allele_lookup <- data.frame(
    probeset_id = c("SNP1", "SNP2", "SNP3"),
    allele_a = c("T", "C", "A"),  # SNP1 has T > A which is wrong order
    allele_b = c("A", "G", "G"),
    stringsAsFactors = FALSE
  )
  
  # Run validation
  result <- validateAlleles(allele_lookup, verbose = FALSE)
  
  # Check rule violations
  violations <- result$rule_violations
  expect_equal(nrow(violations), 1)  # Should find 1 violation
  expect_equal(violations$probeset_id, "SNP1")
  expect_true(grepl("Rule 2 violation", violations$violations[1]))
})

test_that("validateAlleles identifies non-AT/CG SNP violations", {
  # Create test data with non-AT/CG SNP violation
  allele_lookup <- data.frame(
    probeset_id = c("SNP1", "SNP2", "SNP3"),
    allele_a = c("A", "G", "C"),  # SNP2 has G as allele_a, should be A or T
    allele_b = c("C", "A", "T"),  # SNP3 has T as allele_b, should be C or G
    stringsAsFactors = FALSE
  )
  
  # Run validation
  result <- validateAlleles(allele_lookup, verbose = FALSE)
  
  # Check rule violations
  violations <- result$rule_violations
  expect_equal(nrow(violations), 2)  # Should find 2 violations
  expect_equal(sort(violations$probeset_id), c("SNP2", "SNP3"))
  expect_true(all(grepl("Rule 3 violation", violations$violations)))
})

test_that("validateAlleles identifies indel violations", {
  # Create test data with indel violation
  allele_lookup <- data.frame(
    probeset_id = c("INS1", "INS2"),
    allele_a = c("-", "-"),
    allele_b = c("A", "-"),  # INS2 has "-" as allele_b, which is wrong
    stringsAsFactors = FALSE
  )
  
  # Run validation
  result <- validateAlleles(allele_lookup, verbose = FALSE)
  
  # Check rule violations
  violations <- result$rule_violations
  expect_equal(nrow(violations), 1)  # Should find 1 violation
  expect_equal(violations$probeset_id, "INS2")
  expect_true(grepl("Rule 4 violation", violations$violations[1]))
})

test_that("validateAlleles identifies multi-base order violations", {
  # Create test data with multi-base order violation
  allele_lookup <- data.frame(
    probeset_id = c("MB1", "MB2"),
    allele_a = c("ATG", "TCC"),  # MB2 has TCC > GGA which is wrong order
    allele_b = c("CCA", "GGA"),
    stringsAsFactors = FALSE
  )
  
  # Run validation
  result <- validateAlleles(allele_lookup, verbose = FALSE)
  
  # Check rule violations
  violations <- result$rule_violations
  expect_equal(nrow(violations), 1)  # Should find 1 violation
  expect_equal(violations$probeset_id, "MB2")
  expect_true(grepl("Rule 5 violation", violations$violations[1]))
})

test_that("validateAlleles returns empty violations when all rules pass", {
  # Create test data with no violations
  allele_lookup <- data.frame(
    probeset_id = c("SNP1", "SNP2", "INS1", "MB1"),
    allele_a = c("A", "A", "-", "AAG"),
    allele_b = c("G", "C", "T", "TTG"),
    stringsAsFactors = FALSE
  )
  
  # Run validation
  result <- validateAlleles(allele_lookup, verbose = FALSE)
  
  # Check no violations
  violations <- result$rule_violations
  expect_equal(nrow(violations), 0)  # Should find 0 violations
})
