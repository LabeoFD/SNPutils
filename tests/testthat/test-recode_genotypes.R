library(testthat)
library(dplyr)

# Helper function to create test data
create_test_data <- function() {
  data.frame(
    probeset_id = c("SNP1", "SNP2", "SNP3", "SNP4"),
    sample1 = c(-1, 0, 1, 2),
    sample2 = c(0, 1, 2, -1),
    sample3 = c(2, -1, 0, 1),
    stringsAsFactors = FALSE
  )
}

test_that("basic recoding works", {
  test_data <- create_test_data()
  result <- recode_genotypes(test_data)
  
  # Check ID column unchanged
  expect_equal(result$probeset_id, test_data$probeset_id)
  
  # Check recoding
  expect_equal(result$sample1[1], "NoCall")  # -1 -> NoCall
  expect_equal(result$sample1[2], "AA")      # 0 -> AA
  expect_equal(result$sample1[3], "AB")      # 1 -> AB
  expect_equal(result$sample1[4], "BB")      # 2 -> BB
  
  # All genotype columns should be character
  expect_true(all(sapply(result[, -1], is.character)))
})

test_that("handles all standard values", {
  test_data <- data.frame(
    probeset_id = "SNP1",
    genotype = c(-1, 0, 1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- recode_genotypes(test_data)
  expect_equal(result$genotype, c("NoCall", "AA", "AB", "BB"))
})

test_that("works with custom ID column", {
  test_data <- data.frame(
    snp_id = c("SNP1", "SNP2"),
    sample1 = c(-1, 0),
    sample2 = c(1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- recode_genotypes(test_data, id_column = "snp_id")
  
  expect_equal(result$snp_id, test_data$snp_id)
  expect_equal(result$sample1, c("NoCall", "AA"))
  expect_equal(result$sample2, c("AB", "BB"))
})

test_that("handles NA values", {
  test_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    sample1 = c(NA, 0),
    sample2 = c(1, NA),
    stringsAsFactors = FALSE
  )
  
  result <- recode_genotypes(test_data)
  
  expect_true(is.na(result$sample1[1]))
  expect_true(is.na(result$sample2[2]))
  expect_equal(result$sample1[2], "AA")
  expect_equal(result$sample2[1], "AB")
})

test_that("catches invalid inputs", {
  # Invalid numeric values
  bad_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    sample1 = c(0, 3),    # 3 is invalid
    sample2 = c(-1, 5),   # 5 is invalid
    stringsAsFactors = FALSE
  )
  
  expect_error(recode_genotypes(bad_data), 
               "Found invalid genotype values")
  
  # String values
  string_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    sample1 = c("AA", "BB"),
    sample2 = c(0, 1),
    stringsAsFactors = FALSE
  )
  
  expect_error(recode_genotypes(string_data), 
               "Found invalid genotype values")
})


test_that("preserves data structure", {
  test_data <- create_test_data()
  result <- recode_genotypes(test_data)
  
  expect_equal(nrow(result), nrow(test_data))
  expect_equal(ncol(result), ncol(test_data))
  expect_equal(names(result), names(test_data))
  expect_true(is.data.frame(result))
})

test_that("works with single genotype column", {
  test_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    only_sample = c(-1, 2),
    stringsAsFactors = FALSE
  )
  
  result <- recode_genotypes(test_data)
  expect_equal(result$only_sample, c("NoCall", "BB"))
})

test_that("handles mixed data types", {
  test_data <- data.frame(
    probeset_id = c("SNP1", "SNP2"),
    sample1 = c(-1L, 0L),  # Integer
    sample2 = c(1.0, 2.0), # Numeric
    stringsAsFactors = FALSE
  )
  
  result <- recode_genotypes(test_data)
  expect_equal(result$sample1, c("NoCall", "AA"))
  expect_equal(result$sample2, c("AB", "BB"))
})
