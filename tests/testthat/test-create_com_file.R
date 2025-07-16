test_that("create_com_file works with valid input", {
  # Create temporary test file
  test_file <- tempfile(fileext = ".txt")
  
  # Create a simple test file manually
  header <- c(
    "#%chip_type=Axiom_GW_Hu_SNP",
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate\taverage_heterozygosity",
    "Sample_001.CEL\t97.5\t0.25",
    "Sample_002.CEL\t94.2\t0.30",
    "Sample_003.CEL\t98.1\t0.22"
  )
  writeLines(header, test_file)
  
  # Test function
  result <- create_com_file(test_file)
  
  # Test structure
  expect_type(result, "list")
  expect_named(result, c("formatted_date", "data"))
  expect_s3_class(result$data, "tbl_df")
  
  # Test dimensions
  expect_equal(ncol(result$data), 3)
  expect_equal(nrow(result$data), 3)
  
  # Test column names (should have spaces)
  expect_true(all(c("Sample ID", "Call Rate", "Comment") %in% names(result$data)))
  
  # Test data types
  expect_type(result$data$`Sample ID`, "character")
  expect_type(result$data$`Call Rate`, "double")
  expect_type(result$data$Comment, "character")
  
  # Test call rate values are between 0 and 1
  expect_true(all(result$data$`Call Rate` >= 0 & result$data$`Call Rate` <= 1))
  
  # Test comment values are 0 or 1
  expect_true(all(result$data$Comment %in% c("0", "1")))
  
  # Clean up
  unlink(test_file)
})

test_that("create_com_file handles custom thresholds", {
  test_file <- tempfile(fileext = ".txt")
  
  header <- c(
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate",
    "Sample_001.CEL\t97.5",
    "Sample_002.CEL\t94.2",
    "Sample_003.CEL\t93.0"
  )
  writeLines(header, test_file)
  
  # Test with different thresholds
  result_default <- create_com_file(test_file, call_rate_threshold = 0.95)
  result_strict <- create_com_file(test_file, call_rate_threshold = 0.98)
  
  # Should have same structure
  expect_s3_class(result_strict$data, "tbl_df")
  expect_equal(ncol(result_strict$data), 3)
  
  # Stricter threshold should result in more "1" comments (low quality)
  strict_failures <- sum(result_strict$data$Comment == "1")
  default_failures <- sum(result_default$data$Comment == "1")
  expect_gte(strict_failures, default_failures)
  
  unlink(test_file)
})

test_that("create_com_file validates inputs correctly", {
  # Test non-existent file
  expect_error(
    create_com_file("nonexistent_file.txt"),
    "File does not exist"
  )
  
  # Test invalid threshold - too high
  test_file <- tempfile(fileext = ".txt")
  header <- c(
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate",
    "Sample_001.CEL\t97.5"
  )
  writeLines(header, test_file)
  
  expect_error(
    create_com_file(test_file, call_rate_threshold = 1.5),
    "call_rate_threshold must be a numeric value between 0 and 1"
  )
  
  # Test invalid threshold - negative
  expect_error(
    create_com_file(test_file, call_rate_threshold = -0.1),
    "call_rate_threshold must be a numeric value between 0 and 1"
  )
  
  unlink(test_file)
})

test_that("date parsing works correctly", {
  test_file <- tempfile(fileext = ".txt")
  
  # Test with specific date
  header <- c(
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate",
    "Sample_001.CEL\t97.5"
  )
  writeLines(header, test_file)
  
  result <- create_com_file(test_file)
  
  # Check date format (should be DDMMYY)
  expect_type(result$formatted_date, "character")
  expect_equal(nchar(result$formatted_date), 6)
  expect_equal(result$formatted_date, "230525")
  
  unlink(test_file)
})

test_that("call rate conversion is correct", {
  test_file <- tempfile(fileext = ".txt")
  
  header <- c(
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate",
    "Sample_001.CEL\t97.5",
    "Sample_002.CEL\t94.0"
  )
  writeLines(header, test_file)
  
  result <- create_com_file(test_file)
  
  # Check that call rates are properly converted (divided by 100)
  expected_rates <- c(0.975, 0.940)
  expect_equal(result$data$`Call Rate`, expected_rates, tolerance = 1e-10)
  
  unlink(test_file)
})

test_that("comment assignment works correctly", {
  test_file <- tempfile(fileext = ".txt")
  
  header <- c(
    "#%affymetrix-algorithm-param-apt-time-str=May 23 14:30:45 2025",
    "cel_files\ttotal_call_rate",
    "Sample_001.CEL\t97.0",  # Should be "0" (>= 0.95)
    "Sample_002.CEL\t94.0"   # Should be "1" (< 0.95)
  )
  writeLines(header, test_file)
  
  result <- create_com_file(test_file)
  
  # Checcomment assignment
  expect_equal(result$data$Comment, c("0", "1"))
  
  unlink(test_file)
})
