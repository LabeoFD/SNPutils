test_that("create_com_file works with valid input", {
  # Create test data
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "#%other-header-info",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5",
    "sample2.cel\t94.0",
    "sample3.cel\t98.2"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  # Test the function (with save_file = FALSE to avoid file creation)
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # Test structure
  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 4)
  expect_equal(nrow(result), 3)
  
  # Test column names
  expect_true(all(c("Sample ID", "Call Rate", "Comment", "Status") %in% names(result)))
  
  # Test column types
  expect_type(result$`Sample ID`, "character")
  expect_type(result$`Call Rate`, "double")
  expect_type(result$Comment, "character")
  expect_type(result$Status, "character")
  
  # Test data values
  expect_true(all(result$`Call Rate` >= 0 & result$`Call Rate` <= 1))
  expect_true(all(result$Comment %in% c("0", "1")))
  expect_true(all(result$Status %in% c("PASS", "FAIL")))
  
  # Test sample ID extraction (should remove .cel extension)
  # Ordered by call rate desc: sample3 (98.2), sample1 (97.5), sample2 (94.0)
  expect_equal(result$`Sample ID`, c("sample3", "sample1", "sample2"))
  
  # Clean up
  unlink(test_file)
})

test_that("create_com_file handles custom thresholds", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5",
    "sample2.cel\t94.0"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  # Test with strict threshold
  result_strict <- create_com_file(test_file, call_rate_threshold = 0.98, save_file = FALSE, save_log = FALSE)
  result_default <- create_com_file(test_file, call_rate_threshold = 0.95, save_file = FALSE, save_log = FALSE)
  
  expect_s3_class(result_strict, "tbl_df")
  expect_s3_class(result_default, "tbl_df")
  
  # With 0.98 threshold, only no samples should pass (both are below 0.98)
  # With 0.95 threshold, sample1 should pass
  strict_pass_count <- sum(result_strict$Status == "PASS")
  default_pass_count <- sum(result_default$Status == "PASS")
  
  expect_lte(strict_pass_count, default_pass_count)
  
  unlink(test_file)
})

test_that("create_com_file validates inputs correctly", {
  # Test non-existent file
  expect_error(
    create_com_file("nonexistent_file.txt"),
    "File does not exist"
  )
  
  # Test invalid threshold
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5"
  )
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  expect_error(
    create_com_file(test_file, call_rate_threshold = 1.5),
    "Call rate threshold must be between 0 and 1"
  )
  
  expect_error(
    create_com_file(test_file, call_rate_threshold = -0.1),
    "Call rate threshold must be between 0 and 1"
  )
  
  unlink(test_file)
})

test_that("date parsing works correctly", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # The function processes the data but doesn't return the date
  # We can only test that it runs without error
  expect_s3_class(result, "tbl_df")
  
  unlink(test_file)
})

test_that("call rate conversion is correct", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5",    # Should be converted to 0.975
    "sample2.cel\t0.94"     # Should stay as 0.94
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  expected_rates <- c(0.975, 0.94)  # Ordered by desc call rate
  expect_equal(result$`Call Rate`, expected_rates)
  
  unlink(test_file)
})

test_that("comment assignment works correctly", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5",    # Above 0.95 threshold -> Comment "0"
    "sample2.cel\t92.0"     # Below 0.95 threshold -> Comment "1"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # Check comments (ordered by desc call rate: sample1 first, then sample2)
  expect_equal(result$Comment, c("0", "1"))
  expect_equal(result$Status, c("PASS", "FAIL"))
  
  unlink(test_file)
})

test_that("missing required columns triggers error", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "wrong_column\tother_column",
    "sample1.cel\t97.5"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  expect_error(
    create_com_file(test_file, save_file = FALSE, save_log = FALSE),
    "Missing columns"
  )
  
  unlink(test_file)
})

test_that("empty data file triggers error", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "#%only-header-lines"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  expect_error(
    create_com_file(test_file, save_file = FALSE, save_log = FALSE),
    "No data found in file"
  )
  
  unlink(test_file)
})

test_that("sample ID extraction removes file extensions and paths", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "/path/to/sample1.cel\t97.5",
    "C:\\Windows\\Path\\sample2.CEL\t94.0",
    "sample3.cel\t96.0"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # Should extract just the sample names without paths or extensions
  # Ordered by call rate desc: sample1 (97.5), sample3 (96.0), sample2 (94.0) 
  expect_equal(result$`Sample ID`, c("sample1", "sample3", "sample2"))
  
  # Verify no .cel or .CEL extensions remain
  expect_true(all(!grepl("\\.(cel|CEL)$", result$`Sample ID`)))
  
  # Verify no paths remain
  expect_true(all(!grepl("[\\/\\\\]", result$`Sample ID`)))
  
  unlink(test_file)
})

test_that("percentage call rates are converted correctly", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "sample1.cel\t97.5",     # > 1, should be divided by 100
    "sample2.cel\t0.94",     # <= 1, should stay as is
    "sample3.cel\t150.0"     # > 1, should be divided by 100
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # Check conversion: 97.5 -> 0.975, 0.94 -> 0.94, 150.0 -> 1.5 (but capped at reasonable values in practice)
  # Ordered by desc: 1.5, 0.975, 0.94
  expected_call_rates <- c(1.5, 0.975, 0.94)
  expect_equal(result$`Call Rate`, expected_call_rates)
  
  unlink(test_file)
})

test_that("data is ordered by call rate descending", {
  test_content <- c(
    "#%affymetrix-algorithm-param-apt-time-str=Thu May 25 14:30:22 2023",
    "cel_files\ttotal_call_rate",
    "low_sample.cel\t90.0",
    "high_sample.cel\t98.0",
    "mid_sample.cel\t95.0"
  )
  
  test_file <- tempfile(fileext = ".txt")
  writeLines(test_content, test_file)
  
  result <- create_com_file(test_file, save_file = FALSE, save_log = FALSE)
  
  # Should be ordered by call rate descending
  expect_equal(result$`Sample ID`, c("high_sample", "mid_sample", "low_sample"))
  expect_equal(result$`Call Rate`, c(0.98, 0.95, 0.90))
  
  unlink(test_file)
})
