#' Validate TYP and COM Data Consistency
#'
#' This function validates that TYP (genotype) and COM (phenotype) datasets have
#' consistent sample IDs and that the TYP data contains the expected number of SNPs.
#' It performs comprehensive checks and provides detailed reporting of any mismatches.
#'
#' @param typ_data A list containing TYP data with a `data` component that includes
#'   `Sample_ID` and `SNP_Name` columns.
#' @param com_data A data frame containing COM data with a `Sample ID` column.
#' @param expected_snp_count Integer. The expected number of unique SNPs in the TYP data.
#'   Default is 60752.
#' @param show_mismatched_samples Logical. Whether to display samples that don't match
#'   between datasets. Default is TRUE.
#' @param max_show Integer. Maximum number of mismatched samples to display.
#'   Default is 10.
#' @param verbose Logical. Whether to print detailed validation results to console.
#'   Default is FALSE.
#' @param save_log Logical. Whether to save validation results to a log file.
#'   Default is FALSE.
#' @param log_file Character. Name of the log file when `save_log = TRUE`.
#'   Default is "typcom_validation_log.txt".
#' @param output_dir Character. Directory path where the log file should be saved.
#'   Required when `save_log = TRUE`. Default is NULL.
#'
#' @return Invisibly returns a logical value: TRUE if validation passes,
#'   FALSE if validation fails.
#'
#' @details
#' The function performs the following validations:
#' \itemize{
#'   \item Sample count consistency between TYP and COM datasets
#'   \item Complete overlap of sample IDs between datasets
#'   \item Correct number of SNPs in TYP data
#' }
#'
#' When validation fails, the function provides detailed information about:
#' \itemize{
#'   \item Sample count differences
#'   \item Samples present in one dataset but not the other
#'   \item SNP count discrepancies
#' }
#'
#' @examples
#' \dontrun{
#' # Basic validation
#' result <- validate_typcom(typ_data, com_data)
#'
#' # Verbose output with custom SNP count
#' validate_typcom(typ_data, com_data,
#'                 expected_snp_count = 50000,
#'                 verbose = TRUE)
#'
#' # Save validation log
#' validate_typcom(typ_data, com_data,
#'                 save_log = TRUE,
#'                 output_dir = "~/validation_logs")
#' }
#'
#' @importFrom dplyr distinct pull
#' @importFrom glue glue
#'
#' @export
validate_typcom <- function(typ_data, com_data, expected_snp_count = 60752, 
                           show_mismatched_samples = TRUE, max_show = 10, verbose = FALSE,
                           save_log = FALSE, log_file = "typcom_validation_log.txt", output_dir = NULL) {
  
  # Avoid R CMD check NOTEs about global variables
  Sample_ID <- SNP_Name <- `Sample ID` <- NULL
  
  # Get unique samples and SNPs using tidyverse
  typ_samples <- typ_data$data %>% 
    distinct(Sample_ID) %>% 
    pull(Sample_ID)
  
  com_samples <- com_data %>% 
    distinct(`Sample ID`) %>% 
    pull(`Sample ID`)
  
  typ_snps <- typ_data$data %>% 
    distinct(SNP_Name) %>% 
    pull(SNP_Name)
  
  # Calculate counts
  typ_sample_count <- length(typ_samples)
  com_sample_count <- length(com_samples)
  snp_count <- length(typ_snps)
  
  # Check if validation passes
  sample_ids_match <- typ_sample_count == com_sample_count && 
    all(typ_samples %in% com_samples) && 
    all(com_samples %in% typ_samples)
  snp_count_correct <- snp_count == expected_snp_count
  validation_passed <- sample_ids_match && snp_count_correct
  
  # Prepare log content
  log_content <- c()
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_content <- c(log_content, paste("=== TYP/COM Data Validation ==="))
  log_content <- c(log_content, paste("Date:", timestamp))
  log_content <- c(log_content, "")
  
  # Print and log results
  if (validation_passed) {
    main_msg <- "Validation successful PASS [OK]"
    cat(main_msg, "\n")
    log_content <- c(log_content, main_msg)
    
    if (verbose || save_log) {
      details <- c(
        "=== Validation Details ===",
        glue("  Sample count: TYP={typ_sample_count}, COM={com_sample_count} [OK]"),
        glue("  SNP count: {snp_count} (expected: {expected_snp_count}) [OK]"),
        glue("  Sample overlap: 100% ({typ_sample_count}/{typ_sample_count}) [OK]")
      )
      
      if (verbose) {
        cat(paste(details, collapse = "\n"), "\n")
      }
      
      log_content <- c(log_content, details)
    }
    
  } else {
    main_msg <- "Validation FAIL [ERROR]"
    cat(main_msg, "\n")
    log_content <- c(log_content, main_msg)
    
    if (verbose || save_log) {
      log_content <- c(log_content, "=== Validation Details ===")
      if (verbose) {
        cat("=== Validation Details ===\n")
      }
    }
    
    # Check sample ID issues
    if (!sample_ids_match) {
      if (typ_sample_count != com_sample_count) {
        symbol <- if (verbose || save_log) " [ERROR]" else ""
        msg <- glue("  Sample count mismatch: TYP={typ_sample_count}, COM={com_sample_count}, diff={typ_sample_count - com_sample_count}{symbol}")
        cat(msg, "\n")
        log_content <- c(log_content, msg)
      }
      
      if (show_mismatched_samples) {
        # Find samples in typ but not in com
        typ_only <- setdiff(typ_samples, com_samples)
        if (length(typ_only) > 0) {
          symbol <- if (verbose || save_log) " [ERROR]" else ""
          msg <- glue("  Samples in TYP but not in COM ({length(typ_only)} total){symbol}:")
          cat(msg, "\n")
          log_content <- c(log_content, msg)
          
          show_samples <- head(typ_only, max_show)
          sample_msg <- paste("   ", paste(show_samples, collapse = ", "))
          cat(sample_msg, "\n")
          log_content <- c(log_content, sample_msg)
          
          if (length(typ_only) > max_show) {
            more_msg <- glue("    ... and {length(typ_only) - max_show} more")
            cat(more_msg, "\n")
            log_content <- c(log_content, more_msg)
          }
        }
        
        # Find samples in com but not in typ
        com_only <- setdiff(com_samples, typ_samples)
        if (length(com_only) > 0) {
          symbol <- if (verbose || save_log) " [ERROR]" else ""
          msg <- glue("  Samples in COM but not in TYP ({length(com_only)} total){symbol}:")
          cat(msg, "\n")
          log_content <- c(log_content, msg)
          
          show_samples <- head(com_only, max_show)
          sample_msg <- paste("   ", paste(show_samples, collapse = ", "))
          cat(sample_msg, "\n")
          log_content <- c(log_content, sample_msg)
          
          if (length(com_only) > max_show) {
            more_msg <- glue("    ... and {length(com_only) - max_show} more")
            cat(more_msg, "\n")
            log_content <- c(log_content, more_msg)
          }
        }
      }
      
      if ((verbose || save_log) && sample_ids_match) {
        msg <- glue("  Sample overlap: 100% ({typ_sample_count}/{typ_sample_count}) [OK]")
        if (verbose) cat(msg, "\n")
        log_content <- c(log_content, msg)
      }
    }
    
    if (!snp_count_correct) {
      symbol <- if (verbose || save_log) " [ERROR]" else ""
      msg <- glue("  SNP count mismatch: found={snp_count}, expected={expected_snp_count}, diff={snp_count - expected_snp_count}{symbol}")
      cat(msg, "\n")
      log_content <- c(log_content, msg)
    } else if (verbose || save_log) {
      msg <- glue("  SNP count: {snp_count} (expected: {expected_snp_count}) [OK]")
      if (verbose) cat(msg, "\n")
      log_content <- c(log_content, msg)
    }
  }
  
  # Save log file if requested
  if (save_log) {
    # Require user to specify output directory
    if (is.null(output_dir)) {
      stop("output_dir must be specified when save_log = TRUE. Please provide a directory path (e.g., output_dir = '~/Documents')")
    }
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat(glue("Created output directory: {output_dir}"), "\n")
    }
    
    # Construct full log file path
    full_log_path <- file.path(output_dir, log_file)
    
    log_content <- c(log_content, "", paste("Log saved:", timestamp), 
                     paste(rep("=", 50), collapse = ""), "")
    
    # Append to log file
    cat(paste(log_content, collapse = "\n"), "\n", file = full_log_path, append = TRUE)
    cat(glue("Validation log saved to: {full_log_path}"), "\n")
  }
  
  invisible(validation_passed)
}
