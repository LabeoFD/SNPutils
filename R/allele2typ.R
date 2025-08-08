#' Convert Allele Data to TYP Format
#'
#' This function converts allele genotype data from conversion results into TYP format,
#' suitable for genetic analysis workflows. It processes SNP genotype data, extracts
#' processing information from Axiom calls files, and creates formatted output with
#' proper headers and metadata.
#'
#' @param conversion_results A list containing conversion results with a 'results' component.
#'   The 'results' component must be a data frame with columns: Sample_ID, probeset_id, 
#'   and allele_genotype. Sample_ID can contain .CEL extensions which will be removed.
#' @param calls_file_path Character string specifying the path to the AxiomGT1.calls.txt file
#'   from which processing date information will be extracted. This file must contain the
#'   Affymetrix timestamp header line.
#' @param output_dir Character string specifying the output directory path where files will be saved.
#'   Can be absolute or relative path. Required when save_file or save_log is TRUE. Default is NULL.
#' @param save_file Logical indicating whether to save the formatted output as a CSV file.
#'   When TRUE, creates a file named "IFCE_LFR9_TYPEQv1_DDMMYY_XX.csv" where DDMMYY is
#'   the processing date and XX is a sequence number. Default is FALSE.
#' @param save_log Logical indicating whether to save detailed processing logs as a text file.
#'   Creates a corresponding log file with "_log.txt" suffix. Default is FALSE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{data}{Processed tibble with columns: Sample_ID (without .CEL), SNP_Name, 
#'               Allele1_Forward, Allele2_Forward. Ready for further R processing.}
#'   \item{formatted}{Single-column tibble (V1) containing the complete formatted output 
#'                   with header section and data, suitable for file writing.}
#'   \item{header_info}{Tibble containing metadata information with Parameter and Value columns
#'                     including GSGT Version, Processing Date, Content, and sample/SNPounts.}
#'   \item{file_path}{Character string with full path to the saved CSV file, or NULL if 
#'                   save_file was FALSE or saving failed.}
#'   \item{log_content}{Character vector containing detailed processing log messages for
#'                     debugging and verification purposes.}
#' }
#'
#' @details
#' The function performs the following key operations:
#' - Validates input data structure and required columns
#' - Removes .CEL extensions from Sample_ID values
#' - Splits allele_genotype into separate Allele1_Forward and Allele2_Forward columns
#' - Extracts processing date from Axiom calls file header
#' - Creates formatted header with metadata (GSGT version, processing date, counts)
#' - Generates properly formatted CSV output with header section and data section
#' - Implements file versioning (sequence numbers) to prevent overwrites
#' - Provides comprehensive error handling and logging
#'
#' The output CSV file structure includes:
#' \itemize{
#'   \item header section with metadata as Parameter,Value pairs
#'   \item data section with SNP genotype data in proper CSV columns
#'   \item Automatic file versioning with sequence numbers (01, 02, 03)
#' }
#'
#' @section Input Requirements:
#' The conversion_results$results data frame must contain:
#' - Sample_ID: Sample identifiers (can include .CEL extension)
#' - probeset_id: SNP probe set identifiers
#' - allele_genotype: Genotype calls in format "A/T", "G/C", etc.
#'
#' The calls_file_path must point to an Axiom calls filcontaining:
#' - Header line with "#%affymetrix-algorithm-param-apt-time-str" for date extraction
#'
#' @section File Output:
#' When save_file = TRUE, creates files with naming convention:
#' - CSV: "IFCE_LFR9_TYPEQv1_DDMMYY_XX.csv"
#' - Log: "IFCE_LFR9_TYPEQv1_DDMMYY_XX_log.txt"
#' Where DDMMYY is processing date and XX is sequence number (01-03).
#'
#' @examples
#' \dontrun{
#' # Basic usage - process data without saving
#' result <- allele2typ(
#'   conversion_results = my_conversion_data,
#'   calls_file_path = "path/to/AxiomGT1.calls.txt"
#' )
#' 
#' # Access processed data
#' processed_data <- result$data
#' header_info <- result$header_info
#' 
#' # Save files to output directory
#' result <- allele2typ(
#'   conversion_results = my_conversion_data,
#'   calls_file_path = "path/to/AxiomGT1.calls.txt",
#'   output_dir = "output/genotypes",
#'   save_file = TRUE,
#'   save_log = TRUE
#' )
#' 
#' # Check file path and log
#' cat("File saved to:", result$file_path)
#' cat("Processing log:\n", paste(resu$log_content, collapse = "\n"))
#' }
#'
#' @seealso
#' \code{\link[readr]{write_csv}} for CSV writing functionality
#' \code{\link[lubridate]{parse_date_time}} for date parsing
#'
#' @importFrom dplyr relocate rename mutate select
#' @importFrom stringr str_extract str_detect str_split str_remove
#' @importFrom readr write_csv write_lines
#' @importFrom lubridate parse_date_time
#' @importFrom purrr map detect pmap_chr
#' @importFrom tibble tibble
#'
#' @export
allele2typ <- function(conversion_results, calls_file_path, output_dir = NULL, save_file = FALSE, save_log = FALSE) {
  
  # Initialize log content
  log_content <- character(0)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  log_content <- c(log_content, paste("=== allele2typ function started at", timestamp, "==="))
  
  # Error handling wrapper
  tryCatch({
    
    # Input validation
    if (is.null(conversion_results)) {
      stop("conversion_results cannot be NULL")
    }
    
    if (is.null(conversion_results$results)) {
      stop("conversion_results$results is missing or NULL")
    }
    
    if (nrow(conversion_results$results) == 0) {
      stop("conversion_results$results is empty")
    }
    
    if (is.null(calls_file_path) || calls_file_path == "") {
      stop("calls_file_path cannot be NULL or empty")
    }
    
    if (!file.exists(calls_file_path)) {
      stop("calls_file_path does not exist: ", calls_file_path)
    }
    
    if ((save_file || save_log) && is.null(output_dir)) {
      stop("output_dir must be specified when save_file or save_log is TRUE")
    }
    
    log_content <- c(log_content, "Input validation completed")
    
    # Check required columns
    required_cols <- c("Sample_ID", "probeset_id", "allele_genotype")
    missing_cols <- setdiff(required_cols, names(conversion_results$results))
    
    if (length(missing_cols) > 0) {
      stop("Missing required columns in conversion_results$results: ", paste(missing_cols, collapse = ", "))
    }
    
    log_content <- c(log_content, "Column validation completed")
    
    # Original data processing
    tryCatch({
      pre_typ <- conversion_results$results %>% 
        dplyr::relocate(Sample_ID, probeset_id, allele_genotype) %>% 
        dplyr::rename(SNP_Name = probeset_id) %>% 
        dplyr::mutate(
          Sample_ID = stringr::str_remove(Sample_ID, "\\.CEL$"),  # Remove .CEL extension
          Allele1_Forward = stringr::str_extract(allele_genotype, "^[^/]+"),
          Allele2_Forward = stringr::str_extract(allele_genotype, "[^/]+$")
        ) %>% 
        dplyr::select(-allele_genotype)
      
      log_content <- c(log_content, "Data processing completed")
      
    }, error = function(e) {
      log_content <<- c(log_content, paste("ERROR during data processing:", e$message))
      stop("Data processing failed: ", e$message)
    })
    
    # Function to extract date from AxiomGT1.calls.txt
    extract_processing_date <- function(file_path) {
      log_content <<- c(log_content, paste("Extracting date from:", file_path))
      
      tryCatch({
        # Read the file and find the line with the timestamp
        lines <- readLines(file_path, n = 50)  # Read first 50 lines to find the header
        
        date_lines <- lines[str_detect(lines, "#%affymetrix-algorithm-param-apt-time-str")]
        
        if (length(date_lines) == 0) {
          log_content <<- c(log_content, "ERROR: No date found in file")
          stop("No date found in file: ", file_path)
        }
        
        date_string <- str_split(date_lines[1], "\\s*=\\s*", simplify = TRUE)[1, 2]
        
        if (is.na(date_string) || date_string == "") {
          log_content <<- c(log_content, "ERROR: Could not extract date string from line")
          stop("Could not extract date string from line")
        }
        
        date_formats <- c("b! d! H!:M!:S! Y!", "mdy HMS", "ymd HMS", "dmy HMS")
        
        parsed_date <- date_formats %>%
          map(~ parse_date_time(trimws(date_string), orders = .x, quiet = TRUE)) %>%
          detect(~ !is.na(.x))
        
        if (is.null(parsed_date)) {
          log_content <<- c(log_content, paste("ERROR: Could not parse date:", date_string))
          stop("Could not parse date: ", date_string)
        }
        
        # Format for header output (MM/dd/yyyy HH:mm AM/PM)
        formatted_date <- format(parsed_date, "%m/%d/%Y %H:%M %p")
        log_content <<- c(log_content, paste("Date extracted:", formatted_date))
        
        return(list(formatted = formatted_date, parsed = parsed_date))
        
      }, error = function(e) {
        log_content <<- c(log_content, paste("ERROR in date extraction:", e$message))
        stop("Date extraction failed: ", e$message)
      })
    }
    
    # Extract processing date
    date_result <- extract_processing_date(calls_file_path)
    processing_date <- date_result$formatted
    parsed_date <- date_result$parsed
    
    # Get date for filename (ddmmyy format)
    date_for_filename <- format(parsed_date, "%d%m%y")
    
    # Calculate statistics
    tryCatch({
      num_snps <- length(unique(pre_typ$SNP_Name))
      num_samples <- length(unique(pre_typ$Sample_ID))
      
      if (num_snps == 0) {
        warning("No unique SNPs found in data")
      }
      
      if (num_samples == 0) {
        warning("No unique samples found in data")
      }
      
      log_content <- c(log_content, paste("Number of unique SNPs:", num_snps))
      log_content <- c(log_content, paste("Number of unique samples:", num_samples))
      
    }, error = function(e) {
      log_content <<- c(log_content, paste("ERROR calculating number of unique SNPs or samples:", e$message))
      stop("Calculation of unique number of SNPs or samples failed: ", e$message)
    })
    
    # Create header information as a proper tibble
    tryCatch({
      # Create metadata tibble with proper columns
      header_info <- tibble(
        Parameter = c("GSGT Version", "Processing Date", "Content", "Num SNPs", 
                     "Total SNPs", "Num Samples", "Total Samples"),
        Value = c("2.0000", processing_date, "Axiom_Genofil1.r1", 
                 as.character(num_snps), as.character(num_snps), 
                 as.character(num_samples), as.character(num_samples))
      )
      
      log_content <- c(log_content, "Header information created")
      
    }, error = function(e) {
      log_content <<- c(log_content, paste("ERROR creating header:", e$message))
      stop("Header creation failed: ", e$message)
    })
    
    # Helper functis for saving files
    get_sequence_number <- function(date, output_dir) {
      tryCatch({
        dir_path <- normalizePath(output_dir, mustWork = FALSE)
        
        if (!dir.exists(dir_path)) {
          return(1)
        }
        
        pattern <- sprintf("IFCE_LFR9_TYPEQv1_%s_\\d{2}\\.csv", date)
        existing_files <- list.files(dir_path, pattern = pattern)
        
        if (length(existing_files) == 0) {
          return(1)
        }
        
        sequences <- existing_files %>%
          str_extract("_\\d{2}\\.csv$") %>%
          str_extract("\\d{2}") %>%
          as.numeric() %>%
          sort()
        
        next_seq <- max(sequences, na.rm = TRUE) + 1
        if (next_seq > 3) {
          cat("Max sequence reached, using sequence 3\n")
          return(3)
        }
        
        return(next_seq)
        
      }, error = function(e) {
        log_content <<- c(log_content, paste("ERROR getting sequence number:", e$message))
        return(1)
      })
    }
    
    save_typ_file <- function(data, header_info, date, output_dir) {
      tryCatch({
        sequence_num <- get_sequence_number(date, output_dir)
        filename <- sprintf("IFCE_LFR9_TYPEQv1_%s_%02d.csv", date, sequence_num)
        
        dir_path <- normalizePath(output_dir, mustWork = FALSE)
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
          cat("Created directory:", dir_path, "\n")
          log_content <<- c(log_content, paste("Created directory:", dir_path))
        }
        
        full_path <- file.path(dir_path, filename)
        
        # Create a connection to write the file properly
        con <- file(full_path, "w")
        
        tryCatch({
          # Write header section
          writeLines("[Header]", con)
          
          # Write metadata as Parameter,Value pairs using map
          header_info %>%
            purrr::pmap_chr(~ paste(.x, .y, sep = ",")) %>%
            writeLines(con)
          
          # Write data section marker
        writeLines("[Data]", con)
          
          # Close connection to flush header content
          close(con)
          
          # Append the actual data using write_csv (proper CSV formatting)
          write_csv(data, full_path, append = TRUE, col_names = TRUE)
          
          cat("File saved:", filename, "\n")
          cat("Path:", full_path, "\n")
          
          # Verify file exists
          if (file.exists(full_path)) {
            file_info <- file.info(full_path)
            cat("File size:", round(file_info$size / 1024, 2), "KB\n")
            cat("File verification: PASSED [OK]\n")
          } else {
            stop("File was not created")
          }
          
          log_content <<- c(log_content, paste("File saved:", filename))
          log_content <<- c(log_content, paste("Path:", full_path))
          
          return(full_path)
          
        }, error = function(e) {
          # Make sure connection is closed even if error occurs
          if (exists("con") && isOpen(con)) {
            close(con)
          }
          stop("Error writing file: ", e$message)
        })
        
      }, error = function(e) {
        log_content <<- c(log_content, paste("ERROR saving file:", e$message))
        warning("Failed to save file: ", e$message)
        return(NULL)
      })
    }
    
    save_log_file <- function(log_content, date, output_dir) {
      tryCatch({
        sequence_num <- get_sequence_number(date, output_dir)
        log_filename <- sprintf("IFCE_LFR9_TYPEQv1_%s_%02d_log.txt", date, sequence_num)
        
        dir_path <- normalizePath(output_dir, mustWork = FALSE)
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
        }
        
        log_path <- file.path(dir_path, log_filename)
        writeLines(log_content, log_path)
        cat("Log file saved:", log_filename, "\n")
        
      }, error = function(e) {
        cat("Failed to save log file:", e$message, "\n")
        warning("Log file could not be saved: ", e$message)
      })
    }
    
    # Create formatted output for return (combines header info and data)
    # This creates a single-column format similar to your original for compatibility
    formatted_output <- tibble(
      V1 = c(
        "[Header]",
        paste(header_info$Parameter, header_info$Value, sep = ","),
        "[Data]",
        paste(names(pre_typ), collapse = ","),  # Column headers
        map(pre_typ, 1, function(row) paste(row, collapse = ","))  # Data rows
      )
    )
    
    # Save files if requested
    file_path <- NULL
    if (save_file && !is.null(output_dir)) {
      file_path <- save_typ_file(pre_typ, header_info, date_for_filename, output_dir)
    }
    
    # Finalize log
    log_content <- c(log_content, paste("=== allele2typ function completed at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==="))
    
    # Save log file if requested
    if (save_log && !is.null(output_dir)) {
      save_log_file(log_content, date_for_filename, output_dir)
    }
    
    # Return list with both original data and formatted output
    return(list(
      data = pre_typ,              # Original tibble format for further R processing
      formatted = formatted_output, # Formatted output with header (single column format)
      header_info = header_info,   # Header information as separate tibble
      file_path = file_path,       # Path to saved file (if saved)
      log_content = log_content    # Log content for inspection
    ))
    
  }, error = function(e) {
    # Main error handler
    error_msg <- paste("FATAL ERROR in allele2typ:", e$message)
    log_content <- c(log_content, error_msg)
    log_content <- c(log_content, paste("=== Function failed at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "==="))
    
    # Try to save error log if possible
    if (save_log && !is.null(output_dir)) {
      tryCatch({
        timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
        error_dir <- normalizePath(output_dir, mustWork = FALSE)
        error_log_path <- file.path(error_dir, paste0("allele2typ_ERROR_", timestamp_str, ".log"))
        writeLines(log_content, error_log_path)
        cat("Error log saved to:", error_log_path, "\n")
      }, error = function(log_err) {
        cat("Could not save error log:", log_err$message, "\n")
      })
    }
    
    # Print error log to console
    cat("=== ERROR LOG ===\n")
    cat(paste(log_content, collapse = "\n"), "\n")
    cat("================\n")
    
    # Re-throw the error
    stop(error_msg)
  })
}
