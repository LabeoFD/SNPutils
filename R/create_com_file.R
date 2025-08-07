#' Create COM File from Axiom GT Report
#'
#' Reads and processes an Axiom GT report file, extracting date information
#' and sample call rate data with quality comments and pass/fail status.
#'
#' @param file_path Character string. Path to the AxiomGT report file (.txt)
#' @param call_rate_threshold Numeric. Threshold for call rate quality (default: 0.95)
#' @param save_file Logical. Whether to save the COM file with PASS samples only (default: TRUE)
#' @param output_dir Character string. Output directory (default: "Documents/TYPE")
#' @param save_log Logical. Whether to save processing log file (default: TRUE)
#'
#' @return A tibble containing all processed sample data
#'
#' @examples
#' \dontrun{
#' result <- create_com_file("path/to/AxiomGT1.report.txt")
#' result <- create_com_file("file.txt", call_rate_threshold = 0.98, save_file = FALSE)
#' }
#'
#' @export
#' @importFrom dplyr mutate select filter arrange case_when across transmute if_else desc all_of
#' @importFrom readr read_csv write_csv read_delim
#' @importFrom stringr str_detect str_replace str_extract str_split str_trim str_c str_remove str_match
#' @importFrom purrr detect map
#' @importFrom crayon red green blue bold yellow
#' @importFrom here here
#' @importFrom lubridate ymd now today parse_date_time

create_com_file <- function(file_path, 
                            call_rate_threshold = 0.95, 
                            save_file = TRUE,
                            output_dir = "Documents/TYP-COM",
                            save_log = TRUE) {
  
  log_content <- character(0)
  
  if (save_log) {
    log_content <- c(
      "Axiom GT Report Processing Log",
      paste("Generated:", Sys.time()),
      paste("File:", basename(file_path)),
      paste("Threshold:", call_rate_threshold),
      ""
    )
  }
  
  validate_inputs(file_path, call_rate_threshold, output_dir)
  
  cat("Processing file:", green(basename(file_path)), "\n")
  if (save_log) {
    log_content <- c(log_content, paste("Processing file:", basename(file_path)))
  }
  
  lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  
  date_extracted <- extract_date(lines, save_log, log_content)
  if (save_log) {
    log_content <- date_extracted$log_content
  }
  
  processed_data <- process_data(lines, call_rate_threshold)

  summary_info <- display_summary(processed_data, call_rate_threshold, save_log)
  if (save_log) {
    log_content <- c(log_content, summary_info$log_content)
  }
  
  if (save_file && !is.null(date_extracted$date)) {
    file_path_saved <- save_com_file(processed_data, date_extracted$date, output_dir)
    if (save_log) {
      log_content <- c(log_content, paste("COM file saved:", file_path_saved))
    }
  }
  
  if (save_log && !is.null(date_extracted$date)) {
    save_log_file(log_content, date_extracted$date, output_dir)
  }
  
  return(processed_data)
}

validate_inputs <- function(file_path, call_rate_threshold, output_dir) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  if (!is.numeric(call_rate_threshold) || call_rate_threshold < 0 || call_rate_threshold > 1) {
    stop("Call rate threshold must be between 0 and 1")
  }
  
  if (!is.character(output_dir)) {
    stop("Output directory must be a character string")
  }
}

extract_date <- function(lines, save_log = FALSE, log_content = character(0)) {
  date_lines <- lines[str_detect(lines, "#%affymetrix-algorithm-param-apt-time-str")]
  
  if (length(date_lines) == 0) {
    cat(red("No date found in file\n"))
    if (save_log) {
      log_content <- c(log_content, "ERROR: No date found in file")
    }
    return(list(date = NULL, log_content = log_content))
  }
  
  date_string <- str_split(date_lines[1], "\\s*=\\s*", simplify = TRUE)[1, 2]
  
  date_formats <- c("b! d! H!:M!:S! Y!", "mdy HMS", "ymd HMS", "dmy HMS")
  
  parsed_date <- date_formats %>%
    map(~ parse_date_time(trimws(date_string), orders = .x, quiet = TRUE)) %>%
    detect(~ !is.na(.x))
  
  if (is.null(parsed_date)) {
    cat(red("Could not parse date:", date_string, "\n"))
    if (save_log) {
      log_content <- c(log_content, paste("ERROR: Could not parse date:", date_string))
    }
    return(list(date = NULL, log_content = log_content))
  }
  
  formatted_date <- format(parsed_date, "%d%m%y")
  cat("Date extracted:", green(formatted_date), "\n")
  
  if (save_log) {
    log_content <- c(log_content, paste("Date extracted:", formatted_date))
  }
  
  return(list(date = formatted_date, log_content = log_content))
}

process_data <- function(lines, call_rate_threshold) {
  data_lines <- lines[!str_detect(lines, "^#") & nchar(trimws(lines)) > 0]
  
  if (length(data_lines) <= 1) {
    stop("No data found in file")
  }
  
  data_df <- read_delim(
    I(paste(data_lines, collapse = "\n")), 
    delim = "\t", 
    col_names = TRUE, 
    show_col_types = FALSE
  )
  
  required_cols <- c("cel_files", "total_call_rate")
  missing_cols <- setdiff(required_cols, names(data_df))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable: ", paste(names(data_df), collapse = ", "))
  }
  
  data_df %>%
    transmute(
      `Sample ID` = str_remove(cel_files, "\\.(cel|CEL)$") %>% 
        str_remove("^.*[\\/\\\\]"),
      `Call Rate` = case_when(
        total_call_rate > 1 ~ round(total_call_rate / 100, 6),
        TRUE ~ round(as.numeric(total_call_rate), 6)
      ),
      Comment = if_else(`Call Rate` >= call_rate_threshold, "0", "1"),
      Status = if_else(`Call Rate` >= call_rate_threshold, "PASS", "FAIL")
    ) %>%
    arrange(desc(`Call Rate`))
}

display_summary <- function(data, threshold, save_log = FALSE) {
  total <- nrow(data)
  passing <- sum(data$Status == "PASS")
  failing <- total - passing
  
  passing_data <- data %>% filter(Status == "PASS")
  failing_data <- data %>% filter(Status == "FAIL")
  
  passing_avg <- if (nrow(passing_data) > 0) {
    mean(passing_data$`Call Rate`, na.rm = TRUE)
  } else {
    0
  }
  
  qc_threshold <- 0.985
  qc_pass <- passing_avg >= qc_threshold
  
  log_content <- character(0)
  if (save_log) {
    log_content <- c(
      "",
      "QC CALL RATE SUMMARY:",
      "=============================",
      paste("Call rate threshold:", threshold),
      paste("Total samples:", total),
      paste("Passing Threshold:", passing, "(", round(100 * passing / total, 1), "%)"),
      paste("Failing Threshold:", failing, "(", round(100 * failing / total, 1), "%)"),
      paste("Overall mean call rate:", round(mean(data$`Call Rate`, na.rm = TRUE), 3)),
      paste("Overall range:", 
            round(min(data$`Call Rate`, na.rm = TRUE), 3), "-", 
            round(max(data$`Call Rate`, na.rm = TRUE), 3))
    )
  }
  
  cat("\n", bold("QC Call Rate Summary"), "\n")
  cat("===============================================\n")
  cat("Call rate threshold:", bold(threshold), "\n")
  cat("Total samples:", bold(total), "\n")
  cat("Passing Threshold:", green(bold(passing)), "(", green(paste0(round(100 * passing / total, 1), "%")), ")\n")
  cat("Failing Threshold:", red(bold(failing)), "(", red(paste0(round(100 * failing / total, 1), "%")), ")\n")
  cat("\n")
  cat("Overall mean call rate:", round(mean(data$`Call Rate`, na.rm = TRUE), 3), "\n")
  cat("Overall range:", round(min(data$`Call Rate`, na.rm = TRUE), 3), "-", 
      round(max(data$`Call Rate`, na.rm = TRUE), 3), "\n")
  
  if (nrow(passing_data) > 0) {
    passing_avg_percent <- round(passing_avg * 100, 1)
    
    cat("Passing mean:", green(round(passing_avg, 3)), "\n")
    
    if (qc_pass) {
      cat("QC Check (>=98.5%):", green("PASS"), "(", passing_avg_percent, "%)\n")
      if (save_log) {
        log_content <- c(
          log_content, 
          paste("Passing mean:", round(passing_avg, 3)),
          paste("QC Check (>=98.5%): PASS (", passing_avg_percent, "%)")
        )
      }
    } else {
      cat("QC Check (>=98.5%):", red("FAIL"), "(", passing_avg_percent, "%)\n")
      if (save_log) {
        log_content <- c(
          log_content, 
          paste("Passing mean:", round(passing_avg, 3)),
          paste("QC Check (>=98.5%): FAIL (", passing_avg_percent, "%)")
        )
      }
    }
  }
  
  if (nrow(failing_data) > 0) {
    failing_mean <- round(mean(failing_data$`Call Rate`, na.rm = TRUE), 3)
    failing_samples <- failing_data$`Sample ID`
    
    cat("Failing mean:", red(failing_mean), "\n")
    cat("Failing samples:", red(paste(failing_samples, collapse = ", ")), "\n")
    
    if (save_log) {
      log_content <- c(
        log_content, 
        paste("Failing mean:", failing_mean),
        paste("Failing samples:", paste(failing_samples, collapse = ", "))
      )
    }
  }
  
  cat("===============================================\n\n")
  
  return(list(log_content = log_content))
}

save_com_file <- function(data, date, output_dir) {
  sequence_num <- get_sequence_number(date, output_dir)
  
  filename <- sprintf("IFCE_LFR9_COMEQv1_%s_%02d.csv", date, sequence_num)
  
  dir_path <- here::here(output_dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Created directory:", dir_path, "\n")
  }
  
  # Filter to keep only PASS samples and remove Status column
  filtered_data <- data %>%
    filter(Status == "PASS") %>%
    select(-Status)
  
  full_path <- file.path(dir_path, filename)
  write_csv(filtered_data, full_path, na = "")
  
  cat("File saved:", green(filename), "\n")
  cat("Path:", full_path, "\n")
  cat("Samples saved:", green(nrow(filtered_data)), "(PASS only)\n")
  
  return(full_path)
}

get_sequence_number <- function(date, output_dir) {
  dir_path <- here::here(output_dir)
  
  if (!dir.exists(dir_path)) {
    return(1)
  }
  
  pattern <- sprintf("IFCE_LFR9_COMEQv1_%s_\\d{2}\\.csv", date)
  existing_files <- list.files(dir_path, pattern = pattern)
  
  if (length(existing_files) == 0) {
    return(1)
  }
  
  sequences <- existing_files %>%
    str_extract("_\\d{2}\\.csv$") %>%
    str_extract("\\d{2}") %>%
    as.numeric() %>%
    sort()
  
  next_seq <- max(sequences) + 1
  if (next_seq > 3) {
    cat(yellow("Max sequence reached, using sequence 3\n"))
    return(3)
  }
  
  return(next_seq)
}

save_log_file <- function(log_content, date, output_dir) {
  sequence_num <- get_sequence_number(date, output_dir)
  
  log_filename <- sprintf("IFCE_LFR9_COMEQv1_%s_%02d_log.txt", date, sequence_num)
  
  dir_path <- here::here(output_dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Created directory:", dir_path, "\n")
  }
  
  log_path <- file.path(dir_path, log_filename)
  
  tryCatch({
    writeLines(log_content, log_path)
    cat("Log file saved:", green(log_filename), "\n")
  }, error = function(e) {
    cat(red("Failed to save log file:"), e$message, "\n")
  })
}


