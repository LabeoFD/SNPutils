# In your package directory: R/create_com_file.R
#' Create COM File from Axiom GT Report
#'
#' Reads and processes an Axiom GT report file, extracting date information
#' and sample call rate data with quality comments.
#'
#' @param file_path Character string. Path to the AxiomGT report file (.txt)
#' @param call_rate_threshold Numeric. Threshold for call rate quality (default: 0.95)
#'
#' @return A list containing:
#'   \item{formatted_date}{Character string in DDMMYY format}
#'   \item{data}{Tibble with Sample ID, Call Rate, and Comment columns}
#'
#' @examples
#' \dontrun{
#' result <- create_com_file("path/to/AxiomGT1.report.txt")
#' date <- result$formatted_date
#' data <- result$data
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import readr
#' @import lubridate
#' @import purrr
create_com_file <- function(file_path, call_rate_threshold = 0.95) {
  
  # Validate inputs
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  if (!is.numeric(call_rate_threshold) || call_rate_threshold < 0 || call_rate_threshold > 1) {
    stop("call_rate_threshold must be a numeric value between 0 and 1")
  }
  
  # Read all lines from the file
  lines <- readLines(file_path)
  
  # Extract and parse the date
  formatted_date <- lines %>% 
    str_subset("#%affymetrix-algorithm-param-apt-time-str") %>%  
    str_split("\\s*=\\s*", simplify = TRUE) %>%                  
    .[,2] %>%                                                    
    parse_date_time(orders = "b! d! H!:M!:S! Y!") %>%            
    format("%d%m%y")                                             
  
  # Read data table from non-comment lines
  Pre_COM <- lines %>% 
    discard(str_detect, "^#") %>%               
    I() %>%                                     
    read_delim(delim = "\t", col_names = TRUE, show_col_types = FALSE) %>% 
    select(cel_files, total_call_rate) %>% 
    mutate(
      Sample_ID = cel_files,
      Call_Rate = total_call_rate / 100
    ) %>% 
    relocate(Sample_ID, .before = cel_files) %>% 
    select(-cel_files, -total_call_rate) %>% 
    mutate(
      Comment = if_else(
        Call_Rate >= call_rate_threshold,
        "0",
        "1"
      )
    ) %>% 
    rename_with(~ gsub("_", " ", .x))
  
  # Return both date and data
  return(list(
    formatted_date = formatted_date,
    data = Pre_COM
  ))
}
