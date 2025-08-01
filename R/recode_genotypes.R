#' Recode Numeric Genotype Data to Character Labels
#'
#' This function converts numeric genotype codes to their corresponding character labels:
#' -1 -> "NoCall", 0 -> "AA", 1 -> "AB", 2 -> "BB"
#'
#' @param genotype_data A data frame containing genotype data with a probeset_id column
#'   and numeric genotype columns coded as -1, 0, 1, or 2
#' @param id_column Character string specifying the name of the ID column to exclude 
#'   from recoding (default: "probeset_id")
#'
#' @return A data frame with the same structure as input, but with numeric genotype 
#'   values recoded to character labels
#'
#' @examples
#' # Create sample data
#' sample_data <- data.frame(
#'   probeset_id = c("SNP1", "SNP2", "SNP3"),
#'   sample1 = c(-1, 0, 1),
#'   sample2 = c(0, 1, 2),
#'   sample3 = c(2, -1, 0)
#' )
#' 
#' # Recode genotypes
#' recoded_data <- recode_genotypes(sample_data)
#' print(recoded_data)
#'
#' @export
#' @importFrom dplyr mutate across
#' @importFrom dplyr case_when
recode_genotypes <- function(genotype_data, id_column = "probeset_id") {
  # Input validation
  if (!is.data.frame(genotype_data)) {
    stop("genotype_data must be a data frame")
  }
  
  if (!id_column %in% names(genotype_data)) {
    stop(paste("ID column '", id_column, "' not found in genotype_data", sep = ""))
  }
  
  if (ncol(genotype_data) < 2) {
    stop("genotype_data must have at least 2 columns (ID column and at least one genotype column)")
  }
  
  # Check genotype values are valid (-1, 0, 1, 2, or NA)
  genotype_data %>%
    select(-all_of(id_column)) %>%
    unlist() %>%
    .[!is.na(.)] %>%
    {if (any(!. %in% c(-1, 0, 1, 2)))
       stop("Found invalid genotype values. Only -1, 0, 1, 2, and NA are allowed")}

  # Perform the recoding
  recoded_data <- genotype_data %>%
    mutate(across(-all_of(id_column), ~ case_when(
      .x == -1 ~ "NoCall",
      .x == 0 ~ "AA",
      .x == 1 ~ "AB", 
    .x == 2 ~ "BB",
      TRUE ~ as.character(.x)  # Keep any other values as-is
    )))
  
  return(recoded_data)
}
