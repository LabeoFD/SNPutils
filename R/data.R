#' Example Genotype Data
#'
#' A dataset containing example genotype data for 20 SNPs, including 10 indels.
#' This data shows genotype calls for multiple samples using numeric encoding (0, 1, 2).
#'
#' @format A data frame with 20 rows and 4 columns:
#' \describe{
#'   \item{probeset_id}{Character. SNP identifier (SNP1 to SNP20)}
#'   \item{Sample1}{Numeric. Genotype codes (0, 1, 2) for Sample1}
#'   \item{Sample2}{Numeric. Genotype codes (0, 1, 2) for Sample2}
#'   \item{Sample3}{Numeric. Genotype codes (0, 1, 2) for Sample3}
#' }
#' @details 
#' Genotype encoding: 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternative.
#' For indels (SNP11-SNP20): 0 = D/D, 1 = D/I, 2 = I/I.
#' @examples
#' data(example_genotypes)
#' head(example_genotypes)
#' table(example_genotypes$Sample1)
"example_genotypes"

#' Example Annotation Data
#'
#' A dataset containing annotation information for 20 SNPs, including 10 regular SNPs and 10 indels.
#' This data provides genomic positions and allele information.
#'
#' @format A data frame with 20 rows and 5 columns:
#' \describe{
#'   \item{Probe.Set.ID}{Character. SNP identifier matching probeset_id}
#'   \item{Chromosome}{Character. Chromosome location (1-22, X)}
#'   \item{Physical.Position}{Integer. Base pair position on chromosome}
#'   \item{Allele.A}{Character. First allele (deletion marker "-" for indels)}
#'   \item{Allele.B}{Character. Second allele (insertion sequence for indels)}
#' }
#' @details
#' For standard SNPs (SNP1-SNP10): Allele.A and Allele.B contain single nucleotides.
#' For indels (SNP11-SNP20): Allele.A contains "-" and Allele.B contains insertion sequence.
#' @examples
#' data(example_annotations)
#' head(example_annotations)
#' # Check for indels
#' sum(example_annotations$Allele.A == "-")
"example_annotations"

#' Example Axiom Report Data
#'
#' A list containing example Axiom genotyping report data with header information
#' and quality metrics from the Axiom genotyping platform.
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{header}{Character vector with 9 elements containing metadata including
#'     chip type, library set information, processing parameters, and timestamps}
#'   \item{data}{Data frame with 15 rows and 4 variables:
#'     \itemize{
#'       \item \code{cel_files}: Character. CEL file names
#'       \item \code{total_call_rate}: Numeric. Percentage of successful genotype calls
#'       \item \code{average_heterozygosity}: Numeric. Average heterozygosity rate
#'       \item \code{median_intensity}: Numeric. Median signal intensity
#'     }
#'   }
#' }
#' @examples
#' data(example_axiom_report)
#' example_axiom_report$header
#' head(example_axiom_report$data)
#' summary(example_axiom_report$data$total_call_rate)
"example_axiom_report"
