#' Example Genotype Data
#'
#' A dataset containing example genotype data for 20 SNPs, including 10 indels.
#'
#' @format A data frame with 20 rows and 4 columns:
#' \describe{
#'   \item{probeset_id}{SNP identifier}
#'   \item{Sample1}{Genotype codes (0, 1, 2) for Sample1}
#'   \item{Sample2}{Genotype codes (0, 1, 2) for Sample2}
#'   \item{Sample3}{Genotype codes (0, 1, 2) for Sample3}
#' }
#' @examples
#' data(example_genotypes)
#' head(example_genotypes)
"example_genotypes"

#' Example Annotation Data
#'
#' A dataset containing annotation information for 20 SNPs, including 10 regular SNPs and 10 indels.
#'
#' @format A data frame with 20 rows and 5 columns:
#' \describe{
#'   \item{Probe.Set.ID}{SNP identifier}
#'   \item{Chromosome}{Chromosome location}
#'   \item{Physical.Position}{Base pair position}
#'   \item{Allele.A}{First allele (deletion marker "-" for indels)}
#'   \item{Allele.B}{Second allele (insertion sequence for indels)}
#' }
#' @examples
#' data(example_annotations)
#' head(example_annotations)
"example_annotations"
