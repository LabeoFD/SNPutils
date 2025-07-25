% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotype2allelev2.R
\name{genotype2allelev2}
\alias{genotype2allelev2}
\title{Convert Genotype Data to Allele Pairs (Version 2)}
\usage{
genotype2allelev2(
  genotype_data,
  annotation_data,
  validate_alleles = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{genotype_data}{A data frame containing genotype data with probeset_id column}

\item{annotation_data}{A data frame containing annotation data with Probe.Set.ID, Allele.A, and Allele.B columns}

\item{validate_alleles}{Logical indicating whether to validate allele assignments (default: TRUE)}

\item{verbose}{Logical indicating whether to print validation results (default: TRUE)}
}
\value{
A list containing validation results and converted genotypes
}
\description{
Converts numerical genotype data (0, 1, 2) to actual allele pairs using annotation data.
This enhanced version includes improved validation and better handling
of special cases like missing genotypes and multi-character alleles.
}
\details{
This enhanced version provides:
\itemize{
\item Space-to-dash conversion for alleles
\item Handling of missing genotypes (-1 -> "--")
\item Extended indel detection for multi-character alleles
\item Improved validation and error handling
}
}
\examples{
\dontrun{
# Basic usage
results <- genotype2allelev2(genotype_data, annotation_data)

# Silent mode
results <- genotype2allelev2(genotype_data, annotation_data, verbose = FALSE)
}

}
\seealso{
\code{\link{genotype2allele}} for the original version
}
