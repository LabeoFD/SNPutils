% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotype2allele.R
\name{genotype2allele}
\alias{genotype2allele}
\title{Convert Genotype Data to Allele Pairs}
\usage{
genotype2allele(
  genotype_data,
  annotation_data,
  validate_alleles = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{genotype_data}{A data frame containing genotype data with probeset_id column}

\item{annotation_data}{A data frame containing annotation data with Probe.Set.ID, Allele.A, and Allele.B columns}

\item{validate_alleles}{Logical indicating whether to validate allele assignments}

\item{verbose}{Logical indicating whether to print validation results}
}
\value{
A list containing validation results, converted genotypes
}
\description{
Converts numerical genotype data (0, 1, 2) to actual allele pairs using annotation data.
}
\examples{
\dontrun{
results <- genotype2allele(genotype_data, annotation_data)
}
}
