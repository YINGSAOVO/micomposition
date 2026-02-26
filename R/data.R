#' Demo microbiome dataset
#'
#' A simulated microbiome composition dataset for demonstration purposes.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{otu_table}{A numeric matrix (10 samples × 30 OTUs) of relative abundances.}
#'   \item{tax_table}{A data frame with 30 rows and 7 columns: OTU, Kingdom, Phylum, Class, Order, Family, Genus.}
#'   \item{sample_data}{A data frame with 10 rows and 3 columns: SampleID, Group, Subject.}
#' }
#' @examples
#' data(demo_data)
#' str(demo_data)
"demo_data"
