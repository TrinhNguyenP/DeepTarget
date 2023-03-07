\name{OntargetM}
\alias{OntargetM}
\docType{data}
\title{
 an object containing data to be used in DeepTarget pipeline.

}
\description{
an object containing Viability matrix after CRISPR-KO; Viability after Drug Treatment; Drug metadata from Broad, mutation matrix, and expression matrix with common celllines and common drugs.
}
\usage{data("OntargetM")}
\format{
  A list of 5 matrices
  \describe{
    \item{\code{DrugMetadata}}{containing 1448 drugs with 5 columns:  broad_id_trimmed, name, target, drug_category, and moa}
    \item{\code{secondary_prism}}{containing Viability scores after Drug Treatment with 1448 drugs across 319 celllines}
    \item{\code{avana_CRISPR}}{containing Viability scores after CRISPR-KO with 18333 genes across 319 celllines}
    \item{\code{mutations_mat}}{Mutation binary matrix for ~2K cell lines; 0 is WT; 1 is mutated}
    \item{\code{expression_20Q4}}{ Expression profile of 319 celllines for 20K genes}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
DrugMetadata
depmap \url{ https://depmap.org/repurposing/#:~:text=Corsello_supplemental_tables.xlsx}
Secondary prism
depmap \url{https://depmap.org/portal/download/all/?releasename=PRISM+Repurposing+19Q4&filename=secondary-screen-dose-response-curve-parameters.csv}
avana_CRISPR
depmap \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=CRISPRGeneEffect.csv}
mutations_mat
Mutation binary profile for ~2K cell lines; 0 is WT; 1 is mutated

depmap \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=OmicsSomaticMutations.csv}

expression_20Q4
Expression profile of ~2K cell lines for 20K genes; file name: OmicsExpressionProteinCodingGenesTPMLogp1.csv
depmap \url{https://depmap.org/portal/download/all/}

}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(OntargetM)
}
\keyword{datasets}
