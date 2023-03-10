\name{PredTarget}
\alias{PredTarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction the most similarity across the known targeted genes.
}
\description{
Predict the gene having the most similarity associated with drug’s response scores across targeted genes
}

\usage{
PredTarget(Sim.GES.DRS=sim, D.M = Drug.Metadata)
}

\arguments{
  \item{Sim.GES.DRS}{
    similarity between Drug’s response scores and Gene effect scores from Knock-out method such as CRISPR.
}
  \item{D.M}{
   Drug Metadata
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a table of drug information and their most targeted gene with stats of correlation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 sanjusinha7
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
}
