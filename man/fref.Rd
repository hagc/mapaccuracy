\name{fref}
\alias{fref}
\docType{data}
\title{
Reference data to assess the accuracy of a fuzzy classification
}
\description{
Each observation presents single class membership, represented by either TRUE or FALSE, which is equivalent to a matrix of zeros and ones. This is an example of a crisp reference data set.
}
\usage{data("fref")}
\format{
  A matrix with 500 observations identified as belonging to one of the following classes.
  \describe{
    \item{\samp{a}}{membership to class a}
    \item{\samp{b}}{membership to class b}
    \item{\samp{c}}{membership to class c}
    \item{\samp{d}}{membership to class d}
    \item{\samp{e}}{membership to class e}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(fref)
## maybe haed(fref)
}
\keyword{datasets}
