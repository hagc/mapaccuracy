\name{mapaccuracy-package}

\alias{mapaccuracy-package}

\docType{package}

\title{
  Overview of the \pkg{mapaccuracy} package
}

\description{
Unbiased estimators of overall and per-class thematic map accuracy and area for diverse sampling designs (simple, systematic and stratified) and classification types (hard and fuzzy).
}
\details{
Implemented estimators are described in
\itemize{
\item olofsson, P.; Foody, G. M.; Herold, M.; Stehman, S. V.; Woodcock, C. E.; Wulder, M. A. (2014). \href{https://doi.org/10.1016/j.rse.2014.02.015}{Good practices for estimating area and assessing accuracy of land change}. \emph{Remote Sens. Environ.}, 148, 42-57.
\item olofsson, P.; Foody, G. M.; Stehman, S. V.; Woodcock, C. E. (2013). \href{https://doi.org/10.1016/j.rse.2012.10.031}{Making better use of accuracy data in land change studies: Estimating accuracy and area and quantifying uncertainty using stratified estimation}. \emph{Remote Sens. Environ.}, 129, 122-131.
\item Stehman, S. V. (2014).\href{https://doi.org/10.1080/01431161.2014.930207}{Estimating area and map accuracy for stratified random sampling when the strata are different from the map classes}. \emph{Int. J. Remote Sens.}, 35, 4923-4939.
\item Stehman, S. V.; Arora, M. K.; Kasetkasem, T.; Varshney, P. K (2007). \href{https://doi.org/10.14358/PERS.73.2.165}{Estimation of fuzzy error matrix accuracy measures under stratified random sampling}. \emph{Photogramm. Eng. Remote Sensing}, 73, 165–173.
\item Binaghi, E.; Brivio, P. A.; Ghezzi, P.; Rampini, A. (1999). \href{https://doi.org/10.1016/S0167-8655(99)00061-6}{A fuzzy set-based accuracy assessment of soft classification}. \emph{Pattern Recognit. Lett.}, 20, 935-948.
}
.
}
\section{I. Hard (crisp) vs soft (fuzzy) classification}{
Traditional thematic maps (e.g. land cover maps) allocate every single areal unit to a single class. This is often called 'hard' or 'crisp' classification in the remote sensing community. However, areal unit units (e.g. transitional land, mixed pixels, etc.) can present multiple class membership, which can be considered in 'soft' or 'fuzzy' classification. The following estimators are suitable for hard and soft classification:
\tabular{ll}{
\code{hard} \tab \code{\link{olofsson}}, \code{\link{stehman2014}} \cr
\code{soft} \tab To create a RasterStack (multiple layers)\cr
}
}
\section{II. Sampling design}{
\tabular{ll}{
\code{Simple random} \tab \code{\link{olofsson}} \cr
\code{Systematic random} \tab \code{\link{olofsson}} \cr
\code{Stratified random} \tab \code{\link{olofsson}}, \code{\link{stehman2014}} \cr
}
}
\section{III. Areal unit}{
\tabular{ll}{
\code{Point/pixel} \tab \code{\link{olofsson}}, \code{\link{stehman2014}} \cr
}
}

\keyword{package}
\keyword{spatial}
