%%% $Id: Coal.Rd,v 1.1 1999/05/27 19:37:56 bates Exp $
\name{Coal}
\alias{Coal}
\non_function{}
\title{Coal liquefaction data}
\description{
The \code{Coal} data frame has 23 rows and 20 columns of data on a coal
liquefaction process.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Time}{
    the time of the observation, measured from the beginning of the run (min).
}
\item{temp}{
    the temperature of the run (K).
}
\item{pressure}{
    the pressure at which the run was performed (MPa).
}
\item{fe2o3}{
    iron oxide? (Fe2 O3) concentration in the run (wt\% maf).
}
\item{tr}{
    transfer rate of the feed (wt\% maf).
}
\item{x1}{
    inlet composition of unconverted coal (wt \%).
}
\item{x2}{
    inlet composition of thermal residuals (wt \%).
}
\item{x3}{
    inlet composition of C4-822K (wt \%).
}
\item{x4}{
    inlet composition of C1-C3 gases (wt \%).
}
\item{x5}{
    inlet composition of byproduct gases (wt \%).
}
\item{x6}{
    inlet composition of water (wt \%).
}
\item{x7}{
    inlet composition of hydrogen (wt \%).
}
\item{x8}{
    inlet composition of Coal (wt \%).
}
\item{y1}{
    outlet composition of unconverted coal (wt \%).
}
\item{y2}{
    outlet composition of thermal residuals (wt \%).
}
\item{y3}{
    outlet composition of C4-822K (wt \%).
}
\item{y4}{
    outlet composition of C1-C3 gases (wt \%).
}
\item{y5}{
    outlet composition of byproduct gases (wt \%).
}
\item{y6}{
    outlet composition of water (wt \%).
}
\item{y7}{
    outlet composition of hydrogen (wt \%).
}
\item{y8}{
    outlet composition of Coal (wt \%).
}
}
\source{
Bates and Watts (1998), Nonlinear Regression Analysis and Its
Applications, Wiley, Appendix A4.7.  Originally from Lythgoe (1986),
M.Sc. Thesis, Queen's University at Kingston.
}
\examples{
library(nls)
data(Coal)
coplot( y1 ~ x1 | x8 * Time, data = Coal )
coplot( y2 ~ x2 | x8 * Time, data = Coal )
coplot( y3 ~ x3 | x8 * Time, data = Coal )
coplot( y4 ~ x4 | x8 * Time, data = Coal )
coplot( y5 ~ x5 | x8 * Time, data = Coal )
coplot( y6 ~ x6 | x8 * Time, data = Coal )
coplot( y7 ~ x7 | x8 * Time, data = Coal )
}
\keyword(datasets)
