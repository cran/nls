\name{Ethyl}
\alias{Ethyl}
\non_function{}
\title{Ethyl Acrylate data}
\description{
The \code{Ethyl} data frame has 12 rows and 3 columns of ethyl acrylate
concentrations in exhalate of a rat.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{start}{
    start of collection period, measured from the start of the
    experiment (hr).
}
\item{length}{
    duration of the collection period (hr).
}
\item{CO2}{
    exhaled radioactively tagged CO2 (g).
}
}
\source{
    Bates and Watts (1998), Nonlinear Regression Analysis and Its
    Applications, Wiley, Appendix A1.10.  Originally from Watts,
    deBethizy, and Stiratelli (1986), Technical Report, Rohm and Haas Co.
}
\examples{
library(nls)
data(Ethyl)
plot(cumsum(CO2) ~ I(start + length), data = Ethyl,
     type = "S", log = "x")  # compare to Figure 3.10, page 97
}
\keyword(datasets)
