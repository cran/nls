\name{Leaves}
\alias{Leaves}
\non_function{}
\title{Growth of Leaves}
\description{
The \code{Leaves} data frame has 15 rows and 2 columns of leaf length
over time.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Time}{
    time from initial emergence (days).
}
\item{Length}{
    leaf length (cm).
}
}
\source{
Bates and Watts (1998), Nonlinear Regression Analysis and Its
Applications, Wiley, Appendix A4.5.  Originally from Heyes and Brown
(1956) in F.L. Milthorpe (ed), "The Growth of Leaves", London: Butterworth.
}
\examples{
library(nls)
data(Leaves)
plot(Length ~ Time, data = Leaves)
## first fit a logistic model
fm1 <- nls(Length ~ SSlogis(Time, Asym, xmid, scal), data = Leaves)
## compare with Richards growth model
# fm2 <- nls(Length ~ SSRichards(Time, Asym, xmid, scal, pow),
#            data = Leaves)
# anova( fm2, fm1 )
}
\keyword(datasets)
