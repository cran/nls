\name{nls}
\alias{nls}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Nonlinear least squares}
\usage{
nls(formula, data, start, control=nlsControl())
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{formula}{a nonlinear model formula including variables and parameters}
 \item{data}{a data frame in which to evaluate the variables in \code{formula}}
 \item{start}{a named list or named numeric vector of starting estimates}
 \item{control}{an optional list of control settings.  See
   \code{nlsControl} for the names of the settable control values and
   their effect.}
}
\description{
  Determine the nonlinear least squares estimates of the parameters.
}
\details{
  An \code{nls} object is a type of fitted model object.  It has methods
  for the generic functions \code{coef}, \code{formula}, \code{resid},
  \code{print}, \code{summary}, and \code{fitted}.
}
\value{
  A list of
  \item{m}{an \code{nlsModel} object incorporating the model}
  \item{data}{the expression that was passed to \code{nls} as the data
    argument.  The actual data values are present in the environment of
    the \code{m} component.}
}
\references{Bates, D.M. and Watts, D.G. (1988), \emph{Nonlinear Regression Analysis
    and Its Applications}, Wiley}
\author{Jose C. Pinheiro and Douglas M. Bates}
\seealso{\code{\link{nlsModel}}}
\examples{
data( DNase )
DNase1 <- DNase[ DNase$Run == 1, ]
fm1DNase1 <- nls( density ~ SSlogis( log(conc), Asym, xmid, scal ), DNase1,
             list( Asym = 1.8, xmid = 0, scal = 1 ) )
}
\keyword{nonlinear, regression, models}%-- one or more ...
