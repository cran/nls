\name{setClass}
\alias{setClass}
\title{set the class of an object}
\usage{
setClass(object, nm)
}
\arguments{
 \item{object}{an object for which a class can be assigned}
 \item{nm}{a character vector with the name of a class or classes}
}
\description{
  This is a convenience function that sets the class of an object and
  returns the object.  It is most useful at the end of a function
  definition where one is creating the object to be returned and would
  prefer not to store it under a name just so the class can be assigned.
}
\value{
  An object with the contents of \code{object} and the class \code{class}.
}
\author{Jose C. Pinheiro and Douglas M. Bates }
\seealso{\code{\link{setNames}}}
\examples{
setClass( 1:3, c("foo") )
# this is just a short form of
tmp <- 1:3
class(tmp) <-  "foo"
tmp
}
\keyword{data}%-- one or more ...
