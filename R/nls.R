### $Id: nls.R,v 1.8 1999/05/27 19:37:56 bates Exp $
###
###            Nonlinear least squares for R
###
### Copyright 1999-1999 Jose C. Pinheiro <jcp@research.bell-labs.com>,
###                     Douglas M. Bates <bates@stat.wisc.edu>,
###                     Saikat DebRoy <saikat@stat.wisc.edu>
###
### This file is part of the nlme library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

nlsModel.plinear <- function( .form, .data, .start, .temp) {
  on.exit(remove(.i, .data, .start))
  .env <- environment( ) 
  .ind <- as.list( .start )
  .form <- as.formula( .form )
  for( .i in all.vars( .form ) ) {
    if (is.na(match(.i, names(.ind), no = NA))) {
      .temp <- if(is.environment(.data))
        get(.i, env=.data)
      else  .data[[ .i ]]
      if (is.null(.temp)) .temp <- get(.i)  # for cases like "pi"
      storage.mode(.temp) <- "double"
      assign( .i, .temp)
    }
  }
  .p2 <- 0
  for( .i in names( .ind ) ) {
    .temp <- .start[[ .i ]]
    storage.mode(.temp) <- "double"
    assign( .i, .temp )
    .ind[[ .i ]] <- .p2 + seq( along = .start[[ .i ]] )
    .p2 <- .p2 + length( .start[[ .i ]] )
  }
  .lhs <- eval( .form[[2]] )
  storage.mode(.lhs) <- "double"
  .rhs <- eval( .form[[3]] )
  storage.mode(.rhs) <- "double"
  .p1 <- if(is.matrix(.rhs))
    ncol(.rhs)
  else 1
  .p <- .p1 + .p2
  .n <- length(.lhs)
  .fac <- (.n - .p)/.p
  .cc <- .QR.B <- NA
  
  if(is.null(attr( .rhs, "gradient" ))) {
    .getRHS <- function()
      .External("numeric_deriv", .form[[3]], names( .ind), .env)
    .QR.rhs <- qr(.rhs <- .getRHS())
  } else {
    .getRHS <- function()
      eval( .form[[3]])
    .QR.rhs <- qr(.rhs)
  }
  .lin <- qr.coef( .QR.rhs, .lhs )
  .resid <- qr.resid(.QR.rhs, .lhs)
  .topzero <- double(.p1)
  .marg <- length(dim(attr(.rhs, "gradient")))
  .dev <- sum( ( .resid )^2 )
  if(.marg <= 1) {
    .ddot <- function(A,b) A %*% b
    .dtdot <- function(A,b) t(A) %*% b
  } else if(.marg == 2) {
    if(.p1 == 1) {
      .ddot <- function(A,b) A*b
      .dtdot <- function(A,b) t(b) %*% A
    } else if(.p2 == 1) {
      .ddot <- function(A,b) A %*% b
      .dtdot <- function(A,b) t(A) %*% b
    }
  } else {
    .ddot <- function(A,b) apply(A, MARGIN = 3, FUN="%*%", b)
    .dtdot <- function(A,b) apply(A, MARGIN = c(2,3), FUN = "%*%", b)
  }
  .error <- NULL
  .getPars <- function()
    unlist( setNames( lapply( as.list(names(.ind)), get, envir =
                             .env ), names( .ind ) ) )  
  .m <-
    setClass( list(resid = function() .resid,
                   fitted = function() .rhs,
                   formula = function() .form,
                   deviance = function() .dev,
                   gradient = function() attr( .rhs, "gradient" ),
                   conv = function() {
                     assign(".cc", c(.topzero, qr.qty( .QR.rhs,
                                                      .lhs)[-(1:.p1)]),
                            envir = .env)
                     rr <- qr.qy( .QR.rhs, .cc)
                     B <- qr.qty(.QR.rhs, .ddot(attr(.rhs, "gradient"), .lin)) 
                     B[1:.p1, ] <- .dtdot(attr(.rhs, "gradient"), rr)
                     R <- t(qr.R(.QR.rhs)[1:.p1,])
                     if(.p1 == 1) B[1, ] <- B[1, ]/R
                     else B[1:.p1, ] <- solve(R, B[1:.p1, ])
                     assign(".QR.B",qr(B), envir = .env)
                     rr <- qr.qty( .QR.B, .cc)
                     sqrt( .fac*sum( rr[1:.p1]^2 ) / sum( rr[ -(1:.p1) ]^2 ) )
                   },
                   incr = function() {
                     qr.solve( .QR.B, .cc)
                   },
                   setError = function(error) {
                     assign(".error",
                            switch(error,
                                   "step factor reduced below minimum",
                                   "maximum number of iterations exceeded",
                                   "singular design matrix for linear parameters"),
                            envir = .env)
                   },
                   handleAnyError = function(warnOnly = FALSE) {
                     if(!is.null(.error)) {
                       if(warnOnly) warning(.error)
                       else stop(.error)
                     }
                   },
                   setPars = function(newPars) {
                     for(i in names(.ind) ) {
                       assign( i, clearNames(newPars[ .ind[[i]] ]), envir = .env )
                     }
                     assign(".QR.rhs",qr(assign(".rhs",.getRHS(),envir =.env)),
                            envir = .env)
                     assign( ".resid", qr.resid(.QR.rhs, .lhs), envir=.env)
                     assign(".dev",sum( ( .resid )^2), envir = .env)
                     if (.QR.rhs$rank < .p1) {
                       return(1)
                     } else {
                       assign( ".lin", qr.coef(.QR.rhs, .lhs), envir = .env )
                       return(0)
                     }
                   },
                   getPars = .getPars,
                   getAllPars = function() {
                     c(.getPars(), c(.lin = .lin))
                   },
                   getEnv = function() .env,
                   trace = function() cat(format(.dev),":",
                     format(c(.getPars(), .lin)),"\n"),
                   Rmat = function() qr.R( .QR )
                   ),
             "nlsModel.plinear" )
  .m$conv();
  .m
}

nlsModel <- function( .form, .data, .start ) {
  on.exit(remove(.i, .data, .offset, .start))
  .env <- environment( ) 
  .ind <- as.list( .start )
  .form <- as.formula( .form )
  for( .i in all.vars( .form ) ) {
    if (is.na(match(.i, names(.ind), no = NA))) {
      .temp <- if(is.environment(.data))
        get(.i, env=.data)
      else  .data[[ .i ]]
      if (is.null(.temp)) .temp <- get(.i)  # for cases like "pi"
      storage.mode(.temp) <- "double"
      assign( .i, .temp)
    }
  }
  .offset <- 0
  for( .i in names( .ind ) ) {
    .temp <- .start[[ .i ]]
    storage.mode(.temp) <- "double"
    assign( .i, .temp )
    .ind[[ .i ]] <- .offset + seq( along = .start[[ .i ]] )
    .offset <- .offset + length( .start[[ .i ]] )
  }
  .lhs <- eval( .form[[2]] )
  .rhs <- eval( .form[[3]] )
  .resid <- .lhs - .rhs
  .dev <- sum( ( .resid )^2 )
  if(is.null(attr( .rhs, "gradient" ))) {
    .getRHS <- function()
      .External("numeric_deriv", .form[[3]], names( .ind), .env)
    .QR <- qr(attr(.getRHS(), "gradient"))
  } else {
    .getRHS <- function()
      eval( .form[[3]])
    .QR <- qr(attr( .rhs, "gradient" ))
  }
  .dim <- min(dim(.QR$qr))
  if(.QR$rank < .dim)
    stop("singular gradient matrix at initial parameter estimates");
  .getPars <- function()
    unlist( setNames( lapply( as.list(names(.ind)), get, envir =
                             .env ), names( .ind ) ) )  
  .error <- NULL
  
  setClass( list(resid = function() .resid,
                 fitted = function() .rhs,
                 formula = function() .form,
                 deviance = function() sum( ( .resid )^2 ),
                 gradient = function() attr( .rhs, "gradient" ),
                 conv = function() {
                   rr <- qr.qty( .QR, .resid ) # rotated residual vector
                   p <- length( unlist( .ind ) )  # length of parameter vector
                   sqrt( sum( rr[1:p]^2 ) / sum( rr[ -(1:p) ]^2 ) )
                 },
                 incr = function() {
                   qr.coef( .QR, .resid )
                 },
                 setError = function(error) {
                   assign(".error",
                          switch(error,
                                 "step factor reduced below minimum",
                                 "maximum number of iterations exceeded",
                                 "singular gradient matrix"),
                          envir = .env)
                 },
                 handleAnyError = function(warnOnly = FALSE) {
                   if(!is.null(.error)) {
                     if(warnOnly) warning(.error)
                     else stop(.error)
                   }
                 },
                 setPars = function(newPars) {
                   for(i in names(.ind) ) {
                     assign( i, clearNames(newPars[ .ind[[i]] ]), envir = .env )
                   }
                   assign( ".resid",
                          .lhs - assign(".rhs", .getRHS(), envir = .env ),
                          envir = .env )
                   assign(".dev",sum( ( .resid )^2), envir = .env)
                   assign(".QR", qr( attr( .rhs, "gradient")), envir = .env )
                   return(.QR$rank < .dim)  # to catch the singular gradient matrix
                 },
                 getPars = .getPars,
                 getAllPars = .getPars,
                 getEnv = function() .env,                   
                 trace = function() cat(format(.dev),"\t: ",format(.getPars()),
                   "\n"),
                 Rmat = function() qr.R( .QR )
                 ),
           "nlsModel" )
}

nls.control <- function( maxiter = 50, tol = 0.00001, minFactor = 1/1024 ) {
  list( maxiter = maxiter, tol = tol, minFactor = minFactor )
}

"nls" <-
  function (formula, data = sys.frame(sys.parent()),
            start = getInitial(formula, data), control,
            algorithm="default", trace = F, warnOnly = FALSE)
{
  m <- switch(algorithm,
              plinear=nlsModel.plinear(formula, data, start),
              nlsModel(formula, data, start))
  ctrl <- nls.control()
  if(!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  m <- .External("nls_iter", m, ctrl, trace)
  m$handleAnyError(warnOnly)
  setClass(list(m = m, data = substitute(data)), "nls")
}

coef.nls <- function( x, ... ) x$m$getAllPars()
## The next two methods are defined in more generality in the lme library
# residuals.nls <- function(object, ...) as.vector(object$m$resid())
# fitted.nls <- function(object, ...) as.vector(object$m$fitted())
print.nls <- function(x, ...) {
  cat( "Nonlinear regression model\n" )
  cat( "  model: ", deparse( x$m$formula() ), "\n" )
  cat( "   data: ", as.character( x$data ), "\n" )
  print( x$m$getAllPars() )
  cat( " residual sum-of-squares: ", format( x$m$deviance() ), "\n" )
  invisible(x)
}

### Force loading of the nls dynamic library in R
.First.lib <- function(lib, pkg) library.dynam( "nls", pkg, lib )
require(lme)
