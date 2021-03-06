### $Id: nlsFunc.q,v 1.3 1999/10/26 22:41:45 bates Exp $
###
###            Utility functions used with nls
###
### Copyright 1997,1999 Jose C. Pinheiro <jcp$research.bell-labs.com>,
###                     Douglas M. Bates <bates$stat.wisc.edu>
###           1999-1999 Saikat DebRoy <saikat$stat.wisc.edu>
###
### This file is part of the nls library for R and related languages.
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

###
### asOneSidedFormula is extracted from the NLME-3.1 library for S 
###

asOneSidedFormula <-
  ## Converts an expression or a name or a character string
  ## to a one-sided formula
  function(object)
{
  if ((mode(object) == "call") && (object[[1]] == "~")) {
    object <- eval(object)
  }
  if (inherits(object, "formula")) {
    if (length(object) != 2) {
      stop(paste("Formula", deparse(as.vector(object)),
		 "must be of the form \"~expr.\""))
    }
    return(object)
  }
  do.call("~",
	  list(switch(mode(object),
		      name = ,
                      numeric = ,
		      call = object,
		      character = as.name(object),
		      expression = object[[1]],
		      stop(paste(substitute(object), "cannot be of mode",
				 mode(object))))))
}

setNames <- function( object, nm ) {
  names( object ) <- nm
  object
}

clearNames <- function( object ) {
  names( object ) <- NULL
  object
}
