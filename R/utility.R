### $Id: utility.R,v 1.1.1.1 1999/04/01 00:56:43 saikat Exp $
###
###            Utility functions used with nls
###
### Copyright 1999-1999 Jose C. Pinheiro <jcp@research.bell-labs.com>,
###                     Douglas M. Bates <bates@stat.wisc.edu>
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

setClass <- function( object, cl ) {
  class( object ) <- cl
  object
}

setNames <- function( object, nm ) {
  names( object ) <- nm
  object
}

clearNames <- function( object ) {
  names( object ) <- NULL
  object
}
