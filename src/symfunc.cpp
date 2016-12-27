/*
#    Function: Implementation of Andersen's summation algorithm following Verhelst, Glas & van der Sluis (1984);
#    Computes elementary symmetric functions and their derivatives.
#    Uses parts of the Scythe Statistical Library.
#
#    Copyright (C) 2011  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/
*/


#include "matrix.h"
using namespace scythe;

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> symfun2 (const Matrix<T, O, S>&eta){

int k = eta.rows(); 
Matrix <double> M(k, k, false);
M(0, 0) = eta(0,0);

for(int i = 1; i<k; i++){M(0, i) = eta(i, 0) + M(0, i-1); 
			 M(i, i) = eta(i, 0) * M(i-1, i-1);}

for(int j = 1; j<k; j++){
	for(int i = 1; i < j; i++){M(i, j) = M(i, j-1) + eta(j, 0)*M(i-1, j-1);}
}
return(M(_, k-1)); 		     
}		     
		     
		     

