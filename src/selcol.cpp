/*
#    Function: Returns selected columns of a matrix as a new scythe matrix object.
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
#include "la.h"
#include "rng.h"
#include "smath.h"
#include "stat.h"

using namespace scythe;

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> selcol (const Matrix<T, O, S>&M, const Matrix<T, O, S>&W){

int length = W.rows();

Matrix<double> R = M(_ , W(0, 0)); 

for(int i = 1; i < length; i++){R = cbind(R, M(_, W(i, 0)));}

return(R); }

