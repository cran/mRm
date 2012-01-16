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
#include "la.h"
#include "smath.h"
#include "stat.h"

using namespace scythe;

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> symfun2 (const Matrix<T, O, S>&eta); 


template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> selcol (const Matrix<T, O, S>&M, const Matrix<T, O, S>&W); 


template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> symfun (const Matrix<T, O, S>&beta){

int k = beta.rows();
Matrix <double> eta = beta; 
Matrix <double> symfun = symfun2(eta); 
Matrix <double> symfundiff(k, k, true, 0);

for(int i = 0; i < k; i++){
	Matrix<double> intereta = eta; 
	intereta(i, 0) = 0; 
	symfundiff(i, _) = symfun2(intereta);
};
Matrix<double> sel = seqa(0, 1, k-1);
symfundiff = selcol(symfundiff, sel); 

Matrix<double> ones (k, 1, true, 1); 
Matrix<double> Res = cbind(ones, symfundiff);
Res = cbind(symfun, Res);
return(Res); }



