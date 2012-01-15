/*  
#    Function: Core routine of mrm function. Performs iterations of the 
#    EM algorithm (Rost 1991). Uses Parts of the Scythe Statistical Library
#    for matrix manipulations.
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
#include "selcol.cpp"
#include "selrow.cpp"
#include "symfunc.cpp"
#include "symfuncomplete.cpp"
#include <R.h>           
#include <R_ext/Utils.h> 


using namespace std; 
using namespace scythe;

template <typename T, matrix_order O, matrix_style S> 
Matrix<T, O, S> selcol (const Matrix<T, O, S>&M, const Matrix<T, O, S>&W); 

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> selrow (const Matrix<T, O, S>&M, const Matrix<T, O, S>&W); 

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> symfun2 (const Matrix<T, O, S>&beta); 

template <typename T, matrix_order O, matrix_style S>
Matrix<T, O, S> symfun (const Matrix<T, O, S>&beta); 


extern"C"{
void em (const int *Rk, const int *RNwhole, const int *RN, const int *Rdimpirc, 
	 const int *Rcl, const int *Rmaxit, const double *Rconvcrit, 
	 double *Rbeta, double *Rpirc, double *Rresppat, double *Rclasssize, double *Lik, int *Riterations,
	 int *Rnpara, int *Rboundary){

Matrix <double> beta(*Rk, *Rcl, Rbeta); 	
beta = exp(beta); 
Matrix <double> pirc(*Rdimpirc, *Rcl, Rpirc); 			
Matrix <double> resppat(*RN, *Rk + 2, Rresppat); 		
Matrix <double> classsize(*Rcl, 1, Rclasssize); 	
Matrix <double> N(*Rk, *Rcl, true, 0); 
Matrix <double> M = N; 
Matrix <double> M5(*Rk - 1, *Rcl, false);
int iteration = 0; 
int check = 0; 
Matrix<double> numberscore = resppat(_, *Rk);  
Matrix<double> cond(*RN, *Rcl, false);				
Matrix<double> patternsumloc = resppat(_, *Rk +1) - 1;		
Matrix<double> where = seqa(0, 1, *Rk);				
Matrix<double> patterns = selcol(resppat, where); 		
Matrix<double> nxc(*RN, *Rcl, false);
Matrix<double> score = patternsumloc + 1; 
Matrix<double> r = unique(score); 
int nr = r.cols(); 
Matrix <double> symfunscl(*Rk, *Rcl, false);
Matrix <double> symfundiff(*Rk, 1, false);  
Matrix <double>  sel = seqa(1, 1, *Rk);       
double logLik = 0; 
double prelimlogLik = 0;
double logLikdiff = 0;
double AIC = 0; 
double BIC = 0; 
int para = 2*(*Rcl)*(*Rk-1) + 1;
double dRNwhole = (double) *RNwhole; 
*Rnpara = para; 


while(iteration < *Rmaxit){
iteration = iteration + 1; 
R_CheckUserInterrupt();

//1)=== SYMMETRIC FUNCTIONS 

Matrix <double> SF = symfun(beta(_, 0));
symfunscl(_, 0) = SF(_, 0);
Matrix<double> firstdiff = selcol(SF, sel);
symfundiff = firstdiff;

for(int j = 1; j < *Rcl; j++){
	SF = symfun(beta(_, j)); 
	symfunscl(_, j) = SF(_, 0); 
	Matrix<double> S = selcol(SF, sel);
	symfundiff = cbind(symfundiff, S); 
;}; 

//2)=== CONDITIONAL PATTERN FREQUENCIES

Matrix<double> tclasssize = t(classsize); 

for(int j = 0; j < *RN; j++){
	double wo = patternsumloc(j, 0); 
	Matrix<bool> patternloc = t(patterns(j, _)); 
	cond(j, _) = (pirc(wo, _) / symfunscl(wo, _)) % prodc(selif(beta, patternloc)); //cond(i, j) = P(X = x|c) 
	nxc(j, _) =  (resppat(_, *Rk)(j, 0) * (tclasssize % cond(j, _))) / (cond(j,_)*classsize); //nxc(j,i) = \hat(n)(x,c); 
};;


//3)=== LOG-LIKELIHOOD

prelimlogLik = sum(log(cond * classsize) % numberscore); 
logLikdiff = fabs(prelimlogLik - logLik); 
logLik = prelimlogLik; 

AIC = -2*logLik + 2*para;
BIC = -2*logLik + para * log(dRNwhole); 

if(logLikdiff < *Rconvcrit || check == 1 || iteration == *Rmaxit)
	{break; 
	 }

//4)=== M-STEP

for(int j = 0; j < nr; j++){
			Matrix<bool> sel = score == r(0, j); 
			M(r(0, j)-1, _) = sumc(selif(nxc, sel));} 


for(int j = 0; j < *Rcl; j++){
	Matrix <double> selsym = seqa(j* (*Rk), 1, *Rk); 
	Matrix <double> transpSFdiffclass = t(selcol(symfundiff, selsym));
		
	for(int l = 0; l < *Rk; l++){
		N(l, j) = sum(patterns(_, l) % nxc(_, j)); 
		beta(l, j) = N(l, j)/sum((M(_, j)% transpSFdiffclass(_,l))/symfunscl(_,j)); 
;}}

//5) === Product 1 Standardization of exp(beta) and check if abs(beta) exceeds 20

double frac = 1/double(*Rk);
Matrix <double> stand = pow(prodc(beta), frac); 

for(int j = 0; j < *Rcl; j++){
	for(int l = 0; l < *Rk; l++){
	beta(l,j) = beta(l,j) / stand(0, j);
	if(beta(l,j) >= 485165195 || beta(l,j) <= .000000002){check = 1;}
}}


//6)=== NEW VALUES

Matrix <double> nc = sumc(nxc); 
Matrix <double> selprel = seqa(0, 1, (*Rk - 1)); 
for(int j = 0; j < (*Rk - 1); j++){M5(j,_) = nc;}; 
Matrix<double> Q = selrow(M, selprel);
pirc = Q/M5; 
classsize = t(nc/ *RNwhole); 

	
}

//================== TERMINATION =================================

//return estimates; 

Lik[0] = logLik;
Lik[1] = AIC; 
Lik[2] = BIC; 
*Riterations = iteration; 
*Rboundary = check;

beta = log(beta); 
for(int i = 0; i < *Rk; i++){
	for(int j= 0; j < *Rcl; j++){Rbeta[i + j*(*Rk)] = beta(i,j);}}

for(int i = 0; i < *Rdimpirc; i++){
	for(int j= 0; j < *Rcl; j++){Rpirc[i + j*(*Rdimpirc)] = pirc(i,j);}}

int nrow = classsize.rows(); 
int ncol = classsize.cols(); 
classsize = 1/(sum(classsize)) * classsize;

for(int i = 0; i < nrow; i++){
	for(int j= 0; j < ncol; j++){Rclasssize[i + j*nrow] = classsize(i,j);}}


}
}


