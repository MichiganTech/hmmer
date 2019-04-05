#pragma once
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* sre_math.c
 * 
 * Portability for and extensions to C math library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "squid.h"

  
/* Function: Linefit()
 * 
 * Purpose:  Given points x[0..N-1] and y[0..N-1], fit to
 *           a straight line y = a + bx.
 *           a, b, and the linear correlation coefficient r
 *           are filled in for return.
 *           
 * Args:     x     - x values of data
 *           y     - y values of data               
 *           N     - number of data points
 *           ret_a - RETURN: intercept
 *           ret_b - RETURN: slope
 *           ret_r - RETURN: correlation coefficient  
 *           
 * Return:   1 on success, 0 on failure.
 */
int          
Linefit(
  float *x, 
  float *y, 
  int N, 
  float *ret_a, 
  float *ret_b, 
  float *ret_r);


/* Function: WeightedLinefit()
 * 
 * Purpose:  Given points x[0..N-1] and y[0..N-1] with
 *           variances (measurement errors) var[0..N-1],  
 *           fit to a straight line y = mx + b.
 *           
 * Method:   Algorithm from Numerical Recipes in C, [Press88].
 *           
 * Return:   (void)
 *           ret_m contains slope; ret_b contains intercept 
 */                
void
WeightedLinefit(
  float *x, 
  float *y, 
  float *var, 
  int N, 
  float *ret_m, 
  float *ret_b);


/* Function: Gammln()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.  
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
double
Gammln(double x);


/* 2D matrix operations
 */
float **
FMX2Alloc(int rows, int cols);


void
FMX2Free(
  float **mx);


double**
DMX2Alloc(
  int rows, 
  int cols);


void
DMX2Free(
  double **mx);


/* Function: FMX2Multiply()
 * 
 * Purpose:  Matrix multiplication.
 *           Multiply an m x p matrix A by a p x n matrix B,
 *           giving an m x n matrix C.
 *           Matrix C must be a preallocated matrix of the right
 *           size.
 */
void
FMX2Multiply(
  float **A, 
  float **B, 
  float **C, 
  int m, 
  int p, 
  int n);


/* Function: IncompleteGamma()
 * 
 * Purpose:  Returns 1 - P(a,x) where:
 *           P(a,x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt
 *                  = \frac{\gamma(a,x)}{\Gamma(a)}
 *                  = 1 - \frac{\Gamma(a,x)}{\Gamma(a)}
 *                  
 *           Used in a chi-squared test: for a X^2 statistic x
 *           with v degrees of freedom, call:
 *                  p = IncompleteGamma(v/2., x/2.) 
 *           to get the probability p that a chi-squared value
 *           greater than x could be obtained by chance even for
 *           a correct model. (i.e. p should be large, say 
 *           0.95 or more).
 *           
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988. 
 *           
 * Args:     a  - for instance, degrees of freedom / 2     [a > 0]
 *           x  - for instance, chi-squared statistic / 2  [x >= 0] 
 *           
 * Return:   1 - P(a,x).
 */          
double
IncompleteGamma(
  double a, 
  double x);