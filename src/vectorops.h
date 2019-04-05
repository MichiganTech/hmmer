/* vectorops.h
 * Header file for vectorops.c
 */

#pragma once

void 
DSet(
  double *vec, 
  int n, 
  double value);


void 
FSet(
  float *vec, 
  int n, 
  float value);


void 
DScale(
  double *vec, 
  int n,
  double scale);


void 
FScale(
  float *vec, 
  int n, 
  float scale);


double 
DSum(
  double *vec, 
  int n);


float 
FSum(
  float *vec, 
  int n);


void 
DAdd(
  double *vec1, 
  double *vec2, 
  int n);


void 
FAdd(
  float *vec1, 
  float *vec2, 
  int n);


void 
DCopy(
  double *vec1, 
  double *vec2, 
  int n);


void 
FCopy(
  float *vec1, 
  float *vec2, 
  int n);


double 
DDot(
  double *vec1, 
  double *vec2, 
  int n);


float 
FDot(
  float *vec1, 
  float *vec2, 
  int n);


double 
DMax(
  double *vec, 
  int n);


float 
FMax(
  float *vec, 
  int n);


double 
DMin(
  double *vec, 
  int n);


float 
FMin(
  float *vec, 
  int n);


double 
DArgMax(
  double *vec, 
  int n);


float 
FArgMax(
  float *vec, 
  int n);


double 
DArgMin(
  double *vec, 
  int n);


float 
FArgMin(
  float *vec, 
  int n);


void 
DNorm(
  double *vec, 
  int n);


void 
FNorm(
  float *vec, 
  int n);


void 
DLog(
  double *vec, 
  int n);


void 
FLog(
  float *vec, 
  int n);


void 
DExp(
  double *vec, 
  int n);


void 
FExp(
  float *vec, 
  int n);


float 
DLogSum(
  double *vec, 
  int n);


float 
FLogSum(
  float *vec, 
  int n);

int
FChoose(
  float *p, 
  int N);


/* Function:  FEntropy(), DEntropy()
 * Synopsis:  Return Shannon entropy of p-vector, in bits.           
 *
 * Purpose:   Returns the Shannon entropy of a probability vector <p>,
 *            in bits ($\log_2$), defined as \citep{CoverThomas}:
 *            
 *            \[
 *               H = - \sum_x p_x \log_2 p_x.
 *            \]
 */
double
DEntropy(
  const double *p, 
  int n);


float
FEntropy(
  const float *p, 
  int n);