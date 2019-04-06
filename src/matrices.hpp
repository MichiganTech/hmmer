/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1999 Washington University School of Medicine
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* matrices.c
 *
 * for matrix manipulation
 * many functions borrowed from Elena Rivas matrix package (as noted)
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "structs.hpp"
#include "config.hpp"


/*  Function: AssignMatrixNotLog()
 *
 *  Purpose:  assign a matrix to P_emit_ij based on a time and a Qij
 *
 *  Args:     Qij   - rate matrix
 *        matrix_n  - width or height of matrix
 *        time   - time to use as exponentiation factor for Qij
 *        P_emit_ij - conditional matrix
 *
 *  Return:  (void)
 */
void
AssignMatrixNotLog (
  double *Qij,
  int matrix_n,
  double time,
  double *P_emit_ij);



void
ReadAAMatrices (
  double **ret_Sij,
  double **ret_pi,
  char *matrixfile,
  int environ,
  int L);



/*  Function: SymToRateMatrices ()
 *
 *  Purpose:  Convert a symmetric to a rate matrix
 *
 *  Args:     ret_Qij - pointer to a rate matrix
 *        S       - symmetric matrix
 *          p       - residue frequencies
 *        L        - number of elements in a row or col of the matrix
 *
 *  Return:  (void)
 */
void
SymToRateMatrices(
  double *Qij,
  double *Sij,
  double *pi,
  int L,
  int environ);


/*  Function: NormRateMatrices ()
 *
 *  Purpose:  normalize matrices such that columns sum to  1.
 *
 *  Args:     Qij     - pointer to a rate matrix
 *          pi      - residue frequencies
 *        L        - number of elements in a row or col of the matrix
 *        environ - number of different classes being modeled
 *
 *  Return:  (void)
 */
void
NormRateMatrices(
  double *Qij,
  double *pi,
  int L,
  int environ);


/* Function: Check_Accuracy()
   from Elena Rivas' matrix functions package
   */
int
Check_Accuracy(
  double *vec,
  int L);


/* Function: Cal_Id()
 * Purpose:  Creates a Id(LxL) matrix
 *
 * Args:    L - dimension
 *
 * Returns: Id(LxL) = \delta_{ij}
 *          Id is allocated here, freed by caller.
 *
 */
double *
Cal_Id(
  int L);



/* Function: Cal_M_Exp()
 * Purpose:  Given a matrix M, calculate exp{r*M}
 *
 *            exp{rM} = \sum_{n=0}^{\infty} [ r^n * M^n / n!   ]
 *
 * Args:     M  - LL joint prob matrix (prealloc)
 *
 * Returns:  Q(LxL) = exp{rM(LxL)}
 *           Q is alocated here, freed by caller.
 *
 */
double *
Cal_M_Exp(
  double *M,
  int L,
  double r);


/* Function: Cal_M_N_Prod()
 *
 * Purpose:  given M (LixLk) and N (LkxLj) calculate Mprod(LixLk) = M * N.
 *
 * Args:     M  - LixLk (prealloc)
 *           N  - LkxLj (prealloc)
 *
 * Returns:  Mprod(LixLj) = M(LixLk) * N(LkxLj)
 *           Mprod is alocated here, freed by caller.
 */
double *
Cal_M_N_Prod(
  double *M,
  double *N,
  int Li,
  int Lk,
  int Lj);


void
Comp_M_N_Prod(
  double *M,
  double *N,
  int L
);


void
CopyMatrix (
  double *copyQ,
  double *Q,
  int N);


void
AssignWagMatrix (
  double **ret_Sij,
  double **ret_pi);
