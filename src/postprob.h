/************************************************************
 * Copyright (C) 1998 Ian Holmes (ihh@sanger.ac.uk)
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* postprob.h
 *
 * Functions for working with posterior probabilities,
 * including unfussed "backwards" and "optimal accuracy"
 * implementations.
 */

#pragma once

#include "config.h"
#include "structs.h"
#include "funcs.h"

/* Extra algorithms to work with posterior probabilities.
 */

float 
P7OptimalAccuracy(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct p7trace_s **ret_tr);


float 
P7Backward(
  unsigned char *dsq, 
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s **ret_mx);

void  
P7EmitterPosterior(
  int L, 
  struct plan7_s *hmm,
  struct dpmatrix_s *forward,
  struct dpmatrix_s *backward,
  struct dpmatrix_s *mx);

float P7FillOptimalAccuracy(
  int L, 
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);

void  P7OptimalAccuracyTrace(
  int L, 
  int M,
  struct dpmatrix_s *posterior,
  struct dpmatrix_s *mx,
  struct p7trace_s **ret_tr);
