/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* lsjfuncs.h
 * Declarations of external functions used in lsj_eweight.c
 * (Entropy-based sequence weighting)
 *
 * Steve Johnson
 */

#pragma once

#include "structs.hpp"

float
Eweight(
  struct plan7_s *hmm, 
  struct p7prior_s *pri,
  float numb_seqs,
  float entwgt);


void
ModelContent(
	float *ent1,
	float *ent2,
	int M);
