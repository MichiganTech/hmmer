/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* plan9.c
 *
 * alloc, free, and initialization of old Plan9 (HMMER 1.x) functions.
 * Rescued from the wreckage of HMMER 1.9m code.
 */

#pragma once


#include "structs.hpp"

struct plan9_s*
P9AllocHMM(
  int M /* length of model to make */
);


int
P9FreeHMM(
  struct plan9_s *hmm);


/* Function: P9ZeroHMM()
 *
 * Purpose:  Zero emission and transition counts in an HMM.
 */
void
P9ZeroHMM(
  struct plan9_s *hmm);



/* Function: P9Renormalize()
 *
 * Normalize all P distributions so they sum to 1.
 * P distributions that are all 0, or contain negative
 * probabilities, are left untouched.
 *
 * Returns 1 on success, or 0 on failure.
 */
void
P9Renormalize(
  struct plan9_s *hmm);


/* Function: P9DefaultNullModel()
 *
 * Purpose:  Set up a default random sequence model, using
 *           global aafq[]'s for protein or 0.25 for nucleic
 *           acid. randomseq is alloc'ed in caller. Alphabet information
 *           must already be known.
 */
void
P9DefaultNullModel(
  float *null);